"""3D finite-element model of the tube / plate / shield setup (scikit-fem).

Each sheet (shields, plate) is a solid hexahedral mesh, so conduction runs in all three
directions, including through the thickness. The sheets' broad faces exchange heat with:

- the tube, by gray-diffuse radiation; the tube is split into segments along x, each held
  at the measured temperature (astral_temps.csv) at its centre;
- each other, by gray-diffuse radiation and by conduction through the still air between them;
- the room, as a black environment at T_amb, plus Churchill-Chu natural convection on the
  two faces open to room air (the first sheet's front, the last sheet's back).

Every quad on a broad face is one radiation surface. View factors come from
radiation3d.py (tube-to-patch contour integration, Howell C-11 for parallel patches).
Heat carried by hot air rising off the tube, and the thin sheet edges, are ignored.

Sheet mid-planes sit where the thin-sheet models put them: shields evenly spaced between
the tube surface and the plate's front face, which is 10 mm from the tube. With no plate,
the outermost shield's mid-plane takes that 10 mm position.
"""

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import Basis, ElementHex1, MeshHex, asm
from skfem.models.poisson import laplace

import radiation as rad
import radiation3d as r3

C = 273.15


@dataclass
class Sheet:
    y0: float  # front face (tube side), m from the tube axis
    thickness: float
    k: float
    eps: float


@dataclass
class FEAResult:
    n_shields: int
    plate: bool
    x: np.ndarray  # node coordinates, m
    y: list  # per sheet, node y coordinates
    z: np.ndarray
    T: list  # per sheet, nodal temperatures (n_x + 1, n_y + 1, n_z + 1), K
    tube_x: np.ndarray
    tube_T: np.ndarray
    Q_radiation: float  # W, tube's net radiation
    face_mean: list  # per sheet, area-weighted mean of its broad-face temperatures, K

    # "plate" means the outermost sheet: the plate, or the last shield if there is no plate.
    @property
    def T_plate_max(self):
        return self.T[-1].max()

    @property
    def T_plate_mean(self):
        return self.face_mean[-1]


def _sheet_mesh(x, y, z):
    mesh = MeshHex.init_tensor(x, y, z)
    ix = np.abs(mesh.p[0][:, None] - x[None, :]).argmin(axis=1)
    iy = np.abs(mesh.p[1][:, None] - y[None, :]).argmin(axis=1)
    iz = np.abs(mesh.p[2][:, None] - z[None, :]).argmin(axis=1)
    node = np.empty((len(x), len(y), len(z)), dtype=int)
    node[ix, iy, iz] = np.arange(mesh.p.shape[1])
    return mesh, node


def _face_average(node, iy, offset):
    """Sparse (rows, cols) averaging the 4 corner nodes of each quad on face y index iy."""
    nx, nz = node.shape[0] - 1, node.shape[2] - 1
    a, b = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    a, b = a.ravel(), b.ravel()
    cols = np.stack([node[a, iy, b], node[a + 1, iy, b], node[a, iy, b + 1], node[a + 1, iy, b + 1]], axis=1)
    return np.repeat(np.arange(nx * nz), 4), cols.ravel() + offset


def layout(setup: rad.Setup, n_shields, plate=True):
    """Sheets from the tube outward."""
    s = setup
    r = s.tube_od / 2
    n_pos = n_shields + 1 if plate else n_shields
    sheets = []
    for k in range(1, n_shields + 1):
        mid = r + k * s.gap / n_pos
        sheets.append(Sheet(mid - s.shield_thickness / 2, s.shield_thickness, s.k_shield, s.eps_shield))
    if plate:
        sheets.append(Sheet(r + s.gap, s.plate_thickness, s.k_plate, s.eps_plate))
    return sheets


def solve_fea(setup: rad.Setup, profile, n_shields=0, plate=True, tube_x=(0.100, 0.300),
              plate_x0=0.130, n_x=20, n_z=30, n_y=2, n_segments=40, T_outside=None, tol=1e-6) -> FEAResult:
    """Sheets span x in [plate_x0, plate_x0 + setup.length] and the plate's full height.
    The tube spans tube_x with temperatures from `profile` (radiation3d.load_tube_profiles).

    T_outside: if given, the outermost sheet's back face radiates and convects to air and
    surroundings at this temperature instead of setup.T_amb (the sheets form a wall between
    a box interior at setup.T_amb and the outside). Everything else sees setup.T_amb."""
    if not plate and n_shields < 1:
        raise ValueError("need at least one shield when there is no plate")
    s = setup
    r = s.tube_od / 2
    sheets = layout(s, n_shields, plate)
    n_s = len(sheets)

    x = np.linspace(plate_x0, plate_x0 + s.length, n_x + 1)
    z = rad.strip_edges(s.plate_height, n_z)
    n_p = n_x * n_z
    A_p = np.outer(np.diff(x), np.diff(z)).ravel()

    seg_edges = np.linspace(*tube_x, n_segments + 1)
    seg_x = 0.5 * (seg_edges[1:] + seg_edges[:-1])
    seg_T = r3.tube_temperature(profile, seg_x)
    A_seg = np.pi * s.tube_od * np.diff(seg_edges)
    n_f = n_segments

    # --- conduction in each sheet (FEM), and the operator averaging nodes onto face quads
    blocks, rows, cols, ys, nodes, offset = [], [], [], [], [], 0
    for j, sh in enumerate(sheets):
        y = sh.y0 + np.linspace(0, sh.thickness, n_y + 1)
        mesh, node = _sheet_mesh(x, y, z)
        blocks.append(sh.k * asm(laplace, Basis(mesh, ElementHex1())))
        for face, iy in ((2 * j, 0), (2 * j + 1, n_y)):
            rr, cc = _face_average(node, iy, offset)
            rows.append(rr + face * n_p)
            cols.append(cc)
        ys.append(y)
        nodes.append((offset, node))
        offset += mesh.p.shape[1]
    K = sp.block_diag(blocks, format="csr")
    n_faces = 2 * n_s * n_p
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    Mf = sp.csr_matrix((np.full(len(rows), 0.25), (rows, cols)), shape=(n_faces, offset))

    # --- radiation network: tube segments, then the front and back faces of each sheet
    m = n_f + n_faces
    front = lambda j: slice(n_f + 2 * j * n_p, n_f + (2 * j + 1) * n_p)
    back = lambda j: slice(n_f + (2 * j + 1) * n_p, n_f + (2 * j + 2) * n_p)
    F = np.zeros((m, m))
    F_src = r3.tube_segments_to_patches(r, seg_edges, sheets[0].y0, x, z)
    F[:n_f, front(0)] = F_src
    F[front(0), :n_f] = (A_seg[:, None] * F_src).T / A_p[:, None]
    gaps = []
    for j in range(n_s - 1):
        c = sheets[j + 1].y0 - (sheets[j].y0 + sheets[j].thickness)
        F_gap = r3.parallel_patches(x, z, c)
        F[back(j), front(j + 1)] = F_gap
        F[front(j + 1), back(j)] = (A_p[:, None] * F_gap).T / A_p[:, None]
        gaps.append(c)
    F_amb = 1 - F.sum(axis=1)
    eps = np.concatenate([np.full(n_f, s.eps_tube)] + [np.full(2 * n_p, sh.eps) for sh in sheets])
    area = np.concatenate([A_seg, np.tile(A_p, 2 * n_s)])
    T_out = s.T_amb if T_outside is None else T_outside
    T_env = np.full(m, s.T_amb)
    T_env[back(n_s - 1)] = T_out
    Eb_amb = rad.SIGMA * T_env**4
    Mrad = np.eye(m) - (1 - eps)[:, None] * F
    R = area[:, None] * (np.eye(m) - F) @ np.linalg.inv(Mrad)
    q0 = R @ ((1 - eps) * F_amb * Eb_amb) - area * F_amb * Eb_amb
    R = R * eps[None, :]
    R_ss = R[n_f:, n_f:]
    Eb_src = rad.SIGMA * seg_T**4
    q_fixed = R[n_f:, :n_f] @ Eb_src + q0[n_f:]

    # --- face-to-face conduction through the gap air, and which faces see room air
    gr, gc, gv = [], [], []
    for j, c in enumerate(gaps):
        g = rad.k_air(473.15) * A_p / c
        a = np.arange(n_p) + (2 * j + 1) * n_p
        b = np.arange(n_p) + (2 * j + 2) * n_p
        gr += [a, b, a, b]
        gc += [a, b, b, a]
        gv += [g, g, -g, -g]
    G = (sp.csr_matrix((np.concatenate(gv), (np.concatenate(gr), np.concatenate(gc))), shape=(n_faces, n_faces))
         if gaps else sp.csr_matrix((n_faces, n_faces)))
    exposed = np.zeros(n_faces)
    exposed[:n_p] = A_p
    exposed[-n_p:] = A_p
    T_air = T_env[n_f:]
    h_of = lambda T: rad.h_vertical_plate(T, T_air, s.plate_height)

    # --- Newton on K T + Mf^T q_out(Mf T) = 0; GMRES with a sparse preconditioner
    T = np.full(offset, s.T_amb + 100.0)
    for _ in range(100):
        Tf = Mf @ T
        h = h_of(Tf)
        dh = (h_of(Tf + 1e-3) - h) / 1e-3
        q = R_ss @ (rad.SIGMA * Tf**4) + q_fixed + exposed * h * (Tf - T_air) + G @ Tf
        res = K @ T + Mf.T @ q
        d_rad = 4 * rad.SIGMA * Tf**3
        d_conv = exposed * (h + dh * (Tf - T_air))
        lu = spla.splu((K + Mf.T @ (sp.diags(np.diag(R_ss) * d_rad + d_conv) + G) @ Mf).tocsc())

        def matvec(v):
            w = Mf @ v
            return K @ v + Mf.T @ (R_ss @ (d_rad * w) + d_conv * w + G @ w)

        J = spla.LinearOperator((offset, offset), matvec=matvec)
        step, info = spla.gmres(J, -res, M=spla.LinearOperator((offset, offset), matvec=lu.solve),
                                rtol=1e-10, atol=0.0, maxiter=200)
        if info != 0:
            raise RuntimeError(f"GMRES did not converge (info={info})")
        step = np.clip(step, -50, 50)
        T = np.maximum(T + step, 0.5 * min(s.T_amb, T_out))
        if np.max(np.abs(step)) < tol:
            break
    else:
        raise RuntimeError("Newton iteration did not converge")

    Tf = Mf @ T
    Q = (R[:n_f, :n_f] @ Eb_src + R[:n_f, n_f:] @ (rad.SIGMA * Tf**4) + q0[:n_f]).sum()
    T_sheets = [T[o + node] for o, node in nodes]
    face_mean = [np.average(Tf[2 * j * n_p:(2 * j + 2) * n_p], weights=np.tile(A_p, 2)) for j in range(n_s)]
    return FEAResult(n_shields, plate, x, ys, z, T_sheets, seg_x, seg_T, Q, face_mean)


def shields_required_fea(setup, profile, T_limit=150 + C, n_max=6, plate=True, **kw):
    """Smallest number of shields that keeps the outermost sheet's peak below T_limit, or None."""
    for n in range(0 if plate else 1, n_max + 1):
        if solve_fea(setup, profile, n, plate=plate, **kw).T_plate_max < T_limit:
            return n
    return None
