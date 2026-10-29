"""3D version of the tube-beside-a-plate model, for a tube whose temperature varies
along its length (the measured profiles in astral_temps.csv).

The tube is split into segments along x, each held at the measured temperature at its
centre. Every sheet (shields, plate) is split into an x-z grid of patches, so sheet
temperatures can vary both along the tube and up the sheet. Everything else follows
radiation.py: gray diffuse surfaces, a black environment at T_amb, room-air natural
convection on exposed faces, conduction within each sheet and through the still air
between sheets, and no heat from hot air rising off the tube.
"""

import csv
from dataclasses import dataclass

import numpy as np

import radiation as rad

C = 273.15


# --- measured tube temperatures

def load_tube_profiles(path, centre=0.200):
    """Read astral_temps.csv into {column name: (x [m], T [K])}.

    The file covers one half of the tube, x <= centre; the other half is its mirror image.
    Blank cells are skipped.
    """
    with open(path, encoding="utf-8") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    profiles = {}
    for col, name in enumerate(header[1:], start=1):
        pts = [(float(r[0]) / 1e3, float(r[col]) + C) for r in data if col < len(r) and r[col].strip()]
        x, T = np.array(sorted(pts)).T
        mirrored = 2 * centre - x[::-1]
        keep = mirrored > centre + 1e-9
        profiles[name] = (np.r_[x, mirrored[keep]], np.r_[T, T[::-1][keep]])
    return profiles


def tube_temperature(profile, x):
    """Linear interpolation; beyond the last reading the nearest reading is held."""
    xs, Ts = profile
    return np.interp(x, xs, Ts)


# --- view factors on x-z patch grids

def parallel_patches(x_edges, z_edges, c):
    """View factors between the patches of two identical x-z grids on parallel planes a
    distance c apart. Patch index is ix * n_z + iz. Returns F[from, to] (Howell C-11)."""
    def g(x, y):
        rx, ry = np.sqrt(x * x + c * c), np.sqrt(y * y + c * c)
        return (y * rx * np.arctan(y / rx) + x * ry * np.arctan(x / ry)
                - c * c / 2 * np.log(x * x + y * y + c * c)) / (2 * np.pi)

    ex, ez = np.asarray(x_edges), np.asarray(z_edges)
    dx = ex[:, None, None, None] - ex[None, None, :, None]
    dz = ez[None, :, None, None] - ez[None, None, None, :]
    corners = g(dx, dz)  # axes: x edge (from), z edge (from), x edge (to), z edge (to)
    A1F12 = corners
    for axis in range(4):
        A1F12 = np.diff(A1F12, axis=axis)
    nx, nz = len(ex) - 1, len(ez) - 1
    A1F12 = A1F12.reshape(nx * nz, nx * nz)
    area = np.outer(np.diff(ex), np.diff(ez)).ravel()
    return A1F12 / area[:, None]


def tube_segments_to_patches(radius, tube_x_edges, distance, x_edges, z_edges, n_phi=180, per_segment=5):
    """View factors from tube segments (axis along x, segment i spans tube_x_edges[i:i+2])
    to the patches of an x-z grid in the plane y = distance. Returns (n_segments, n_patches)."""
    n_seg = len(tube_x_edges) - 1
    phi = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    frac = (np.arange(per_segment) + 0.5) / per_segment
    x = (tube_x_edges[:-1, None] + np.diff(tube_x_edges)[:, None] * frac[None, :]).ravel()
    P, X = np.meshgrid(phi, x, indexing="ij")  # (n_phi, n_seg * per_segment)
    P, X = P.ravel(), X.ravel()
    cos_p, sin_p = np.cos(P), np.sin(P)
    points = np.stack([X, radius * cos_p, radius * sin_p], axis=1)
    normals = np.stack([np.zeros_like(P), cos_p, sin_p], axis=1)

    # Each element sees only the part of the plane where cos_p*distance + sin_p*z > radius.
    lo = np.full_like(P, -np.inf)
    hi = np.full_like(P, np.inf)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_cut = (radius - distance * cos_p) / sin_p
    up, down = sin_p > 1e-12, sin_p < -1e-12
    lo[up], hi[down] = z_cut[up], z_cut[down]
    faces_plane = np.where(~(up | down), distance * cos_p > radius, True)

    nx, nz = len(x_edges) - 1, len(z_edges) - 1
    F = np.zeros((len(P), nx * nz))
    for ix in range(nx):
        for iz in range(nz):
            z0 = np.maximum(z_edges[iz], lo)
            z1 = np.minimum(z_edges[iz + 1], hi)
            v = faces_plane & (z1 > z0)
            if not v.any():
                continue
            y = np.full(v.sum(), distance)
            xa, xb = np.full_like(y, x_edges[ix]), np.full_like(y, x_edges[ix + 1])
            verts = np.stack([
                np.stack([xa, y, z0[v]], axis=1),
                np.stack([xb, y, z0[v]], axis=1),
                np.stack([xb, y, z1[v]], axis=1),
                np.stack([xa, y, z1[v]], axis=1),
            ], axis=1)
            F[v, ix * nz + iz] = rad.differential_to_polygon(points[v], normals[v], verts)
    # Average over the elements of each segment (all elements have equal area).
    return F.reshape(n_phi, n_seg, per_segment, -1).mean(axis=(0, 2))


def grid_conductance(k, thickness, x_edges, z_edges):
    """(n_p, n_p) W/K conductance matrix between neighbouring patches of one sheet."""
    dx, dz = np.diff(x_edges), np.diff(z_edges)
    xc, zc = 0.5 * (x_edges[1:] + x_edges[:-1]), 0.5 * (z_edges[1:] + z_edges[:-1])
    nx, nz = len(dx), len(dz)
    idx = np.arange(nx * nz).reshape(nx, nz)
    G = np.zeros((nx * nz, nx * nz))
    g_x = k * thickness * dz[None, :] / np.diff(xc)[:, None]  # (nx - 1, nz)
    g_z = k * thickness * dx[:, None] / np.diff(zc)[None, :]  # (nx, nz - 1)
    G[idx[:-1].ravel(), idx[1:].ravel()] = g_x.ravel()
    G[idx[:, :-1].ravel(), idx[:, 1:].ravel()] = g_z.ravel()
    return G + G.T


# --- the 3D problem

@dataclass
class Result3D:
    n_shields: int
    x_edges: np.ndarray  # m
    z_edges: np.ndarray  # m
    T_sheets: np.ndarray  # (n_sheets, n_x, n_z) K, from the tube outward
    tube_x: np.ndarray  # tube segment centres, m
    tube_T: np.ndarray  # K
    Q_radiation: float  # W, tube's net radiation

    @property
    def x(self):
        return 0.5 * (self.x_edges[1:] + self.x_edges[:-1])

    @property
    def z(self):
        return 0.5 * (self.z_edges[1:] + self.z_edges[:-1])

    # "plate" means the outermost sheet: the plate, or the last shield if there is no plate.
    @property
    def T_plate_max(self):
        return self.T_sheets[-1].max()

    @property
    def T_plate_mean(self):
        w = np.outer(np.diff(self.x_edges), np.diff(self.z_edges))
        return np.average(self.T_sheets[-1], weights=w)


def solve3d(setup: rad.Setup, profile, n_shields=0, plate=True, tube_x=(0.100, 0.300),
            plate_x0=0.130, n_x=20, n_z=30, n_segments=40) -> Result3D:
    """Sheets span x in [plate_x0, plate_x0 + setup.length]. The tube spans tube_x, with
    temperatures from `profile` (see load_tube_profiles); setup.T_tube is not used."""
    if not plate and n_shields < 1:
        raise ValueError("need at least one shield when there is no plate")
    s = setup
    r = s.tube_od / 2
    spacing = s.gap / (n_shields + 1 if plate else n_shields)

    x_edges = np.linspace(plate_x0, plate_x0 + s.length, n_x + 1)
    z_edges = rad.strip_edges(s.plate_height, n_z)
    A_patch = np.outer(np.diff(x_edges), np.diff(z_edges)).ravel()
    seg_edges = np.linspace(*tube_x, n_segments + 1)
    seg_x = 0.5 * (seg_edges[1:] + seg_edges[:-1])
    seg_T = tube_temperature(profile, seg_x)
    A_seg = np.pi * s.tube_od * np.diff(seg_edges)

    F_src = tube_segments_to_patches(r, seg_edges, r + spacing, x_edges, z_edges)
    F_gap = parallel_patches(x_edges, z_edges, spacing)
    sheet_eps = [(s.eps_shield, s.eps_shield)] * n_shields + [(s.eps_plate, s.eps_plate)] * plate
    G_shield = grid_conductance(s.k_shield, s.shield_thickness, x_edges, z_edges)
    G_plate = grid_conductance(s.k_plate, s.plate_thickness, x_edges, z_edges)

    T, Q = rad.solve_stack(
        seg_T, A_seg, s.eps_tube, F_src, A_patch, sheet_eps, F_gap, s.T_amb,
        G_lateral=[G_shield] * n_shields + [G_plate] * plate,
        G_gap=rad.k_air(473.15) * A_patch / spacing,
        h_exposed=lambda T: rad.h_vertical_plate(T, s.T_amb, s.plate_height),
    )
    return Result3D(n_shields, x_edges, z_edges, T.reshape(len(sheet_eps), n_x, n_z), seg_x, seg_T, Q)


def shields_required3d(setup, profile, T_limit=150 + C, n_max=6, plate=True, **kw):
    """Smallest number of shields that keeps the outermost sheet's peak below T_limit, or None."""
    for n in range(0 if plate else 1, n_max + 1):
        if solve3d(setup, profile, n, plate=plate, **kw).T_plate_max < T_limit:
            return n
    return None
