"""Steady-state model of a hot horizontal tube beside a vertical plate, with optional
radiation shields in the gap between them (issue #1).

Geometry (z up):
    tube    - axis along x, outer radius r, centred at y = 0, z = 0
    plate   - thin vertical sheet in the plane y = r + gap, spanning x in [0, L]
              and z in [-H/2, H/2], i.e. centred on the tube's height
    shields - thin sheets the same size as the plate, evenly spaced in the gap

Each sheet (shield or plate) is split into horizontal strips, so temperature can
vary with height; with one strip per sheet the model is lumped. Temperature is
uniform along x and through each sheet's thickness. Surfaces are gray and
diffuse; everything not seen by a surface is a black environment at T_amb.

Heat carried by the tube's rising hot air onto the sheets is ignored. Sheets lose
heat to the surrounding air (at T_air, by default T_amb) by natural convection on
their exposed faces, conduct heat along their height, and exchange heat through the
still air in the gaps between them.

Inside a closed box of still air (issue #4), T_amb is the box's wall temperature and
T_air its air temperature. The walls are taken as black: a box much larger than the
sheets acts as a cavity, so its effective emissivity seen from the sheets is close
to 1 whatever the wall finish.
All units SI, temperatures in kelvin.
"""

from dataclasses import dataclass

import numpy as np

SIGMA = 5.670374419e-8
G = 9.81
PR_AIR = 0.70


# --- air properties: power-law fits to Incropera Table A.4 (within ~2 % for 300-900 K)

def k_air(T):
    return 0.0263 * (T / 300.0) ** 0.80


def nu_air(T):
    return 1.589e-5 * (T / 300.0) ** 1.73


def h_vertical_plate(T_s, T_inf, height):
    """Average h from the Churchill-Chu correlation (full Ra range)."""
    T_f = 0.5 * (T_s + T_inf)
    nu = nu_air(T_f)
    Ra = G / T_f * np.abs(T_s - T_inf) * height**3 / (nu * nu / PR_AIR)
    Nu = (0.825 + 0.387 * Ra ** (1 / 6) / (1 + (0.492 / PR_AIR) ** (9 / 16)) ** (8 / 27)) ** 2
    return Nu * k_air(T_f) / height


# --- view factors

def parallel_rectangles(a, b, c):
    """Two directly opposed, aligned a x b rectangles a distance c apart."""
    X, Y = a / c, b / c
    x1, y1 = np.sqrt(1 + X**2), np.sqrt(1 + Y**2)
    return 2 / (np.pi * X * Y) * (
        np.log(x1 * y1 / np.sqrt(1 + X**2 + Y**2))
        + X * y1 * np.arctan(X / y1)
        + Y * x1 * np.arctan(Y / x1)
        - X * np.arctan(X)
        - Y * np.arctan(Y)
    )


def parallel_strips(length, z_edges_1, z_edges_2, c):
    """View factors between horizontal strips on two parallel planes a distance c apart.

    Both planes span x in [0, length]; strip i of plane 1 covers z_edges_1[i:i+2].
    Returns F[i, j] from strip i on plane 1 to strip j on plane 2 (Howell catalog C-11).
    """
    def g(x, y):
        rx, ry = np.sqrt(x * x + c * c), np.sqrt(y * y + c * c)
        return (y * rx * np.arctan(y / rx) + x * ry * np.arctan(x / ry)
                - c * c / 2 * np.log(x * x + y * y + c * c)) / (2 * np.pi)

    z1 = np.asarray(z_edges_1)[:, None]
    z2 = np.asarray(z_edges_2)[None, :]
    # Both strips share x in [0, length], so the sum over x corners collapses.
    corners = 2 * g(length, z1 - z2) - 2 * g(0.0, z1 - z2)
    A1F12 = -(corners[1:, 1:] - corners[1:, :-1] - corners[:-1, 1:] + corners[:-1, :-1])
    return A1F12 / (length * np.diff(z_edges_1)[:, None])


def differential_to_polygon(points, normals, vertices):
    """View factor from differential areas to planar polygons by contour integration.

    points, normals: (n, 3). vertices: (n, m, 3), polygon corners in order. The
    polygon must lie entirely in front of each differential area.
    """
    R = vertices - points[:, None, :]
    R_next = np.roll(R, -1, axis=1)
    cross = np.cross(R, R_next)
    cross_norm = np.linalg.norm(cross, axis=2)
    cos_t = np.einsum("nmk,nmk->nm", R, R_next) / (
        np.linalg.norm(R, axis=2) * np.linalg.norm(R_next, axis=2)
    )
    theta = np.arccos(np.clip(cos_t, -1.0, 1.0))
    with np.errstate(invalid="ignore", divide="ignore"):
        g = np.where(cross_norm > 0, theta / cross_norm, 0.0)
    total = np.einsum("nk,nmk,nm->n", normals, cross, g)
    return np.abs(total) / (2 * np.pi)


def tube_to_rectangle(radius, length, distance, z_lo, z_hi, n_phi=360, n_x=80):
    """View factor from a tube (axis along x, x in [0, length]) to the rectangle in the
    plane y = distance (from the axis), x in [0, length], z in [z_lo, z_hi]."""
    phi = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    x = (np.arange(n_x) + 0.5) * length / n_x
    P, X = np.meshgrid(phi, x, indexing="ij")
    P, X = P.ravel(), X.ravel()
    cos_p, sin_p = np.cos(P), np.sin(P)
    points = np.stack([X, radius * cos_p, radius * sin_p], axis=1)
    normals = np.stack([np.zeros_like(P), cos_p, sin_p], axis=1)

    # Part of the rectangle in front of each element: cos_p*distance + sin_p*z > radius.
    lo = np.full_like(P, z_lo)
    hi = np.full_like(P, z_hi)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_cut = (radius - distance * cos_p) / sin_p
    up, down = sin_p > 1e-12, sin_p < -1e-12
    lo[up] = np.maximum(lo[up], z_cut[up])
    hi[down] = np.minimum(hi[down], z_cut[down])
    side = ~(up | down)
    visible = np.where(side, distance * cos_p > radius, hi > lo)

    F = np.zeros_like(P)
    v = visible
    y = np.full(v.sum(), distance)
    verts = np.stack([
        np.stack([np.zeros_like(y), y, lo[v]], axis=1),
        np.stack([np.full_like(y, length), y, lo[v]], axis=1),
        np.stack([np.full_like(y, length), y, hi[v]], axis=1),
        np.stack([np.zeros_like(y), y, hi[v]], axis=1),
    ], axis=1)
    F[v] = differential_to_polygon(points[v], normals[v], verts)
    return F.mean()


# --- radiation/conduction network for a source facing a stack of parallel sheets

def solve_stack(T_src, A_src, eps_src, F_src, A_strip, sheet_eps, F_gap, T_amb,
                G_lateral=None, G_gap=None, h_exposed=None, T_last=None, T_air=None,
                tol=1e-8):
    """Solve for strip temperatures on a stack of sheets facing a hot source.

    F_src:     (n_k,) view factors from the source to the first sheet's strips.
    A_strip:   (n_k,) strip areas, the same on every sheet.
    sheet_eps: [(eps_front, eps_back), ...] from the source outward; the last is the plate.
    F_gap:     (n_k, n_k) strip-to-strip view factors from one sheet's back to the next
               sheet's front (equal spacing assumed).
    G_lateral: per sheet, (n_k - 1,) W/K conductances between neighbouring strips.
    G_gap:     (n_k,) W/K conductance through the air between facing strips.
    h_exposed: h(T) for faces open to the air (first sheet's front, last sheet's back).
    T_last:    fix the last sheet's temperature instead of solving for it.
    T_air:     air temperature for h_exposed; defaults to T_amb.

    Returns (T, Q_src): T is (n_sheets, n_k) in K, Q_src the source's net radiation in W.
    """
    A_strip = np.atleast_1d(np.asarray(A_strip, dtype=float))
    T_air = T_amb if T_air is None else T_air
    F_src = np.atleast_1d(F_src)
    F_gap = np.atleast_2d(F_gap)
    n_k, n_s = len(A_strip), len(sheet_eps)
    n_u = n_s * n_k
    m = 1 + 2 * n_u
    front = lambda j: slice(1 + 2 * j * n_k, 1 + (2 * j + 1) * n_k)
    back = lambda j: slice(1 + (2 * j + 1) * n_k, 1 + (2 * j + 2) * n_k)

    F = np.zeros((m, m))
    F[0, front(0)] = F_src
    F[front(0), 0] = A_src * F_src / A_strip
    for j in range(n_s - 1):
        F[back(j), front(j + 1)] = F_gap
        F[front(j + 1), back(j)] = (A_strip[:, None] * F_gap).T / A_strip[:, None]
    F_amb = 1 - F.sum(axis=1)
    eps = np.concatenate([[eps_src]] + [np.r_[np.full(n_k, a), np.full(n_k, b)] for a, b in sheet_eps])
    area = np.concatenate([[A_src], np.tile(A_strip, 2 * n_s)])
    Eb_amb = SIGMA * T_amb**4

    # Net radiation leaving each surface is linear in E_b = sigma T^4: q = R E_b + q0.
    M = np.eye(m) - (1 - eps)[:, None] * F
    R = area[:, None] * (np.eye(m) - F) @ np.linalg.inv(M)
    q0 = R @ ((1 - eps) * F_amb * Eb_amb) - area * F_amb * Eb_amb
    R = R * eps[None, :]

    # P maps unknown strip temperatures onto the front and back surfaces they set.
    P = np.zeros((m - 1, n_u))
    for j in range(n_s):
        block = slice(j * n_k, (j + 1) * n_k)
        P[2 * j * n_k:(2 * j + 1) * n_k, block] = np.eye(n_k)
        P[(2 * j + 1) * n_k:(2 * j + 2) * n_k, block] = np.eye(n_k)

    # Linear conduction network (W/K): along each sheet and across the air gaps.
    K = np.zeros((n_u, n_u))

    def link(a, b, g):
        K[a, a] -= g
        K[b, b] -= g
        K[a, b] += g
        K[b, a] += g

    for j in range(n_s):
        if G_lateral is not None:
            for k in range(n_k - 1):
                link(j * n_k + k, j * n_k + k + 1, G_lateral[j][k])
        if j < n_s - 1 and G_gap is not None:
            for k in range(n_k):
                link(j * n_k + k, (j + 1) * n_k + k, G_gap[k])

    exposed = np.zeros(n_u)
    exposed[:n_k] += A_strip
    exposed[-n_k:] += A_strip

    def evaluate(T):
        Eb = SIGMA * np.concatenate([[T_src], P @ T]) ** 4
        q = R @ Eb + q0
        res = -P.T @ q[1:] + K @ T
        dEb = 4 * SIGMA * (P @ T) ** 3
        jac = -P.T @ (R[1:, 1:] * dEb[None, :]) @ P + K
        if h_exposed is not None:
            h = h_exposed(T)
            dh = (h_exposed(T + 1e-3) - h) / 1e-3
            res += exposed * h * (T_air - T)
            jac -= np.diag(exposed * (h + dh * (T - T_air)))
        return res, jac, q[0]

    T = np.full(n_u, 0.5 * (T_src + max(T_amb, T_air)))
    free = np.ones(n_u, dtype=bool)
    if T_last is not None:
        free[-n_k:] = False
        T[~free] = T_last
    for _ in range(200 if free.any() else 0):
        res, jac, _ = evaluate(T)
        step = np.clip(np.linalg.solve(jac[np.ix_(free, free)], -res[free]), -50, 50)
        T[free] = np.maximum(T[free] + step, 0.5 * min(T_amb, T_air))
        if np.max(np.abs(step)) < tol:
            break
    else:
        if free.any():
            raise RuntimeError("sheet temperatures did not converge")

    return T.reshape(n_s, n_k), evaluate(T)[2]


# --- the issue #1 problem

@dataclass
class Setup:
    tube_od: float = 0.75 * 0.0254
    length: float = 0.140
    plate_height: float = 0.220
    gap: float = 0.010
    T_tube: float = 650 + 273.15
    T_amb: float = 25 + 273.15  # surroundings seen by radiation (box walls)
    T_air: float | None = None  # air around the sheets; None means T_amb
    eps_tube: float = 0.8
    eps_plate: float = 0.8
    eps_shield: float = 0.2
    plate_thickness: float = 0.001
    k_plate: float = 45.0  # carbon steel
    shield_thickness: float = 0.0005
    k_shield: float = 16.0  # stainless steel


@dataclass
class Result:
    n_shields: int
    z_edges: np.ndarray  # strip edges, m
    T_sheets: np.ndarray  # (n_shields + 1, n_strips) K, shields then plate
    Q_radiation: float  # W, tube's net radiation

    @property
    def z(self):
        return 0.5 * (self.z_edges[1:] + self.z_edges[:-1])

    @property
    def T_plate_max(self):
        return self.T_sheets[-1].max()

    @property
    def T_shield_max(self):
        """Hottest shield temperature, or nan with no shields."""
        return self.T_sheets[:-1].max() if self.n_shields else np.nan

    @property
    def T_plate_mean(self):
        return np.average(self.T_sheets[-1], weights=np.diff(self.z_edges))


def strip_edges(height, n_strips):
    """Strip edges graded finer near z = 0, where the tube's radiation peaks."""
    u = np.linspace(-1, 1, n_strips + 1)
    return height / 2 * np.sign(u) * np.abs(u) ** 1.5


def solve(setup: Setup, n_shields: int = 0, n_strips: int = 60) -> Result:
    """n_strips = 1 gives the lumped model (each sheet at one temperature)."""
    s = setup
    T_air = s.T_amb if s.T_air is None else s.T_air
    r = s.tube_od / 2
    spacing = s.gap / (n_shields + 1)
    A_tube = np.pi * s.tube_od * s.length

    edges = strip_edges(s.plate_height, n_strips)
    z = 0.5 * (edges[1:] + edges[:-1])
    A_strip = s.length * np.diff(edges)

    F_src = np.array([tube_to_rectangle(r, s.length, r + spacing, lo, hi)
                      for lo, hi in zip(edges[:-1], edges[1:])])
    F_gap = parallel_strips(s.length, edges, edges, spacing)
    sheet_eps = [(s.eps_shield, s.eps_shield)] * n_shields + [(s.eps_plate, s.eps_plate)]

    dz = np.diff(z)
    g_shield = s.k_shield * s.shield_thickness * s.length / dz
    g_plate = s.k_plate * s.plate_thickness * s.length / dz
    # Air in the narrow gaps between sheets is treated as still (conduction, Nu = 1),
    # evaluated at a 200 C film temperature, or at the air temperature if hotter.
    G_gap = k_air(max(473.15, T_air)) * A_strip / spacing

    T, Q = solve_stack(
        s.T_tube, A_tube, s.eps_tube, F_src, A_strip, sheet_eps, F_gap, s.T_amb,
        G_lateral=[g_shield] * n_shields + [g_plate],
        G_gap=G_gap,
        h_exposed=lambda T: h_vertical_plate(T, T_air, s.plate_height),
        T_air=T_air,
    )
    return Result(n_shields, edges, T, Q)


def shields_required(setup: Setup, T_limit=150 + 273.15, n_max=10, use_peak=True, **kw):
    """Smallest number of shields that keeps the plate below T_limit, or None."""
    for n in range(n_max + 1):
        res = solve(setup, n, **kw)
        if (res.T_plate_max if use_peak else res.T_plate_mean) < T_limit:
            return n
    return None


def box_sweep(T_box, n_shields=(0, 1, 2, 3), eps_shield=(0.2, 0.45), setup=None, **kw):
    """Peak plate and shield temperatures inside a box of still air (issue #4).

    The box's air and walls are both at each temperature in T_box (K). Returns a list
    of dicts with keys T_box, eps_shield, n_shields, T_plate_max and T_shield_max (K).
    """
    from dataclasses import replace

    base = setup or Setup()
    rows = []
    for Tb in T_box:
        for eps in eps_shield:
            for n in n_shields:
                if n == 0 and eps != eps_shield[0]:
                    continue
                res = solve(replace(base, T_amb=Tb, T_air=Tb, eps_shield=eps), n, **kw)
                rows.append(dict(T_box=Tb, eps_shield=eps if n else np.nan, n_shields=n,
                                 T_plate_max=res.T_plate_max, T_shield_max=res.T_shield_max))
    return rows
