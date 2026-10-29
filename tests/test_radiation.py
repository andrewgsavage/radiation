import numpy as np
import pytest

import radiation as rad


def test_differential_to_parallel_rectangle_corner():
    # Element parallel to an a x b rectangle, facing one corner, distance c (Incropera Table 13.1 form).
    a, b, c = 0.3, 0.5, 0.2
    A, B = a / c, b / c
    expected = (A / np.sqrt(1 + A**2) * np.arctan(B / np.sqrt(1 + A**2))
                + B / np.sqrt(1 + B**2) * np.arctan(A / np.sqrt(1 + B**2))) / (2 * np.pi)
    verts = np.array([[[0, 0, c], [a, 0, c], [a, b, c], [0, b, c]]], dtype=float)
    F = rad.differential_to_polygon(np.zeros((1, 3)), np.array([[0, 0, 1.0]]), verts)
    assert F[0] == pytest.approx(expected, rel=1e-10)


def test_parallel_rectangles_limits():
    assert rad.parallel_rectangles(1e3, 1e3, 1.0) == pytest.approx(1.0, abs=2e-3)
    # Incropera Fig 13.4: X = Y = 1 gives about 0.2.
    assert rad.parallel_rectangles(1.0, 1.0, 1.0) == pytest.approx(0.1998, abs=1e-3)


def test_parallel_strips_match_aligned_rectangles():
    length, height, c = 0.14, 0.22, 0.004
    whole = rad.parallel_strips(length, [-height / 2, height / 2], [-height / 2, height / 2], c)
    assert whole[0, 0] == pytest.approx(rad.parallel_rectangles(length, height, c), rel=1e-9)

    # Splitting into strips must conserve A1 F12 for the whole pair.
    edges = rad.strip_edges(height, 25)
    F = rad.parallel_strips(length, edges, edges, c)
    A = length * np.diff(edges)
    assert (A[:, None] * F).sum() == pytest.approx(length * height * whole[0, 0], rel=1e-9)
    assert np.allclose(A[:, None] * F, (A[:, None] * F).T)


def test_tube_strips_sum_to_whole_rectangle():
    r, length, c, height = 0.0095, 0.14, 0.0195, 0.22
    edges = rad.strip_edges(height, 12)
    parts = sum(rad.tube_to_rectangle(r, length, c, lo, hi) for lo, hi in zip(edges[:-1], edges[1:]))
    whole = rad.tube_to_rectangle(r, length, c, -height / 2, height / 2)
    assert parts == pytest.approx(whole, rel=1e-3)


@pytest.mark.parametrize("half_width", [0.02, 0.11, 0.5])
def test_long_tube_matches_2d_strip(half_width):
    # Infinite cylinder to a symmetric strip: F = atan(w / c) / pi.
    r, c = 0.0095, 0.0195
    F = rad.tube_to_rectangle(r, 200.0, c, -half_width, half_width, n_phi=720, n_x=400)
    assert F == pytest.approx(np.arctan(half_width / c) / np.pi, rel=5e-3)


@pytest.mark.parametrize("n", [1, 2, 4])
def test_shields_between_large_plates_cut_heat_to_one_over_n_plus_one(n):
    side, gap, eps = 1e4, 0.01, 0.5
    T1, T2 = 900.0, 300.0
    A = side**2
    F = rad.parallel_rectangles(side, side, gap)

    def heat(n_shields):
        _, Q = rad.solve_stack(T1, A, eps, F, A, [(eps, eps)] * (n_shields + 1),
                               F, T2, T_last=T2)
        return Q

    assert heat(n) / heat(0) == pytest.approx(1 / (n + 1), rel=5e-3)


def test_single_shield_temperature():
    side, gap = 100.0, 0.01
    A = side**2
    F = rad.parallel_rectangles(side, side, gap)
    T, _ = rad.solve_stack(900.0, A, 0.5, F, A, [(0.5, 0.5)] * 2, F, 300.0, T_last=300.0)
    assert T[0] == pytest.approx(((900.0**4 + 300.0**4) / 2) ** 0.25, rel=2e-3)


def test_more_shields_cool_the_plate():
    s = rad.Setup()
    temps = [rad.solve(s, n, n_strips=20).T_plate_max for n in range(4)]
    assert all(a > b for a, b in zip(temps, temps[1:]))


def test_lumped_plate_energy_balance():
    s = rad.Setup()
    res = rad.solve(s, 0, n_strips=1)
    assert res.T_sheets.shape == (1, 1)
    assert s.T_amb < res.T_plate_max < s.T_tube


def test_air_defaults_to_surroundings():
    a = rad.solve(rad.Setup(), 1, n_strips=10)
    b = rad.solve(rad.Setup(T_air=rad.Setup().T_amb), 1, n_strips=10)
    assert np.allclose(a.T_sheets, b.T_sheets)


def test_hotter_box_air_warms_the_plate():
    s = rad.Setup()
    cool = rad.solve(s, 1, n_strips=10)
    warm_air = rad.solve(rad.Setup(T_air=s.T_amb + 100), 1, n_strips=10)
    warm_box = rad.solve(rad.Setup(T_amb=s.T_amb + 100, T_air=s.T_amb + 100), 1, n_strips=10)
    assert cool.T_plate_max < warm_air.T_plate_max < warm_box.T_plate_max


def test_cold_tube_in_hot_box_reaches_box_temperature():
    # With no heat source the sheets must settle at the box temperature.
    T_box = 200 + 273.15
    s = rad.Setup(T_tube=T_box, T_amb=T_box, T_air=T_box)
    res = rad.solve(s, 2, n_strips=8)
    assert np.allclose(res.T_sheets, T_box, atol=1e-6)


def test_plate_runs_above_a_hot_box():
    T_box = 150 + 273.15
    s = rad.Setup(T_amb=T_box, T_air=T_box)
    res = rad.solve(s, 3, n_strips=20)
    assert T_box < res.T_plate_max < res.T_shield_max < s.T_tube
    assert rad.shields_required(s, n_max=3, n_strips=20) is None


def test_box_sweep_rows():
    rows = rad.box_sweep([400.0, 450.0], n_shields=(0, 1), eps_shield=(0.2, 0.45), n_strips=8)
    # One unshielded row plus one per emissivity, for each box temperature.
    assert len(rows) == 6
    assert np.isnan(rows[0]["T_shield_max"])
    by_key = {(r["T_box"], r["n_shields"], r["eps_shield"]): r for r in rows if r["n_shields"]}
    assert by_key[(450.0, 1, 0.2)]["T_plate_max"] > by_key[(400.0, 1, 0.2)]["T_plate_max"]
    assert by_key[(400.0, 1, 0.45)]["T_shield_max"] > by_key[(400.0, 1, 0.2)]["T_shield_max"]


def test_many_shields_leave_the_plate_above_a_hot_box():
    # Air conduction across the gap sets a floor that extra shields can't remove.
    T_box = 150 + 273.15
    s = rad.Setup(T_amb=T_box, T_air=T_box, shield_thickness=0.1e-3)
    p0, p4, p10 = (rad.solve(s, n, n_strips=20).T_plate_max for n in (0, 4, 10))
    assert p10 < p4
    assert p4 - p10 < 0.1 * (p0 - p4)
    assert p10 > T_box
