from dataclasses import replace

import numpy as np
import pytest

import fea3d
import radiation as rad
import radiation3d as r3

C = 273.15
S = rad.Setup()
UNIFORM = (np.array([0.0, 0.14]), np.array([S.T_tube, S.T_tube]))
COARSE = dict(tube_x=(0.0, 0.14), plate_x0=0.0, n_x=10, n_z=16, n_segments=14)


@pytest.mark.parametrize("n_shields", [0, 1])
def test_matches_patch_model_for_uniform_tube(n_shields):
    f = fea3d.solve_fea(S, UNIFORM, n_shields, **COARSE)
    p = r3.solve3d(S, UNIFORM, n_shields, n_z=16, **{k: v for k, v in COARSE.items() if k != "n_z"})
    assert f.T_plate_max == pytest.approx(p.T_plate_max, abs=3.0)
    assert f.T_plate_mean == pytest.approx(p.T_plate_mean, abs=2.0)


def test_isothermal_when_tube_is_at_room_temperature():
    cold = (np.array([0.0, 0.14]), np.array([S.T_amb, S.T_amb]))
    f = fea3d.solve_fea(S, cold, 1, **COARSE)
    assert all(np.allclose(T, S.T_amb, atol=1e-3) for T in f.T)


def test_sheets_cool_outward():
    f = fea3d.solve_fea(S, UNIFORM, 2, **COARSE)
    peaks = [T.max() for T in f.T]
    assert S.T_amb < peaks[-1] < peaks[0] < S.T_tube


def test_outside_at_box_temperature_matches_default():
    a = fea3d.solve_fea(S, UNIFORM, 1, **COARSE)
    b = fea3d.solve_fea(S, UNIFORM, 1, T_outside=S.T_amb, **COARSE)
    assert b.T_plate_max == pytest.approx(a.T_plate_max, abs=1e-6)


def test_cooler_outside_cools_the_outer_sheet():
    hot_box = replace(S, T_amb=150 + C)
    a = fea3d.solve_fea(hot_box, UNIFORM, 1, **COARSE)
    b = fea3d.solve_fea(hot_box, UNIFORM, 1, T_outside=30 + C, **COARSE)
    assert b.T_plate_max < a.T_plate_max


def test_layout_spacing():
    sheets = fea3d.layout(S, 1, plate=True)
    r = S.tube_od / 2
    assert sheets[-1].y0 == pytest.approx(r + S.gap)
    assert sheets[0].y0 + S.shield_thickness / 2 == pytest.approx(r + S.gap / 2)
    no_plate = fea3d.layout(S, 1, plate=False)
    assert no_plate[0].y0 + S.shield_thickness / 2 == pytest.approx(r + S.gap)
