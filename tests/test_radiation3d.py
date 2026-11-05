from pathlib import Path

import numpy as np
import pytest

import radiation as rad
import radiation3d as r3

CSV = Path(__file__).resolve().parents[1] / "astral_temps.csv"


def test_parallel_patches_conserve_whole_rectangle():
    length, height, c = 0.14, 0.22, 0.004
    xe, ze = np.linspace(0, length, 6), rad.strip_edges(height, 8)
    F = r3.parallel_patches(xe, ze, c)
    A = np.outer(np.diff(xe), np.diff(ze)).ravel()
    AF = A[:, None] * F
    assert AF.sum() == pytest.approx(length * height * rad.parallel_rectangles(length, height, c), rel=1e-9)
    assert np.allclose(AF, AF.T)


def test_parallel_patches_reduce_to_strips():
    length, c = 0.14, 0.005
    ze = rad.strip_edges(0.22, 6)
    F = r3.parallel_patches(np.array([0, length]), ze, c)
    assert np.allclose(F, rad.parallel_strips(length, ze, ze, c))


def test_tube_segments_sum_to_whole_tube():
    r, c, length, height = 0.0095, 0.0195, 0.14, 0.22
    seg = np.linspace(0, length, 5)
    xe, ze = np.linspace(0, length, 4), np.array([-height / 2, 0, height / 2])
    F = r3.tube_segments_to_patches(r, seg, c, xe, ze, n_phi=360, per_segment=20)
    whole = rad.tube_to_rectangle(r, length, c, -height / 2, height / 2)
    # Equal-length segments: the tube's view factor is the mean over segments.
    assert F.sum(axis=1).mean() == pytest.approx(whole, rel=2e-3)


def test_uniform_tube_matches_strip_model_on_average():
    s = rad.Setup()
    flat = (np.array([0.0, s.length]), np.full(2, s.T_tube))
    res3 = r3.solve3d(s, flat, 0, tube_x=(0.0, s.length), plate_x0=0.0, n_x=6, n_z=16, n_segments=12)
    res1 = rad.solve(s, 0, n_strips=16)
    assert res3.T_plate_mean == pytest.approx(res1.T_plate_mean, abs=1.0)
    # x-averaging in the strip model can only lower the peak.
    assert res3.T_plate_max >= res1.T_plate_max - 0.5


def test_profiles_are_mirrored_about_the_centre():
    for name, (x, T) in r3.load_tube_profiles(CSV).items():
        assert np.allclose(x + x[::-1], 0.4), name
        assert np.allclose(T, T[::-1]), name
        assert T.argmax() == len(T) // 2, name
