import numpy as np

from thznerve.physics.tissue_change import delay_phase, required_change


def test_matches_sim43_phase_at_2thz():
    # Sim 43 (0.2-1.4 THz band): arg T at 2.0 THz = -7.90 rad
    assert abs(-delay_phase(2.0e12) - (-7.90)) < 0.01


def test_required_change_at_0p6():
    loss, delay = required_change(np.array([0.6e12]), 1.09)
    assert abs(loss[0] - 1.089) < 1e-3
    assert abs(delay[0] - 0.227) < 1e-3          # 0.09 * (2 pi 0.6 THz 0.611 ps + 0.22)
    _, d0 = required_change(np.array([0.6e12]), 1.09, offset=False)
    assert d0[0] < delay[0]


def test_min_required_delay_matches_sim47():
    import json
    from pathlib import Path
    from thznerve.physics.tissue_change import SHIFT_PS, min_required_delay
    assert abs(min_required_delay(0.6e12) - 0.0735) < 1e-3          # 5-min lower bound, 0.0195 ps
    assert abs(min_required_delay(0.6e12, ("2 min", "4 min")) - 0.0784) < 1e-3
    res = json.loads((Path(__file__).resolve().parents[1] / "results" / "sim47" / "sim47.json").read_text())
    for scan, key in (("2 min", "2 min (H2)"), ("4 min", "4 min (H1)"), ("5 min", "5 min (H3)")):
        assert abs(SHIFT_PS[scan][0] - res[key]["s"]) < 1e-3
        assert abs(SHIFT_PS[scan][1] - res[key]["prof_lo95"]) < 1e-3


def test_uneven_thickness_matches_gaussian_phase_average():
    from thznerve.physics.slice_change import uneven
    vg = np.arange(-0.95, 2.0, 0.002); a = np.array([1.0, 4.0, 8.0])
    T = np.exp(1j * np.outer(vg, a))  # a pure phase per unit thickness change
    for m, s in ((0.1, 0.1), (0.0, 0.2)):
        got = uneven(T, vg, m, s)
        assert np.allclose(got, np.exp(1j * a * m - (a * s) ** 2 / 2), atol=1e-6)
    assert np.allclose(uneven(T, vg, 0.3, 0), np.exp(1j * a * vg[np.argmin(abs(vg - 0.3))]))
