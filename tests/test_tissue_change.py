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


def test_min_required_delay_4_5_min():
    from thznerve.physics.tissue_change import min_required_delay
    assert abs(min_required_delay(0.6e12) - 0.147) < 1e-3
