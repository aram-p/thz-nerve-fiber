import numpy as np
import pytest

from thznerve.physics import mie_cyl, tmm
from thznerve.physics.materials import C0, WATER_KS1996, alpha_per_cm, refractive_index


def test_water_ks1996_literature_values():
    # Kindt & Schmuttenmaer 1996: n ~ 2.1, alpha ~ 220 cm^-1 at 1 THz
    e = WATER_KS1996(1e12)
    assert e.imag < 0
    n = refractive_index(e)
    assert 2.0 < n.real < 2.25
    assert 190 < alpha_per_cm(e, 1e12) < 250


def test_tmm_lossless_slab_conserves_energy_and_matches_airy():
    f, n, d = 1e12, 1.5, 100e-6
    R, T, A, r, t = tmm.slab_rta(f, 1.0, [(n ** 2, d)], 1.0)
    assert abs(A) < 1e-12
    k0 = 2 * np.pi * f / C0
    r12 = (1 - n) / (1 + n)
    phase = np.exp(-2j * k0 * n * d)
    r_airy = (r12 - r12 * phase) / (1 - r12 ** 2 * phase)
    assert abs(R - abs(r_airy) ** 2) < 1e-12


def test_tmm_lossy_slab_beer_lambert_limit():
    f = 1e12
    e = WATER_KS1996(f)
    d = 500e-6
    R, T, A, *_ = tmm.slab_rta(f, 1.0, [(e, d)], 1.0)
    n = refractive_index(e)
    t_in = 4 * n.real / abs(1 + n) ** 2  # rough single-interface power transmission
    T_bl = (abs(2 / (1 + n)) ** 2 * n.real) * (abs(2 * n / (1 + n)) ** 2 / n.real) * np.exp(-alpha_per_cm(e, f) * d * 100)
    assert T == pytest.approx(T_bl, rel=0.02)
    assert 0 < A < 1 and t_in > 0


@pytest.mark.parametrize("pol", ["TM", "TE"])
def test_mie_lossless_no_absorption(pol):
    s = mie_cyl.solve(1e12, [5e-6, 7e-6], [2.0, 4.0, 1.7 ** 2], pol)
    assert np.allclose(s.absorbed_per_layer(), 0, atol=1e-12 * s.I0 * 1e-6)


def test_mie_small_core_quasistatic_tm():
    f = 0.5e12
    e = WATER_KS1996(f)
    a = 0.2e-6
    s = mie_cyl.solve(f, [a], [e, 1.0 + 0j], "TM")
    k0 = 2 * np.pi * f / C0
    c_qs = k0 * (-e.imag) * np.pi * a ** 2 / 1.0
    assert s.c_abs_per_layer()[0] == pytest.approx(c_qs, rel=1e-3)


@pytest.mark.parametrize("pol", ["TM", "TE"])
def test_mie_flux_matches_volume_integral(pol):
    f = 1.5e12
    ew = WATER_KS1996(f)
    s = mie_cyl.solve(f, [5e-6, 7e-6], [ew, 2.5 - 0.3j, ew], pol)
    w = 2 * np.pi * f
    from thznerve.physics.materials import EPS0
    # myelin shell absorbed power by quadrature
    rs = np.linspace(5e-6, 7e-6, 81)[1:-1]
    ph = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    integ = np.array([[s.field_E2(r, p) for p in ph] for r in rs]).mean(axis=1) * 2 * np.pi * rs
    p_vol = 0.5 * w * EPS0 * 0.3 * np.trapezoid(integ, rs)
    p_vol += 0.5 * w * EPS0 * 0.3 * 0  # endpoints negligible
    assert s.absorbed_per_layer()[1] == pytest.approx(p_vol, rel=0.03)
