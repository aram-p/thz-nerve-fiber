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


def test_mixtures_and_literature_values():
    from thznerve.physics.materials import (WATER_KS1996, ConstantEps, Looyenga, axoplasm, myelin,
                                            refractive_index)
    f = 1e12
    w = complex(WATER_KS1996(f))
    assert abs(w - (4.14 - 2.14j)) < 0.02                      # docs/literature.md 2a table
    assert abs(complex(myelin()(f)) - (3.01 - 0.70j)) < 0.02    # 2b table, 40 % water
    d = ConstantEps(2.3 - 0.05j)
    assert abs(complex(Looyenga(WATER_KS1996, d, 1.0)(f)) - w) < 1e-12
    assert abs(complex(Looyenga(WATER_KS1996, d, 0.0)(f)) - d.eps) < 1e-12
    # myelin index stays below water's over the whole band (anti-guiding)
    fs = np.linspace(0.1e12, 3e12, 30)
    assert np.all(refractive_index(myelin(phi=0.5)(fs)).real < refractive_index(WATER_KS1996(fs)).real)
    assert np.all(refractive_index(axoplasm()(fs)).real < refractive_index(WATER_KS1996(fs)).real)


def test_quasistatic_transverse_vs_2d_coated_cylinder():
    """E perpendicular, node filled with myelin = infinite coated cylinder (exact 2D solution)."""
    from thznerve.physics.quasistatic import solve_node_qs
    e1, e2, e3, a, b = 4 - 2j, 2.5 - 0.3j, 4 - 2j, 3.5, 5.0
    r = solve_node_qs(a, b, 1.0, e1, e2, e2, e3, field="transverse")
    M = np.array([[a, -a, -1 / a, 0], [e1, -e2, e2 / a ** 2, 0], [0, b, 1 / b, -1 / b],
                  [0, e2, -e2 / b ** 2, e3 / b ** 2]], complex)
    A, B, C, _ = np.linalg.solve(M, np.array([0, 0, -b, -e3]))
    rr = np.linspace(a, b, 2001)
    shell = np.trapezoid((abs(B) ** 2 + abs(C) ** 2 / rr ** 4) * rr, rr) / np.trapezoid(rr, rr)
    assert abs(r.F_axon_at_node / abs(A) ** 2 - 1) < 0.01
    assert abs(r.F_node / shell - 1) < 0.03
    u = solve_node_qs(a, b, 1.0, e3, e3, e3, e3, field="transverse")
    assert abs(u.F_node - 1) < 1e-9


def test_liebe_matches_tabulated_37c():
    from thznerve.physics.materials import WATER_LIEBE37, water_liebe
    w = water_liebe(37.0)
    for f in (0.1e12, 0.5e12, 1e12):
        assert abs(complex(w(f)) / complex(WATER_LIEBE37(f)) - 1) < 5e-3


def test_node_factor_extrapolation_converges():
    from thznerve.physics.quasistatic import node_factor, solve_node_qs
    e = [6.97 - 10.4j, 4.53 - 3.16j, 7.55 - 13.8j, 7.55 - 13.8j]     # 0.1 THz, strongest contrast
    F = node_factor(3.5, 5.0, 0.5, *e)                       # from h, h/2
    F_fine = node_factor(3.5, 5.0, 0.5, *e, h_fine=0.5 / 80)  # from h/2, h/4
    single = solve_node_qs(3.5, 5.0, 0.5, *e).F_node
    assert abs(F / F_fine - 1) < 1e-3          # extrapolation is grid-independent
    assert abs(single / F_fine - 1) > 0.02     # the single default grid is not
