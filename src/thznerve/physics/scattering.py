"""Coherent wave in a slab holding spheres: Mie forward amplitude and the Foldy effective index.

Bohren-Huffman conventions inside (exp(-i w t), Im n > 0 for loss); `foldy_eps` takes and returns the
engineering permittivity used elsewhere in the package (exp(+j w t), eps = eps' - j eps'').

Foldy / van de Hulst: m_eff = 1 + i 2 pi N S(0) / k^3, N the number density, S(0) the forward amplitude
of one sphere in the host, k the host wavenumber. Valid for a dilute suspension (volume fraction up to
~0.1-0.2); for small spheres it reduces to Maxwell Garnett. The host's absorption enters through its
complex index; the size parameter uses the real part of the host wavenumber (an approximation for an
absorbing host).
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn

C0 = 299792458.0


def mie_ab(m, x, nmax=None):
    """Mie coefficients a_n, b_n (n = 1..nmax) for relative index m and size parameter x (BH 4.53)."""
    x = float(x)
    if nmax is None:
        nmax = int(np.ceil(x + 4 * x ** (1 / 3) + 2)) + 2
    n = np.arange(1, nmax + 1)
    mx = m * x
    jx, djx = spherical_jn(n, x), spherical_jn(n, x, derivative=True)
    yx, dyx = spherical_yn(n, x), spherical_yn(n, x, derivative=True)
    jm, djm = spherical_jn(n, mx), spherical_jn(n, mx, derivative=True)
    psi_x, dpsi_x = x * jx, jx + x * djx
    xi_x, dxi_x = x * (jx + 1j * yx), (jx + 1j * yx) + x * (djx + 1j * dyx)
    psi_m, dpsi_m = mx * jm, jm + mx * djm
    a = (m * psi_m * dpsi_x - psi_x * dpsi_m) / (m * psi_m * dxi_x - xi_x * dpsi_m)
    b = (psi_m * dpsi_x - m * psi_x * dpsi_m) / (psi_m * dxi_x - m * xi_x * dpsi_m)
    return n, a, b


def forward_amp(m, x):
    """S(0) = 1/2 sum (2n + 1)(a_n + b_n)."""
    n, a, b = mie_ab(m, x)
    return 0.5 * np.sum((2 * n + 1) * (a + b))


def q_ext_q_sca(m, x):
    n, a, b = mie_ab(m, x)
    q_ext = 2 / x ** 2 * np.sum((2 * n + 1) * (a + b).real)
    q_sca = 2 / x ** 2 * np.sum((2 * n + 1) * (abs(a) ** 2 + abs(b) ** 2))
    return q_ext, q_sca


def foldy_eps(f_hz, eps_host, eps_sph, r_um, fv):
    """Effective engineering permittivity of a host holding spheres of radius r_um at volume fraction fv."""
    if fv == 0:
        return complex(eps_host)
    nh = np.conj(np.sqrt(complex(eps_host)))            # BH: Im > 0
    ns = np.conj(np.sqrt(complex(eps_sph)))
    k0 = 2 * np.pi * f_hz / C0
    k = k0 * nh
    r = r_um * 1e-6
    s0 = forward_amp(ns / nh, (k.real) * r)
    N = fv / (4 / 3 * np.pi * r ** 3)
    n_eff = nh * (1 + 1j * 2 * np.pi * N * s0 / k.real ** 3)       # S(0) taken at x = Re(k) r
    return complex(np.conj(n_eff) ** 2)
