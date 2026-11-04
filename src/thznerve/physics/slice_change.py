"""Changes of the measured slice under bias, as a slab model (Sims 46 and 48).

The slice is an effective-medium white-matter slab (E-par, packing 0.6, axon fraction 0.7) on PTFE, 220 um
optically thick (Sim 45), with Kindt-Schmuttenmaer water as base. Each candidate change returns a function
f -> (eps, thickness in um); `ratio` gives the change of the slab transfer t1/t0 relative to the unbiased
slab, every transfer taken against the same thickness of air (a thicker slice displaces air). The sign is
the numpy FFT one (exp(+jwt)): a delay makes arg negative, so multiplying a pulse's rfft by `ratio`
applies the change to the pulse.
"""
import numpy as np

from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import EPS0, PROTEIN, WATER_KS1996, axoplasm, myelin, water_liebe
from thznerve.physics.tmm import slab_rta

D_UM, P, N_PTFE, WF = 220.0, 0.6, 1.43, 0.77
W0 = WATER_KS1996


class HeatedWater:
    """Base water plus k times Liebe's change of eps from 20 C to 20 + dT C."""

    def __init__(self, dT, k=1.0):
        self.dT, self.k = dT, k

    def __call__(self, f):
        return W0(f) + self.k * (water_liebe(20.0 + self.dT)(f) - water_liebe(20.0)(f))


def eps_tissue(f, water=W0, sig_ax=0.0, sig_w=0.0):
    ax = axoplasm(water)(f) - 1j * sig_ax / (2 * np.pi * f * EPS0)
    e, _ = tissue_eps(ax, myelin(water)(f), water(f), 0.7, P)
    return e - 1j * WF * sig_w / (2 * np.pi * f * EPS0)


def t_slab(f, e, d_um):
    """Transfer relative to the same thickness of air on PTFE: a thicker slice displaces air, as a
    cuvette window moving outwards would (extra phase k0 (n - 1) dd)."""
    return slab_rta(f, 1.0, [(complex(e), d_um * 1e-6)], N_PTFE)[4] / slab_rta(f, 1.0, [(1.0 + 0j, d_um * 1e-6)], N_PTFE)[4]


def _mix_delta(f, x):
    """Change of tissue eps when its water content changes by the fraction x (dry matter fixed)."""
    dry = PROTEIN(f)
    mix = lambda wf: (wf * W0(f) ** (1 / 3) + (1 - wf) * dry ** (1 / 3)) ** 3
    wf1 = WF * (1 + x) / (1 + WF * x)
    return mix(wf1) - mix(WF)


def water_frac(x, mode="thickness"):
    """Water content changed by the fraction x (x < 0: loss). mode: 'thickness' (added volume into
    thickness), 'eps' (permittivity only, thickness fixed), 'geom' (thickness only, eps fixed),
    'isotropic' (volume spread in 3D: thickness x cube root)."""
    def fn(f):
        de = 0.0 if mode == "geom" else _mix_delta(f, x)
        d = {"thickness": D_UM * (1 + WF * x), "geom": D_UM * (1 + WF * x), "eps": D_UM,
             "isotropic": D_UM * (1 + WF * x) ** (1 / 3)}[mode]
        return eps_tissue(f) + de, d
    return fn


def mechanisms(k_heat):
    """Candidate changes: name -> (value -> (f -> (eps, d_um)), value range). k_heat scales Liebe's
    temperature change (Sim 46 fits it to the measured 4-21 C absorption change at 1 THz)."""
    return {
        "heating, Liebe (dT, K)": (lambda v: (lambda f: (eps_tissue(f, HeatedWater(v)), D_UM)), (0.0, 80.0)),
        f"heating, Liebe x{k_heat:.2f} (dT, K)": (lambda v: (lambda f: (eps_tissue(f, HeatedWater(v, k_heat)), D_UM)), (0.0, 60.0)),
        "water content (fractional change)": (water_frac, (-0.5, 0.5)),
        "ions: all tissue water (S/m)": (lambda v: (lambda f: (eps_tissue(f, sig_w=v), D_UM)), (0.0, 500.0)),
        "axoplasm conductivity (S/m)": (lambda v: (lambda f: (eps_tissue(f, sig_ax=v), D_UM)), (0.0, 2000.0)),
        "thicker path (fraction)": (lambda v: (lambda f: (eps_tissue(f), D_UM * (1 + v))), (0.0, 0.5)),
    }


def ratio(f_hz, fn):
    """t1/t0 at the frequencies f_hz (> 0) for a change f -> (eps, d_um)."""
    return np.array([t_slab(f, *fn(f)) / t_slab(f, eps_tissue(f), D_UM) for f in np.atleast_1d(f_hz)])


def k_millon():
    """Factor that scales Liebe's temperature change of water to the measured 4 -> 21 C absorption change at
    1 THz (46 /cm; Millon et al. 2022)."""
    from thznerve.physics.materials import alpha_per_cm
    d = alpha_per_cm(water_liebe(21.0)(1e12), 1e12) - alpha_per_cm(water_liebe(4.0)(1e12), 1e12)
    return 46.0 / d


def uneven(T, vg, m, s, vmin=-0.95):
    """Beam-averaged transfer of a slice whose relative thickness change v varies across the beam as a Gaussian
    of mean m and rms s, cut at v > vmin (thickness >= 5 % of D_UM). T: transfers per v (len(vg) x n), each
    taken at thickness D_UM (1 + v). The detector sees the coherent mean, so a spread of thickness loses
    amplitude as exp(-(phase spread)^2 / 2), a loss that rises with frequency. s = 0 gives the nearest v."""
    vg = np.asarray(vg, float); T = np.asarray(T)
    if s == 0:
        return T[int(np.argmin(abs(vg - m)))]
    w = np.exp(-0.5 * ((vg - m) / s) ** 2) * (vg > vmin)
    return (w / w.sum()) @ T
