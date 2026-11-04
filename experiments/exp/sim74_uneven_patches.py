"""Sim 74 - Review round 34: does one slice under the beam see the ensemble mean that Sim 72's uneven-thickness
model assumes?

`uneven` averages the uniform-slab transfer over a Gaussian distribution of thickness: the mean over many
independent patches. A real slice under a beam ~5 mm across holds a finite number of patches. Here random
thickness fields v(x, y) (Gaussian, mean m, rms s, correlation exp(-r^2 / 2 L^2); Gaussian-filtered white noise) are laid
over a 12 x 12 mm slice on a 0.04 mm grid; the detected field for mode-matched detection is the transfer
weighted by the product of input and detector modes, exp(-2 r^2 / w^2), with w = 2.5 mm (5 mm across), fixed
or scaled with wavelength from 0.6 THz. Each slab transfer is local (valid for L above both the wavelength and
the 0.22 mm optical thickness). For L = 0.25, 0.5, 1, 2 mm and the sizes Sim 72 fitted, 200 realisations give the
rms departure of |T| (dB) and of arg T (rad) from the ensemble mean at 0.3, 0.6, 1.0, 1.5 THz.
Writes results/sim74/sim74.{txt,json}.
"""
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

from thznerve.physics.slice_change import uneven
from thznerve.runlog import RunLog

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim74"


def load(name, file):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / "exp" / file)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


S64 = load("s64", "sim64_size_fitted_tests.py")
VG = np.round(np.arange(-0.95, 1.801, 0.01), 3)
FQ = np.array([0.3, 0.6, 1.0, 1.5])
SIZES = [(0.12, 0.2), (0.18, 0.25), (-0.15, 0.4)]
LS = [0.25, 0.5, 1.0, 2.0]
DX, HALF, W0 = 0.04, 6.0, 2.5
NR = 200


def field(rng, L, m, s, x):
    g = gaussian_filter(rng.standard_normal((len(x), len(x))), sigma=L / DX / np.sqrt(2), mode="wrap")
    return m + s * (g - g.mean()) / g.std()


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fr = S64.fr
    idx = [int(np.argmin(abs(fr - f))) for f in FQ]
    T = np.array([S64.transfer("thicker path (fraction)", float(v))[idx] for v in VG])  # (len VG, 4)
    x = np.arange(-HALF, HALF + 1e-9, DX); X, Y = np.meshgrid(x, x); R2 = X ** 2 + Y ** 2
    lines, out = [], {}
    say = lambda s="": (print(s, flush=True), lines.append(s))
    with RunLog("sim74", "round 34: finite patches under a 5 mm beam vs the ensemble mean of uneven thickness (Sim 72)") as log:
        say("rms departure of one slice from the ensemble mean, over 200 realisations; |T| in dB / arg T in rad, at " +
            " / ".join(f"{f:g}" for f in FQ) + " THz")
        for m, s in SIZES:
            ens = uneven(T, VG, m, s)
            for scale in ("fixed", "with wavelength"):
                w = np.full(len(FQ), W0) if scale == "fixed" else W0 * 0.6 / FQ
                wt = [np.exp(-2 * R2 / wi ** 2) for wi in w]; wt = [a / a.sum() for a in wt]
                for L in LS:
                    rng = np.random.default_rng(int(1e3 * L) + int(1e2 * s))
                    dA, dP = [], []
                    for _ in range(NR):
                        v = np.clip(field(rng, L, m, s, x), VG[0], VG[-1])
                        k = np.clip(np.round((v - VG[0]) / 0.01).astype(int), 0, len(VG) - 1)
                        t = np.array([np.sum(wt[j] * T[k, j]) for j in range(len(FQ))])
                        dA.append(20 * np.log10(abs(t) / abs(ens))); dP.append(np.angle(t / ens))
                    a = np.sqrt(np.mean(np.square(dA), axis=0)); p = np.sqrt(np.mean(np.square(dP), axis=0))
                    out[f"m {m} s {s} | {scale} | L {L}"] = dict(dB=a.tolist(), rad=p.tolist())
                    say(f"  m {m:+.2f} rms {s:.2f}, waist {scale:15s} L {L:4.2f} mm: " + " / ".join(f"{u:4.1f}" for u in a) +
                        " dB;  " + " / ".join(f"{u:4.2f}" for u in p) + " rad")
        say("ensemble mean itself (0.3 / 0.6 / 1.0 / 1.5 THz), relative to the unchanged slab: " + "; ".join(
            f"m {m:+.2f} rms {s:.2f}: " + " / ".join(f"{20 * np.log10(abs(u)):.1f} dB" for u in uneven(T, VG, m, s)) for m, s in SIZES))
        log.outcome = "; ".join(f"{k}: {max(v['dB']):.1f} dB" for k, v in out.items() if "fixed" in k)[:250]
    (OUT / "sim74.json").write_text(json.dumps(out))
    (OUT / "sim74.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
