"""Sim 54b - How much gas would a swelling slice need at 4 and 5 min?

Sim 54 put the best gas fraction at the edge of its grid (0.2) at 4-5 min. Here gas in a swelling slice
(r 30 and 50 um) is tested on a grid to 0.4, beyond the dilute (Foldy) range: a fraction that passes only
above ~0.2 says what the data would ask for, not that the model holds there. Same test as Sim 54
(`waveform_fit.supply_test`), references with the increment held above 2 THz. Writes results/sim54/sim54b.txt.
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
from thznerve.physics.scattering import foldy_eps  # noqa: E402
from thznerve.physics.slice_change import D_UM, eps_tissue, ratio  # noqa: E402
from thznerve.waveform_fit import TwoPulse, f_test, load_scans, n_eff, supply_test  # noqa: E402

tt, K, SCANS = load_scans(ROOT)
GRID = np.linspace(0.0, 0.4, 21)
tp = TwoPulse(tt, K, f_loss=2.0)
KAPPA0 = float(TwoPulse(tt, K).fit(*SCANS["30 s (H1)"], loss="zero", a_max=1.0)[0][4])
fr = tp.fr; band = (fr > 0) & (fr <= 4.0)
FLOOR = tp.loss_transfer(KAPPA0, "zero")
lines = []


def tr_gas(rad):
    def tr(v):
        r = np.ones(len(fr), complex)
        if v:
            r[band] = ratio(fr[band] * 1e12, lambda f: (foldy_eps(f, eps_tissue(f), 1.0, rad, v), D_UM / (1 - v)))
            r[fr > 4.0] = r[band][-1]
        return r
    return tr


for ab in (None, 1.0):
    fit = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
    for scan in ("2 min (H2)", "4 min (H1)", "5 min (H3)"):
        y, x0 = SCANS[scan]
        xf = tp.apply_transfer(x0, FLOOR)
        pr, s_ref, res = fit(y, xf, loss="minphase")
        ne, _ = n_eff(res)
        q0, _, _ = fit(y, xf, s=0.0)
        for rad in (30, 50):
            d = supply_test(tp, fit, y, xf, tr_gas(rad), GRID, q0, pr, s_ref, ne)
            s = (f"{'a + b <= 1' if ab else 'a <= 1':10s} {scan}: gas r {rad} um, {d['verdict']:4s} p = {d['p']:.3f} "
                 f"({d['p_own_neff']:.3f}); best v {d['best']:.2f}; v passing "
                 + (f"{d['pass_range'][0]:.2f}-{d['pass_range'][1]:.2f}" if d['pass_range'] else "none"))
            print(s, flush=True); lines.append(s)
(ROOT / "results" / "sim54" / "sim54b.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
