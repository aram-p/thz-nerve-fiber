"""Sim 65 - Is the phase-randomised null of Sim 64 valid? Test it on data with no change.

The three reports each show a no-voltage pulse (H1, H2, H3), traced from three figures. One fitted with another is
a comparison with no tissue change, only tracing error and figure-to-figure differences, so a valid test should
reject 'no change' in about 5 % of such pairs, and a null that rejects most of them is too strict. For each of
the six ordered pairs, the no-change test (s = 0) and the shift-alone test (s free) of Sim 64
(thznerve.waveform_fit.size_fitted_test, grid = no candidate, phase-randomised surrogates, up to 400 draws),
with a <= 1 and with a free. Writes results/sim65/sim65.txt and results/sim65/sim65.json.
"""
import json
import zlib
from itertools import permutations
from pathlib import Path

import numpy as np

from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, multistart, size_fitted_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim65"; OUT.mkdir(parents=True, exist_ok=True)
tt, K, SCANS = load_scans(ROOT)
tp = TwoPulse(tt, K, f_loss=2.0)
REF = {"H1": SCANS["30 s (H1)"][1], "H2": SCANS["2 min (H2)"][1], "H3": SCANS["5 min (H3)"][1]}
lines, out = [], {}
with RunLog("sim65", "phase-randomised null of Sim 64 tested on the six unbiased tracing pairs (no change)") as log:
    for a, b in permutations(REF, 2):
        for amax in (1.0, np.inf):
            one = lambda y, x, **kw: tp.fit(y, x, a_max=amax, **kw)
            for mode in ("alone", "drift"):
                key = f"{a} from {b} | a <= {amax:g} | {'no change' if mode == 'alone' else 'shift alone'}"
                d = size_fitted_test(multistart(one), one, REF[a], [REF[b]], 0.0 if mode == "alone" else None,
                                     2 if mode == "alone" else 1, lambda q, x, s: tp.model(q, x, s=s),
                                     rng=np.random.default_rng(zlib.crc32(key.encode())), nb_max=400)
                out[key] = d
                s = (f"{key:40s}: {d['verdict']:9s} p_cal {d['p_cal']:.3f} (95 % {d['ci'][0]:.3f}-{d['ci'][1]:.3f}, "
                     f"{d['nb']} draws); nominal p {d['p_nominal']:.4f}; N_eff data {d['ne_obs']:.1f}, draws {d['ne_draws_median']:.1f}")
                print(s, flush=True); lines.append(s)
    n_rej = sum(d["verdict"] == "fail" for d in out.values())
    lines.append(f"\nrejected (no change or shift alone) in {n_rej} of {len(out)} no-change comparisons")
    log.outcome = lines[-1].strip()
(OUT / "sim65.json").write_text(json.dumps(out, indent=1, default=float))
(OUT / "sim65.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
print(lines[-1])
