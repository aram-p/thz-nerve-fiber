"""Sim 64 - Review round 28: the waveform tests with the candidate's size fitted inside the statistic and a null
whose noise keeps the data's correlation.

Round 28 found that the block bootstrap of Sims 61 and 63 gave synthetic scans about twice the data's N_eff, that
holding the size fixed in the synthetic scans favours passing, and that 60-100 draws cannot resolve verdicts near
0.05. Here (thznerve.waveform_fit.size_fitted_test):
  statistic  null: the candidate at the best size on a grid (grid[0] = no candidate), no loss term, s = 0 (alone)
             or s free (with a drift); alternative: the best size on the same grid with the causal loss increment,
             s free; F(2) alone, F(1) with a drift; N_eff from the best alternative's residual. The same fits
             (warm starts from multi-start fits on the data) on the data and on every synthetic scan.
  null       the null model at its best size plus a phase-randomised surrogate of the alternative's residual
             (same amplitude spectrum, so the same autocorrelation). The F statistic uses the data's N_eff on
             every synthetic scan too, so it is one monotone function of the misfit ratio; the draws' own
             N_eff (refitting raises it 1.1-1.4x) is reported as a diagnostic.
  draws      batches of 50 until the 95 % Clopper-Pearson interval of the fraction p_b <= p excludes 0.05
             (pass above, fail below), at most 400; else 'undecided'. Seed per reference, scan and mode.
Each candidate is tested over its plausible range and over a wide range:
  thicker path 0-30 %; ions in all tissue water 0-5 S/m added (wide 0-60); axoplasm 0-5 S/m added (wide 0-150);
  warming Liebe x2.31 0-24 K; gas in a swelling slice r 30 and 50 um, v 0-0.2 (wide 0-0.4); node states (the
  state Sim 61 used for that reference, scan and mode) lambda 0-1 of the every-extreme change (wide 0-6).
References: floors 1.0 and 2.5 dB/THz x (a <= 1, a + b <= 1), increment held above 2 THz; scans 2, 4, 5 min.
Also the no-change and shift-alone tests (grid = no candidate) on all eight references of Sim 59 and all four
scans. Four worker processes; resumes from a partial results/sim64/sim64.json. Writes results/sim64/sim64.txt.
"""
import json
import sys
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from thznerve.physics.scattering import foldy_eps
from thznerve.physics.slice_change import D_UM, eps_tissue, k_millon, mechanisms, ratio
from thznerve.runlog import RunLog
from thznerve.waveform_fit import TwoPulse, load_scans, multistart, node_transfer, size_fitted_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim64"
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 400
tt, K, SCANS = load_scans(ROOT)
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
MECH = mechanisms(k_millon())
HX = [n for n in MECH if n.startswith("heating, Liebe x")][0]
NODES = {x["state"]: x for x in json.loads((ROOT / "results" / "sim30" / "sim36_signedband.json").read_text())["states"]}
R61 = json.loads((ROOT / "results" / "sim61" / "sim61.json").read_text())["cases"]
tp = TwoPulse(tt, K, f_loss=2.0)
fr = tp.fr; band = (fr > 0) & (fr <= 4.0)
REFS8 = [f"floor {k:g}, " + ab for k in (1.0, 1.5, 2.07, 2.5) for ab in ("a <= 1", "a + b <= 1")]
REFS4 = [f"floor {k:g}, " + ab for k in (1.0, 2.5) for ab in ("a <= 1", "a + b <= 1")]
GRIDS = {  # name: (plausible grid, wide grid)
    "thicker path (fraction)": (np.linspace(0, 0.3, 11), None),
    "ions: all tissue water (S/m)": (np.array([0, 1, 2, 3, 4, 5.0]), np.array([0, 2.5, 5, 10, 15, 20, 30, 40, 60.0])),
    "axoplasm conductivity (S/m)": (np.array([0, 1, 2, 3, 4, 5.0]), np.array([0, 5, 10, 20, 30, 50, 75, 100, 150.0])),
    HX: (np.array([0, 1, 2, 3, 4, 6, 8, 12, 16, 24.0]), None),
    "gas, slice swells, r 30 um": (np.linspace(0, 0.2, 6), np.linspace(0, 0.4, 11)),
    "gas, slice swells, r 50 um": (np.linspace(0, 0.2, 6), np.linspace(0, 0.4, 11)),
    "node": (np.linspace(0, 1, 5), np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6.0])),
}


def transfer(name, v):
    if name.startswith("node "):
        return node_transfer(NODES[name[5:]], fr, lam=v)
    r = np.ones(len(fr), complex)
    if v != 0:
        if name.startswith("gas"):
            rad = int(name.split("r ")[1].split()[0])
            r[band] = ratio(fr[band] * 1e12, (lambda f: (foldy_eps(f, eps_tissue(f), 1.0, rad, v), D_UM / (1 - v))))
        else:
            r[band] = ratio(fr[band] * 1e12, MECH[name][0](v))
        r[fr > 4.0] = r[band][-1]
    return r


def node_state(vn, sc, mode):
    """The node state Sim 61 used (best of three screened at lambda = 1, chosen on the data)."""
    return next(k.split(" | ")[2] for k in R61 if k.startswith(f"{vn} | {sc} | node ") and k.endswith(f"| {mode}"))


def cases():
    out = [(vn, sc, "none", t, "-") for t in ("alone", "drift") for sc in SCANS for vn in REFS8]
    for vn in REFS4:
        for sc in LATE:
            for mode in ("alone", "drift"):
                for c, (g, w) in GRIDS.items():
                    out.append((vn, sc, c, mode, "plausible"))
                    if w is not None:
                        out.append((vn, sc, c, mode, "wide"))
    return out


def run(case):
    vn, sc, c, mode, rng_name = case
    k0 = float(vn.split()[1].rstrip(",")); ab = 1.0 if "a + b" in vn else None
    one = lambda y, x, **kw: tp.fit(y, x, a_max=1.0, ab_max=ab, **kw)
    y, x0 = SCANS[sc]
    xf = tp.apply_transfer(x0, tp.loss_transfer(k0, "zero"))
    if c == "none":
        name, grid = "none", np.array([0.0])
    else:
        name = node_state(vn, sc, mode) if c == "node" else c
        grid = GRIDS[c][0 if rng_name == "plausible" else 1]
    xs = [xf] + [tp.apply_transfer(xf, transfer(name, v)) for v in grid[1:]]
    seed = zlib.crc32(f"{vn}|{sc}|{mode}".encode())
    d = size_fitted_test(multistart(one), one, y, xs, 0.0 if mode == "alone" else None, 2 if mode == "alone" else 1,
                         lambda q, x, s: tp.model(q, x, s=s), rng=np.random.default_rng(seed), nb_max=NB)
    return case, dict(name=name, grid=[float(v) for v in grid], size_null=float(grid[d["size_index_null"]]),
                      size_alt=float(grid[d["size_index_alt"]]), **d)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    part = OUT / "sim64.json"
    res = json.loads(part.read_text()).get("cases", {}) if part.exists() else {}
    with RunLog("sim64", f"round 28: size fitted inside the statistic, phase-randomised null, up to {NB} draws, plausible and wide ranges") as log:
        cs = [c for c in cases() if " | ".join(c) not in res]
        print(f"{len(res)} cases done, {len(cs)} to run", flush=True)
        with Pool(4) as pool:
            for case, d in pool.imap_unordered(run, cs):
                res[" | ".join(case)] = d
                print(f"{' | '.join(case)}: {d['verdict']} (p_cal {d['p_cal']:.3f}, 95 % {d['ci'][0]:.3f}-{d['ci'][1]:.3f}, "
                      f"{d['nb']} draws); nominal p {d['p_nominal']:.4f}; best size {d['size_null']:g}; N_eff data "
                      f"{d['ne_obs']:.1f}, draws {d['ne_draws_median']:.1f}", flush=True)
                part.write_text(json.dumps({"cases": res}, indent=1, default=float))
        lines = [f"{k}: {d['verdict']} (p_cal {d['p_cal']:.3f}, 95 % {d['ci'][0]:.3f}-{d['ci'][1]:.3f}, {d['nb']} draws); "
                 f"nominal p {d['p_nominal']:.4f}; best size {d['size_null']:g} ({d['name'][:30]}); N_eff data "
                 f"{d['ne_obs']:.1f}, draws {d['ne_draws_median']:.1f}" for k, d in sorted(res.items())] + [""]
        summ = {}
        for k, d in res.items():
            vn, sc, c, mode, rn = k.split(" | ")
            row = summ.setdefault(f"{c} | {mode} | {rn}", {s: {"pass": 0, "fail": 0, "undecided": 0} for s in SCANS})
            row[sc][d["verdict"]] += 1
        for k, row in sorted(summ.items()):
            lines.append(f"{k:58s} " + "; ".join(f"{sc[:5]} {v['pass']}/{v['fail']}/{v['undecided']}" for sc, v in row.items()
                                                  if sum(v.values())) + "  (pass/fail/undecided)")
        ratio_ne = [d["ne_draws_median"] / d["ne_obs"] for d in res.values()]
        lines.append(f"N_eff of draws / data: median {np.median(ratio_ne):.2f}, range {min(ratio_ne):.2f}-{max(ratio_ne):.2f}")
        print("\n".join(lines[-len(summ) - 1:]), flush=True)
        log.outcome = f"{len(res)} cases; N_eff draws/data {np.median(ratio_ne):.2f}"
    part.write_text(json.dumps({"cases": res, "summary": summ}, indent=1, default=float))
    (OUT / "sim64.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
