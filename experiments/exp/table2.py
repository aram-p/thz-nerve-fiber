"""Table 2 of the manuscript, built from the Sim 69-72 result files with one rule for every row.

Each candidate is judged in each reference by the better of its two fits, alone or with a drift (the lower
largest ratio; round 33: the drift fit contains the alone fit, so a cause is excluded only if both are).
"No change" and "shift only" stay separate rows: they are the two claims that the change is real and that it
changes shape. Per reference and fit: threshold per scan = hypot(quantile(noise draws, 1 - 0.05/3), tracing D); excluded if the largest
D / threshold over 2, 4, 5 min exceeds 1 (thznerve.waveform_fit.threshold_verdict, 200 bootstrap resamples of the
draws); borderline if the verdict flips in > 10 % of resamples or the largest ratio lies within 0.87-1.15.
Floors 1.0-2.5 dB/THz: eight references; floors 0-0.5: four. Writes results/table2.md and results/table2.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

from thznerve.waveform_fit import threshold_verdict

ROOT = Path(__file__).resolve().parents[2]
LATE = ("2 min (H2)", "4 min (H1)", "5 min (H3)")
LEVEL = 1 - 0.05 / 3
J = {n: json.loads((ROOT / "results" / f"sim{n}" / f"sim{n}.json").read_text()) for n in (69, 70, 71, 72)}
TR = J[70]["tracing"]
REFS8 = [f"floor {k:g}, " + ab for k in (1.0, 1.5, 2.07, 2.5) for ab in ("a <= 1", "a + b <= 1")]
REFS4 = [f"floor {k:g}, " + ab for k in (0.0, 0.5) for ab in ("a <= 1", "a + b <= 1")]
HX = "heating, Liebe x2.31 (dT, K)"

ROWS = [  # label, (sim, key) for floors 1.0-2.5, (sim, key) for floors 0-0.5, size formatter, fits
    ("No change", (69, "none"), (70, "none"), None, ("alone",)),
    ("Shift only", (69, "none"), (70, "none"), None, ("drift",)),
    ("Thicker tissue path, uniform (0–51 %)", (71, "thicker"), (70, "thicker path (fraction)"), "thick", None),
    ("Uneven thickness across the beam (mean −30 to +51 %, rms 0–40 %)", (72, "uneven"), (72, "uneven"), "uneven", None),
    ("Gas in a slice that swells (r 20–100 µm, v 0–0.2)", (71, "gas"), (71, "gas"), "gas", None),
    ("Node states, all 80 (λ 0.5–6 of every extreme)", (69, "node"), None, "node", None),
    ("Ions in all tissue water (0–60 S/m)", (69, "ions: all tissue water (S/m)"), (70, "ions: all tissue water (S/m)"), "S/m", None),
    ("Axoplasm conduction (0–150 S/m)", (69, "axoplasm conductivity (S/m)"), (70, "axoplasm conductivity (S/m)"), "S/m", None),
    ("Warming A, Liebe ×2.31", (69, HX), (70, "warming A"), "K", None),
    ("Warming B, slow relaxation (×1.04)", (70, "warming B"), (70, "warming B"), "K", None),
    ("Warming C, fast relaxation (×3.26)", (70, "warming C"), (70, "warming C"), "K", None),
    ("Warming, slow and fast shared freely (w B + (1 − w) C)", (71, "warming"), (71, "warming"), "mix", None),
]


def verdicts(sim, key, refs, rng, modes=None):
    """Per reference, the verdict of the better fit among modes (default alone and drift)."""
    per = {m: _verdicts(sim, f"{key} | {m}", refs, rng) for m in (modes or ("alone", "drift"))}
    return {vn: dict(min((p[vn] for p in per.values()), key=lambda v: v["largest"])) for vn in refs}


def _verdicts(sim, key, refs, rng):
    out = {}
    for vn in refs:
        rows = [J[sim]["cases"][f"{vn} | {sc}"][key] for sc in LATE]
        v = threshold_verdict([r["D"] for r in rows], [r["draws"] for r in rows], TR[vn][key.split(" | ")[1]], LEVEL, 200, rng)
        v["borderline"] = v["flip"] > 0.1 or 0.87 < v["largest"] < 1.15
        v["mode"] = key.split(" | ")[1]; v["D"] = [r["D"] for r in rows]; v["size"] = [r.get("size") for r in rows]; v["state"] = [r.get("state") for r in rows]
        out[vn] = v
    return out


def rng_(a):
    return f"{min(a):.2f}–{max(a):.2f}"


def sizes(fmt, per):
    s = [x for v in per.values() for x in v["size"]]
    if fmt is None:
        return "–"
    if fmt == "thick":
        f = [x[1] if isinstance(x, list) else x for x in s]; return f"{100 * min(f):.0f}–{100 * max(f):.0f} %"
    if fmt == "gas":
        r = sorted({x[0] for x in s}); v = [x[1] for x in s]
        return f"r {r[0]}–{r[-1]} µm, v {min(v):.2f}–{max(v):.2f}"
    if fmt == "uneven":
        m = [x[0] for x in s]; r = [x[1] for x in s]
        return f"mean {100 * min(m):.0f} to {100 * max(m):.0f} %, rms {100 * min(r):.0f}–{100 * max(r):.0f} %"
    if fmt == "node":
        return f"λ {min(s):g}–{max(s):g}"
    if fmt == "mix":
        w = [x[0] for x in s]; t = [x[1] for x in s]; return f"w {min(w):.1f}–{max(w):.1f}, {min(t):g}–{max(t):g} K"
    return f"{min(s):g}–{max(s):g} {fmt}"


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    rng = np.random.default_rng(2)
    out, lines = {}, ["| Candidate (sizes searched) | Best sizes | Fits with a drift (of 8) | D (%), 2 / 4 / 5 min | Largest D ÷ threshold | Excluded in (of 8): all (borderline) | Floors 0–0.5: D ÷ threshold; excluded (of 4), borderline |",
                      "|---|---|---|---|---|---|---|"]
    for label, a, b, fmt, modes in ROWS:
        per = verdicts(*a, REFS8, rng, modes)
        D = [[v["D"][i] for v in per.values()] for i in range(3)]
        L = [v["largest"] for v in per.values()]
        ne = sum(v["excluded"] for v in per.values()); nb = sum(v["borderline"] for v in per.values())
        nd = sum(v["mode"] == "drift" for v in per.values())
        if b:
            p4 = verdicts(*b, REFS4, rng, modes); L4 = [v["largest"] for v in p4.values()]
            c4 = f"{rng_(L4)}; {sum(v['excluded'] for v in p4.values())}, {sum(v['borderline'] for v in p4.values())}"
        else:
            p4, c4 = None, "not refitted"
        lines.append(f"| {label} | {sizes(fmt, per)} | {nd} | " + " / ".join(f"{min(d):.1f}–{max(d):.1f}" for d in D) +
                     f" | {rng_(L)} | {ne} ({nb}) | {c4} |")
        out[label] = dict(source=a, source_low=b, per=per, per_low=p4)
    (ROOT / "results" / "table2.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (ROOT / "results" / "table2.json").write_text(json.dumps(out, default=float))
    print("\n".join(lines))
