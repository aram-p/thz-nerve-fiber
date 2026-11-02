"""Sim 44 - Review round 8 follow-up: is paper 1's 2 THz dip (30 s at 100 V) a feature of the sample,
or one of the narrow dips (water-vapour lines, noise) that the published spectra show everywhere
at this signal level?

Paper 1 page 3, Figs. 3b (empty cuvette, sample) and 4b (no voltage, 30 s, 4 min). Each curve is
traced per pixel column as its top (level) and bottom (narrow dips draw vertical strokes) with
thznerve.digitise.column_extent; dip depth is the running maximum of the level within +-0.05 THz
minus the bottom (thznerve.digitise.narrow_dips). Dips of every curve over 1.5-2.9 THz are listed
and matched in frequency across curves: a dip at the same frequency in the empty cuvette is a gas
line, not the sample.
Writes results/sim44/sim44.txt and results/sim44/dips.png.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.digitise import axis_map, column_extent, copies_disagree, narrow_dips, render_page

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "sim44"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


p1 = render_page(ROOT / "docs" / "references" / "paper 1.pdf", 3, 300)
SPEC = {
    "3b": dict(box=(930, 1400, 1750, 2030), gx=([121.0, 238.5, 355.5, 473.5, 590.5, 708.5], [0, 0.5, 1, 1.5, 2, 2.5]),
               gy=([56.5, 111.0, 164.5, 218.0, 271.5, 325.0, 378.0, 431.0, 485.0, 538.0, 592.5], list(range(0, -101, -10))),
               curves={"empty": "blue", "sample": "black"}),
    "4b": dict(box=(930, 2370, 1750, 3010), gx=([116.5, 235.0, 353.5, 472.0, 590.5, 709.5], [0, 0.5, 1, 1.5, 2, 2.5]),
               gy=([61.5, 122.5, 182.0, 243.0, 302.5, 362.0, 422.0, 482.5, 542.5, 603.5], list(range(-10, -101, -10))),
               curves={"none": "black", "30s": "magenta", "4min": "red"}),
}
CM = {"blue": lambda R, G, B: (B > 170) & (R < 120) & (G > 90) & (G < 200),
      "black": lambda R, G, B: (R < 90) & (G < 90) & (B < 90),
      "magenta": lambda R, G, B: (R > 170) & (B > 170) & (G < 120),
      "red": lambda R, G, B: (R > 180) & (G < 110) & (B < 110)}
F_LO, F_HI, HALF, MIN_DEPTH = 1.5, 2.9, 0.05, 6.0
dips, traces = {}, {}
for name, p in SPEC.items():
    l, t, r, b = p["box"]; rgb = p1[t:b, l:r]; R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    fxs, _ = axis_map(*p["gx"]); fys, _ = axis_map(*p["gy"])
    for key, col in p["curves"].items():
        m = CM[col](R, G, B)
        m[:, :int(p["gx"][0][0]) + 3] = False; m[:, 815:] = False
        m[:int(p["gy"][0][0]) - 3, :] = False
        m[int(p["gy"][0][-1]) - 2:, :] = False           # the frame line at -100 dB
        if name == "3b" and col == "black":
            m[:118, :] = False                             # label "1"
        if name == "4b" and col == "black":
            m[395:, 530:600] = False                       # arrow and label "2" (1.75-2.0 THz, below -65 dB)
        xs, top, bot = column_extent(m)
        f = fxs(xs); s = (f >= F_LO) & (f <= F_HI)
        traces[key] = (name, xs[s], top[s], bot[s], rgb)
        dips[key] = narrow_dips(f[s], fys(top[s]), fys(bot[s]), HALF, MIN_DEPTH)

def to_db(key):
    name, xs, top, bot, _ = traces[key]
    fxs, _ = axis_map(*SPEC[name]["gx"]); fys, _ = axis_map(*SPEC[name]["gy"])
    return fxs(xs), fys(top), fys(bot)


fa, ta, ba = to_db("sample"); fb_, tb, bb = to_db("none")
say("Copy check: Fig. 3b sample and Fig. 4b no voltage are one scan; bands differing by > 6 dB:")
say(f"  level: {copies_disagree(fa, ta, fb_, tb, 6.0) or 'none'}; bottom: {copies_disagree(fa, ba, fb_, bb, 6.0) or 'none'}")
say(f"Narrow dips >= {MIN_DEPTH:.0f} dB below the level within +-{HALF} THz, {F_LO}-{F_HI} THz (paper 1 page 3)")
for k, d in dips.items():
    say(f"  {k:6s}: {len(d)} dips; >= 10 dB: {sum(x[1] >= 10 for x in d)}; >= 14 dB: {sum(x[1] >= 14 for x in d)}; "
        + ", ".join(f"{fi:.3f} ({di:.0f})" for fi, di in d))
d30 = [x for x in dips["30s"] if 1.95 <= x[0] <= 2.05]
say()
if d30:
    f0, dep0 = min(d30, key=lambda x: abs(x[0] - 2.0))
    _, xs, top, bot, _ = traces["30s"]
    fxs, _ = axis_map(*SPEC["4b"]["gx"]); fys, _ = axis_map(*SPEC["4b"]["gy"])
    j = int(np.argmin(np.abs(fxs(xs) - f0)))
    say(f"30 s dip near 2 THz: {f0:.3f} THz, {dep0:.1f} dB below its own level; bottom at {fys(bot[j]):.1f} dB absolute")
    for k in ("empty", "sample", "none", "4min"):
        near = [x for x in dips[k] if abs(x[0] - f0) <= 0.015]
        deeper = [x for x in dips[k] if x[1] >= dep0]
        say(f"  {k:6s}: dip within 0.015 THz: {', '.join(f'{a:.3f} ({b:.0f} dB)' for a, b in near) or 'none'}; "
            f"dips at least as deep anywhere in {F_LO}-{F_HI} THz: {len(deeper)}")
# noise floor: the sample curve of Fig. 3b above 3.6 THz (the axis runs to 4.5 THz)
l, t = 930, 1400; rgb = p1[t:2030, l:2130]; R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
fxs, _ = axis_map(*SPEC["3b"]["gx"]); fys, _ = axis_map(*SPEC["3b"]["gy"])
for key, col in (("sample", "black"), ("empty", "blue")):
    m = CM[col](R, G, B); m[:, :int(SPEC["3b"]["gx"][0][0]) + 3] = False; m[:118, :] = False
    m[int(SPEC["3b"]["gy"][0][-1]) - 2:, :] = False
    xs, top, bot = column_extent(m); f = fxs(xs); s = (f >= 3.6) & (f <= 4.4)
    say(f"Fig. 3b {key} at 3.6-4.4 THz (noise): level median {np.median(fys(top[s])):.1f} dB, bottom median {np.median(fys(bot[s])):.1f} dB")
(OUT / "sim44.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
for i, (fig_name, keys) in enumerate((("3b", ("empty", "sample")), ("4b", ("none", "30s", "4min")))):
    ax[i].imshow(traces[keys[0]][4].astype(np.uint8))
    for k, c in zip(keys, ("yellow", "cyan", "lime")):
        _, xs, top, bot, _ = traces[k]
        ax[i].vlines(xs, top, bot, color=c, lw=0.4)
    ax[i].set_title(f"paper 1 Fig. {fig_name}: traced level-to-bottom per column, {F_LO}-{F_HI} THz")
fig.tight_layout(); fig.savefig(OUT / "dips.png", dpi=130)
