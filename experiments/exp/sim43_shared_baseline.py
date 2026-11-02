"""Sim 43 - Review round 7: (a) are the three reports one measurement session? (b) does the 2-min
waveform of paper 2 fall where the undelayed component of Sims 39/42 says? (c) does the component
explain paper 1's published spectra at 1.5-3 THz, including the 2 THz feature?

(a) The empty-cuvette and no-voltage waveforms of paper 1 (Figs. 3a/4a, stated 70 um) and paper 3
    (Fig. 1a, stated 100 um) are compared after tracing (results/sim42, results/sim39); paper 2
    Fig. 2 (stated 100 um) no-voltage trace is traced here and compared too.
(b) Paper 2 Fig. 2: no voltage (blue) and 2 min at 100 V (red), same axes as paper 1 Fig. 4a.
    Fitted with paper 1's empty-cuvette trace: biased = a * no voltage + b * empty(t - tau).
(c) Paper 1 Figs. 3b (empty, sample) and 4b (no voltage, 30 s, 4 min): published power spectra in
    dB, traced as the upper envelope per pixel column (water-vapour lines are narrow downward
    spikes) and binned to 0.05 THz (median). For a biased spectrum S = |a T + b e^{i phi}| |K|, any
    phase gives |S/K| >= |b - a |T||; the observed ratio is compared with that floor.
Writes results/sim43/sim43.txt, results/sim43/spectra.npz, results/sim43/trace_check.png.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from thznerve.digitise import axis_map, render_page, trace_track

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "docs" / "references"
OUT = ROOT / "results" / "sim43"; OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


rms = lambda r: float(np.sqrt(np.mean(np.square(r))))
H1 = np.load(ROOT / "results" / "sim42" / "digitised.npz")
H3 = np.load(ROOT / "results" / "sim39" / "digitised.npz")
tt = H1["t_ps"]

# ---------------------------------------------------------------- (a) shared baseline
say("(a) Are the reports one session? Traced waveforms on paper 1's time grid")


def compare(x, y, name):
    best = min((rms(np.interp(tt - s, tt, x) - y) / rms(y), s) for s in np.arange(-0.1, 0.1001, 0.005))
    say(f"  {name}: rms difference {rms(x - y) / rms(y):.1%} of signal (best shift {best[1]:+.3f} ps: {best[0]:.1%}); "
        f"minima {x.min():.3f} at {tt[np.argmin(x)]:.3f} ps vs {y.min():.3f} at {tt[np.argmin(y)]:.3f} ps")
    return best[0]


h3 = {k: np.interp(tt, H3["t_ps"], H3["wave_" + k]) for k in ("black", "blue", "red")}
compare(H1["empty"], h3["black"], "empty cuvette, paper 1 vs paper 3")
compare(H1["none"], h3["blue"], "no voltage, paper 1 vs paper 3")
compare(H1["4min"], h3["red"], "biased, paper 1 4 min vs paper 3 5 min (different scans expected)")

# ---------------------------------------------------------------- (b) paper 2, 2 min
p2 = render_page(REF / "paper 2.pdf", 3, 300)
l, t, r, b = 480, 1060, 1240, 1570
rgb = p2[t:b, l:r]; R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
fx, ex = axis_map([120.0, 222.0, 324.5, 426.5, 528.5, 630.5, 732.5], [-363, -362, -361, -360, -359, -358, -357])
fy, ey = axis_map([27.5, 72.5, 117.5, 163.0, 207.5, 252.5, 297.5, 342.0, 387.5, 432.0, 478.0],
                  [2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2, -2.5, -3])
masks = {"blue": (B > 150) & (R < 120) & (G < 170) & (G > 60), "red": (R > 200) & (G < 100) & (B < 100)}
tr2 = {}
for c, m in masks.items():
    m = m.copy(); m[:, :125] = False; m[:, 728:] = False
    xs, ys = trace_track(m, y_start=207.5, max_jump=25)
    tr2[c] = (xs, ys)
w2 = {c: np.interp(tt, fx(xs), fy(ys)) for c, (xs, ys) in tr2.items()}
say()
say(f"(b) Paper 2 Fig. 2 (page 3), axis misfit {ex:.3f} ps / {ey:.3f}")
compare(w2["blue"], H1["none"], "no voltage, paper 2 vs paper 1")


def shift(y, tau):
    return np.interp(tt - tau, tt, y, left=0.0, right=0.0)


def fit_td(y, x0, e0):
    res = lambda p: p[0] * x0 + p[1] * shift(e0, p[2]) - y
    best = min((least_squares(res, [0.6, b0, t0]) for b0 in (0.05, 0.2) for t0 in (-0.1, 0.0, 0.1)), key=lambda s: s.cost)
    one = min((least_squares(lambda p: p[0] * shift(x0, p[1]) - y, [0.6, t0]) for t0 in np.arange(-0.8, 0.81, 0.2)),
              key=lambda s: s.cost)
    return best.x, rms(best.fun) / rms(y), rms(one.fun) / rms(y)


say("  biased = a * no voltage + b * empty cuvette(t - tau), empty cuvette from paper 1 Fig. 3a")
rows = {}
for name, y, x0 in (("30 s (paper 1)", H1["30s"], H1["none"]), ("2 min (paper 2)", w2["red"], w2["blue"]),
                    ("4 min (paper 1)", H1["4min"], H1["none"]), ("5 min (paper 3)", h3["red"], h3["blue"])):
    p, res, one = fit_td(y, x0, H1["empty"])
    spread = []
    for d in (-1, 1):                                   # +-1 pixel (0.0133 ps) shift of the biased trace
        spread.append(fit_td(np.interp(tt - d * 0.0133, tt, y), x0, H1["empty"])[0][1])
    rows[name] = p
    say(f"  {name:16s}: a = {p[0]:.3f}, b = {p[1]:.3f} ({min(spread + [p[1]]):.3f}-{max(spread + [p[1]]):.3f}), "
        f"tau = {p[2]:+.3f} ps, residual {res:.0%}; one path {one:.0%}")
say("  delay-stage drift moves both terms: biased = a * no voltage(t - s) + b * empty(t - s - tau)")
for name, y, x0 in (("30 s (paper 1)", H1["30s"], H1["none"]), ("2 min (paper 2)", w2["red"], w2["blue"]),
                    ("4 min (paper 1)", H1["4min"], H1["none"]), ("5 min (paper 3)", h3["red"], h3["blue"])):
    r0 = lambda p: p[0] * shift(x0, p[2]) + p[1] * shift(H1["empty"], p[2]) - y
    s0 = min((least_squares(r0, [0.6, b0, t0]) for b0 in (0.05, 0.2) for t0 in (-0.1, 0.0, 0.1)), key=lambda s: s.cost)
    r1 = lambda p: p[0] * shift(x0, p[2]) + p[1] * shift(H1["empty"], p[2] + p[3]) - y
    s1 = min((least_squares(r1, [0.6, b0, t0, 0.0]) for b0 in (0.05, 0.2) for t0 in (-0.1, 0.0, 0.1)), key=lambda s: s.cost)
    say(f"  {name:16s}: tau = 0: b = {s0.x[1]:.3f}, s = {s0.x[2]:+.3f} ps, residual {rms(s0.fun) / rms(y):.0%}; "
        f"tau free: b = {s1.x[1]:.3f}, s = {s1.x[2]:+.3f} ps, tau = {s1.x[3]:+.3f} ps, residual {rms(s1.fun) / rms(y):.0%}")
b30, b4 = rows["30 s (paper 1)"][1], rows["4 min (paper 1)"][1]
b2 = rows["2 min (paper 2)"][1]
lin = b30 + (b4 - b30) * (2.0 - 0.5) / (4.0 - 0.5)
say(f"  check: 2 min lies between 30 s and 4 min if b30 < b2 < b4: {b30:.3f} < {b2:.3f} < {b4:.3f} -> "
    f"{'yes' if b30 < b2 < b4 else 'NO'}; linear-in-time interpolation gives {lin:.3f}")

# ---------------------------------------------------------------- (c) paper 1 spectra
p1 = render_page(REF / "paper 1.pdf", 3, 300)
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
fb = np.arange(0.20, 2.951, 0.05)
spec = {}
check = []
for name, p in SPEC.items():
    l, t, r, b = p["box"]; rgb = p1[t:b, l:r]; R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    fxs, exs = axis_map(*p["gx"]); fys, eys = axis_map(*p["gy"])
    say() if name == "3b" else None
    say(f"(c) paper 1 Fig. {name}: axis misfit {exs:.4f} THz / {eys:.3f} dB")
    for key, col in p["curves"].items():
        m = CM[col](R, G, B); m[:, :int(p["gx"][0][0]) + 2] = False; m[:, 815:] = False
        m[:int(p["gy"][0][0]) - 3, :] = False                 # drop the labels above the plot
        if col == "black":                               # numbers and arrows are black: keep the plot box only
            m[:, :int(p["gx"][0][0]) + 3] = False
            if name == "3b":                             # label "1" of the empty-cuvette curve sits at -7 dB, above the sample curve (max -12.5 dB)
                m[:118, :] = False
        xs = np.flatnonzero(m.any(axis=0))
        top = np.array([np.flatnonzero(m[:, x])[0] for x in xs], float)    # upper envelope (dB maximum)
        f, db = fxs(xs), fys(top)
        spec[key] = np.array([np.median(db[(f >= c - 0.025) & (f < c + 0.025)]) if np.any((f >= c - 0.025) & (f < c + 0.025))
                              else np.nan for c in fb])
        check.append((name, xs, top, col))
np.savez(OUT / "spectra.npz", f_thz=fb, **spec)
say(f"  same normalisation: sample in Fig. 3b vs no voltage in Fig. 4b, median |difference| "
    f"{np.nanmedian(np.abs(spec['sample'] - spec['none'])):.1f} dB over 0.2-2.95 THz")
rel = {k: spec[k] - spec["empty"] for k in ("30s", "4min")}
rel["none"] = spec["sample"] - spec["empty"]   # Fig. 3b copy of the no-voltage scan: in Fig. 4b the black arrow "1" lies on it at 0.8-1.25 THz
bad = (fb > 0.7) & (fb < 1.3)
say(f"  Fig. 4b black vs Fig. 3b sample: median |difference| {np.nanmedian(np.abs(spec['sample'] - spec['none'])[~bad]):.1f} dB outside 0.7-1.3 THz, "
    f"{np.nanmax(np.abs(spec['sample'] - spec['none'])[bad]):.1f} dB max inside (arrow)")
a30, b30_, a4, b4_ = rows["30 s (paper 1)"][0], rows["30 s (paper 1)"][1], rows["4 min (paper 1)"][0], rows["4 min (paper 1)"][1]
floor30 = 20 * np.log10(np.maximum(np.abs(b30_ - a30 * 10 ** (rel["none"] / 20)), 1e-6))
floor4 = 20 * np.log10(np.maximum(np.abs(b4_ - a4 * 10 ** (rel["none"] / 20)), 1e-6))
say("  ratio to the empty cuvette (dB), band medians; floor = |b - a|T|| from the time-domain fit (any phase)")
say("   f THz   no volt   30 s  floor30   4 min  floor4")
for i, c in enumerate(fb):
    if abs(c * 20 - round(c * 20)) < 1e-6 and (round(c * 20) % 2 == 0 or 1.85 <= c <= 2.2):
        say(f"   {c:4.2f}  {rel['none'][i]:8.1f} {rel['30s'][i]:7.1f} {floor30[i]:7.1f} {rel['4min'][i]:7.1f} {floor4[i]:7.1f}")
for lo, hi in ((1.5, 1.9), (1.95, 2.1), (2.3, 3.0)):
    s = (fb >= lo - 1e-9) & (fb <= hi + 1e-9)
    say(f"  {lo:.2f}-{hi:.2f} THz: 30 s {np.nanmean(rel['30s'][s]):.1f} dB vs floor {np.nanmean(floor30[s]):.1f} dB "
        f"(deficit {np.nanmean(floor30[s] - rel['30s'][s]):+.1f}); 4 min {np.nanmean(rel['4min'][s]):.1f} vs floor "
        f"{np.nanmean(floor4[s]):.1f} (field fraction {10 ** (np.nanmean(rel['4min'][s]) / 20):.2f})")
say("Reading: a positive deficit means the published spectrum lies below anything a constant undelayed")
say("component of the fitted size allows, whatever its phase.")

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
ax[0].imshow(rgb.astype(np.uint8))
for name, xs, top, col in check:
    if name == "4b":
        ax[0].plot(xs, top, ".", ms=0.6, color="yellow")
ax[0].set_title("paper 1 Fig. 4b, traced upper envelope")
ax[1].imshow(p2[1060:1570, 480:1240].astype(np.uint8))
for c, (xs, ys) in tr2.items():
    ax[1].plot(xs, ys, "-", lw=0.6, color="yellow")
ax[1].set_title("paper 2 Fig. 2, traced")
for k, c in (("none", "k"), ("30s", "m"), ("4min", "r")):
    ax[2].plot(fb, rel[k], "-", color=c, label=k)
ax[2].plot(fb, floor30, "m--", lw=0.8, label="30 s floor"); ax[2].plot(fb, floor4, "r--", lw=0.8, label="4 min floor")
ax[2].set_xlabel("f (THz)"); ax[2].set_ylabel("dB re empty cuvette"); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(OUT / "trace_check.png", dpi=120)
(OUT / "sim43.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

# ---------------------------------------------------------------- (d) where the components cancel
# Model ratio S/K = a T^q + b exp(-i w tau). |T| from the published spectra (band medians, valid to 3 THz);
# the phase of T from the traced waveforms, which are valid only below ~1.5 THz: fitted there as
# phi(f) = phi0 - 2 pi f dt and extrapolated. Dips sit where the two terms are in antiphase; depth needs
# |a T^q| ~ b. No parameter is fitted to the spectra.
say()
say("(d) Antiphase of the two terms, |T| from the published spectra, phase of T from the traced waveforms")
dt_s = (tt[1] - tt[0]) * 1e-12
win = (tt > -362.6) & (tt < -357.4)
fq = np.fft.rfftfreq(4096, dt_s) / 1e12
Kf = np.fft.rfft(np.where(win, H1["empty"], 0.0), 4096); Bf = np.fft.rfft(np.where(win, H1["none"], 0.0), 4096)
ph = np.unwrap(np.angle(Bf / Kf))
s = (fq > 0.2) & (fq < 1.4)
c1, c0 = np.polyfit(fq[s], ph[s], 1)
dt_T = -c1 / (2 * np.pi)
say(f"  tissue phase over 0.2-1.4 THz: linear, delay {dt_T:.3f} ps relative to the empty cuvette, offset {c0:+.2f} rad, "
    f"rms deviation {np.std(ph[s] - (c1 * fq[s] + c0)):.2f} rad")
magT = 10 ** (rel["none"] / 20)
phT = c1 * fb + c0
for name, key, q in (("30 s", "30 s (paper 1)", 1.0), ("4 min", "4 min (paper 1)", 1.0), ("4 min, q 1.09", "4 min (paper 1)", 1.09)):
    a_, b_, tau_ = rows[key]
    model = 20 * np.log10(np.abs(a_ * magT ** q * np.exp(1j * q * phT) + b_ * np.exp(-1j * 2 * np.pi * fb * tau_)))
    spec_key = "30s" if name == "30 s" else "4min"
    ap = [(2 * k + 1) * np.pi for k in range(3)]
    rel_phase = q * phT + 2 * np.pi * fb * tau_                     # phase of a T^q relative to b e^{-i w tau}
    fa = [float(np.interp(-x, rel_phase[::-1], fb[::-1])) for x in ap if -x >= rel_phase.min()]
    say(f"  {name:13s}: antiphase at {', '.join(f'{x:.2f}' for x in fa)} THz; model minima "
        + ", ".join(f"{fb[i]:.2f} THz {model[i]:.1f} dB" for i in range(1, len(fb) - 1)
                    if model[i] < model[i - 1] and model[i] <= model[i + 1] and model[i] < -15)
        + "; published minima " + ", ".join(f"{fb[i]:.2f} THz {rel[spec_key][i]:.1f} dB" for i in range(1, len(fb) - 1)
                                           if rel[spec_key][i] < rel[spec_key][i - 1] and rel[spec_key][i] <= rel[spec_key][i + 1]
                                           and rel[spec_key][i] < np.nanmedian(rel[spec_key]) - 5))
    for lo, hi in ((0.5, 0.8), (1.85, 2.15), (2.3, 3.0)):
        m_ = (fb >= lo - 1e-9) & (fb <= hi + 1e-9)
        say(f"      {lo:.2f}-{hi:.2f} THz: model {np.nanmean(model[m_]):6.1f} dB, published {np.nanmean(rel[spec_key][m_]):6.1f} dB")
(OUT / "sim43.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

# ---------------------------------------------------------------- figure for the paper
FIG = ROOT / "paper" / "figures"
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for axi, (name, key, spec_key, col, q) in zip(ax, (("30 s at 100 V", "30 s (paper 1)", "30s", "m", 1.0),
                                                  ("4 min at 100 V", "4 min (paper 1)", "4min", "r", 1.09))):
    a_, b_, tau_ = rows[key]
    model = 20 * np.log10(np.abs(a_ * magT ** q * np.exp(1j * q * phT) + b_ * np.exp(-1j * 2 * np.pi * fb * tau_)))
    axi.plot(fb, rel["none"], "k-", lw=1, label="no voltage (published)")
    axi.plot(fb, rel[spec_key], "-", color=col, lw=1.4, label=f"{name} (published)")
    axi.plot(fb, model, "--", color=col, lw=1, label=f"{a_:.2f}·T{'' if q == 1 else '^%.2f' % q} + {b_:.3f}·empty")
    axi.axvspan(1.95, 2.05, color="0.9", zorder=0)
    axi.set_title(f"paper 1, {name}", fontsize=10); axi.set_xlabel("frequency (THz)")
    axi.legend(fontsize=7, loc="lower left")
ax[0].set_ylabel("power re empty cuvette (dB)")
ax[0].text(2.0, -6, "30-s dip", ha="center", fontsize=7, color="0.4")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_paper1_spectra.{ext}", dpi=200)
