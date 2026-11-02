"""Sim 42 - Does paper 1 (Armenian J. Phys. 2024, 70 um slice, 18 mm gap) show the same undelayed
component as paper 3 (Sim 39), and how does it change between 30 s and 4 min of bias?

Review round 6 asked for this replication: paper 1 is the source of the drift (0.7 -> 0.5 THz) and of
the 2 THz feature. Paper 1 page 3 has Fig. 3a (empty cuvette, blue; sample, black) and Fig. 4a (no
voltage, black; 30 s at 100 V, magenta; 4 min, red) on the same -363..-357 ps axis and the same a.u.
Traces are followed with thznerve.digitise.trace_track (labels and arrows cross the curves).
Tracing noise is estimated by redoing every fit with each trace shifted by +-1 pixel in time.
Writes results/sim42/sim42.txt, results/sim42/digitised.npz, paper/figures/fig_paper1_bypass.png/pdf
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from thznerve.digitise import axis_map, render_page, trace_track

ROOT = Path(__file__).resolve().parents[2]
PDF = ROOT / "docs" / "references" / "paper 1.pdf"
OUT = ROOT / "results" / "sim42"; OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "paper" / "figures"
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


page = render_page(PDF, 3, 300)
# crops (left, top, right, bottom) on the 300 dpi page; gridline pixels (crop coords) -> values,
# read off the rendered crops (vertical lines every 2 ps plus the right frame at -357 ps)
PANELS = {
    "3a": dict(box=(330, 1420, 960, 2010), gx=([205.5, 356.5, 508.0, 583.0], [-362, -360, -358, -357]),
               gy=([37.5, 104.5, 171.0, 305.0, 371.5, 438.5, 505.5, 574.5], [6, 4, 2, -2, -4, -6, -8, -10]),
               xlim=(140, 580)),
    "4a": dict(box=(330, 2390, 960, 2990), gx=([201.5, 351.5, 502.0, 577.0], [-362, -360, -358, -357]),
               gy=([42.5, 96.5, 150.5, 204.5, 259.0, 312.5, 366.0, 420.0, 474.0, 528.5, 585.0],
                   [2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2, -2.5, -3]),
               xlim=(130, 575)),
}
MASKS = {
    "blue": lambda R, G, B: (B > 170) & (R < 120) & (G > 90) & (G < 200),
    "black": lambda R, G, B: (R < 90) & (G < 90) & (B < 90),
    "magenta": lambda R, G, B: (R > 170) & (B > 170) & (G < 120),
    "red": lambda R, G, B: (R > 180) & (G < 110) & (B < 110),
}


def panel(name, curves):
    p = PANELS[name]
    l, t, r, b = p["box"]
    rgb = page[t:b, l:r]
    fx, ex = axis_map(*p["gx"]); fy, ey = axis_map(*p["gy"])
    R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    y0 = float(np.interp(0.0, p["gy"][1][::-1], p["gy"][0][::-1]))
    out = {}
    for c in curves:
        m = MASKS[c](R, G, B)
        m[:, :p["xlim"][0]] = False; m[:, p["xlim"][1]:] = False
        xs, ys = trace_track(m, y_start=y0, max_jump=25)
        out[c] = (xs, ys, fx, fy)
    return out, ex, ey, rgb


f3, ex3, ey3, rgb3 = panel("3a", ("blue", "black"))
f4, ex4, ey4, rgb4 = panel("4a", ("black", "magenta", "red"))
say("Paper 1 Fig. 3a / 4a traced (page 3, 300 dpi)")
say(f"  axis maps: Fig. 3a misfit {ex3:.3f} ps / {ey3:.3f}; Fig. 4a {ex4:.3f} ps / {ey4:.3f}; "
    f"1 pixel = {abs(np.diff(PANELS['4a']['gx'][1])[-1] / np.diff(PANELS['4a']['gx'][0])[-1]):.4f} ps")

tt = np.arange(-362.8, -357.2, 0.005)


def series(d, c, px_shift=0.0):
    xs, ys, fx, fy = d[c]
    # gaps (curve hidden under another) interpolated; outside the traced extent NaN, then 0 with the
    # fit restricted to where every curve was traced (round 9: held ends had entered the fit)
    return np.interp(tt, fx(xs + px_shift), fy(ys), left=np.nan, right=np.nan)


def finish(Wn):
    valid = np.all([np.isfinite(v) for v in Wn.values()], axis=0)
    return {k: np.nan_to_num(v) for k, v in Wn.items()}, valid


W, VALID = finish({"empty": series(f3, "blue"), "sample3a": series(f3, "black"), "none": series(f4, "black"),
                   "30s": series(f4, "magenta"), "4min": series(f4, "red")})
say(f"  fit window (every curve traced): {tt[VALID][0]:.2f} to {tt[VALID][-1]:.2f} ps")
np.savez(OUT / "digitised.npz", t_ps=tt, **W)
say(f"  same scale check, no-voltage trace in Fig. 3a vs Fig. 4a: minima {W['sample3a'].min():.2f} vs {W['none'].min():.2f} "
    f"at {tt[np.argmin(W['sample3a'])]:.2f} vs {tt[np.argmin(W['none'])]:.2f} ps; rms difference "
    f"{np.sqrt(np.mean((W['sample3a'] - W['none']) ** 2)):.3f}")

rms = lambda r: float(np.sqrt(np.mean(np.square(r))))


def shift(y, tau):
    return np.interp(tt - tau, tt, y, left=0.0, right=0.0)


def fits(Wd):
    """Return dict of (params, residual fraction) for the time-domain models."""
    out = {}
    for key in ("30s", "4min"):
        y = Wd[key]; s = rms(y[VALID])
        r = least_squares(lambda q: (q[0] * Wd["none"] + q[1] * shift(Wd["empty"], q[2]) - y)[VALID], [0.8, 0.1, 0.0])
        out[(key, "bypass")] = (r.x, rms(r.fun) / s)
        best = min((rms((np.dot(shift(Wd["none"], t), y) / np.dot(shift(Wd["none"], t), shift(Wd["none"], t))
                         * shift(Wd["none"], t) - y)[VALID]), t) for t in np.arange(-1.0, 1.0, 0.005))
        out[(key, "one path")] = ((best[1],), best[0] / s)
    return out


base = fits(W)
# tracing noise: every trace shifted by -1, 0, +1 pixel independently (27 combinations)
var = []
for s_e in (-1, 0, 1):
    for s_n in (-1, 0, 1):
        for s_v in (-1, 0, 1):
            Wv = dict(W, empty=series(f3, "blue", s_e), none=series(f4, "black", s_n),
                      **{"30s": series(f4, "magenta", s_v), "4min": series(f4, "red", s_v)})
            var.append(fits({k: np.nan_to_num(v) for k, v in Wv.items()}))
say("\nTime domain: biased = a * no voltage + b * empty cuvette(t - tau); spread over +-1 pixel shifts of each trace")
for key in ("30s", "4min"):
    (a_, b_, t_), res = base[(key, "bypass")]
    bs = [v[(key, "bypass")][0][1] for v in var]; ts = [v[(key, "bypass")][0][2] for v in var]
    rs = [v[(key, "bypass")][1] for v in var]; os_ = [v[(key, "one path")][1] for v in var]
    say(f"  {key:4s}: a = {a_:.3f}, b = {b_:.3f} ({min(bs):.3f}-{max(bs):.3f}), tau = {t_:+.3f} ps ({min(ts):+.3f} to {max(ts):+.3f}), "
        f"residual {res:.0%} ({min(rs):.0%}-{max(rs):.0%}); one path (no voltage scaled+shifted) {base[(key, 'one path')][1]:.0%} "
        f"({min(os_):.0%}-{max(os_):.0%})")

# frequency domain: tissue transfer T = none/empty, tissue path raised to power q
dt = (tt[1] - tt[0]) * 1e-12; N = 16384
fr = np.fft.rfftfreq(N, dt) / 1e12; w = 2 * np.pi * fr * 1e12
F = {k: np.fft.rfft(v, N) for k, v in W.items()}
vb = (fr >= 0.15) & (fr <= 1.5)
T = F["none"][vb] / F["empty"][vb]
phi = np.unwrap(np.angle(T))
Tp = lambda q: np.abs(T) ** q * np.exp(1j * q * phi)
ph = lambda tau: np.exp(-1j * w[vb] * tau * 1e-12)
say(f"  tissue index from the traced phase at the stated 70 um: "
    + ", ".join(f"{v:.1f} THz {np.interp(v, fr[vb], 1 - 299792458.0 * phi / (w[vb] * 70e-6)):.2f}" for v in (0.3, 0.6, 1.0)))
say("\nFrequency domain, 0.15-1.5 THz: biased = A * tissue path^q + b * empty(tau); dip of model / no voltage (0.35-1.2 THz)")
sel = (fr[vb] > 0.35) & (fr[vb] < 1.2)
models = {}
for key in ("30s", "4min"):
    FR = F[key][vb]; E = float(np.sum(np.abs(FR) ** 2))
    meas = 20 * np.log10(np.abs(FR) / np.abs(F["none"][vb]))
    say(f"  {key}: traced spectrum biased/no voltage: dip {fr[vb][sel][np.argmin(meas[sel])]:.2f} THz {meas[sel].min():.1f} dB")
    for name, model, x0 in (
            ("air component only, tissue path unchanged", lambda q: q[0] * F["none"][vb] + q[1] * F["empty"][vb] * ph(q[2]), [0.8, 0.1, 0.0]),
            ("air component + tissue path q x thicker", lambda q: q[0] * F["empty"][vb] * Tp(q[3]) + q[1] * F["empty"][vb] * ph(q[2]), [0.8, 0.1, 0.0, 1.0]),
            ("tissue path change only (no air component)", lambda q: q[0] * F["empty"][vb] * Tp(q[1]), [0.8, 1.0])):
        r = least_squares(lambda q: np.concatenate([(model(q) - FR).real, (model(q) - FR).imag]), x0)
        m = model(r.x); res = np.sqrt(np.sum(np.abs(m - FR) ** 2) / E)
        g = 20 * np.log10(np.abs(m) / np.abs(F["none"][vb]))
        models[(key, name)] = (r.x, g)
        say(f"    {name:44s} residual {res:.0%}  dip {fr[vb][sel][np.argmin(g[sel])]:.2f} THz {g[sel].min():6.1f} dB  "
            + ", ".join(f"{v:.3f}" for v in r.x))

say("\nReading")
say("  Paper 1 shows the same signature as paper 3 if b > 0 at tau ~ 0 with one-path fits much worse; how b and")
say("  q change from 30 s to 4 min tests whether the dip moves because the undelayed component grows or because")
say("  the tissue path changes. Above ~1.5 THz the 75 px/ps tracing cannot test the 2 THz feature.")

# figure: traces and fits
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for i, key in enumerate(("30s", "4min")):
    (a_, b_, t_), res = base[(key, "bypass")]
    ax[i].plot(tt, W["empty"] / 4, "C0", lw=0.8, label="empty cuvette / 4 (Fig. 3a)")
    ax[i].plot(tt, W["none"], "k", lw=0.9, label="no voltage (Fig. 4a)")
    ax[i].plot(tt, W[key], "C3" if key == "4min" else "m", lw=1.3, label=f"{key} at 100 V (Fig. 4a)")
    ax[i].plot(tt, a_ * W["none"] + b_ * shift(W["empty"], t_), "k--", lw=1.0,
               label=f"{a_:.2f} × no voltage + {b_:.3f} × empty")
    ax[i].set_xlim(-362.5, -357.5); ax[i].set_ylim(-2.2, 1.5)
    ax[i].set_xlabel("Time (ps)"); ax[i].set_title(f"paper 1, {key} at 100 V", fontsize=9)
    ax[i].legend(fontsize=7, loc="lower right")
ax[0].set_ylabel("THz field (a.u.)")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_paper1_bypass.{ext}", dpi=200)
# tracing check overlay (not a paper figure)
fig2, ax2 = plt.subplots(1, 2, figsize=(10, 5))
for a2, (rgb, d) in zip(ax2, ((rgb3, f3), (rgb4, f4))):
    a2.imshow(rgb.astype(np.uint8))
    for c, (xs, ys, _, _) in d.items():
        a2.plot(xs, ys, "y-", lw=0.6)
fig2.savefig(OUT / "trace_check.png", dpi=150)
(OUT / "sim42.txt").write_text("\n".join(lines), encoding="utf-8")
