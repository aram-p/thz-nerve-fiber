"""Sim 39 - Is the voltage-induced change in paper 3 an air bypass? Test it on the published waveforms,
with no tissue model and no thickness.

Sim 38 fitted a gas column (x, t) through modelled tissue to hand-read points; review round 4 found
the points biased (up to 5 dB at 2 THz) and t = 225-240 um implausible for a 100 um cuvette. Here:
  1. Render paper 3 page 2 (Fig. 1a waveforms, Fig. 1b spectra) at 300 dpi, trace the three curves by
     colour, calibrate both axes on the gridlines (helpers in thznerve.digitise, tested).
  2. Time domain: red (5 min at 100 V) = a * blue (no voltage) + b * black (empty cuvette)(t - tau).
     A bypass of fraction x at the empty-cuvette delay gives a = 1 - x, b = x, tau = 0. Controls:
     single-path changes (tissue only rescaled and shifted; tissue that lost water, i.e. blue moved
     towards black) must fit worse if the bypass reading is right.
  3. Frequency domain: spectrum of the fitted sum over spectrum of blue, against the Fig. 1b
     red - blue difference read by pixel (not by eye).
  4. Sim 38's tissue gas-column model refitted to the pixel points with chi^2, free and t <= 100 um.
Writes results/sim39/sim39.txt, results/sim39/digitised.npz and paper/figures/fig_bypass.png/pdf
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from thznerve.digitise import axis_map, gridlines, render_page, trace
from thznerve.physics.effmedium import tissue_eps
from thznerve.physics.materials import WATER_BROADBAND, WATER_KS1996, axoplasm, myelin, refractive_index

ROOT = Path(__file__).resolve().parents[2]
PDF = ROOT / "docs" / "references" / "paper 3.pdf"
OUT = ROOT / "results" / "sim39"; OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "paper" / "figures"
C0 = 299792458.0
lines = []


def say(s=""):
    print(s, flush=True); lines.append(s)


# ---- 1. digitise -----------------------------------------------------------------------------
page = render_page(PDF, 2, 300)
# crops (row, col, height, width) located on the 300 dpi page; legend boxes and text masked (crop coords)
PANELS = {
    "a": dict(box=(420, 198, 1054, 1043), legend=(515, 650, 950, 830), rows=(40, 842),
              xtick=(-364.0, 1.0), ytick=(6.0, -3.0), xlim=(125, 945)),
    "b": dict(box=(1508, 198, 1018, 1043), legend=(175, 630, 645, 785), rows=(40, 812),
              xtick=(0.0, 0.5), ytick=(0.0, -10.0), xlim=(150, 968)),
}


def ticks(px, v0, step):
    d = np.min(np.diff(px))
    return v0 + np.round((px - px[0]) / d) * step


def digitise(p):
    r0, c0, h, w = p["box"]
    rgb = page[r0:r0 + h, c0:c0 + w]
    gray = rgb.mean(axis=2)
    cx = gridlines(gray, 0, frac=0.7); cx = cx[(cx > 100) & (cx < 980)]
    ry = gridlines(gray, 1, frac=0.7); ry = ry[(ry > 20) & (ry < 900)]
    fx, ex = axis_map(cx, ticks(cx, *p["xtick"]))
    fy, ey = axis_map(ry, ticks(ry, *p["ytick"]))
    R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    masks = {"red": (R > 180) & (G < 120) & (B < 120),
             "blue": (B > 170) & (R < 120) & (G < 170),
             "black": (R < 70) & (G < 70) & (B < 70)}
    out = {}
    for k, m in masks.items():
        m = m.copy()
        x0, y0, x1, y1 = p["legend"]; m[y0:y1, x0:x1] = False
        m[:p["rows"][0]] = False; m[p["rows"][1]:] = False
        m[:, :p["xlim"][0]] = False; m[:, p["xlim"][1]:] = False
        xs, ys = trace(m)
        out[k] = (fx(xs), fy(ys))
    return out, (len(cx), ex, len(ry), ey)


wave, cal_a = digitise(PANELS["a"])
spec, cal_b = digitise(PANELS["b"])
say("Digitising paper 3 Fig. 1 (page 2, 300 dpi)")
for name, c in (("1a", cal_a), ("1b", cal_b)):
    say(f"  Fig. {name}: {c[0]} x gridlines (max misfit {c[1]:.3f}), {c[2]} y gridlines (max misfit {c[3]:.3f}) in data units")

tt = np.arange(-363.8, -357.2, 0.01)
W = {k: np.interp(tt, *v) for k, v in wave.items()}
ff = np.arange(0.15, 3.0001, 0.01)
S = {k: np.interp(ff, *v) for k, v in spec.items()}
np.savez(OUT / "digitised.npz", t_ps=tt, f_thz=ff, **{f"wave_{k}": v for k, v in W.items()},
         **{f"spec_{k}": v for k, v in S.items()})


# ---- 2. time-domain fits ----------------------------------------------------------------------
def shift(y, tau):
    return np.interp(tt - tau, tt, y, left=0.0, right=0.0)


rms = lambda r: float(np.sqrt(np.mean(np.square(r))))
sig = rms(W["red"])
say(f"\nTime domain, {tt[0]:.1f} to {tt[-1]:.1f} ps, signal rms of the 5-min trace {sig:.3f}")
fits = {}


def fit(name, fun, x0, show):
    r = least_squares(fun, x0)
    fits[name] = r
    say(f"  {name:58s} rms residual {rms(r.fun):.3f} ({rms(r.fun)/sig:.0%} of signal)  {show(r.x)}")
    return r


r2 = fit("bypass: a*blue + b*black(t - tau)",
         lambda q: q[0] * W["blue"] + q[1] * shift(W["black"], q[2]) - W["red"], [0.8, 0.2, 0.0],
         lambda q: f"a = {q[0]:.3f}, b = {q[1]:.3f}, tau = {q[2]:+.3f} ps")
fit("bypass, a + b = 1: (1 - x)*blue + x*black(t - tau)",
    lambda q: (1 - q[0]) * W["blue"] + q[0] * shift(W["black"], q[1]) - W["red"], [0.2, 0.0],
    lambda q: f"x = {q[0]:.3f}, tau = {q[1]:+.3f} ps")
fit("bypass, blue also shifted: a*blue(t - t1) + b*black(t - t2)",
    lambda q: q[0] * shift(W["blue"], q[2]) + q[1] * shift(W["black"], q[3]) - W["red"], [0.8, 0.2, 0.0, 0.0],
    lambda q: f"a = {q[0]:.3f}, b = {q[1]:.3f}, t1 = {q[2]:+.3f}, t2 = {q[3]:+.3f} ps")
fit("one path: a*blue(t - tau)",
    lambda q: q[0] * shift(W["blue"], q[1]) - W["red"], [1.0, 0.0],
    lambda q: f"a = {q[0]:.3f}, tau = {q[1]:+.3f} ps")
# best single-path shift of blue over a grid (least squares has local minima for a pure shift)
best = min((rms(np.dot(shift(W["blue"], s), W["red"]) / np.dot(shift(W["blue"], s), shift(W["blue"], s))
               * shift(W["blue"], s) - W["red"]), s) for s in np.arange(-1.5, 1.5, 0.01))
say(f"  {'one path, global search over the shift of blue':58s} rms residual {best[0]:.3f} ({best[0]/sig:.0%} of signal)  tau = {best[1]:+.2f} ps")
best = min((rms(np.dot(shift(W["black"], s), W["red"]) / np.dot(shift(W["black"], s), shift(W["black"], s))
               * shift(W["black"], s) - W["red"]), s) for s in np.arange(-1.5, 1.5, 0.01))
say(f"  {'one path, scaled and shifted empty-cuvette pulse':58s} rms residual {best[0]:.3f} ({best[0]/sig:.0%} of signal)  tau = {best[1]:+.2f} ps")

# ---- 3. frequency domain ----------------------------------------------------------------------
dt = (tt[1] - tt[0]) * 1e-12
N = 8192
fr = np.fft.rfftfreq(N, dt) / 1e12
win = np.hanning(len(tt)) ** 0.25           # soft edges only; the traces are ~0 at both ends


def spectrum_db(y):
    return 20 * np.log10(np.abs(np.fft.rfft(y * win, N)) + 1e-30)


a, b, tau = r2.x
model = a * W["blue"] + b * shift(W["black"], tau)
F_SEL = (fr > 0.2) & (fr < 3.0)
d_fft = {k: np.interp(ff, fr[F_SEL], spectrum_db(W[k])[F_SEL]) for k in W}
d_model = np.interp(ff, fr[F_SEL], spectrum_db(model)[F_SEL])
ref_fig = S["black"] - S["blue"]
ref_fft = d_fft["black"] - d_fft["blue"]
band = (ff >= 0.3) & (ff <= 2.0)
say("\nFrequency domain (dB). Check of the waveform digitising: tissue insertion loss black - blue from the")
say(f"  spectra of the traced 7-ps waveforms vs Fig. 1b: rms difference {rms((ref_fft - ref_fig)[band]):.1f} dB over 0.3-2 THz")
say("  (Fig. 1b was computed from a longer record; outside 0.3-2 THz the traced window's noise floor dominates)")
red_blue = S["red"] - S["blue"]
pred = d_model - d_fft["blue"]
say("  f (THz)     " + " ".join(f"{v:6.2f}" for v in (0.3, 0.45, 0.56, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)))
pick = lambda y: " ".join(f"{np.interp(v, ff, y):6.1f}" for v in (0.3, 0.45, 0.56, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0))
say("  Fig.1b red-blue (pixel) " + pick(red_blue))
say("  FFT red-blue (traces)   " + pick(d_fft["red"] - d_fft["blue"]))
say("  bypass fit / blue (FFT) " + pick(pred))
fft_rb = d_fft["red"] - d_fft["blue"]
bad = np.flatnonzero((ff > 1.0) & (np.abs(fft_rb - red_blue) > 1.5))
F_OK = ff[bad[0]] if len(bad) else ff[-1]
say(f"  the traced 7-ps window reproduces the published red - blue within 1.5 dB up to {F_OK:.2f} THz; above,")
say("  the window and tracing noise dominate, so the fit is judged only below that frequency")
ok = ff <= F_OK
say(f"  bypass fit vs same-window spectrum (like with like), 0.3-{F_OK:.2f} THz: rms {rms((pred - fft_rb)[ok & (ff >= 0.3)]):.1f} dB; "
    f"vs Fig. 1b: rms {rms((pred - red_blue)[ok & (ff >= 0.3)]):.1f} dB")
lo = (ff > 0.4) & (ff < 0.9)
say(f"  dip: Fig. 1b {ff[lo][np.argmin(red_blue[lo])]:.2f} THz ({red_blue[lo].min():.1f} dB); "
    f"bypass fit {ff[lo][np.argmin(pred[lo])]:.2f} THz ({pred[lo].min():.1f} dB)")
hi = (ff > 1.5) & (ff < 2.6)
d2 = np.diff(red_blue[hi])
say(f"  1.5-2.6 THz on Fig. 1b: red - blue rises from {red_blue[hi][0]:.1f} to {red_blue[hi].max():.1f} dB; "
    f"largest local drop {max(0.0, -np.min(np.cumsum(d2) - np.maximum.accumulate(np.cumsum(d2)))):.1f} dB (no dip near 2 THz)")


# ---- 3b. review round 5: what else fits, how much is noise, where the dip goes ---------------
pre = tt < -362.8
jit = max(float(np.std(np.polyval(np.polyfit(tt[pre], W[k][pre], 1), tt[pre]) - W[k][pre])) for k in W)
say(f"\nTracing jitter before the pulse (detrended rms, worst trace) {jit:.4f}: the 28 % residual is "
    f"{rms(r2.fun)/jit:.0f}x that, so the two-path sum is rejected as an exact description")
FR, FB, FK = (np.fft.rfft(W[k], N) for k in ("red", "blue", "black"))
w_rad = 2 * np.pi * fr * 1e12
vb = (fr >= 0.15) & (fr <= F_OK)
T = FB / FK                                         # measured tissue transfer (cuvette with tissue / empty)
phi = np.unwrap(np.angle(T[vb]))
n_tr = 1 - C0 * phi / (w_rad[vb] * 100e-6)          # numpy FFT: a delay multiplies by exp(-i w t)
say("  index from the traced tissue phase at the stated d = 100 um: " + ", ".join(
    f"{v:.1f} THz {np.interp(v, fr[vb], n_tr):.2f}" for v in (0.3, 0.6, 1.0, 1.5))
    + " (water-based tissue model 1.95-2.2): the slice delays like ~2x its stated thickness")
Tp = lambda q: np.abs(T[vb]) ** q * np.exp(1j * q * phi)   # tissue path scaled to q x its optical thickness
E_R = float(np.sum(np.abs(FR[vb]) ** 2))


def band_fit(name, model, x0):
    r = least_squares(lambda q: np.concatenate([(model(q) - FR[vb]).real, (model(q) - FR[vb]).imag]), x0)
    m = model(r.x)
    res = np.sqrt(np.sum(np.abs(m - FR[vb]) ** 2) / E_R)
    g = 20 * np.log10(np.abs(m) / np.abs(FB[vb]))
    sel = (fr[vb] > 0.35) & (fr[vb] < 1.0)
    say(f"  {name:62s} residual {res:.0%}  dip {fr[vb][sel][np.argmin(g[sel])]:.2f} THz {g[sel].min():5.1f} dB  "
        + ", ".join(f"{v:.3f}" for v in r.x))
    return r.x


say(f"Frequency-domain fits over 0.15-{F_OK:.2f} THz (residual = band-limited time-domain rms / signal), "
    "dip of model/no-voltage; reported dip 0.56 THz, -11.4 dB")
ph = lambda tau: np.exp(-1j * w_rad[vb] * tau * 1e-12)
band_fit("air bypass: a*B + b*K(tau)                         [a b tau]", lambda q: q[0] * FB[vb] + q[1] * FK[vb] * ph(q[2]), [0.6, 0.2, 0.0])
band_fit("air bypass, a + b = 1                              [x tau]", lambda q: (1 - q[0]) * FB[vb] + q[0] * FK[vb] * ph(q[1]), [0.2, 0.0])
m4 = band_fit("tissue path q x thicker + air bypass               [a b tau q]", lambda q: q[0] * FK[vb] * Tp(q[3]) + q[1] * FK[vb] * ph(q[2]), [0.6, 0.2, 0.0, 1.1])
for qf in (1.1, 1.2):
    band_fit(f"tissue path {qf} x thicker + air bypass               [a b tau]", lambda q, qf=qf: q[0] * FK[vb] * Tp(qf) + q[1] * FK[vb] * ph(q[2]), [0.6, 0.2, 0.0])
for s_ in (0.1, 0.2):
    band_fit(f"bypass crosses {s_:.0%} of the tissue                   [a b tau]", lambda q, s_=s_: q[0] * FB[vb] + q[1] * FK[vb] * Tp(s_) * ph(q[2]), [0.6, 0.2, 0.0])
band_fit("one path: tissue q x thicker, no bypass            [a q]", lambda q: q[0] * FK[vb] * Tp(q[1]), [0.6, 0.9])
say("Dip of a pure bypass, (1 - x)*B + x*K, as the bypass grows (model/no-voltage, 0.35-1 THz):")
sel = (fr[vb] > 0.35) & (fr[vb] < 1.0)
for qf in (1.0, 1.1, 1.2):
    row = []
    for x in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35):
        g = 20 * np.log10(np.abs((1 - x) * FK[vb] * Tp(qf) + x * FK[vb]) / np.abs(FB[vb]))
        row.append(f"x {x:.2f}: {fr[vb][sel][np.argmin(g[sel])]:.2f} THz {g[sel].min():5.1f} dB")
    say(f"  tissue path {qf} x: " + "; ".join(row))
hb = (ff >= 1.5) & (ff <= 3.0)
rk = S["red"] - S["black"]
say("Fig. 1b alone, 1.5-3 THz, 5 min at 100 V minus empty cuvette (dB): "
    + ", ".join(f"{np.interp(v, ff, rk):.1f}" for v in (1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0))
    + f" -> flat, as an undelayed component that dominates where tissue loss is 30-50 dB; amplitude "
    f"{10 ** (rk[hb].max() / 20):.2f}-{10 ** (rk[hb].min() / 20):.2f} of the empty-cuvette field")


# ---- 4. Sim 38 gas column through modelled tissue, refitted to the pixel points ----------------
def n_tissue(f, water):
    e_par, _ = tissue_eps(axoplasm(water)(f), myelin(water)(f), water(f), 0.7, 0.6)
    return refractive_index(e_par)


def ratio_db(f_thz, x, t_um, water):
    f = np.asarray(f_thz) * 1e12
    n = n_tissue(f, water); k0 = 2 * np.pi * f / C0
    g = (1 - x) + x * np.exp(1j * k0 * (n.real - 1) * t_um * 1e-6) * np.exp(-k0 * n.imag * t_um * 1e-6)
    return 20 * np.log10(np.abs(g))


FP = np.array([0.3, 0.45, 0.56, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0])
DP = np.interp(FP, ff, red_blue)
say("\nSim 38 gas column through modelled tissue (E || fibres), refitted to pixel points at "
    + ", ".join(f"{v:g}" for v in FP) + " THz, sigma 2 dB")
for w, name in ((WATER_KS1996, "KS1996"), (WATER_BROADBAND, "broadband")):
    for tmax in (600.0, 100.0):
        r = least_squares(lambda q: (ratio_db(FP, q[0], q[1], w) - DP) / 2.0, x0=[0.15, min(200.0, tmax - 1)],
                          bounds=([0.0, 10.0], [1.0, tmax]))
        chi2 = float(np.sum(r.fun ** 2))
        say(f"  {name:9s} t <= {tmax:3.0f} um: x = {r.x[0]:.3f}, t = {r.x[1]:.0f} um, chi2 = {chi2:.0f} on {len(FP) - 2} dof, "
            f"rms {rms(2 * r.fun):.1f} dB")
        if name == "KS1996" and tmax == 600.0:
            gas_fit = r.x

say("\nReading")
say("  The 5-min trace contains ~21 % of the empty-cuvette pulse at the empty-cuvette delay (+-0.02 ps): a")
say("  component that crosses at most ~20 % of the tissue's optical thickness. Without it nothing fits (94 %).")
say("  With it, the tissue path must also change: 0.59 x in amplitude, or ~10 % more optical thickness, which")
say("  gives the dip at 0.57 THz with 10.7 of 11.4 dB. Neither the cause of the undelayed component (gas,")
say("  slice pulling away from a window, beam or sample movement) nor of the tissue change is identified;")
say("  the sum is rejected as an exact description (residual 8x tracing jitter).")

# ---- figure -----------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
ax[0].plot(tt, W["black"], "k", lw=0.8, label="empty cuvette (traced)")
ax[0].plot(tt, W["blue"], "C0", lw=0.9, label="no voltage (traced)")
ax[0].plot(tt, W["red"], "C3", lw=1.2, label="5 min at 100 V (traced)")
ax[0].plot(tt, model, "k--", lw=1.0, label=f"{a:.3f} × no voltage + {b:.3f} × empty")
ax[0].set_xlim(-363, -357.2); ax[0].set_ylim(-2.5, 2.0)
ax[0].set_xlabel("Time (ps)"); ax[0].set_ylabel("THz field (a.u., paper 3 Fig. 1a)")
ax[0].legend(fontsize=7, loc="lower right")
ax[1].plot(ff, red_blue, "C3", lw=1.4, label="reported: 5 min at 100 V − no voltage (Fig. 1b, traced)")
ax[1].plot(ff[ok], fft_rb[ok], "C3:", lw=0.9, label="same, from the traced 7-ps waveforms")
ax[1].plot(ff[ok], pred[ok], "k--", lw=1.0, label=f"{a:.3f} × no voltage + {b:.3f} × empty cuvette (waveform fit)")
g4 = 20 * np.log10(np.abs(m4[0] * FK[vb] * Tp(m4[3]) + m4[1] * FK[vb] * ph(m4[2])) / np.abs(FB[vb]))
ax[1].plot(fr[vb], g4, "C2-.", lw=1.0, label=f"{m4[0]:.2f} × tissue path {m4[3]:.2f}× thicker + {m4[1]:.2f} × empty cuvette (fit)")
ax[1].axvspan(F_OK, 3.0, color="0.92", lw=0, zorder=0)
ax[1].plot(ff, ratio_db(ff, *gas_fit, WATER_KS1996), "C7", lw=0.8,
           label=f"withdrawn: gas column in modelled tissue, x = {gas_fit[0]:.2f}, t = {gas_fit[1]:.0f} µm")
ax[1].axhline(0, color="0.6", lw=0.6)
ax[1].set_xlabel("Frequency (THz)"); ax[1].set_ylabel("Δ power spectrum (dB)")
ax[1].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, frameon=False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_bypass.{ext}", dpi=200)
(OUT / "sim39.txt").write_text("\n".join(lines), encoding="utf-8")
