"""Sim 73a - Trace the waveforms of H4 (Hovhannisyan, Makaryan, Martirosyan, Reports NAS RA 124(2):16-23, 2024),
obtained in review round 33 (docs/references/paper 4.pdf, open on arar.sci.am).

H4 has Fig. 3a (page 4: empty 0.1 mm PTFE cuvette, blue; cuvette with a spinal-cord sample, red) and Fig. 4a / 4v
(page 5: two samples, no voltage, black; 6 min at 100 V, red), all on the same -363..-357 ps axis as H1.
Traced with thznerve.digitise.trace_track on the 300 dpi page, gridlines found with thznerve.digitise.gridlines.
Writes results/sim73/digitised.npz (on the H1 time grid of load_scans), results/sim73/trace_check.png and
results/sim73/trace.txt (axis misfits, pulse minima, arrival times, spectral ratios).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from thznerve.digitise import axis_map, render_page, trace_track
from thznerve.waveform_fit import load_scans

ROOT = Path(__file__).resolve().parents[2]
PDF = ROOT / "docs" / "references" / "paper 4.pdf"
OUT = ROOT / "results" / "sim73"
# crops (left, top, right, bottom) on the 300 dpi page; gridline pixels (crop coords, from gridlines()) -> values
PANELS = {
    "3a": dict(page=4, box=(450, 1600, 1040, 2250), gx=([225.0, 363.0, 500.5, 569.0], [-362, -360, -358, -357]),
               gy=([97.5, 158.0, 219.0, 279.2, 341.0, 401.5, 462.5, 523.5, 587.0], [6, 4, 2, 0, -2, -4, -6, -8, -10]),
               xlim=(166, 566), curves=("blue", "red")),
    "4a": dict(page=5, box=(600, 1680, 1060, 2150), gx=([72.0, 264.5, 392.0], [-363, -360, -358]),
               gy=([141.5, 180.0, 218.5, 256.0, 294.5, 333.0, 411.0], [0.5, 0, -0.5, -1, -1.5, -2, -3]),
               xlim=(76, 455), curves=("black", "red"), ymax=405, erase=[(170, 270, 200, 320), (255, 290, 280, 340)]),
    "4v": dict(page=5, box=(600, 2200, 1060, 2690), gx=([69.5, 265.0, 395.0], [-363, -360, -358]),
               gy=([33.0, 71.0, 110.0, 149.0, 188.0, 226.0, 304.0, 342.0, 383.5], [1.5, 1, 0.5, 0, -0.5, -1, -2, -2.5, -3]),
               xlim=(74, 455), curves=("black", "red"), ymax=395),
}
MASKS = {
    "blue": lambda R, G, B: (B > 170) & (R < 140) & (G > 90),
    "black": lambda R, G, B: (R < 110) & (G < 110) & (B < 110),
    "red": lambda R, G, B: (R > 200) & (G < 160) & (B < 160),
}
NAMES = {("3a", "blue"): "empty", ("3a", "red"): "sample3a", ("4a", "black"): "s1_none", ("4a", "red"): "s1_6min",
         ("4v", "black"): "s2_none", ("4v", "red"): "s2_6min"}


def panel(p, pages):
    l, t, r, b = p["box"]
    rgb = pages[p["page"]][t:b, l:r].astype(float)
    fx, ex = axis_map(*p["gx"]); fy, ey = axis_map(*p["gy"])
    R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    y0 = float(np.interp(0.0, p["gy"][1][::-1], p["gy"][0][::-1]))
    out = {}
    for c in p["curves"]:
        m = MASKS[c](R, G, B)
        m[:, :p["xlim"][0]] = False; m[:, p["xlim"][1]:] = False; m[p.get("ymax", m.shape[0]):] = False
        for x0, y0_, x1, y1 in p.get("erase", ()):  # number labels (x0, y0, x1, y1), crop pixels
            m[y0_:y1, x0:x1] = False
        xs, ys = trace_track(m, y_start=y0, max_jump=25)
        out[c] = (xs, ys, fx, fy)
    return out, ex, ey, rgb


def traces(px_shift=None):
    """All six traces on the H1 grid; px_shift maps a name to a time shift in pixels (tracing noise)."""
    W = {}
    for name, p in PANELS.items():
        for c, (xs, ys, fx, fy) in TR[name][0].items():
            k = NAMES[(name, c)]
            dx = (px_shift or {}).get(k, 0.0)
            W[k] = np.interp(TT, fx(xs + dx), fy(ys), left=np.nan, right=np.nan)
    valid = np.all([np.isfinite(v) for v in W.values()], axis=0)
    return {k: np.nan_to_num(v) for k, v in W.items()}, valid


TT = load_scans(ROOT)[0]
PAGES = {pg: render_page(PDF, pg, 300) for pg in (4, 5)}
TR = {name: panel(p, PAGES) for name, p in PANELS.items()}


def spectrum_db(y, valid, n=8192):
    dt = TT[1] - TT[0]
    Y = np.fft.rfft(np.where(valid, y, 0.0), n); f = np.fft.rfftfreq(n, dt)
    return f, 20 * np.log10(np.abs(Y) + 1e-12)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    say = lambda s="": (print(s, flush=True), lines.append(s))
    W, VALID = traces()
    say(f"H4 traced (300 dpi); fit window where all six curves are traced: {TT[VALID][0]:.2f} to {TT[VALID][-1]:.2f} ps")
    for name, (d, ex, ey, _) in TR.items():
        px = abs(np.diff(PANELS[name]["gx"][1])[-1] / np.diff(PANELS[name]["gx"][0])[-1])
        say(f"  Fig. {name}: axis misfit {ex:.3f} ps / {ey:.3f} a.u.; 1 pixel = {px:.4f} ps; traced columns " +
            ", ".join(f"{c} {len(v[0])}" for c, v in d.items()))
    say("")
    say("pulse minima (a.u. at ps) and rms over the window")
    for k, v in W.items():
        say(f"  {k:9s} min {v.min():6.2f} at {TT[np.argmin(v)]:.3f} ps; max {v.max():5.2f} at {TT[np.argmax(v)]:.3f}; rms {np.sqrt(np.mean(v[VALID] ** 2)):.3f}")
    say("")
    for s in ("s1", "s2"):
        a, b = W[f"{s}_none"], W[f"{s}_6min"]
        lags = np.arange(-400, 401); c = [np.dot(np.roll(a, k)[VALID], b[VALID]) for k in lags]
        say(f"  {s}: 6 min pulse vs no voltage, best cross-correlation shift {lags[int(np.argmax(c))] * (TT[1] - TT[0]):+.3f} ps (negative = earlier)")
        f, A = spectrum_db(a, VALID); _, B = spectrum_db(b, VALID)
        say(f"  {s}: 6 min / no voltage at 0.3 / 0.5 / 1.0 / 1.5 / 2.0 THz: " +
            " / ".join(f"{np.interp(x, f, B - A):+.1f}" for x in (0.3, 0.5, 1.0, 1.5, 2.0)) + " dB")
    f, E = spectrum_db(W["empty"], VALID); _, S = spectrum_db(W["sample3a"], VALID)
    say(f"  Fig. 3a sample / empty at 0.25 / 0.5 / 1.0 / 2.0 THz: " + " / ".join(f"{np.interp(x, f, S - E):+.1f}" for x in (0.25, 0.5, 1.0, 2.0)) +
        " dB (H4 text: about -5 dB at 0.25 THz, below -25 dB at 2 THz)")
    np.savez(OUT / "digitised.npz", t_ps=TT, valid=VALID, **W)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    for i, (name, (d, _, _, rgb)) in enumerate(TR.items()):
        ax[i].imshow(rgb.astype(np.uint8))
        for c, (xs, ys, _, _) in d.items():
            ax[i].plot(xs, ys, "y-", lw=0.6)
        ax[i].set_title(f"H4 Fig. {name}: traces in yellow")
    fig.savefig(OUT / "trace_check.png", dpi=150); plt.close(fig)
    (OUT / "trace.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
