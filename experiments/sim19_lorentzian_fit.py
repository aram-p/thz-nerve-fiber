"""Sim 19 — Lorentzian fit of the resonance peak.

Fits a Lorentzian + linear background to the sim 1 spectrum
(peak |E|-at-node vs frequency) to quantify the 0.6 THz feature:
centre frequency f₀, FWHM γ, quality factor Q = f₀ / γ, peak height.

Also produces the same fit on sim 18's data if available (the E ∥ fibre
configuration).
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from thznerve.plots.style import DOUBLE_COL_MM, OKABE_ITO, apply_thesis_style, save_figure

REPO_ROOT = Path(__file__).resolve().parent.parent


def lorentzian(f, A, f0, gamma, b, m):
    """Lorentzian + linear baseline."""
    return A * (gamma / 2) ** 2 / ((f - f0) ** 2 + (gamma / 2) ** 2) + b + m * f


def _load_spectrum(csv_path: Path, freq_col: str, peak_col: str):
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    f_thz = np.array([float(r[freq_col]) for r in rows])
    e_peak = np.array([float(r[peak_col]) for r in rows])
    return f_thz, e_peak


def _fit_lorentzian(f, y, *, f0_guess: float, window: float = 0.4,
                     smooth_window: int = 3):
    """Fit Lorentzian + linear baseline to data in a window around f0_guess.

    Data is first smoothed with a centred moving average to reduce
    point-to-point noise so the fit converges to a meaningful Q.
    """
    mask = (f >= f0_guess - window) & (f <= f0_guess + window)
    f_fit = f[mask]
    y_fit = y[mask]
    if smooth_window > 1:
        k = np.ones(smooth_window) / smooth_window
        pad = smooth_window // 2
        y_fit = np.convolve(np.pad(y_fit, pad, mode="edge"), k, mode="valid")
    A0 = float(np.max(y_fit) - np.min(y_fit))
    p0 = [A0, f0_guess, 0.15, float(np.min(y_fit)), 0.0]
    bounds = (
        [1e-4,             f0_guess - 0.20, 0.03,        -10.0, -50.0],
        [10 * (A0 + 0.5),  f0_guess + 0.20, window * 0.8, 10.0,  50.0],
    )
    try:
        popt, pcov = curve_fit(lorentzian, f_fit, y_fit, p0=p0,
                               bounds=bounds, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
    except Exception as e:
        print(f"  fit failed: {e}")
        return None, None, None
    return popt, perr, (f_fit, y_fit)


def _format_fit(popt, perr):
    A, f0, gamma, b, m = popt
    sA, sf, sg, sb, sm = perr
    Q = f0 / abs(gamma) if abs(gamma) > 0 else float("inf")
    return (f"  f0 = {f0:.4f} ± {sf:.4f} THz\n"
            f"  γ (FWHM) = {abs(gamma):.4f} ± {sg:.4f} THz\n"
            f"  Q = f0/γ = {Q:.2f}\n"
            f"  peak amplitude A = {A:.3f} ± {sA:.3f}")


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.5)
    fig, axes = plt.subplots(1, 2)

    # Sim 1 — peak |E| on axis vs frequency
    sim1_csv = REPO_ROOT / "results" / "sim1" / "spectrum.csv"
    f1, e1 = _load_spectrum(sim1_csv, "frequency_THz", "peak_E_node")
    print("Sim 1 (E _perp_ fibre, axial sampling) — Lorentzian fit around 0.6 THz")
    popt, perr, win1 = _fit_lorentzian(f1, e1, f0_guess=0.63, window=0.35)
    if popt is not None:
        print(_format_fit(popt, perr))
        f_grid = np.linspace(f1.min(), f1.max(), 400)
        axes[0].plot(f1, e1, "o", color=OKABE_ITO[0], ms=4, label="Sim 1 data")
        axes[0].plot(f_grid, lorentzian(f_grid, *popt), "-", color=OKABE_ITO[5],
                     lw=1.4, label=f"Lorentzian fit: f₀={popt[1]:.3f} THz, Q={popt[1]/abs(popt[2]):.1f}")
        if win1 is not None:
            axes[0].axvspan(win1[0].min(), win1[0].max(), color="grey", alpha=0.08)
    axes[0].axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.65)
    axes[0].set_xlabel("Frequency (THz)")
    axes[0].set_ylabel("peak |E| in node region")
    axes[0].set_title("Sim 1 — E _perp_ fibre, axial sampling")
    axes[0].legend(loc="best", fontsize=7, framealpha=0.92)
    axes[0].grid(True, alpha=0.3)

    # Sim 18 if available
    sim18_csv = REPO_ROOT / "results" / "sim18" / "spectrum.csv"
    if sim18_csv.exists():
        f18, e18 = _load_spectrum(sim18_csv, "frequency_THz", "peak_E_node_annulus")
        print("\nSim 18 (E ∥ fibre, annulus sampling) — Lorentzian fit search")
        # Pick peak location from the data
        peak_idx = int(np.argmax(e18))
        f0_guess = float(f18[peak_idx])
        print(f"  data peak at f = {f0_guess:.3f} THz, |E| = {e18[peak_idx]:.3f}")
        popt18, perr18, win18 = _fit_lorentzian(f18, e18, f0_guess=f0_guess, window=0.35)
        if popt18 is not None:
            print(_format_fit(popt18, perr18))
            f_grid = np.linspace(f18.min(), f18.max(), 400)
            axes[1].plot(f18, e18, "o", color=OKABE_ITO[0], ms=4, label="Sim 18 data")
            axes[1].plot(f_grid, lorentzian(f_grid, *popt18), "-",
                         color=OKABE_ITO[5], lw=1.4,
                         label=f"f₀={popt18[1]:.3f} THz, Q={popt18[1]/abs(popt18[2]):.1f}")
        axes[1].set_xlabel("Frequency (THz)")
        axes[1].set_ylabel("peak |E| in node annulus")
        axes[1].set_title("Sim 18 — E ∥ fibre, annular sampling")
        axes[1].axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.65)
        axes[1].axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.65)
        axes[1].legend(loc="best", fontsize=7, framealpha=0.92)
        axes[1].grid(True, alpha=0.3)
    else:
        axes[1].text(0.5, 0.5, "Sim 18 data not yet available",
                     transform=axes[1].transAxes, ha="center", va="center",
                     fontsize=10, color="0.45")
        axes[1].set_axis_off()

    fig.suptitle("Sim 19 — Lorentzian fits to the resonance peaks")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim19_lorentzian_fit")
    print(f"\nwrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
