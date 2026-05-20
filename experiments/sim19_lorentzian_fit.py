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


def _plot_one(ax, f, y, *, title, f0_guess, label_data="data", window=0.35,
              colour=OKABE_ITO[0]):
    popt, perr, win = _fit_lorentzian(f, y, f0_guess=f0_guess, window=window)
    ax.plot(f, y, "o", color=colour, ms=4, label=label_data)
    if popt is not None:
        f_grid = np.linspace(f.min(), f.max(), 400)
        ax.plot(f_grid, lorentzian(f_grid, *popt), "-", color=OKABE_ITO[5], lw=1.4,
                label=fr"$f_0$={popt[1]:.3f} THz, Q={popt[1]/abs(popt[2]):.1f}")
        if win is not None:
            ax.axvspan(win[0].min(), win[0].max(), color="grey", alpha=0.08)
    ax.axvline(0.6, color="red", ls=":", lw=0.7, alpha=0.65)
    ax.axvline(2.0, color="red", ls=":", lw=0.7, alpha=0.65)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("peak |E|")
    ax.set_title(title, fontsize=8)
    ax.legend(loc="best", fontsize=6.5, framealpha=0.92)
    ax.grid(True, alpha=0.3)
    return popt, perr


def main() -> None:
    apply_thesis_style(width_mm=DOUBLE_COL_MM, aspect=0.55)
    fig, axes = plt.subplots(2, 2)

    sim1_csv = REPO_ROOT / "results" / "sim1" / "spectrum.csv"
    sim18_csv = REPO_ROOT / "results" / "sim18" / "spectrum.csv"

    # Sim 1 — E perp fibre, axial sampling
    f1, e1_axial = _load_spectrum(sim1_csv, "frequency_THz", "peak_E_node")
    print("Sim 1 (E perp fibre, axial sampling) — fit around 0.6 THz")
    p, e = _plot_one(
        axes[0, 0], f1, e1_axial,
        title="Sim 1 — E ⊥ fibre, axial (z) sampling",
        f0_guess=0.63, window=0.35, colour=OKABE_ITO[2],
    )
    if p is not None:
        print(_format_fit(p, e))

    # Sim 18 — E parallel fibre, axial sampling
    if sim18_csv.exists():
        f18, e18_axial = _load_spectrum(sim18_csv, "frequency_THz", "peak_E_node_axis")
        print("\nSim 18 (E parallel fibre, axial sampling) — fit around 0.6 THz")
        p, e = _plot_one(
            axes[0, 1], f18, e18_axial,
            title="Sim 18 — E ∥ fibre, axial (z) sampling",
            f0_guess=0.6, window=0.35, colour=OKABE_ITO[1],
        )
        if p is not None:
            print(_format_fit(p, e))

        # Sim 18 — annular
        _, e18_ann = _load_spectrum(sim18_csv, "frequency_THz", "peak_E_node_annulus")
        peak_idx = int(np.argmax(e18_ann))
        print(f"\nSim 18 (E parallel fibre, annulus sampling) — empirical peak at "
              f"f = {f18[peak_idx]:.3f}, |E| = {e18_ann[peak_idx]:.3f}")
        p, e = _plot_one(
            axes[1, 0], f18, e18_ann,
            title="Sim 18 — E ∥ fibre, node-annulus sampling",
            f0_guess=float(f18[peak_idx]), window=0.4, colour=OKABE_ITO[6],
        )
        if p is not None:
            print(_format_fit(p, e))

        # Sim 18 — mean annular
        _, e18_mean = _load_spectrum(sim18_csv, "frequency_THz", "mean_E_node_annulus")
        peak_idx2 = int(np.argmax(e18_mean))
        print(f"\nSim 18 mean-annulus — empirical peak at "
              f"f = {f18[peak_idx2]:.3f}, |E| = {e18_mean[peak_idx2]:.3f}")
        p, e = _plot_one(
            axes[1, 1], f18, e18_mean,
            title="Sim 18 — E ∥ fibre, mean |E| in node annulus",
            f0_guess=float(f18[peak_idx2]), window=0.4, colour=OKABE_ITO[3],
        )
        if p is not None:
            print(_format_fit(p, e))
    else:
        for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
            ax.text(0.5, 0.5, "Sim 18 data not yet available",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="0.45")
            ax.set_axis_off()

    fig.suptitle("Sim 19 — Lorentzian fits: ⊥ vs ∥ fibre × axial vs annular sampling")
    fig.tight_layout()
    pdf, png = save_figure(fig, "sim19_lorentzian_fit")
    print(f"\nwrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
