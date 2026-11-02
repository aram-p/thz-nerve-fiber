"""Figure for Sim 27: legacy 'node field' spectra depend on box length, not on the fibre."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sim27"
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def load(name):
    p = R / f"summary_{name}.csv"
    a = np.genfromtxt(p, delimiter=",", names=True, dtype=None, encoding="utf-8")
    return a


A = load("A")
D = load("D")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(A["f_THz"], A["ax_peak_node"], "o-", label="legacy, box 240 µm (Sims 1-26)")
ax[0].plot(D["f_THz"], D["ax_peak_node"], "s-", label="legacy, box 300 µm")
ax[0].plot(A["f_THz"], A["true_planewave_node"], "k--", label="plane wave in water (no fibre)")
ax[0].set_xlabel("Frequency (THz)"); ax[0].set_ylabel("peak |E| on axis at node (V/m)")
ax[0].legend(fontsize=7); ax[0].set_title("(a) legacy 'resonance' metric")
z = np.load(R / "A_fibre_L240_0.60_axial.npy")
ax[1].plot(z[0], z[1], label="legacy FEM |E|(z), f = 0.6 THz")
ax[1].plot(z[0], z[3], "k--", label="plane wave exp(-j k_w z)")
ax[1].axvspan(100, 140, color="0.85", label="node")
ax[1].set_xlabel("z (µm)"); ax[1].set_ylabel("|E| on axis (V/m)")
ax[1].legend(fontsize=7); ax[1].set_title("(b) vacuum-k background in water")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_legacy_artefact.{ext}", dpi=200)
print("saved")
