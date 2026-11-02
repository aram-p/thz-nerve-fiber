"""Figure: node absorption spectra (Sim 30) vs the quasi-static law."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sim30"
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

par = np.genfromtxt(R / "base_par.csv", delimiter=",", names=True)
perp = np.genfromtxt(R / "base_perp.csv", delimiter=",", names=True)
qs = np.loadtxt(R / "qs.csv", delimiter=",")
f = par["f_THz"]

fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
ax[0].plot(f, par["eta"], "o", ms=4, color="C0", label="full wave, E ∥ fibre")
ax[0].plot(qs[:, 0], qs[:, 1], "-", color="C0", lw=1, label="quasi-static (extrapolated), E ∥")
ax[0].plot(f, perp["eta"], "s", ms=4, color="C1", label="full wave, E ⊥ fibre")
ax[0].plot(qs[:, 0], qs[:, 2], "-", color="C1", lw=1, label="quasi-static (extrapolated), E ⊥")
ax[0].axhline(1, color="0.6", lw=0.7)
ax[0].set_xlabel("Frequency (THz)"); ax[0].set_ylabel("η = ⟨|E|²⟩_node / |E0|²")
ax[0].set_ylim(0, 2); ax[0].legend(fontsize=7); ax[0].set_title("(a) node field factor")

ax[1].plot(f, par["sigma_node"] * 1e12, "o-", ms=3, label="E ∥")
ax[1].plot(f, perp["sigma_node"] * 1e12, "s-", ms=3, label="E ⊥")
ax[1].set_xlabel("Frequency (THz)"); ax[1].set_ylabel("node absorption cross-section (µm²)")
ax[1].legend(fontsize=7); ax[1].set_title("(b) absolute node absorption")

ax[2].semilogy(f, par["share"] * 100, "o-", ms=3, label="E ∥")
ax[2].semilogy(f, perp["share"] * 100, "s-", ms=3, label="E ⊥")
ax[2].set_xlabel("Frequency (THz)"); ax[2].set_ylabel("node share of fibre absorption (%)")
ax[2].legend(fontsize=7); ax[2].set_title("(c) node share (1 mm internode)")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_node_spectrum.{ext}", dpi=200)
print("saved", FIG / "fig_node_spectrum.png")
