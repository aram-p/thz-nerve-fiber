"""Figure 1: model geometry (axisymmetric unit cell, not to scale) and THz material data."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from thznerve.physics.materials import (WATER_JEPSEN2007, WATER_KS1996, alpha_per_cm, axoplasm, myelin,
                                        refractive_index)

FIG = Path(__file__).resolve().parents[2] / "paper" / "figures"
f = np.linspace(0.1, 3.0, 300) * 1e12

fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), gridspec_kw=dict(width_ratios=[1.1, 1, 1]))

# (a) r-z half plane, one period, z compressed away from the node
a, b, Ln, R, pml = 3.5, 5.0, 1.0, 14.0, 4.0
H = 12.0                                                   # drawn half-period (not to scale)
A = ax[0]
A.add_patch(Rectangle((0, -H), a, 2 * H, fc="#cfe3f5", ec="k", lw=0.6))
A.add_patch(Rectangle((a, -H), b - a, H - Ln / 2, fc="#f4d9a7", ec="k", lw=0.6))
A.add_patch(Rectangle((a, Ln / 2), b - a, H - Ln / 2, fc="#f4d9a7", ec="k", lw=0.6))
A.add_patch(Rectangle((a, -Ln / 2), b - a, Ln, fc="#3b7dc4", ec="k", lw=0.6))
A.add_patch(Rectangle((b, -H), R - b, 2 * H, fc="#e8f1fa", ec="k", lw=0.6))
A.add_patch(Rectangle((R, -H), pml, 2 * H, fc="0.85", ec="k", lw=0.6, hatch="//"))
A.plot([R - 2, R - 2], [-H, H], "r--", lw=0.8)
A.text(a / 2, 6, "axon\n(axoplasm)", ha="center", fontsize=7)
A.text((a + b) / 2, 6, "myelin", ha="center", rotation=90, fontsize=7)
A.annotate("node gap\n(water, Ln)", xy=(b - 0.2, 0), xytext=(7.5, 3.0), fontsize=7,
           arrowprops=dict(arrowstyle="->", lw=0.6))
A.text(9.5, -6, "water", fontsize=7)
A.text(R + pml / 2, 0, "PML", rotation=90, ha="center", va="center", fontsize=7)
A.text(R - 2.3, -H + 1, "flux contour", rotation=90, color="r", fontsize=6, ha="right")
A.text(0.3, H + 0.6, "periodic in z, period Lint + Ln", fontsize=7)
A.annotate("", xy=(0, -H - 1.2), xytext=(a, -H - 1.2), arrowprops=dict(arrowstyle="<->", lw=0.6))
A.text(a / 2, -H - 2.6, "a", ha="center", fontsize=7)
A.annotate("", xy=(0, -H - 3.4), xytext=(b, -H - 3.4), arrowprops=dict(arrowstyle="<->", lw=0.6))
A.text(b / 2, -H - 4.8, "b", ha="center", fontsize=7)
A.set_xlim(-0.5, R + pml + 0.5); A.set_ylim(-H - 5.5, H + 2)
A.set_xlabel("r (schematic)"); A.set_ylabel("z (schematic)")
A.set_xticks([]); A.set_yticks([])
A.set_title("(a) axisymmetric unit cell, orders m = 0..3", fontsize=9)

# (b) permittivity
mats = [("water (Kindt & Schmuttenmaer)", WATER_KS1996(f), "C0", "-"),
        ("water (Jepsen et al.)", WATER_JEPSEN2007(f), "C0", ":"),
        ("axoplasm (85 % water)", axoplasm()(f), "C2", "-"),
        ("myelin (40 % water)", myelin()(f), "C3", "-")]
for lab, e, c, ls in mats:
    ax[1].plot(f / 1e12, e.real, ls, color=c, label=lab)
    ax[1].plot(f / 1e12, -e.imag, ls, color=c, lw=0.8, alpha=0.6)
ax[1].set_xlabel("Frequency (THz)"); ax[1].set_ylabel("ε′ (thick), ε″ (thin)")
ax[1].set_title("(b) permittivity, 20–25 °C", fontsize=9); ax[1].legend(fontsize=6)

# (c) index, attenuation length and electrical size of the node
nw = refractive_index(WATER_KS1996(f))
nm = refractive_index(myelin()(f))
ax[2].plot(f / 1e12, nw.real, "C0", label="n water")
ax[2].plot(f / 1e12, nm.real, "C3", label="n myelin")
ax[2].plot(f / 1e12, refractive_index(axoplasm()(f)).real, "C2", label="n axoplasm (core V ≤ 0.4: leaky)")
ax[2].set_xlabel("Frequency (THz)"); ax[2].set_ylabel("refractive index")
t2 = ax[2].twinx()
t2.plot(f / 1e12, 1e4 / alpha_per_cm(WATER_KS1996(f), f), "k--", lw=0.8, label="water power 1/α (µm)")
t2.plot(f / 1e12, 2 * np.pi * nw.real * f / 2.99792458e8 * 1e-6 * 1000, "k:", lw=0.8,
        label="k_w·Ln ×1000 (Ln = 1 µm)")
t2.set_ylabel("µm  /  dimensionless ×1000")
h1, l1 = ax[2].get_legend_handles_labels(); h2, l2 = t2.get_legend_handles_labels()
ax[2].legend(h1 + h2, l1 + l2, fontsize=6, loc="upper center")
ax[2].set_title("(c) why the node is quasi-static", fontsize=9)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_model.{ext}", dpi=200)
print("saved", FIG / "fig_model.png")
print("k_w Ln max", (2 * np.pi * nw.real * f / 2.99792458e8 * 1e-6).max(),
      "1/alpha um range", (1e4 / alpha_per_cm(WATER_KS1996(f), f)).min(), (1e4 / alpha_per_cm(WATER_KS1996(f), f)).max())
