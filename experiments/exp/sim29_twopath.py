"""Sim 29 - Does a partial air path reproduce the reported voltage-induced features?

Inputs are only the authors' own numbers (Hovhannisyan, Makaryan, Martirosyan,
IRPhE 2024, Figs 1-2): d = 100 um, n = 3.9 (low f) -> 3.6 (2 THz), power loss of
the sample relative to the empty cuvette read off Fig. 1b. One free parameter:
the field fraction p of the beam that crosses air. Prediction without fitting:
notches at f_k = (2k-1) c / (2 (n-1) d).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.physics.materials import C0
from thznerve.physics.twopath import notch_frequencies, retrieve_n, two_path

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

f = np.linspace(0.15, 3.0, 2000) * 1e12
d = 100e-6
# authors' sample index and loss (approximate reading of Figs 1b and 2, no voltage)
n_t = np.interp(f / 1e12, [0.15, 1.0, 2.0, 3.0], [3.9, 3.85, 3.6, 3.5])
loss_db = np.interp(f / 1e12, [0.15, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0], [4, 5, 13, 16, 25, 33, 33])
T_t = 10 ** (-loss_db / 20) * np.exp(-1j * (n_t - 1) * 2 * np.pi * f / C0 * d)

fn = notch_frequencies(3.9, d), notch_frequencies(3.6, d)
fn70 = notch_frequencies(3.9, 70e-6), notch_frequencies(3.6, 70e-6)
print("notches d=100um n=3.9:", np.round(fn[0] / 1e12, 3), " n=3.6:", np.round(fn[1] / 1e12, 3))
print("notches d=70um  n=3.9:", np.round(fn70[0] / 1e12, 3), " n=3.6:", np.round(fn70[1] / 1e12, 3))

# p chosen so the first notch is deep: p/(1-p) = |T_t(f1)|
i1 = np.argmin(np.abs(f - fn[0][0]))
p = abs(T_t[i1]) / (1 + abs(T_t[i1]))
T_on = two_path(T_t, p)
print(f"p (field fraction through air) = {p:.3f}  -> area fraction ~ {p:.2f} of detected field")
n_off = retrieve_n(f, T_t, d)
n_on = retrieve_n(f, T_on, d)

# time domain with a generic TDS reference pulse
t = np.linspace(-5e-12, 15e-12, 4000)
spec = (f / 1e12) ** 1.5 * np.exp(-((f / 1e12) / 0.9) ** 2)
def td(T):
    return np.array([np.sum(spec * T * np.exp(2j * np.pi * f * tt)).real for tt in t]) / len(f)

fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
ax[0].plot(f / 1e12, 20 * np.log10(np.abs(T_t)), label="tissue (authors' loss, no voltage)")
ax[0].plot(f / 1e12, 20 * np.log10(np.abs(T_on)), "--", label=f"tissue + {p:.2f} air path")
for x in fn[0][:2]:
    ax[0].axvline(x / 1e12, color="0.6", lw=0.7)
ax[0].plot([0.6, 1.5, 2.0, 3.0], [-11.5 - 13, -30 + 8, -33 + 11, -31 + 18], "kx", label="authors, voltage on (Fig. 1b)")
ax[0].set_xlabel("Frequency (THz)"); ax[0].set_ylabel("|T|² rel. empty cuvette (dB)")
ax[0].legend(fontsize=7); ax[0].set_title("(a) transmission")
ax[1].plot(f / 1e12, n_off, label="retrieved n, tissue")
ax[1].plot(f / 1e12, n_on, "--", label="retrieved n, two-path")
alt = np.where(f > fn[0][0], n_on + C0 / (f * d), np.nan)
ax[1].plot(f / 1e12, alt, ":", color="C1", label="two-path, other 2π branch")
ax[1].plot([0.6, 1.2, 2.1], [7.8, 4.2, 2.9], "kx", label="authors, voltage on (Fig. 2)")
ax[1].plot([0.58], [0.7], "kx")
ax[1].axhline(1, color="0.6", lw=0.7)
ax[1].set_ylim(-1, 9); ax[1].set_xlim(0.15, 2.2)
ax[1].set_xlabel("Frequency (THz)"); ax[1].set_ylabel("n = 1 + c·φ/(ωd)")
ax[1].legend(fontsize=7); ax[1].set_title("(b) phase retrieval artefact")
ax[2].plot(t * 1e12, td(np.ones_like(f)), "k", lw=0.8, label="empty cuvette")
ax[2].plot(t * 1e12, td(T_t), label="tissue")
ax[2].plot(t * 1e12, td(T_on), "--", label="two-path")
ax[2].set_xlim(-2, 6); ax[2].set_xlabel("t (ps)"); ax[2].legend(fontsize=7)
ax[2].set_title("(c) early air-path pulse")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_twopath.{ext}", dpi=200)
k = (f > 0.4e12) & (f < 0.8e12)
print("two-path n range near notch:", n_on[k].min().round(2), n_on[k].max().round(2))
for fx in (1.5, 2.0, 2.5, 3.0):
    j = np.argmin(np.abs(f / 1e12 - fx))
    print(f"gain at {fx} THz: {20*np.log10(abs(T_on[j])/abs(T_t[j])):.1f} dB")
