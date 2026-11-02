"""Figure + table: Sim 32 scaling, full wave vs extrapolated quasi-static vs closed-form law."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thznerve.physics.materials import WATER_KS1996, axoplasm, myelin
from thznerve.physics.nodelaw import node_factor_closed
from thznerve.physics.quasistatic import node_factor

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sim30"
FIG = ROOT / "paper" / "figures"
runs = {}
for fn in ("sim32_all.json", "sim32_a.json", "sim32_Lint.json"):   # sweeps run as separate jobs
    if (R / fn).exists():
        runs |= json.load(open(R / fn))
F_LIST = np.array([0.1, 0.6, 1.0, 2.0, 3.0]) * 1e12
EPS = {f: [complex(axoplasm()(f)), complex(myelin()(f)), complex(WATER_KS1996(f))] for f in F_LIST}

rows = []
for tag, r in runs.items():
    p = r["params"]
    for pol, field in (("par", "axial"), ("perp", "transverse")):
        for i, f in enumerate(F_LIST):
            e = EPS[f]
            Fq = node_factor(p["a"], p["b"], p["Ln"], e[0], e[1], e[2], e[2], field=field)
            Fc = node_factor_closed(p["a"], p["b"], p["Ln"], *e, pol=pol)
            rows.append([tag, pol, f / 1e12, p["a"], p["b"], p["Ln"], p["Lint"], r[pol]["eta"][i], Fq, Fc,
                         r[pol]["P_node"][i], r[pol]["share"][i]])
cols = ["tag", "pol", "f_THz", "a", "b", "Ln", "Lint", "eta_fem", "F_qs", "F_closed", "P_node", "share"]
with open(R / "sim32_table.csv", "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in rows:
        fh.write(",".join(str(x) for x in r) + "\n")
eta = np.array([r[7] for r in rows]); Fq = np.array([r[8] for r in rows]); Fc = np.array([r[9] for r in rows])
print(f"full wave vs QS: max |ratio-1| = {np.max(abs(eta / Fq - 1)):.3f}; vs interpolation law: {np.max(abs(eta / Fc - 1)):.3f}")
for tag in runs:
    k = [i for i, r in enumerate(rows) if r[0] == tag]
    print(tag, runs[tag]["params"], "eta/F_qs", np.round(eta[k] / Fq[k], 3).min(), np.round(eta[k] / Fq[k], 3).max())

fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
for pol, mk, col in (("par", "o", "C0"), ("perp", "s", "C1")):
    sel = [r for r in rows if r[1] == pol and abs(r[2] - 1.0) < 1e-9 and r[3] == 3.5 and r[4] == 5.0]
    x = np.array([r[5] / (r[4] - r[3]) for r in sel])
    ax[0].plot(x, [r[7] for r in sel], mk, color=col, label=f"full wave, E {'∥' if pol == 'par' else '⊥'}")
    xs = np.logspace(-1, 1, 60)
    e = EPS[1e12]
    ax[0].plot(xs, [node_factor_closed(3.5, 5.0, xx * 1.5, *e, pol=pol) for xx in xs], "-", color=col, lw=1)
ax[0].set_xscale("log"); ax[0].set_xlabel("node length / myelin thickness"); ax[0].set_ylabel("η at 1 THz")
ax[0].set_title("(a) node length, a = 3.5, b = 5 µm\n(line: interpolation law)", fontsize=10); ax[0].legend(fontsize=7)
ax[1].plot(Fq, eta, ".", ms=5)
lim = [0, max(eta.max(), Fq.max()) * 1.05]
ax[1].plot(lim, lim, "k-", lw=0.7); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
ax[1].set_xlabel("F quasi-static (Laplace)"); ax[1].set_ylabel("η full wave"); ax[1].set_title("(b) all cases, 0.1–3 THz")
for pol, mk in (("par", "o"), ("perp", "s")):
    sel = [r for r in rows if r[1] == pol and abs(r[2] - 1.0) < 1e-9 and r[0].startswith(("Lint",))]
    sel = sorted(sel, key=lambda r: r[6])
    ax[2].plot([r[6] for r in sel], [r[10] * 1e15 for r in sel], mk + "-", label=f"E {'∥' if pol == 'par' else '⊥'}")
ax[2].set_xlabel("internode length (µm)"); ax[2].set_ylabel("P_node at 1 THz (fW per node, |E0| = 1 V/m)")
ax[2].set_title("(c) no internode dependence"); ax[2].legend(fontsize=7)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig_scaling.{ext}", dpi=200)
