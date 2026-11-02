#!/bin/sh
# Run after chain A finishes: SBC-vs-R check, then Sim 33 (node states, robustness).
cd /c/Users/aramp/dev/thz-nerve-fiber
uv run python -u experiments/prod/sim30_node_spectrum.py sbcR tag=bc_sbc_R300 outer=sbc R=300 fmin=0.1 fmax=0.3 nf=2 > results/sim30/sbcR.log 2>&1
uv run python -u experiments/prod/sim33_robust_nodestate.py state > results/sim30/sim33_state.log 2>&1
uv run python -u experiments/prod/sim33_robust_nodestate.py robust > results/sim30/sim33_robust.log 2>&1
echo CHAIN_B_DONE
