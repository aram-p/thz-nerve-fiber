#!/bin/sh
# Sim 32 remaining sweeps (axon radius at Lint = 200, internode length), then chain B.
cd /c/Users/aramp/dev/thz-nerve-fiber
uv run python -u experiments/prod/sim32_scaling.py a > results/sim30/sim32_a.log 2>&1
uv run python -u experiments/prod/sim32_scaling.py Lint > results/sim30/sim32_Lint.log 2>&1
sh scripts/chain_b.sh
