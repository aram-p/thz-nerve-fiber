#!/bin/sh
# Sim 36: one COMSOL process at a time.
cd "$(dirname "$0")/.."
for m in sigma aniso dense epsnode; do
  uv run python -u experiments/prod/sim36_review.py $m > results/sim30/sim36_$m.log 2>&1 || echo "FAILED $m" >> results/sim30/sim36_$m.log
done
