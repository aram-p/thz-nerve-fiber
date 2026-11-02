#!/bin/sh
cd "$(dirname "$0")/.."
for m in aniso axon epsfine; do
  uv run python -u experiments/prod/sim36_review.py $m > results/sim30/sim36_$m.log 2>&1 || echo "FAILED $m" >> results/sim30/sim36_$m.log
done
