"""Dispatcher CLI for sweep configs.

Usage:
    uv run python scripts/run.py <sweep.yaml> [--dry-run] [--force]
"""

import argparse
import sys
from pathlib import Path

import yaml

from thznerve.sweep.expand import expand_grid
from thznerve.sweep.schema import SweepConfig

# Placeholder per-sim runtime (minutes). Phase 2 replaces this with a real estimate.
_EST_MIN_PER_SIM = 2


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a YAML-defined sweep.")
    parser.add_argument("config", type=Path, help="Path to sweep YAML")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the parameter grid and an estimated runtime, then exit.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run sims that have already completed (skips idempotency check).",
    )
    args = parser.parse_args()

    with args.config.open() as f:
        raw = yaml.safe_load(f)
    cfg = SweepConfig.model_validate(raw)

    grid = expand_grid(cfg.sweep)
    n = len(grid)

    print(f"Sweep: {cfg.name}")
    if cfg.description:
        print(f"  {cfg.description}")
    print(f"Grid size: {n} simulation{'s' if n != 1 else ''}")

    if args.dry_run:
        for i, row in enumerate(grid, start=1):
            params_str = ", ".join(f"{k}={v:g}" for k, v in row.items())
            print(f"  [{i:>3}] {params_str}")
        est = n * _EST_MIN_PER_SIM
        print(f"Estimated runtime: ~{est} min ({est / 60:.1f} h) at {_EST_MIN_PER_SIM} min/sim")
        return

    # Without --dry-run: Phase 2 wires this to the COMSOL runner. For now, log.
    for i, row in enumerate(grid, start=1):
        print(f"  [{i}/{n}] would run params={row} force={args.force}")


if __name__ == "__main__":
    sys.exit(main())
