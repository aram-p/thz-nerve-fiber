"""Dispatcher CLI for sweep configs.

Usage:
    uv run python scripts/run.py <sweep.yaml> [--dry-run] [--force]

Loads a YAML sweep config, expands the parameter grid, and drives the
per-simulation runner. Body is filled in by Phase 1.4 (issue #5).
"""


def main() -> None:
    raise NotImplementedError("Phase 1.4 — see issue #5")


if __name__ == "__main__":
    main()
