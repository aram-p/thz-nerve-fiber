"""Smoke test: confirm `mph` can find and drive the local COMSOL install.

Run via: `uv run python scripts/smoke_test.py`

Captured invocation (Windows, 2026-05-19):
    - mph 1.3.1 + COMSOL 6.3 auto-detected; no env vars needed.
    - `mph.start()` cold-starts COMSOL in ~25-30s.
    - `client.version` returns the string `"6.3"`.

If `mph.start()` cannot locate COMSOL, set `COMSOL_HOME` to the install
root (e.g. `C:\\Program Files\\COMSOL\\COMSOL63\\Multiphysics`).
"""

import mph


def main() -> None:
    client = mph.start()
    print(f"COMSOL version: {client.version}")
    client.clear()


if __name__ == "__main__":
    main()
