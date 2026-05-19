"""Top-level CSV manifest of simulations (ADR-0004).

One row per simulation, appended after the per-sim HDF5 is written. Uses a
sidecar lock file (`<manifest>.lock`) for cross-process safety on the
Windows box.
"""

from __future__ import annotations

import csv
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

import pandas as pd

MANIFEST_COLUMNS: list[str] = [
    "sweep_name",
    "sim_id",
    "params_json",
    "output_path",
    "comsol_version",
    "git_sha",
    "config_hash",
    "timestamp",
    "status",
]


@contextmanager
def _file_lock(lock_path: Path, timeout: float = 30.0, poll: float = 0.05) -> Iterator[None]:
    """Advisory lock via O_CREAT|O_EXCL. Sufficient for single-machine, low-concurrency use."""

    start = time.monotonic()
    while True:
        try:
            fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            break
        except FileExistsError:
            if time.monotonic() - start > timeout:
                raise TimeoutError(f"Could not acquire lock: {lock_path}") from None
            time.sleep(poll)
    try:
        yield
    finally:
        try:
            os.remove(str(lock_path))
        except FileNotFoundError:
            pass


def append(manifest_path: Path | str, row: dict[str, Any]) -> None:
    """Append a single row to the manifest. Creates the file (with header) if missing."""

    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = manifest_path.with_suffix(manifest_path.suffix + ".lock")

    missing_cols = set(MANIFEST_COLUMNS) - set(row)
    if missing_cols:
        raise ValueError(f"manifest row missing columns: {sorted(missing_cols)}")
    extra_cols = set(row) - set(MANIFEST_COLUMNS)
    if extra_cols:
        raise ValueError(f"manifest row has unknown columns: {sorted(extra_cols)}")

    with _file_lock(lock_path):
        write_header = (not manifest_path.exists()) or manifest_path.stat().st_size == 0
        with manifest_path.open("a", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=MANIFEST_COLUMNS)
            if write_header:
                writer.writeheader()
            writer.writerow({k: row[k] for k in MANIFEST_COLUMNS})


def query(manifest_path: Path | str, **filters: Any) -> pd.DataFrame:
    """Read the manifest and apply equality filters. Returns the filtered DataFrame.

    All columns read as strings (sim_ids like ``"001"`` must not be coerced to int).
    """

    df = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    for col, val in filters.items():
        if col not in df.columns:
            raise KeyError(f"unknown manifest column: {col!r}")
        df = df[df[col] == val]
    return df.reset_index(drop=True)
