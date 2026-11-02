"""Per-run bookkeeping: wall time, peak memory, DOF -> results/run_log.csv.

Usage::

    with RunLog("sim27", "legacy box-length artefact test") as log:
        ...solve...
        log.dof = n_dof
        log.outcome = "peaks move with L: artefact"
        log.wasted_min, log.why = 0, ""

Peak memory is the peak private commit of this process (Windows
``peak_pagefile``). In mph stand-alone mode COMSOL lives in this process,
so that is the solver's peak.
"""

from __future__ import annotations

import csv
import datetime as _dt
import time
from dataclasses import dataclass, field
from pathlib import Path

import psutil

REPO_ROOT = Path(__file__).resolve().parents[2]
LOG_PATH = REPO_ROOT / "results" / "run_log.csv"
COLUMNS = [
    "date", "sim", "description", "dof", "peak_mem_gb", "wall_s",
    "outcome", "time_wasted_min", "why_wasted",
]


def peak_private_gb() -> float:
    info = psutil.Process().memory_info()
    peak = getattr(info, "peak_pagefile", None) or getattr(info, "rss", 0)
    return peak / 1e9


def append_row(row: dict, path: Path = LOG_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists()
    with path.open("a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in COLUMNS})


@dataclass
class RunLog:
    sim: str
    description: str
    dof: int | str = ""
    outcome: str = ""
    wasted_min: float = 0.0
    why: str = ""
    path: Path = LOG_PATH
    _t0: float = field(default=0.0, repr=False)

    def __enter__(self) -> "RunLog":
        self._t0 = time.monotonic()
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc is not None and not self.outcome:
            msg = " ".join(str(exc).split())          # one line: COMSOL messages are multi-line
            self.outcome = f"FAILED: {exc_type.__name__}: {msg[-200:]}"
        append_row(
            {
                "date": _dt.datetime.now().isoformat(timespec="seconds"),
                "sim": self.sim,
                "description": self.description,
                "dof": self.dof,
                "peak_mem_gb": f"{peak_private_gb():.2f}",
                "wall_s": f"{time.monotonic() - self._t0:.1f}",
                "outcome": self.outcome,
                "time_wasted_min": self.wasted_min,
                "why_wasted": self.why,
            },
            self.path,
        )
        return False
