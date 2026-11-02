"""Read curves off published raster figures, reproducibly.

render_page   PDF page -> RGB array (poppler pdftoppm)
gridlines     pixel centres of light-grey gridlines along one axis
axis_map      least-squares pixel -> data map from gridline pixels and their tick values
trace         one y pixel per x column from a boolean colour mask (mean of the hit rows)
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

_POPPLER_GLOB = "AppData/Local/Microsoft/WinGet/Packages/oschwartz10612.Poppler*/poppler-*/Library/bin/pdftoppm.exe"


def _pdftoppm() -> str:
    exe = os.environ.get("PDFTOPPM") or shutil.which("pdftoppm")
    if exe:
        return exe
    hits = sorted(Path.home().glob(_POPPLER_GLOB))
    if hits:
        return str(hits[-1])
    raise FileNotFoundError("pdftoppm not found: install poppler or set PDFTOPPM")


def render_page(pdf: str | Path, page: int, dpi: int = 300) -> np.ndarray:
    from PIL import Image
    with tempfile.TemporaryDirectory() as td:
        stem = Path(td) / "p"
        subprocess.run([_pdftoppm(), "-r", str(dpi), "-f", str(page), "-l", str(page), "-png",
                        "-singlefile", str(pdf), str(stem)], check=True)
        return np.asarray(Image.open(f"{stem}.png").convert("RGB")).astype(int)


def gridlines(gray: np.ndarray, axis: int, lo: int = 190, hi: int = 240, frac: float = 0.3) -> np.ndarray:
    """Centres of grey lines spanning > frac of the image. axis=0: vertical lines (x pixels), 1: rows."""
    m = (gray > lo) & (gray < hi)
    count = m.sum(axis=0 if axis == 0 else 1)
    span = gray.shape[0 if axis == 0 else 1]
    idx = np.where(count > frac * span)[0]
    groups = np.split(idx, np.where(np.diff(idx) > 2)[0] + 1) if len(idx) else []
    return np.array([g.mean() for g in groups if len(g)])


def axis_map(pixels, values):
    """Linear pixel -> value map; returns (function, max residual in data units)."""
    p, v = np.asarray(pixels, float), np.asarray(values, float)
    k, c = np.polyfit(p, v, 1)
    return (lambda x: k * np.asarray(x, float) + c), float(np.max(np.abs(k * p + c - v)))


def trace(mask: np.ndarray, min_hits: int = 1):
    """x pixels and mean y pixel of the masked curve in each column that has >= min_hits pixels."""
    xs, ys = [], []
    for x in range(mask.shape[1]):
        yy = np.flatnonzero(mask[:, x])
        if len(yy) >= min_hits:
            xs.append(x); ys.append(yy.mean())
    return np.array(xs, float), np.array(ys, float)
