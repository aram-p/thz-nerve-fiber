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


def trace_track(mask: np.ndarray, y_start: float, max_jump: float = 40.0):
    """Follow one curve column by column, keeping the run of mask pixels that touches (or is
    nearest) the curve's extent in the last kept column. Survives labels, arrows and stray pixels
    of the same colour; where the curve is hidden under another curve the allowed jump grows by
    max_jump per missing column, so tracking resumes where it reappears (gap left for the caller
    to interpolate)."""
    xs, ys = [], []
    lo = hi = float(y_start)                 # pixel extent of the curve in the last column kept
    last = None
    for x in range(mask.shape[1]):
        yy = np.flatnonzero(mask[:, x])
        if not len(yy):
            continue
        runs = np.split(yy, np.where(np.diff(yy) > 1)[0] + 1)
        # gap between each run and the previous extent (0 if they overlap): steep edges are long
        # vertical runs whose centres jump, but consecutive runs of one line always touch
        gap = np.array([max(0.0, r[0] - hi, lo - r[-1]) for r in runs])
        if gap.min() == 0:                   # several touch: the one sharing most pixels
            overlap = np.array([min(hi, r[-1]) - max(lo, r[0]) + 1 for r in runs])
            k = int(np.argmax(np.where(gap == 0, overlap, -np.inf)))
        else:
            k = int(np.argmin(gap))
        allowed = max_jump * (1 if last is None else x - last)
        if gap[k] <= allowed:
            xs.append(x); ys.append(runs[k].mean()); lo, hi = float(runs[k][0]), float(runs[k][-1])
            last = x
    return np.array(xs, float), np.array(ys, float)


def trace(mask: np.ndarray, min_hits: int = 1):
    """x pixels and mean y pixel of the masked curve in each column that has >= min_hits pixels."""
    xs, ys = [], []
    for x in range(mask.shape[1]):
        yy = np.flatnonzero(mask[:, x])
        if len(yy) >= min_hits:
            xs.append(x); ys.append(yy.mean())
    return np.array(xs, float), np.array(ys, float)


def column_extent(mask: np.ndarray):
    """x pixels with the top and bottom y pixel of the masked curve in each column. A plotted spectrum
    with a narrow dip draws a vertical stroke, so the bottom pixel holds the dip and the top the level."""
    xs = np.flatnonzero(mask.any(axis=0))
    top = np.array([np.flatnonzero(mask[:, x])[0] for x in xs], float)
    bot = np.array([np.flatnonzero(mask[:, x])[-1] for x in xs], float)
    return xs.astype(float), top, bot


def narrow_dips(f, top_db, bottom_db, half_window: float, min_depth: float):
    """Narrow dips of a traced spectrum: depth = running maximum of the upper envelope within
    +-half_window minus the curve's lowest value in that column. Contiguous columns deeper than
    min_depth form one dip; returns (frequency, depth) of the deepest column of each."""
    f = np.asarray(f, float); top_db = np.asarray(top_db, float); bottom_db = np.asarray(bottom_db, float)
    env = np.array([top_db[np.abs(f - fi) <= half_window].max() for fi in f])
    depth = env - bottom_db
    deep = depth >= min_depth
    out = []
    i = 0
    while i < len(f):
        if deep[i]:
            j = i
            while j + 1 < len(f) and deep[j + 1]:
                j += 1
            k = i + int(np.argmax(depth[i:j + 1]))
            out.append((float(f[k]), float(depth[k])))
            i = j + 1
        else:
            i += 1
    return out


def copies_disagree(f1, y1, f2, y2, tol: float, step: float = 0.01):
    """Bands where two traced copies of the same published curve differ by more than tol.
    Plot annotations (arrows, labels) in the curve's own colour corrupt a trace silently; a second
    copy of the same scan, from another figure, exposes them. Returns [(f_lo, f_hi, max |diff|)]."""
    lo, hi = max(np.min(f1), np.min(f2)), min(np.max(f1), np.max(f2))
    g = np.arange(lo, hi + step / 2, step)
    d = np.abs(np.interp(g, f1, y1) - np.interp(g, f2, y2))
    bad = d > tol
    out, i = [], 0
    while i < len(g):
        if bad[i]:
            j = i
            while j + 1 < len(g) and bad[j + 1]:
                j += 1
            out.append((float(g[i]), float(g[j]), float(d[i:j + 1].max())))
            i = j + 1
        else:
            i += 1
    return out
