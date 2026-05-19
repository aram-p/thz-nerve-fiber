"""Expand a sweep spec into a flat list of parameter dicts.

Top-level keys in `sweep` produce a Cartesian product. A key literally named
`zip` (whose value is a `ZipGroup`) ties its inner parameters together — they
vary in lockstep and contribute one factor to the Cartesian product.
"""

from itertools import product

import numpy as np

from thznerve.sweep.schema import SweepExpr, SweepValue, ZipGroup


def _expand_expr(expr: SweepExpr) -> list[float]:
    if expr.linspace is not None:
        ls = expr.linspace
        return np.linspace(ls.start, ls.stop, ls.n).tolist()
    if expr.logspace is not None:
        lg = expr.logspace
        return np.logspace(lg.start, lg.stop, lg.n, base=lg.base).tolist()
    if expr.list_ is not None:
        return list(expr.list_)
    raise ValueError("empty SweepExpr — schema validation should have caught this")


def expand_grid(sweep: dict[str, SweepValue]) -> list[dict[str, float]]:
    """Return the full parameter grid as a list of dicts."""

    groups: list[list[dict[str, float]]] = []

    for key, val in sweep.items():
        if isinstance(val, ZipGroup):
            param_lists = {name: _expand_expr(e) for name, e in val.zip.items()}
            lengths = {len(v) for v in param_lists.values()}
            if len(lengths) != 1:
                raise ValueError(
                    f"zip block has mismatched lengths: "
                    f"{ {name: len(v) for name, v in param_lists.items()} }"
                )
            n = lengths.pop()
            zipped = [
                {name: param_lists[name][i] for name in param_lists} for i in range(n)
            ]
            groups.append(zipped)
        else:
            values = _expand_expr(val)
            groups.append([{key: v} for v in values])

    rows: list[dict[str, float]] = []
    for combo in product(*groups):
        row: dict[str, float] = {}
        for d in combo:
            row.update(d)
        rows.append(row)
    return rows
