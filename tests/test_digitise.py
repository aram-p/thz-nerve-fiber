import numpy as np

from thznerve.digitise import axis_map, gridlines, trace


def _synthetic():
    h, w = 200, 300
    g = np.full((h, w), 255)
    for x in (20, 120, 220):          # vertical gridlines at data x = 0, 1, 2
        g[:, x] = 220
    for y in (30, 130):               # rows at data y = 1, 0
        g[y, :] = 220
    xs = np.arange(20, 221)
    ys = np.round(130 - 100 * np.sin(np.pi * (xs - 20) / 200)).astype(int)
    mask = np.zeros((h, w), bool)
    mask[ys, xs] = True
    mask[ys + 1, xs] = True           # 2-pixel-thick line: mean is half a pixel low
    return g, mask


def test_gridlines_and_axis_map():
    g, _ = _synthetic()
    cols, rows = gridlines(g, 0), gridlines(g, 1)
    assert np.allclose(cols, [20, 120, 220]) and np.allclose(rows, [30, 130])
    fx, err = axis_map(cols, [0, 1, 2])
    assert err < 1e-9 and abs(fx(170) - 1.5) < 1e-9


def test_trace_recovers_curve():
    g, mask = _synthetic()
    fx, _ = axis_map(gridlines(g, 0), [0, 1, 2])
    fy, _ = axis_map(gridlines(g, 1), [1, 0])
    x, y = trace(mask)
    xd, yd = fx(x), fy(y - 0.5)
    assert np.max(np.abs(yd - np.sin(np.pi * xd / 2))) < 0.01


def test_trace_track_ignores_crossing_arrow():
    from thznerve.digitise import trace_track
    g, mask = _synthetic()
    for k in range(60):                      # a diagonal "arrow" of the same colour crossing the curve
        mask[20 + k, 100 + k] = True
    x, y = trace_track(mask, y_start=130.0)
    ref_x, ref_y = trace(_synthetic()[1])
    common = np.isin(x, ref_x)
    assert np.max(np.abs(y[common] - ref_y[np.isin(ref_x, x)])) < 1.0


def test_trace_track_prefers_overlapping_run_and_resumes_after_occlusion():
    from thznerve.digitise import trace_track
    _, mask = _synthetic()
    ref_x, ref_y = trace(mask)
    mask[ref_y[50].astype(int) - 1, 70] = False
    mask[int(ref_y[50]) + 8, 70] = True      # stray pixel touching only the edge of the extent
    mask[:, 150:175] = False                 # curve hidden under another one for 25 columns
    x, y = trace_track(mask, y_start=130.0, max_jump=5.0)
    assert x.max() >= 220 and not np.any((x >= 150) & (x < 175))
    common = np.isin(x, ref_x)
    assert np.max(np.abs(y[common] - ref_y[np.isin(ref_x, x)])) < 1.0
