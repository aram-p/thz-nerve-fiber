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


def test_narrow_dips_finds_each_dip_once_with_its_depth():
    from thznerve.digitise import narrow_dips
    f = np.linspace(0, 3, 301)
    level = -10 - 10 * f                      # sloping baseline
    bottom = level.copy()
    bottom[100] -= 20                          # 20 dB dip at 1.0 THz, one column wide
    bottom[200:203] -= np.array([8, 15, 8])    # 15 dB dip at 2.01 THz, three columns wide
    bottom[250] -= 5                           # 5 dB dip: below threshold
    dips = narrow_dips(f, level, bottom, half_window=0.05, min_depth=10)
    assert len(dips) == 2
    (f1, d1), (f2, d2) = dips
    assert abs(f1 - 1.0) < 1e-9 and abs(d1 - (20 + 10 * 0.05)) < 1e-6   # envelope is the max within the window
    assert abs(f2 - 2.01) < 1e-9 and 15 < d2 < 16


def test_column_extent_top_and_bottom():
    from thznerve.digitise import column_extent
    m = np.zeros((50, 10), bool)
    m[10:12, 2] = True
    m[10:40, 3] = True                          # vertical stroke: a narrow dip
    xs, top, bot = column_extent(m)
    assert list(xs) == [2, 3] and list(top) == [10, 10] and list(bot) == [11, 39]


def test_copies_disagree_flags_only_the_corrupted_band():
    from thznerve.digitise import copies_disagree
    f = np.linspace(0, 3, 301)
    y = -10 - 10 * f
    y2 = y + 0.3                                   # tracing noise
    y2[(f > 0.8) & (f < 1.2)] -= 11                # an arrow drawn over the curve
    bands = copies_disagree(f, y, f, y2, tol=3.0)
    assert len(bands) == 1
    lo, hi, mx = bands[0]
    assert 0.8 <= lo <= 0.82 and 1.18 <= hi <= 1.2 and abs(mx - 10.7) < 1e-6
