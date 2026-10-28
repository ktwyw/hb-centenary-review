"""Export the data drawn in a Matplotlib figure to a tidy CSV file.

Called by hbstyle.save() for every figure, so that the numbers behind each curve, marker,
bar, band and reference line of the article and of its supplementary material are available
as plain text in data/figure_data/<figure name>.csv.

One row per plotted point. Columns:
  panel     panel letter as printed, e.g. "(a)"; "(b) right axis" for a twin axis
  x_label, y_label   axis labels as typeset (Matplotlib mathtext), with the axis scale
  series    legend label of the artist, or a description such as "line 3 (#0072B2, --)"
  kind      line | markers | errorbar | bar | band | hline | vline | segment
  index     point number within the series (or segment/bar number)
  x, y      data coordinates; for bars and bands x/y are one corner and x_end/y_end the other
  x_end, y_end, category   used by bars, segments and reference lines where relevant
Schematic elements (shapes without data coordinates) are not exported.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from matplotlib.collections import LineCollection, PathCollection, PolyCollection
from matplotlib.colors import to_hex
from matplotlib.container import BarContainer, ErrorbarContainer
from matplotlib.lines import Line2D

FIELDS = ["panel", "x_label", "y_label", "series", "kind", "index", "x", "y",
          "x_end", "y_end", "category"]


def _panel_names(fig):
    """Map every axes to the panel letter printed next to it."""
    names = {}
    labelled = []
    for ax in fig.axes:
        for t in ax.texts:
            s = t.get_text().strip()
            if len(s) == 3 and s.startswith("(") and s.endswith(")") and s[1].isalpha():
                names[ax] = s
                labelled.append(ax)
                break
    for ax in fig.axes:
        if ax in names:
            continue
        # a twin axis shares the position of a labelled axis
        for other in labelled:
            if np.allclose(ax.get_position().bounds, other.get_position().bounds):
                names[ax] = names[other] + " right axis"
                break
    count = 0
    for ax in fig.axes:
        if ax not in names:
            names[ax] = "(single panel)" if len(fig.axes) == 1 else f"axes {count + 1}"
        count += 1
    return names


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    v = float(v)
    if not np.isfinite(v):
        return ""
    return f"{v:.10g}"


def _label(artist, kind, i):
    gid = artist.get_gid() if hasattr(artist, "get_gid") else None
    if gid:
        return gid
    lab = artist.get_label() if hasattr(artist, "get_label") else ""
    if lab and not lab.startswith("_"):
        return lab
    desc = kind
    try:
        if isinstance(artist, Line2D):
            desc = f"{kind} {i} ({to_hex(artist.get_color())}, {artist.get_linestyle()})"
        else:
            fc = artist.get_facecolor()
            if len(fc):
                desc = f"{kind} {i} ({to_hex(fc[0])})"
    except Exception:
        pass
    return desc


def _axis_label(ax, which):
    axis = ax.xaxis if which == "x" else ax.yaxis
    text = axis.get_label().get_text()
    scale = axis.get_scale()
    return f"{text} [{scale}]" if scale != "linear" else text


def _write_table(path, table, title):
    fields = list(table[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# {title}\n")
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in table:
            w.writerow({k: _fmt(r.get(k)) if not isinstance(r.get(k), str) else r.get(k) for k in fields})


def export(fig, name, outdir):
    """Write data/figure_data/<name>.csv.

    A script can attach fig.hb_table (list of dicts) when a table describes the figure better than
    its drawn coordinates (e.g. a matrix of findings), and fig.hb_extra = {suffix: list of dicts}
    for quantities that are drawn as images or contours (e.g. a residual map on a grid)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for suffix, table in (getattr(fig, "hb_extra", None) or {}).items():
        _write_table(outdir / f"{name}_{suffix}.csv", table,
                     f"Data drawn in {name}.pdf: {suffix.replace('_', ' ')} (see data/figure_data/README.md).")
    if getattr(fig, "hb_table", None):
        path = outdir / f"{name}.csv"
        _write_table(path, fig.hb_table, f"Data shown in {name}.pdf (see data/figure_data/README.md).")
        return path, len(fig.hb_table)
    panels = _panel_names(fig)
    rows = []
    for ax in fig.axes:
        if not ax.get_visible() or not ax.axison:
            # axes switched off (schematics, legend-only axes): nothing quantitative
            continue
        base = {"panel": panels[ax], "x_label": _axis_label(ax, "x"), "y_label": _axis_label(ax, "y")}
        data_tr = ax.transData
        err_lines, err_cols = set(), set()
        for c in ax.containers:
            if isinstance(c, ErrorbarContainer):
                line, caps, cols = c.lines
                if line is not None:
                    err_lines.add(line)
                for col in cols:
                    err_cols.add(col)
                    for k, seg in enumerate(col.get_segments()):
                        if len(seg) == 2:
                            rows.append({**base, "series": _label(c, "errorbar", 0), "kind": "errorbar",
                                         "index": k, "x": seg[0][0], "y": seg[0][1],
                                         "x_end": seg[1][0], "y_end": seg[1][1]})
                for cap in caps:
                    err_lines.add(cap)
            elif isinstance(c, BarContainer):
                horizontal = getattr(c, "orientation", None) == "horizontal"
                ticks = {round(float(t), 6): lab.get_text() for t, lab in
                         zip(ax.get_yticks() if horizontal else ax.get_xticks(),
                             ax.get_yticklabels() if horizontal else ax.get_xticklabels())}
                for k, rect in enumerate(c.patches):
                    x0, y0 = rect.get_x(), rect.get_y()
                    x1, y1 = x0 + rect.get_width(), y0 + rect.get_height()
                    centre = (y0 + y1) / 2 if horizontal else (x0 + x1) / 2
                    rows.append({**base, "series": rect.get_gid() or _label(c, "bar", 0), "kind": "bar", "index": k,
                                 "x": x0, "y": y0, "x_end": x1, "y_end": y1,
                                 "category": ticks.get(round(float(centre), 6), "")})
        for i, line in enumerate(ax.get_lines()):
            if line in err_lines or not line.get_visible():
                continue
            x = np.asarray(line.get_xdata(orig=False), dtype=float)
            y = np.asarray(line.get_ydata(orig=False), dtype=float)
            tr = line.get_transform()
            has_line = line.get_linestyle() not in ("None", "", " ")
            kind = "line" if has_line else "markers"
            if tr == ax.get_yaxis_transform():      # axhline: x in axes fraction
                rows.append({**base, "series": _label(line, "hline", i), "kind": "hline", "index": 0,
                             "y": y[0]})
                continue
            if tr == ax.get_xaxis_transform():      # axvline
                rows.append({**base, "series": _label(line, "vline", i), "kind": "vline", "index": 0,
                             "x": x[0]})
                continue
            if not tr.contains_branch(data_tr):
                continue                            # drawn in axes/figure coordinates: decoration
            lab = _label(line, kind, i)
            for k, (a, b) in enumerate(zip(x, y)):
                rows.append({**base, "series": lab, "kind": kind, "index": k, "x": a, "y": b})
        for i, col in enumerate(ax.collections):
            if col in err_cols or not col.get_visible():
                continue
            if isinstance(col, PathCollection):        # scatter
                if not col.get_offset_transform().contains_branch(data_tr):
                    continue
                lab = _label(col, "markers", i)
                for k, (a, b) in enumerate(np.asarray(col.get_offsets(), dtype=float)):
                    rows.append({**base, "series": lab, "kind": "markers", "index": k, "x": a, "y": b})
            elif isinstance(col, PolyCollection):      # fill_between / axvspan bands
                tr = col.get_transform()
                for p, path in enumerate(col.get_paths()):
                    v = path.vertices
                    if tr.contains_branch(data_tr):
                        lab = _label(col, "band", i)
                        for k, (a, b) in enumerate(v):
                            rows.append({**base, "series": lab, "kind": "band", "index": k, "x": a, "y": b})
                    elif tr == ax.get_xaxis_transform():   # axvspan: x range only
                        rows.append({**base, "series": _label(col, "band", i), "kind": "band", "index": p,
                                     "x": v[:, 0].min(), "x_end": v[:, 0].max()})
            elif isinstance(col, LineCollection):
                if not col.get_transform().contains_branch(data_tr):
                    continue
                lab = _label(col, "segment", i)
                for k, seg in enumerate(col.get_segments()):
                    if len(seg) >= 2:
                        rows.append({**base, "series": lab, "kind": "segment", "index": k,
                                     "x": seg[0][0], "y": seg[0][1], "x_end": seg[-1][0], "y_end": seg[-1][1]})
        for i, patch in enumerate(ax.patches):     # axvspan in recent Matplotlib is a Rectangle
            if any(patch in getattr(c, "patches", []) for c in ax.containers):
                continue
            tr = patch.get_transform()
            if hasattr(patch, "get_x") and tr.contains_branch(ax.get_xaxis_transform()) and \
                    not tr.contains_branch(data_tr):
                x0 = patch.get_x()
                rows.append({**base, "series": _label(patch, "band", i), "kind": "band", "index": 0,
                             "x": x0, "x_end": x0 + patch.get_width()})
    path = outdir / f"{name}.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# Data drawn in {name}.pdf (exported by scripts/figdata.py when the figure is made).\n")
        fh.write("# Columns are described in data/figure_data/README.md.\n")
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in FIELDS})
    return path, len(rows)
