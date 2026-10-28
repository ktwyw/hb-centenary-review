"""Shared plotting style and Herschel-Bulkley helper functions.

Style modelled on Journal of Fluid Mechanics figures (cf. Boyko & Stone,
JFM 923, R5, 2021): Times/STIX text, full box axes with inward ticks on all
four sides, minor ticks, italic panel labels "(a)", "(b)" outside the axes,
compact frameless legends, the colour-blind-safe Okabe-Ito palette with dashed/dotted lines
for asymptotes and markers for numerical results.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "axes.linewidth": 0.6,
    "axes.spines.top": True, "axes.spines.right": True,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
    "xtick.major.pad": 3, "ytick.major.pad": 3,
    "lines.linewidth": 1.1,
    "lines.markersize": 3.5,
    "legend.frameon": False,
    "legend.handlelength": 2.0,
    "legend.handletextpad": 0.5,
    "legend.labelspacing": 0.25,
    "legend.borderaxespad": 0.3,
    "axes.formatter.use_mathtext": True,
    "figure.dpi": 150,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

# Okabe-Ito colour-blind-safe palette (Okabe & Ito 2008), as in the colour cycle of
# github.com/ktwyw/publication-figures; keys keep their semantic roles so that every
# figure of the review uses the same colour for the same kind of curve.
C = {
    "black": "#000000", "red": "#D55E00", "blue": "#0072B2", "cyan": "#56B4E9",
    "purple": "#CC79A7", "green": "#009E73", "orange": "#E69F00", "yellow": "#F0E442",
    "grey": "#808080", "lightgrey": "#c8c8c8",
}
PALETTE = [C["black"], C["red"], C["blue"], C["green"], C["purple"], C["cyan"]]
LS = ["-", "--", ":", "-.", (0, (5, 1, 1, 1)), (0, (3, 1))]

SINGLE = (3.3, 2.5)      # single-panel figure
DOUBLE = (6.6, 3.05)     # two panels side by side; included at \\textwidth (164 mm in REVTeX preprint, 170 mm in print)


def hb_stress(gdot, tau_y, K, n):
    gdot = np.asarray(gdot, dtype=float)
    return tau_y + K * gdot**n


def hb_viscosity(gdot, tau_y, K, n):
    gdot = np.asarray(gdot, dtype=float)
    return tau_y / gdot + K * gdot**(n - 1.0)


def papanastasiou_viscosity(gdot, tau_y, K, n, m):
    gdot = np.asarray(gdot, dtype=float)
    return tau_y * (-np.expm1(-m * gdot)) / gdot + K * gdot**(n - 1.0)


def biviscosity(gdot, tau_y, K, n, eta_r):
    gdot = np.asarray(gdot, dtype=float)
    lo, hi = 1e-12, 1e6
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if eta_r * mid - (tau_y + K * mid**n) > 0:
            hi = mid
        else:
            lo = mid
    gc = np.sqrt(lo * hi)
    return np.where(gdot < gc, eta_r, hb_viscosity(gdot, tau_y, K, n)), gc


def label_panels(axes, labels=None, x=-0.22, y=1.0):
    """JFM-style italic panel labels outside the top-left corner."""
    if labels is None:
        labels = [f"({chr(97+i)})" for i in range(len(axes))]
    for ax, lab in zip(axes, labels):
        ax.text(x, y, lab, transform=ax.transAxes, fontsize=9.5, style="italic",
                va="bottom", ha="left")


def slope_triangle(ax, x0, y0, dx_dec, slope, label, color="k", side="below", fontsize=7.5, run_label="1"):
    """Log-log slope triangle in JFM style. Hypotenuse from (x0, y0) to (x1, y1) with
    x1 = x0*10**dx_dec, y1 = y0*10**(slope*dx_dec); the right-angle corner is placed
    below or above the hypotenuse. The run is labelled `run_label` (default '1') and the rise `label`;
    for a negative slope the convention 'up is positive, left is negative' uses run_label='$-1$'
    and a positive rise label."""
    x1 = x0 * 10**dx_dec; y1 = y0 * 10**(slope * dx_dec)
    if side == "below":
        cx, cy = (x1, y0) if slope > 0 else (x0, y1)
    else:
        cx, cy = (x0, y1) if slope > 0 else (x1, y0)
    ax.plot([x0, cx, x1, x0], [y0, cy, y1, y0], color=color, lw=0.7, solid_joinstyle="miter",
            gid=f"slope triangle, slope {slope:g} (guide)")
    xm = np.sqrt(x0 * x1); ym = np.sqrt(y0 * y1)
    run_is_bottom = cy <= min(y0, y1) * 1.0000001
    ax.text(xm, cy / 1.3 if run_is_bottom else cy * 1.3, run_label, ha="center",
            va="top" if run_is_bottom else "bottom", fontsize=fontsize, color=color)
    rise_is_right = cx >= max(x0, x1) * 0.9999999
    ax.text(cx * 1.3 if rise_is_right else cx / 1.3, ym, label, ha="left" if rise_is_right else "right",
            va="center", fontsize=fontsize, color=color)


def zero_as_integer(ax, xdec=None, ydec=None):
    """Print the shared origin as a bare "0" so the x and y zero labels do not collide.

    xdec/ydec: number of decimals for the remaining linear tick labels (None leaves that axis unchanged)."""
    from matplotlib.ticker import FuncFormatter
    def fmt(dec):
        return FuncFormatter(lambda v, _p: "0" if abs(v) < 1e-12 else f"{v:.{dec}f}")
    if xdec is not None:
        ax.xaxis.set_major_formatter(fmt(xdec))
    if ydec is not None:
        ax.yaxis.set_major_formatter(fmt(ydec))


def schematic_axes(ax):
    """Turn an axes into a clean canvas for a schematic."""
    ax.set_aspect("equal"); ax.axis("off")


SINGLE_COLUMN = {"figS7_friction_factor"}


def save(fig, name, outdir=None):
    from pathlib import Path
    from io import BytesIO
    import os
    import tempfile
    from matplotlib.text import Text
    # Journal of Rheology: figures are included at \textwidth (164 mm in the REVTeX
    # preprint; 170 mm double column in print) or 0.5\textwidth (82 mm; 85 mm single).
    # Text must be >= 8 pt and lines >= 0.5 pt at that final size.
    target_mm = 82.0 if name in SINGLE_COLUMN else 164.0
    MIN_PT, MIN_LW = 8.25, 0.5
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.collections import Collection
    for _ in range(3):
        fig.canvas.draw()
        width_inches = fig.get_tightbbox(fig.canvas.get_renderer()).width + 0.04
        scale = (target_mm / 25.4) / width_inches   # final size / drawn size
        minimum = MIN_PT / scale
        changed = False
        for item in fig.findobj(match=Text):
            if item.get_text() and item.get_fontsize() < minimum - 1e-9:
                item.set_fontsize(minimum)
                changed = True
        if not changed:
            break
        if fig.get_layout_engine() is None and name not in {'fig04_timeline', 'fig08_roadmap', 'fig05_microstructure', 'fig08_pipe_profiles'}:
            fig.tight_layout(pad=0.8, w_pad=1.5, h_pad=1.5)
    fig.canvas.draw()
    width_inches = fig.get_tightbbox(fig.canvas.get_renderer()).width + 0.04
    scale = (target_mm / 25.4) / width_inches
    lw_min = MIN_LW / scale
    for item in fig.findobj(match=Line2D):
        if item.get_linestyle() not in ('None', '', ' ') and 0 < item.get_linewidth() < lw_min:
            item.set_linewidth(lw_min)
        if item.get_marker() not in (None, 'None', '', ' ') and 0 < item.get_markeredgewidth() < lw_min:
            item.set_markeredgewidth(lw_min)
    for item in fig.findobj(match=Patch):
        if 0 < item.get_linewidth() < lw_min and item.get_edgecolor()[3] > 0:
            item.set_linewidth(lw_min)
    for item in fig.findobj(match=Collection):
        lws = item.get_linewidths()
        if len(lws) and 0 < min(lws) < lw_min:
            item.set_linewidths([max(w, lw_min) if w > 0 else w for w in lws])
    for ax in fig.axes:
        for sp in ax.spines.values():
            sp.set_linewidth(max(sp.get_linewidth(), lw_min))
        ax.tick_params(which='major', width=max(0.6, lw_min))
        ax.tick_params(which='minor', width=max(0.5, lw_min))
    globals()['LAST_SCALE'] = scale
    destination = Path(outdir) if outdir is not None else Path(__file__).resolve().parent.parent / 'figures' 
    destination.mkdir(parents=True, exist_ok=True)
    for suffix in ['pdf', 'png']:
        buffer = BytesIO()
        fig.savefig(buffer, format=suffix)
        data = buffer.getvalue()
        if len(data) < 1000:
            raise RuntimeError(f'Incomplete figure output: {name}.{suffix}')
        fd, temporary = tempfile.mkstemp(prefix=f'.{name}_', suffix=f'.{suffix}', dir=destination)
        with os.fdopen(fd, 'wb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        Path(temporary).replace(destination / f'{name}.{suffix}')
        if (destination / f'{name}.{suffix}').stat().st_size != len(data):
            raise RuntimeError(f'Figure size mismatch: {name}.{suffix}')
    if outdir is None:
        from figdata import export
        export(fig, name, Path(__file__).resolve().parent.parent / 'data' / 'figure_data')
    plt.close(fig)
    print("saved", name)
