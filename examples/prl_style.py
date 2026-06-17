"""prl_style - APS/PRL publication style for LSQUANT matplotlib figures.

Publication figures, not slides. The rules below are split into HARD specs (APS
"Information for Contributors", invariant across PRL/PRB/PRA/PRE) and the one
NEGOTIABLE knob (aspect ratio / target).

HARD specs
  * Width is fixed by destination, never by eye:
        single-column = 3.375 in (8.6 cm),  double-column = 7.0 in,
        web (e.g. a GitHub page) ~ 6.8 in.
    Only the aspect ratio is free; default to the golden ratio, trim height if it
    leaves empty margins. Never rescale the figure after export (the 2 mm / ~5.7 pt
    minimum-lettering rule binds at final size).
  * No title -- the scientific message lives in the manuscript caption. Multi-panel
    figures carry (a), (b), ... tags INSIDE the axes (use panel()).
  * Math in LaTeX/mathtext: axis labels, legends, annotations use math mode,
    e.g. r"$\\sigma_{xy}\\,(e^2/h)$". usetex=True is used when latex+dvipng are
    present; otherwise mathtext with the Computer-Modern fontset gives the same
    look with no TeX dependency (set by use()).
  * Redundant encoding for greyscale + colour-vision-deficient readers: every trace
    is separated THREE ways at once -- an Okabe-Ito colour (ordered by luminance so
    it also separates in greyscale), a distinct linestyle, and a distinct marker.
    Get all three from trace(i).
  * Weights above the floor: line >= 0.5 pt (we use ~1.4), markers ~4.5 pt; inward
    major+minor ticks on all four sides; frameless legend; no heavy grid or
    decorative background. Raster export at 600 dpi; prefer vector (PDF/SVG) for
    articles.

NEGOTIABLE
  Only figsize, base font size, and dpi differ per target -- the PRESETS table.

Usage
    import prl_style as prl
    prl.use("single")                       # or "double" / "web"
    fig, ax = plt.subplots()
    for i, (x, y, lbl) in enumerate(series):
        ax.plot(x, y, label=lbl, **prl.trace(i))
    ax.set_xlabel(r"$E/t$"); ax.set_ylabel(r"$\\rho(E)\\,(1/t)$")
    prl.save(fig, "fig_dos")                 # writes fig_dos.pdf (+ .png preview)
"""

GOLDEN = (5.0 ** 0.5 - 1.0) / 2.0          # 0.618

# Okabe-Ito, ordered by luminance (dark -> light) so a greyscale print also
# separates the traces in this order.
COLORS  = ["#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
LINES   = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (1, 1))]
MARKERS = ["o", "s", "^", "D", "v", ">", "P"]

PRESETS = {
    "single": dict(width=3.375, base=9.0,  dpi=600),   # APS single column, 8.6 cm
    "double": dict(width=7.0,   base=9.5,  dpi=600),   # APS double column, ~17.2 cm
    "web":    dict(width=6.8,   base=13.0, dpi=150),   # GitHub page / slide, bigger type
}


def trace(i, line=True, marker=False):
    """Return {color, linestyle, marker} for trace i -- three redundant channels.

    line/marker toggle which channels are active (a smooth curve usually wants
    color+linestyle; sampled data wants color+marker)."""
    s = {"color": COLORS[i % len(COLORS)]}
    if line:   s["linestyle"] = LINES[i % len(LINES)]
    if marker: s["marker"]    = MARKERS[i % len(MARKERS)]
    return s


def use(target="single", aspect=GOLDEN, height=None, usetex=None):
    """Install the rcParams for `target`. Returns the (width, height) in inches."""
    import shutil
    import matplotlib as mpl

    if target not in PRESETS:
        raise ValueError("target must be one of %s" % list(PRESETS))
    p = PRESETS[target]
    w = p["width"]
    h = height if height is not None else w * aspect
    base = p["base"]

    if usetex is None:                       # autodetect a usable TeX stack
        usetex = bool(shutil.which("latex") and shutil.which("dvipng"))

    rc = {
        "figure.figsize": (w, h),
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": p["dpi"],
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,

        "text.usetex": usetex,
        "font.family": "serif",
        "mathtext.fontset": "cm",            # Computer Modern look without TeX
        "axes.unicode_minus": False,

        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": base - 1.0,
        "xtick.labelsize": base - 1.0,
        "ytick.labelsize": base - 1.0,

        "axes.linewidth": 0.6,
        "lines.linewidth": 1.4,
        "lines.markersize": 4.5,
        "lines.markeredgewidth": 0.7,

        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,

        "axes.grid": False,
        "legend.frameon": False,
        "legend.handlelength": 2.2,
        "legend.labelspacing": 0.3,
    }
    if usetex:
        rc["font.serif"] = ["Computer Modern Roman"]
    mpl.rcParams.update(rc)
    return (w, h)


def panel(ax, tag, loc="upper left", dx=0.04, dy=0.04, **kw):
    """Place a panel tag like '(a)' INSIDE the axes (not a title)."""
    xa = dx if "left" in loc else 1.0 - dx
    ya = 1.0 - dy if "upper" in loc else dy
    ha = "left" if "left" in loc else "right"
    va = "top" if "upper" in loc else "bottom"
    ax.text(xa, ya, tag, transform=ax.transAxes, ha=ha, va=va,
            fontweight="bold", **kw)


def save(fig, path_noext, formats=("png",), dpi=None):
    """tight_layout + save in each format. PNG is the committed/web copy (renders on
    GitHub); pass formats=("pdf","png") to also emit the vector article copy. Pass a
    bare stem, e.g. save(fig, 'fig_dos'). Use dpi to override the per-target dpi for a
    single dense figure (keeps a committed PNG under the repo's 200 KB limit)."""
    fig.tight_layout(pad=0.3)
    written = []
    for ext in formats:
        out = "%s.%s" % (path_noext, ext)
        fig.savefig(out, dpi=dpi) if dpi else fig.savefig(out)
        written.append(out)
    print("wrote", ", ".join(written))
    return written
