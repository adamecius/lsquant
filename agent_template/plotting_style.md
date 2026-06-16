# PLOTTING STYLE — AGENT DIRECTIVE (PRL / APS)

**Status: MANDATORY.** Any agent producing a figure for this project MUST comply
with every rule below and MUST run the self-check in §8 before declaring a figure
done. Rules marked **[HARD]** come from APS *Information for Contributors* and are
not negotiable. Rules marked **[HOUSE]** are project conventions.

The reference implementation is `prl_style.py` (presets `single` / `double` /
`web`). Agents SHOULD import it rather than re-deriving rcParams. Re-deriving is
allowed only if every value below is reproduced exactly.

---

## 1. Geometry — fix width by destination, never by eye

| Target  | Width        | Source |
|---------|--------------|--------|
| single-column article | **3.375 in (8.6 cm)** | [HARD] APS column width |
| double-column article | **7.0 in** (REVTeX full text width) | [HARD] |
| 1.5 column (only if detail demands) | ~5.0 in | [HARD] allowed exception |
| web / GitHub Pages | ~6.5 in at screen dpi | [HOUSE] |

- Width is **locked** to the target. Only the **aspect ratio** is free.
- Default aspect = golden ratio; shorten height to remove empty top/bottom
  margins (APS: avoid wasted vertical space). **[HARD]**
- **NEVER** rescale a figure after export. Do **not** write
  `\includegraphics[width=\columnwidth]{...}` on a figure exported at a
  different width — that breaks the 2 mm lettering rule (§4). Export at the final
  size and include at natural size. **[HARD]**

## 2. No title — information lives in the caption

- **NEVER** call `set_title()` / `suptitle()` on a publication figure. **[HARD/HOUSE]**
- Multi-panel: tag panels `(a) (b) (c) …` **inside** the axes
  (`ax.text(0.04, 0.90, "(a)", transform=ax.transAxes, fontweight="bold")`),
  never as titles. **[HOUSE]**
- All descriptive content goes in the caption (§7).

## 3. Math in LaTeX

- All axis labels, tick formats, legend entries, and annotations use math mode:
  `r"$E/t$"`, `r"$\rho(E)\,[t^{-1}]$"`, `r"$\eta = 0.25\,t$"`. **[HOUSE]**
- **Preferred:** `text.usetex=True` with a REVTeX-matching preamble when the
  environment has `latex` **and** `dvipng` (or `dvisvgm`).
- **Fallback (no dvipng):** `text.usetex=False`, `mathtext.fontset="cm"`,
  `font.family="serif"`. Visually equivalent Computer-Modern glyphs.
- The agent MUST detect which path is available and report which it used.

## 4. Lettering and line/marker weights — the 2 mm rule

- Lettering height on the **final page** ≥ **2 mm** (~5.7 pt). In practice:
  tick labels ≥ 8 pt, axis labels ≈ 9–10 pt, at the export width. **[HARD]**
- Line weight ≥ **0.5 pt** final; project default **1.4 pt**. **[HARD]**
- Marker (symbol) size ≥ **2 mm**; project default **5 pt**. Avoid small open
  symbols that fill in on reproduction. **[HARD]**
- Ticks: inward, on all four sides; minor ticks on. **[HOUSE]**

## 5. Grayscale + colour-vision-deficiency: REDUNDANT ENCODING

Every distinct trace MUST be separable by **at least two** of the following, and
in practice all three: **[HARD principle, HOUSE implementation]**

1. **Colour** — Okabe-Ito palette only, ordered light→dark so it also separates
   in grayscale:
   `#0072B2 #D55E00 #009E73 #CC79A7 #E69F00 #56B4E9 #000000`
2. **Linestyle** — `- -- -. :` then custom dashes.
3. **Marker** — `o s ^ D v P X`, with **`markevery`** set so markers are sparse
   (mandatory on dense/analytic curves; markers may be dropped entirely if
   linestyle alone disambiguates).

Forbidden: encoding traces by colour **only**; rainbow/`jet`; red+green as the
sole contrast; fills/hatching too fine to survive 600-dpi reproduction.

## 6. Output

- Raster export at **≥ 600 dpi**; vector (PDF/EPS) preferred for line art. **[HARD]**
- `bbox_inches="tight"`, small `pad_inches`. White/clear background. **[HARD]**
- File naming: `fig_<topic>_<target>.{pdf,png}`.

## 7. CAPTION — required, and the figure is incomplete without it

There is no title, so the caption carries **all** interpretive information and
MUST be self-contained (intelligible without the body text). The agent MUST emit
the caption alongside the figure (in the `.tex`, the PR description, or a sidecar
`.caption.txt`). **[HARD + HOUSE]**

A compliant caption MUST:

1. **Begin** with `FIG. N.` — `FIG.` in all caps, the arabic numeral, a period.
   (REVTeX: use `\caption{...}`; the `FIG. N.` prefix is auto-generated, so write
   the text after it.) **[HARD]**
2. Be **one concise paragraph**; a sentence fragment is acceptable.
3. **Define every visual encoding actually used**: what each colour / linestyle /
   marker represents (e.g. "solid blue: η = 0.25 t; dashed orange: η = 0.33 t").
   Never leave a curve unexplained.
4. State **all units and axes** if not fully obvious from the labels.
5. State the **physical/numerical parameters** needed to reproduce the panel
   (system size, broadening, KPM moment count, k-mesh, temperature, disorder
   realization count, etc.).
6. Reference **panels by tag**: "(a) …; (b) …".
7. Put any **cited work** as a normal reference (APS: figure-caption citations go
   into the reference list, not as inline footnotes).
8. **No information appears only in a title** (there is none) or only in the body
   text that the reader needs to read the figure.

### Caption template (fill every slot, delete none silently)
```
FIG. N. <one-line statement of what is plotted> as a function of <x quantity>.
<Per-trace legend in prose: encoding -> meaning, for each trace.>
<Parameters: system, size, broadening/eta, KPM moments M, k-mesh, T, disorder>.
(a) <panel a>; (b) <panel b>; ... Units: <if not on axes>. <Method note / Ref [n]>.
```

### Worked example
```
FIG. 2. Local density of states rho(E) of the linear-chain tight-binding model
versus energy E/t for four Lorentzian broadenings: solid blue eta = 0.25 t,
dashed orange 0.33 t, dash-dotted green 0.41 t, dotted purple 0.49 t. Curves are
KPM reconstructions with M = 3000 Chebyshev moments and Jackson kernel on a chain
of N = 2^20 sites; markers sample every 45th point and carry no extra data.
```

## 8. SELF-CHECK — run before declaring any figure done

The agent MUST verify and report PASS/FAIL on each:

- [ ] Figure width equals the target column width; not rescaled on inclusion.
- [ ] No `set_title`/`suptitle`; panels tagged inside axes.
- [ ] All text is LaTeX/mathtext; usetex-vs-mathtext path reported.
- [ ] Tick labels ≥ 8 pt and lettering ≥ 2 mm at export size.
- [ ] Lines ≥ 0.5 pt; markers ≥ 2 mm; `markevery` set on dense curves.
- [ ] Every trace separable in grayscale AND under deuter/protanopia
      (Okabe-Ito + linestyle + marker). No colour-only encoding.
- [ ] Export ≥ 600 dpi (raster) or vector; tight bbox; white background.
- [ ] A `FIG. N.` caption exists, is self-contained, and defines every encoding
      and every reproduction parameter (§7 checklist all satisfied).

A figure that fails any box is **not done**. Fix or escalate; do not ship.

## 9. Minimal compliant code pattern
```python
from prl_style import prl_context   # presets: 'single' | 'double' | 'web'

with prl_context("single", usetex=False):   # usetex=True if dvipng present
    fig, ax = plt.subplots()
    ax.plot(E, rho, label=r"$\eta=0.25\,t$", markevery=55)  # cycle gives style+marker
    ax.set_xlabel(r"$E/t$")
    ax.set_ylabel(r"$\rho(E)\,[t^{-1}]$")
    ax.legend(loc="upper right")          # frameless, set in rcParams
    fig.savefig("fig_dos_single.pdf")     # 600 dpi / vector, tight bbox
# Caption emitted separately per §7.
```
