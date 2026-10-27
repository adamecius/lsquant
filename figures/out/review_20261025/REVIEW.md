# Bi and PdSe2 material review figures

Generated on Aron, 2026-10-25, from the frozen SIESTA inputs selected for
the saved stage-6b references. Production library source tree:
bf70c9b8c8b8ec66df329265d6c170a1cb315e0e.
Numerical and rendering script hashes are recorded separately in each
material JSON. The frozen input hashes match the projection campaign.

Files: fig_bi_material_review_double.png and
fig_pdse2_material_review_double.png; SVG versions contain vector line
art and rasterized spin-color fields. The corresponding caption.txt
files define all panels and conventions. The NPZ files contain the
plotted arrays, including the full band range.

Each figure has bands, DOS, longitudinal conductivity and spin Hall
response across the top, with original and projected upper-valence
spin textures underneath. Gray solid lines are the full model; blue
dashed lines are the retained shift-projected subspace. Orange dotted
lines in the spin Hall inset are the existing reduced M=1024, R=16 KPM
calculation, not a newly accepted result.

Settings: clean W=0; 24x24 mesh; covariant physical velocities; exact
Lowdin metric; shift projection in absolute stored [-5,5] eV; no radius
truncation. The band/DOS display is centred on mu_ref=0 eV for Bi and
mu_ref=1.052 eV for PdSe2, relative to stored SIESTA Fermi zero. The
latter selects the previously studied upper gap, not the native Fermi
level. Gaussian widths are 50 meV for DOS and charge conductivity and
10 meV for spin Hall occupations. The finite-width charge response is
a clean finite-mesh diagnostic, not a converged DC conductivity.

| Quantity | Bi | PdSe2 |
|---|---:|---:|
| Retained-band path MAE, meV | 2.3435e-12 | 2.5550e-12 |
| Maximum spin-texture difference, Pauli units | 9.1038e-15 | 8.6597e-15 |
| Full exact SHC at reference | 1.14415384 | 1.97909882 |
| Projected exact SHC at reference | 1.07975850 | 1.98650293 |
| Saved projected KPM M=1024, R=16 SHC | 1.13534094 | 2.52038425 |

Spin textures use the actual nondegenerate upper valence state wherever
available. At energy-degenerate points, the uniquely maximal eigenstate
of physical sigma_z within the valence multiplet is used. This convention
is stated explicitly; arbitrary degenerate eigenvectors are not plotted.
Bi's sampled valence multiplicity is 1; PdSe2 has multiplicities 1 and 2.
No unresolved spin points occur. The projected state and its spin
expectations are evaluated explicitly, not copied from the original.

## Figure style and visual checks

- PASS: exported width exactly 4200 pixels at 600 dpi, seven inches.
- PASS: no axes titles or figure title; panel tags inside the axes.
- PASS: Computer Modern mathtext fallback; no usable TeX stack selected.
- PASS: ticks at least 8 pt; default ticks 8.5 pt.
- PASS: line weights at least 0.5 pt; no point-marker traces.
- PASS: full/projected/KPM curves have redundant solid/dashed/dotted
  encoding. Original gray is explicitly requested by the user; other
  curve colors use Okabe-Ito. Spin color is supplemented by the numeric
  color bar, and arrows represent in-plane spin.
- PASS: 600 dpi PNG and vector SVG; white background and tight export.
- PASS: separate self-contained FIG. 1 / FIG. 2 captions.
- PASS: both final PNGs visually inspected; band and SHC insets, axis
  labels, legends and color bars are legible and within the image.

The initial exploratory rendering showed separate spin sectors; it was
replaced after observing a nondegenerate Bi valence band. Final panels
compare original and explicitly projected textures. Original high-energy
bands remain visible in the band inset. The KPM trace is confined to the
spin Hall gap inset so its noisy off-gap tails do not set the overview
scale. No curve is rescaled or normalized to force agreement.

These figures are review artifacts. They do not close the material,
disorder, performance or final-validation gates.
