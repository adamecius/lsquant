# Original orthogonal downfolding tools

Unmodified snapshot supplied by the maintainer from `/root/scratch/gb` on
Ionos, copied into this checkout before stage 6d development on 2026-10-23.
No source revision was supplied. `SHA256SUMS` records the exact supplied bytes;
the origin Python files are kept unchanged, including their original language
and conventions. The other scripts are retained as provenance, not imported by
LSQUANT. No external project history or licence is inferred.

`downfold.py` implements shell/polar and frozen-window reduction;
`downfold_scdm.py` implements SCDM and imports the former. Both validate bands
on deterministic off-grid points and a symmetry path. `siesta2w90_sisl.py`
and `study.py` record the input and basis-study conventions. No separate notes
or paper PDFs were included in the supplied snapshot. Literature citations
are in the source docstrings; paper-fetch scripts are archival utilities.

Standalone usage (NumPy 2.2, SciPy 1.15, sisl 0.16):

```sh
python tools/downfold_origin/downfold.py RUN_DIRECTORY LABEL --out OUTPUT
python tools/downfold_origin/downfold_scdm.py RUN_DIRECTORY LABEL --out OUTPUT
```

The LSQUANT wrapper and FFT/Wigner-Seitz reuse live in
`python/lsquant/projection/`; changes belong there, never in this snapshot.
