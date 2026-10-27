# Legacy defect D2: energy axis of the uniform-grid rank-2 reconstruction (evidence only)

`KuboBastin_NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexactJACKSON.dat` is the legacy
`inline_kuboBastinFromChebmom <mom> 10` output for the G4c moments
(`goldens/stage0/graphene/offcentre_enclosing/`, graphene VX-VY, `BOUNDS` = `-8.5 7.5`, exact
trace, M = 128), produced by this repository's legacy tree (`build_legacy`, see
`PROVENANCE.md`, G4c). It is **not a reference**; no harness compares against it.

The legacy `src/reconstruction.cpp` uniform-grid path writes the energy column with
`e / sf - shf` (`sf = ScaleFactor() = 1/HalfWidth`, `shf = ShiftFactor() = -BandCenter/HalfWidth`),
i.e. `E = x*HalfWidth + BandCenter/HalfWidth` instead of `x*HalfWidth + BandCenter`: wrong for
every off-centre window (it affects Bastin, Sea, Surf and Greenwood on the uniform grid,
`reconstruction.cpp:147-148, 199, 223`). The FFT-grid path (`e * hw + bc`,
`reconstruction.cpp:150`) is correct. For centred windows (`BandCenter = 0`, all other stage-0
goldens) both forms agree.

Compared with the new code (`python -m lsquant.apps.stage3_kubo`, review 3, `aba179a`), the
value column of this file is byte-identical as text on all 3839 lines; the energy column is
shifted by `BandCenter/HalfWidth - BandCenter = -0.5/8 + 0.5 = +0.4375` eV on every line.
The G4c curve golden was therefore blessed from the new code (see `PROVENANCE.md`, G4c).

`BOUNDS` is the window file the legacy run used.
