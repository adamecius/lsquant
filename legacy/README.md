# `legacy/` — quarantined previous-working code

This directory holds source files that were **previously working but are not built
and not tested**: orphaned single-purpose drivers not yet ported to the unified
`lsquant` entrypoint, and the standalone non-orthogonal ("Santiago") drivers kept as
a comparison corpus. Their build entries remain in `src/CMakeLists.txt`, commented
out and prefixed `# [LEGACY]`, as a record of the old target → source mapping.

`legacy/` has no `add_subdirectory` and is never linked. The commented `LINQT_APPS`
lines are its only build reference.

## Rules

- **Do not compile, link, modify, or "fix" anything in `legacy/`.** It is not part of
  the build graph and is not covered by the golden suite.
- **Read `legacy/` only when you need a previously-working function** — e.g. to port
  its physics into the `lsquant` path, or to diff non-orthogonal logic against the
  preserved implementation at merge-base `5cc0924`.
- Any capability revived from `legacy/` re-enters through the `HamiltonianOp` /
  reconstruction interfaces and must arrive with **golden coverage before** it is
  wired into the build. Goldens first, always.
- The supported, tested surface is `lsquant <command>` (`compute`, `reconstruct`,
  `inspect`) plus the oracle drivers still in `src/`. Prefer these; treat `legacy/`
  as read-only reference, never as the current API.

## What is here (Phase D, L1)

Eighteen drivers: fifteen orphaned/unported single-purpose drivers, plus the three
standalone non-orthogonal drivers (`*-nonOrth`, `kpm-nonOrth-precision-test`) kept as
the non-orthogonal comparison corpus. The core non-orthogonal logic itself stays in
the library (`chebyshev_moments`, `chebyshev_solver`, `metric_policy.hpp`), where it
is an intentional throwing placeholder with the working implementation preserved at
`5cc0924`.
