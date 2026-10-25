# Profiling and scaling results

Status: diagnostic measurements complete; optimization and project closure pending.
Measurement date: 2026-10-25. Host: Aron, two AMD EPYC 7552 sockets,
48 physical cores per socket, 96 total. This report consolidates the evidence in
[tests.md](tests.md). The user retains the final closure decision.

## Conclusion

The observed plateau is **not an established fundamental hardware limit**.
There is a specific implementation limit in the tested native
`linalg::MomentTable`: a 128-by-128 output split into fixed 16-by-16
tiles produces only 64 independent multiplication tasks. That multiplication
phase cannot keep all 96 workers busy. This limit is changeable and applies
to this batch size, not to every moment calculation.

This explains a parallelism ceiling, but does not by itself establish the
cause of the measured slowdown. Packing, allocation, memory placement,
inter-socket traffic, scheduling and small matrix multiplication costs have
not been timed separately. Hardware bandwidth saturation has not been measured.

Conventional changes are worth testing before any specialized implementation:
better thread and tile selection, reusable buffers, and reuse of coefficient
tables. No custom SIMD, assembly, GPU port or engine rewrite is justified by
the present evidence. No gain from these proposals has yet been measured.

## Material profiles: where the time goes

The existing reduced Bi and PdSe2 workloads completed with exit code 0.
Settings: 24x24 mesh, covariant velocity, exact Lowdin metric, shift projection
[-5,5] eV, W=0, M=256 followed by M=1024, R=16, seed=12345, original table
reconstruction. Eight OpenMP threads were requested on CPUs 0-7; BLAS was
limited to one thread.

These are **cProfile timings with the golden build**, not unbiased benchmarks.
They include both moment orders, setup, oracle evaluation and output.

| Material | Total | Moments including preparation | SHC reconstruction | Other | Peak RSS |
|---|---:|---:|---:|---:|---:|
| Bi | 304.753 s | 207.748 s (68.2%) | 85.793 s (28.2%) | 11.212 s | 2.79 GiB |
| PdSe2 | 1070.146 s | 859.114 s (80.3%) | 84.509 s (7.9%) | 126.523 s | 17.00 GiB |

The relevant function chain is
`compute_projection_moments -> compute_moments -> engine.run`.
cProfile does not separate the internal C++ engine kernels. It therefore
does not establish how much of that material cost comes from recursion,
operator applications or accumulation.

Within reconstruction:

| Function | Bi | PdSe2 | Interpretation |
|---|---:|---:|---|
| `_core.ContractRank2Imag` | 66.844 s | 65.343 s | Main measured contraction cost |
| `legacy_delta_table` | 15.206 s | 15.399 s | Python coefficient generation, including cosine calls |

The first two rows are directly supported by the retrieved profiles.
Nested function times are already included in the reconstruction totals.
The material CPU-time/wall-time ratios were 7.26 for Bi and 6.82 for PdSe2.
They do not support a whole-run single-core bottleneck. The recorded affinity
[0] is the calling thread after OpenMP binding; it is not a worker-affinity
inventory.

## Unprofiled native kernel scaling

The existing stage-8 driver was rebuilt using candidate-v5 headers and the
integrated native static library. The candidate-v1 and candidate-v5
`lsquant/core` trees compared identical, permitting reuse of the native
library. These measurements are separate from the golden material profiles.

Each point is the median of three fresh processes. Each process excludes
setup and one warm-up, then times three operations. Runs were sequential;
CPU affinity was 0 through T-1, with close binding, core places, dynamic
OpenMP teams disabled and BLAS limited to one thread. Ninety-six cores span
both sockets.

| Function / workload | 1 core | 8 cores | 24 cores | 48 cores | 96 cores |
|---|---:|---:|---:|---:|---:|
| `linalg::MomentTable`, dimension 165888, batch 128 | 1835.867 ms | 290.141 ms | 147.750 ms | 133.418 ms | 174.241 ms |
| `ContractRank2Imag`, M=256, 7680 points | 181.227 ms | 23.438 ms | 9.050 ms | 6.475 ms | 4.606 ms |

| Function | Speedup at 24 | Speedup at 48 | Speedup at 96 |
|---|---:|---:|---:|
| `MomentTable` | 12.43x | 13.76x | 10.54x |
| `ContractRank2Imag` | 20.03x | 27.99x | 39.35x |

`MomentTable` is the confirmed poor-scaling function in these measurements:
48 cores are only 1.11x faster than 24, and 96 cores are **30.6% slower**
than 48. Its trial ranges at 48 and 96 cores were 131.809-137.722 ms and
170.266-176.364 ms, respectively. The slowdown exceeds the observed run spread.

`ContractRank2Imag` continues improving: 96 cores are 1.41x faster than 48.
Neither result predicts the speedup of a complete material run.

Output checksums agreed at every thread count and repetition:
`8221316592892497927` for MomentTable and
`11479475180656984162` for contraction. This does not replace golden gates.

## Straightforward improvements to evaluate

1. **Choose threads per workload.** Forty-eight cores are the fastest tested
   choice for this MomentTable case. This is a configuration finding, not a
   global recommendation to cap every function at 48: contraction benefits
   from 96. Independent calculations can also be benchmarked one per socket
   for throughput; that does not reduce a single calculation's latency.

2. **Measure packing and multiplication separately, then tune tiles.**
   Inspect [MomentTable](../../lsquant/core/algebra/linalg.cpp).
   An 8-by-8 tile would provide 256 tasks for batch 128, but smaller products
   can be less efficient and increase repeated reads. Compare a small fixed
   set of tile choices with the existing harness; do not assume smaller is
   faster. With batch 1024 the current tile already produces 4096 tasks,
   so the same shortage-of-tasks explanation does not apply.

3. **Reuse allocated and packed buffers where lifetimes permit.**
   MomentTable currently constructs and fills two dense matrices per call.
   Workspace reuse may avoid repeated allocation/packing overhead. Preserve
   initialization and NUMA placement semantics, and measure before changing
   data layouts or the engine interface.

4. **Reuse coefficient tables for identical reconstruction settings.**
   `legacy_delta_table` is a visible serial Python cost. A bounded cache
   keyed by every relevant grid/order/convention parameter can avoid repeated
   generation where inputs match. Vectorization is another candidate, but
   trigonometric rounding can change legacy byte identity.

5. **Compare the existing folded reconstruction on the same saved moments.**
   This option already exists; it requires no new contraction algorithm.
   Record numerical residuals and end-to-end time. Keep the legacy-compatible
   route until the appropriate acceptance is obtained.

These are proposals, not implemented optimizations. First isolate the native
material engine cost so that kernel improvements target the actual production
workload. Do not multiply a native microbenchmark gain by the golden profile's
68-80% engine share.

## Numerical and measurement limits

Candidate-v5 M1024 reference-point SHC values differ from the earlier saved
reduced runs:

| Material / point | Earlier saved run | Current profiled run |
|---|---:|---:|
| Bi reference | 1.1353409434012425 | 1.1550715643600988 |
| PdSe2 upper gap | 2.520384253866164 | 2.3441496803075585 |

These differences have not been explained or accepted. Verify the source,
operator preparation and actual trace vectors, not just the declared seed,
before an end-to-end before/after comparison. Profiler overhead alone is not
an explanation for changed numerical output.

Hardware counters are unavailable: perf is absent and perf_event_paranoid=4.
No measured DRAM bandwidth, cache-miss rate or hardware roof is claimed.
A different tile or reduction order may affect floating-point results; no
golden, seed or tolerance may be relaxed to accommodate an optimization.

## Evidence and next acceptance step

Remote evidence:
`~/lsquant-wrapup-20261024/profiling-20261025/`

- `bi/` and `pdse2/`: profile.pstats, profile.txt, receipt.json,
  run.log, exit.txt and per-case numerical outputs.
- `scaling/`: 30 per-process JSON receipts, summary.json,
  provenance.sha256, core-source-diff.txt, launcher and completion status.
- Candidate source tree: `bf70c9b8c8b8ec66df329265d6c170a1cb315e0e`.
- Material Python source: candidate-v5; compiled golden engine: build-v1-golden.
- Native scaling library: build-v1-native, verified unchanged core sources.

Both profiling controllers and the targeted scaling campaign have finished.
No new optimization was applied for this report. Next: resolve the numerical
comparison, attribute native engine phases, then test the smallest proposed
change with the existing harness, unchanged numerical checks and unprofiled
whole-run timing. W08/W09 and project closure remain incomplete.

