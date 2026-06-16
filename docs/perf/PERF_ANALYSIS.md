# Performance & scaling — threaded SpMV (graphene, dual EPYC 9754)

Measured at HEAD `86d6946` + the Phase-A/B threaded-SpMV commits, on a dual
**AMD EPYC 9754** (Bergamo): 2 sockets × 128 cores = 256 cores, SMT off, 2 NUMA
nodes, 512 MiB aggregate L3 (32 × 16 MiB CCX), 1.5 TiB RAM, gcc 11.5
`-O3 -march=native` (AVX-512). Reproduce with
[`bench/run_graphene_scaling.sh`](../../bench/run_graphene_scaling.sh) +
[`bench/plot_graphene_scaling.py`](../../bench/plot_graphene_scaling.py).

## The change

The KPM moment recursion's one O(nnz) hot loop is the leaf SpMV
`SparseMatrixType::Multiply` (`y = a·M·x + b·y`). It was a single-threaded Eigen
product; **Phase A** replaced it with a row-parallel OpenMP fused kernel (bit-identical
on every recursion path — see the golden-safety note), and **Phase B** implemented
the `BatchMultiply` SpMM block kernel. The moment dot-products / AXPY are left
**serial on purpose** (threading a cross-vector reduction changes summation order and
would move the byte-exact goldens).

## System

Graphene (honeycomb, 2-atom basis), nearest-neighbour `t1 = −2.7 eV` plus a complex
Haldane-type next-nearest hopping `t2 = 0.10 eV, φ = π/2` → a genuinely **complex
Hermitian** operator, ~9 nonzeros/row. Generated with
[`bench/make_graphene.py`](../../bench/make_graphene.py) (self-contained; the
`examples/utilities/kwant2wannier` path needs kwant, absent here). Thread sweep at
`N = 2,000,000` sites (18 M nonzeros), DOS, `M = 512` Chebyshev moments.

## Result 1 — thread scaling (fixed N = 2 M)

| threads | SpMV ms/call | SpMV GB/s | % socket roofline | full recursion (ms) |
|---:|---:|---:|---:|---:|
| 1   | 33.17 | 11.3  | —     | 17930 |
| 2   | 16.92 | 22.2  | —     | 9595 |
| 4   | 11.82 | 31.8  | —     | 6969 |
| 8   | 11.67 | 32.2  | 9%    | 6864 |
| 16  | 11.81 | 31.8  | 9%    | 7053 |
| 32  | 5.64  | 66.7  | 19%   | 3924 |
| 64  | 2.53  | 148.8 | 43%   | 2301 |
| 128 | 1.67  | **225.7** | **65%** | 1786 |
| 256¹| 1.08  | **347.6** | 50%²  | 1767 |

¹ 256 threads spans both sockets (`numactl --interleave=all`).
² of the 700 GB/s 2-socket read roofline (65% of the per-socket roofline 348 is the
clean single-socket number at 128).

- **SpMV speedup: 19.9× at 128 threads (one socket), 30.7× at 256** (33.17 → 1.08 ms).
- The **4–16 thread plateau (~32 GB/s)** is the Bergamo few-CCD bandwidth ceiling;
  bandwidth resumes climbing as more CCDs engage (matches the STREAM behaviour and the
  synthetic-matrix study in `perf/`).
- **Full DOS recursion: 17.9 s → 1.77 s (10.1×)**, plateauing earlier than the SpMV
  because the serial vdot/AXPY (deliberately unthreaded) is the Amdahl floor — visible
  as the gap between the two traces in `fig_scaling_threads`. Threading those (behind an
  opt-in flag, off the golden path) is the next lever.

`perf stat` on the kernel: IPC ≈ 0.8, dTLB-miss ≈ 99% → decisively memory-bandwidth
bound (the SpMV streams the matrix with no reuse). So the ceiling is the DRAM read
roofline, and ~65% of it at one socket is the achievable envelope for a generic CSR
SpMV without matrix reordering.

## Result 2 — size scaling (fixed 128 threads, M = 256)

`N` swept 12.8 k → 2 M. The SpMV cost is **linear in N** in the DRAM-bound regime
(large N). Mid-size points (N ≈ 0.5–0.8 M) report >100 % of the DRAM roofline because
the ~150–250 MB matrix is partly **L3-cache-served** there (served faster than DRAM) —
the same cache-vs-DRAM crossover the in-cache control showed in `perf/`. Beyond L3 the
curve is firmly linear and DRAM-bound. See `fig_scaling_size`.

## NUMA policy (Phase C)

The SpMV/SpMM use `#pragma omp parallel for schedule(static)`, so each thread owns a
fixed, contiguous row range across every one of the M iterations. With kernel NUMA
balancing on (`/proc/sys/kernel/numa_balancing = 1` on this host) the OS migrates each
thread's matrix/vector pages to its local node and they stay there — the recursion's
hundreds of identical-pattern passes converge to NUMA-local placement automatically.
A clean in-process first-touch is not cleanly available (Eigen owns the matrix storage;
`std::vector` zero-initialises on the master node), and would buy only the cold-start
transient given autonuma. **Operational guidance:** for large dual-socket runs with
numa_balancing OFF, launch with `numactl --interleave=all` (measured ~440 GB/s) or pin
to one socket with `numactl --cpunodebind=0 --membind=0` (clean 65 % of roofline);
first-touch-local two-socket placement reaches ~630 GB/s (≈90 % of the dual roofline).

## Golden-safety (why threading did not move a single byte)

Each output `y[i]` is an independent reduction over row `i` in stored order, so the
result is independent of the row→thread partition (thread-invariant for *any* scalars).
The recursion only ever calls `Multiply` with `a ∈ {1, 2}`, `b ∈ {0, −1}`; multiplying
a double by an exact power-of-two / unit / zero is exact in IEEE-754, so
`a·(Σ vₚxₚ) = Σ(a·vₚ)xₚ` bit-for-bit on every recursion path. The full 29-golden suite
(DOS, MSD, Hall, spin precession, Kubo-Greenwood, Kubo-Bastin sea/surf, FFT-grid,
correlations) stays **byte-identical**, and the gates `spmv_threaded_bitexact` /
`spmm_batchmultiply_equiv` assert the two properties (thread-invariance; Eigen
bit-identity on the recursion scalars) with `==`, not a tolerance.
