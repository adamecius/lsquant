# Phase-F figure captions (PRL/prl_style; usetex=False)

FIG. 1. Strong-thread scaling of the threaded lsquant Chebyshev DOS recursion on a
graphene (complex Hermitian) Hamiltonian of N=2000000 sites (honeycomb, 2-atom basis,
nearest-neighbour t1=-2.7 eV plus complex Haldane next-nearest hopping t2=0.10 eV,
phi=pi/2; ~9 nonzeros/row), M=512 Chebyshev moments, vs thread count p (log-log).
Solid blue circles: the SpMV kernel region ("spmv" timer); dashed orange squares: the
full DOS recursion ("compute_spectral/SpectralMoments"); dotted grey: ideal 1/p. The
SpMV scales toward the bandwidth knee; the full recursion saturates earlier because the
moment dot-products / AXPY are deliberately left serial (golden-safe). Pinned
OMP_PROC_BIND=close OMP_PLACES=cores, numactl per socket. dual-EPYC-9754.

FIG. 2. Size scaling: SpMV total time for M=256 moments at 128 threads vs system size N
(log-log); dotted grey is a proportional-to-N reference. The KPM cost is linear in N as
expected. Same graphene model and pinning as FIG. 1.

FIG. 3. SpMV matrix-stream bandwidth (GB/s) vs thread count for the FIG.-1 graphene
system; dashed black and dotted grey are the measured STREAM read-only roofline for one
socket (348 GB/s) and two sockets (700 GB/s). Peak achieved 348 GB/s
(~50% of the 2-socket roofline). Bytes counted: values(16 B)+index(4 B) per
nonzero + outer(8 B) per row, times the number of SpMV calls.
