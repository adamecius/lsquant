"""Stage-6c adapters: FFT disorder with block metric, and direct CSR reference.

All physical matrices enter the existing word compiler before the metric.
Only LinearOperators and StateGenerators enter the C++ engine. In particular
no per-k eigentrace or per-k oracle is used for a disordered calculation.
"""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np
import scipy.sparse as sparse

from . import _core
from .bases import BlochMetricBasis, make_anderson_disorder
from .models.disorder import anderson_potential, checked_csr, supercell_csr
from .models.nonorthogonal import applied_spin_operators, exact_block_metric


def _transform(model, supercell):
    points = model.k_mesh(supercell)
    phases = np.ascontiguousarray(np.exp(-1j*(points @ model.tau.T))).ravel()
    return _core.BlochTransform(*supercell, model.num_wann, phases,
                                _core.BlochDirection.RealToBloch)


def _values(dimension, width, seed, values):
    if not np.isfinite(width) or width < 0:
        raise ValueError("disorder W must be nonnegative and finite")
    w = make_anderson_disorder(dimension, width, seed) if values is None else np.asarray(values)
    if w.shape != (dimension,) or not np.all(np.isfinite(w)) or np.any(w.imag != 0):
        raise ValueError("disorder must be a finite real vector with one value per coefficient")
    return np.ascontiguousarray(w, dtype=np.complex128)


def _triple(first, middle, last, coefficient=0.5):
    """metric_forms: c*A*X*B+conj(c)*B*X*A with Hermitian factors."""
    return _core.HermitianTripleProduct(first, middle, last, coefficient)


def _velocity(H, X, dH, dS, A=None):
    """metric_forms.velocity_forms, or disorder_velocity_forms after H+=V.

    Derivatives are periodic position derivations, including dV={D,dS}/2.
    X applications stay matrix-free, even when a velocity is not sparse.
    """
    result = _core.OperatorSum(dH, _triple(dS, X, H, -0.5))
    if A is not None:
        result = _core.OperatorSum(result, _triple(H, X, A, 1j))
    return result


def _bind_metric_trace(basis, plan, operators, trace_mode):
    """Same dual bra as 6b review 2, for FULL-system seeds on both routes.

    K=Y on the generalised route, K=I on Lowdin; the copied original r is
    prepared on BOTH sides. The compiler's asymmetric bra then carries
    S*T_m*X*A*K*r. Never use independent bra phases or block-local seeds.
    """
    if trace_mode not in ("random", "exact"):
        raise ValueError("trace must be random or exact")
    operators = list(operators)
    K = basis.metric.inverse_sqrt_overlap if basis.route == "generalised" else basis.identity
    dual = basis.overlap() if basis.route == "generalised" else basis.identity
    ki, bi = len(operators), len(operators)+1
    operators.extend([K, dual])
    code, mirrored = [_core.ApplyOperator(ki)], False
    for instruction in plan.code:
        if isinstance(instruction, _core.MirrorToBra):
            code.extend([instruction, _core.ApplyOperator(ki), _core.ApplyOperator(bi)])
            mirrored = True
        elif isinstance(instruction, _core.Trace) and not mirrored:
            code.extend([_core.MirrorToBra(), _core.ApplyOperator(ki), _core.ApplyOperator(bi), instruction])
        else:
            code.append(instruction)
    return _core.Plan(code, plan.word), operators


def _states(basis, trace_mode, num_states, seed):
    if trace_mode == "random":
        count = 1 if num_states is None else int(num_states)
        if count < 1:
            raise ValueError("random trace needs a positive state count")
        return _core.RandomPhaseStates(basis.dimension, count, basis.seed if seed is None else int(seed))
    if trace_mode == "exact":
        if num_states is not None and int(num_states) != basis.dimension:
            raise ValueError("a disordered exact trace requires the full dimension, not W block states")
        return _core.LocalStates(basis.dimension, list(range(basis.dimension)))
    raise ValueError("trace must be random or exact")


def _window(clean_bounds, overlap_bounds, values, pad):
    """Weyl enclosure: ||Y V Y|| <= sqrt(cond(S))*max|w|.

    Y V Y=(Y D S^1/2+S^1/2 D Y)/2. W/2 alone is not an enclosure
    for noncommuting S and D. Both adapters use the same finite-mesh
    clean spectrum and overlap extrema, independent of stochastic seed.
    """
    if not np.isfinite(pad) or pad < 0:
        raise ValueError("window padding must be nonnegative and finite")
    smin, smax = overlap_bounds
    delta = np.sqrt(smax/smin)*float(np.max(np.abs(values)))
    lo, hi = clean_bounds[0]-delta, clean_bounds[1]+delta
    margin = pad*max(hi-lo, 1.)
    return _core.SpectralWindow(lo-margin, hi+margin)


class DisorderedBlochMetricBasis(BlochMetricBasis):
    """H(k)+{B D B^+,S(k)}/2 with the unchanged exact per-block metric.

    Here B is stage 5's RealToBloch transform. The brief writes the
    real-directed transform as B and thus calls the same map B^+ D B.
    No N_supercell by N_supercell sparse or dense matrix is allocated.
    """

    def __init__(self, model, supercell, *, disorder_w=0., disorder_seed=7,
                 disorder=None, velocity_mode="covariant", route="lowdin", seed=12345):
        if route not in ("lowdin", "generalised"):
            raise ValueError("stage6c requires lowdin or generalised metric")
        super().__init__(model, supercell, velocity_mode=velocity_mode, route=route, seed=seed)
        self.disorder_w, self.disorder_seed = float(disorder_w), int(disorder_seed)
        self.disorder_values = _values(self.dimension, self.disorder_w, self.disorder_seed, disorder)
        self.has_disorder = bool(np.any(self.disorder_values))
        self.transform = _transform(model, self.supercell)
        self.identity = _core.IdentityOperator(self.dimension)
        self._D = _core.SandwichOperator(self.transform, _core.DiagonalOperator(self.disorder_values))
        # Exact identity only: preserve stage 5's operation order in its
        # orthogonal limit, rather than accumulating two half-potentials.
        self._orthogonal_overlap = np.array_equal(
            self.sk_array, np.broadcast_to(np.eye(self.num_wann), self.sk_array.shape))
        self._potential = (self._D if self._orthogonal_overlap else
                           _triple(self._D, self.identity, self._s_op))
        self._clean_h_op = self._h_op
        self._clean_bounds = (float(self.eigenvalues.min()), float(self.eigenvalues.max()))
        overlap_eigenvalues = np.linalg.eigvalsh(self.sk_array)
        self.overlap_bounds = (float(overlap_eigenvalues.min()), float(overlap_eigenvalues.max()))
        if self.has_disorder:
            self._h_op = _core.OperatorSum(self._clean_h_op, self._potential)
            if self._orthogonal_overlap:
                self._applied_h_op = self._h_op
            elif route == "lowdin":
                self._applied_h_op = _core.LowdinOperator(self.metric.inverse_sqrt_overlap, self._h_op)
            else:
                self._applied_h_op = _core.MetricInverseProduct(self.metric.inverse_overlap, self._h_op)
        self._disordered_ops = {}

    def potential(self):
        """Return the Hermitian physical potential; S is left unchanged."""
        return self._potential

    def operator(self, name):
        """Physical disorder-dependent velocity/spin current; metric added by compiler."""
        key = name.upper()
        if key in ("H", "HAM"):
            return self.hamiltonian()
        if not self.has_disorder or key in ("1", "I", "IDENTITY", "S", "SX", "SY", "SZ"):
            return super().operator(key)
        if key not in self._disordered_ops:
            if key in ("VX", "VY", "VZ"):
                axis = "XYZ".index(key[1])
                _, ds = self.model.derivatives(self.k_points, axis, cell=self.velocity_mode == "bare")
                dS = _core.BlockDiagonalOperator(ds)
                # metric_forms.disorder_velocity_forms: clean velocity plus
                # delta(V)-{dS X V+V X dS}/2+i(V X A-A X V).
                dV = _triple(self._D, self.identity, dS)
                correction = _core.OperatorSum(dV, _triple(dS, self.metric.inverse_overlap,
                                                                          self._potential, -0.5))
                if self.velocity_mode == "covariant":
                    A = _core.BlockDiagonalOperator(self.model.Ak(self.k_points, axis))
                    correction = _core.OperatorSum(correction, _triple(
                        self._potential, self.metric.inverse_overlap, A, 1j))
                clean_velocity = _core.BlockDiagonalOperator(BlochMetricBasis.operator_blocks(self, key))
                self._disordered_ops[key] = _core.OperatorSum(clean_velocity, correction)
            elif len(key) == 4 and key[0] == "J" and key[1] in "XYZ" and key[2] == "S" and key[3] in "XYZ":
                self._disordered_ops[key] = _triple(self.operator("V"+key[1]), self.metric.inverse_overlap,
                                                    self.operator("S"+key[3]))
            else:
                raise ValueError(f"unknown stage6c operator {name!r}")
        return self._disordered_ops[key]

    def operator_blocks(self, name):
        """Only S/spin stay block diagonal; never return clean disorder operators."""
        key = name.upper()
        if self.has_disorder and key not in ("S", "SX", "SY", "SZ"):
            raise ValueError("disordered H, velocity and spin current have no per-k blocks; use operator()")
        return super().operator_blocks(key)

    def states(self, trace_mode="random", num_states=None, seed=None):
        """Whole-layout orbital Bloch phases, as in 6b review 2."""
        return _states(self, trace_mode, num_states, seed)

    def bind_eigen_trace(self, plan, operators, trace_mode="random"):
        """Keep the 6b prepare_moments interface; there is no per-k eigentrace."""
        if self._orthogonal_overlap and self.route == "lowdin":
            return plan, operators
        return _bind_metric_trace(self, plan, operators, trace_mode)

    def window(self, pad=0.10):
        return _window(self._clean_bounds, self.overlap_bounds, self.disorder_values, pad)

    def eigenbasis(self):
        raise ValueError("stage6c uses the real-space reference; a disordered eigenbasis is not per-k")


class RealSpaceMetricBasis:
    """Direct SIESTA supercell CSR reference with stage-6 overlap appliers.

    H,S,dH,dS,A and spin are assembled from real-space blocks. Velocities
    and spin currents remain compositions: S^-1 is generally dense, so
    materialising those products as CSR would destroy linear storage.
    Exact inverse/square root is for small systems only; Chebyshev uses
    the exact finite-mesh overlap bounds, NOT a guessed positive bound.
    """
    name = "real"

    def __init__(self, model, supercell, *, disorder_w=0., disorder_seed=7, disorder=None,
                 velocity_mode="covariant", route="lowdin", seed=12345,
                 applier="exact", order=256, coeff_tol=1e-14, max_dim=4096, tol=1e-11):
        if route not in ("lowdin", "generalised") or applier not in ("exact", "chebyshev"):
            raise ValueError("stage6c needs lowdin/generalised and exact/chebyshev")
        if velocity_mode not in ("bare", "berry_connection", "covariant"):
            raise ValueError("unknown velocity mode")
        if velocity_mode == "covariant" and model.dipoles is None:
            raise ValueError("covariant velocity requires SIESTA dipoles")
        self.model, self.supercell = model, tuple(supercell)
        self.k_points = model.k_mesh(self.supercell)  # validates dimensions
        self.num_cells, self.num_wann = len(self.k_points), model.num_wann
        self.dimension = self.num_cells*self.num_wann
        if applier == "exact" and self.dimension > max_dim:
            raise ValueError("exact supercell metric exceeds max_dim; select chebyshev")
        self.velocity_mode, self.route, self.seed = velocity_mode, route, int(seed)
        self.applier = applier
        self.disorder_w, self.disorder_seed = float(disorder_w), int(disorder_seed)
        self.disorder_values = _values(self.dimension, self.disorder_w, self.disorder_seed, disorder)
        self.has_disorder = bool(np.any(self.disorder_values))
        self.transform = _transform(model, self.supercell)
        self.identity = _core.IdentityOperator(self.dimension)
        self._matrices, self._owners, self._op_cache = {}, {}, {}
        self.overlap_csr = supercell_csr(model, self.supercell, model.overlaps)
        self.potential_csr = anderson_potential(self.overlap_csr, self.disorder_values)
        self.hamiltonian_csr = checked_csr(supercell_csr(model, self.supercell, model.hoppings)+self.potential_csr,
                                         name="disordered Hamiltonian", tolerance=model.input_tolerance)
        self._s_op = self._csr_operator("S", self.overlap_csr)
        self._h_op = self._csr_operator("H", self.hamiltonian_csr)
        self._D = _core.DiagonalOperator(self.disorder_values)
        self._potential = self._csr_operator("potential", self.potential_csr)
        # Only bounds are borrowed from translation symmetry of the CLEAN
        # problem. Every numerical metric application uses the real-space S.
        sk = model.Sk(self.k_points)
        overlap_eigenvalues = np.linalg.eigvalsh(sk)
        self.overlap_bounds = (float(overlap_eigenvalues.min()), float(overlap_eigenvalues.max()))
        if self.overlap_bounds[0] <= 0:
            raise ValueError("overlap is not positive definite")
        _, Y = exact_block_metric(sk)
        clean_values = np.linalg.eigvalsh(Y @ model.Hk(self.k_points) @ Y)
        self._clean_bounds = (float(clean_values.min()), float(clean_values.max()))
        s_mat = self._matrices["S"]
        if applier == "exact":
            X = _core.ExactInverseOverlap(s_mat, max_dim)
            Y = _core.ExactInverseSqrtOverlap(s_mat, max_dim)
        else:
            if order < 1 or not np.isfinite(coeff_tol) or coeff_tol <= 0:
                raise ValueError("Chebyshev metric needs positive order and coefficient tolerance")
            # Pad outward for CSR/Fourier last-bit differences while keeping
            # the certified strictly positive lower bound positive.
            lo, hi = self.overlap_bounds
            bounds = _core.SpectralBounds(lo*(1-1e-12), hi*(1+1e-12))
            X = _core.ChebyshevInverseOverlap(s_mat, bounds, order, coeff_tol)
            Y = _core.ChebyshevInverseSqrtOverlap(s_mat, bounds, order, coeff_tol)
        self.metric = SimpleNamespace(is_orthogonal=False, overlap=self._s_op,
                                      inverse_overlap=X, inverse_sqrt_overlap=Y)
        self._applied_h_op = (_core.LowdinOperator(Y, self._h_op) if route == "lowdin" else
                              _core.MetricInverseProduct(X, self._h_op))
        self.metric_residuals = self._metric_diagnostics()
        if not np.isfinite(tol) or tol <= 0:
            raise ValueError("metric probe tolerance must be positive and finite")
        if max(self.metric_residuals["inverse_probe"], self.metric_residuals["inverse_sqrt_probe"]) > tol:
            raise ValueError(f"metric probes exceed {tol:g}; increase Chebyshev order or reduce coeff_tol: "
                             f"{self.metric_residuals}")

    def _csr_operator(self, name, matrix):
        matrix = sparse.csr_matrix(matrix, dtype=np.complex128)
        self._owners[name] = matrix
        mat = _core.SparseMatrix.from_csr(self.dimension, np.asarray(matrix.indptr, np.int64),
                                         np.asarray(matrix.indices, np.int64),
                                         np.ascontiguousarray(matrix.data, dtype=np.complex128))
        self._matrices[name] = mat
        return _core.SparseMatrixOperator(mat)

    def _metric_diagnostics(self):
        # Deterministic probes diagnose truncation but do not certify an
        # operator norm. The tester must compare exact/Chebyshev moments.
        r = np.random.default_rng(12345).normal(size=self.dimension).astype(complex)
        r /= np.linalg.norm(r)
        X, Y = self.metric.inverse_overlap, self.metric.inverse_sqrt_overlap
        return dict(inverse_probe=float(np.linalg.norm(self._s_op.multiply(X.multiply(r))-r)),
                    inverse_sqrt_probe=float(np.linalg.norm(Y.multiply(self._s_op.multiply(Y.multiply(r)))-r)),
                    applier=self.applier, overlap_bounds=list(self.overlap_bounds),
                    inverse_effective_order=getattr(X, "effective_order", None),
                    inverse_sqrt_effective_order=getattr(Y, "effective_order", None))

    def hamiltonian(self):
        return self._h_op

    def overlap(self):
        return self._s_op

    def potential(self):
        return self._potential

    def rescaled_hamiltonian(self, window):
        """a X(H+V)+b I = X(a(H+V)+b S); no shift by b X."""
        return _core.RescaledOperator(self._applied_h_op, window)

    def operator(self, name):
        key = name.upper()
        if key in ("H", "HAM"):
            return self._h_op
        if key in ("1", "I", "IDENTITY", "S"):
            return self._s_op
        if key not in self._op_cache:
            if key in ("VX", "VY", "VZ"):
                axis = "XYZ".index(key[1])
                options = dict(derivative=axis, cell=self.velocity_mode == "bare")
                dh = supercell_csr(self.model, self.supercell, self.model.hoppings, **options)
                ds = supercell_csr(self.model, self.supercell, self.model.overlaps, **options)
                dS = self._csr_operator("dS"+key, ds)
                dH = self._csr_operator("dH"+key, dh)
                if self.has_disorder:
                    dH = _core.OperatorSum(dH, _triple(self._D, self.identity, dS))
                A = None
                if self.velocity_mode == "covariant":
                    A = self._csr_operator("A"+key, supercell_csr(
                        self.model, self.supercell, self.model.dipoles[axis]))
                self._op_cache[key] = _velocity(self._h_op, self.metric.inverse_overlap, dH, dS, A)
            elif key in ("SX", "SY", "SZ"):
                if self.model.physical_spin is not None:
                    raise ValueError("real-space transformed spin requires real-space spin blocks")
                q = applied_spin_operators(self.num_wann)[key]
                qR = sparse.kron(sparse.eye(self.num_cells), sparse.csr_matrix(q), format="csr")
                # Refuse a non spin-scalar overlap before using Sigma=S*q.
                sigma = checked_csr(self.overlap_csr @ qR, name=key)
                self._op_cache[key] = self._csr_operator(key, sigma)
            elif len(key) == 4 and key[0] == "J" and key[1] in "XYZ" and key[2] == "S" and key[3] in "XYZ":
                self._op_cache[key] = _triple(self.operator("V"+key[1]), self.metric.inverse_overlap,
                                             self.operator("S"+key[3]))
            else:
                raise ValueError(f"unknown stage6c operator {name!r}")
        return self._op_cache[key]

    def states(self, trace_mode="random", num_states=None, seed=None):
        """Map the SAME orthonormal orbital Bloch phases by B^+ to real space."""
        phases = _states(self, trace_mode, num_states, seed)
        return _core.TransformedStates(self.transform.adjoint(), phases)

    def bind_eigen_trace(self, plan, operators, trace_mode="random"):
        return _bind_metric_trace(self, plan, operators, trace_mode)

    def window(self, pad=0.10):
        return _window(self._clean_bounds, self.overlap_bounds, self.disorder_values, pad)
