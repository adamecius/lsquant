"""Spectral-projection basis adapters for Bloch, Eigen, and real-space layouts.

Operators are projected from the ORIGINAL physical velocity ladder, never
recomputed by differentiating H_P/H_s. The cached metric factors are reused
for hard/shift; the optional reduced tool uses identity overlap.
"""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np
from scipy.linalg import eigvalsh

from .. import _core
from ..bases import BlochMetricBasis, EigenMetricBasis
from .spectral import spectral_projection, resolve_projection, adjoint, physical_blocks
from .orthogonal import downfold_orthogonal, validate_bands
from .real_space import LocalMatrix


class ProjectedEigenBasis(EigenMetricBasis):
    """Eigen adapter whose explicit identity slot means the selected identity.

    An identity-free DOS word still traces ALL modes, including zero/edge
    complement states. This behaviour is recorded and never silently masked.
    """

    def operator(self, name):
        if name.upper() in ("1", "I", "IDENTITY", "WINDOW"):
            blocks = adjoint(self.eigenvectors) @ self.bloch.operator_blocks("1") @ self.eigenvectors
            return _core.BlockDiagonalOperator(physical_blocks(blocks, name="selected identity"))
        return super().operator(name)

    def states(self, trace_mode="random", num_states=None, seed=None):
        """Return random-phase or exact states in the projected eigenbasis."""
        if trace_mode == "random":
            return _core.RandomPhaseStates(self.dimension, 16 if num_states is None else int(num_states),
                                            self.bloch.seed if seed is None else int(seed))
        return super().states(trace_mode, num_states, seed)


class ProjectedBlochBasis(BlochMetricBasis):
    """Build a projected pencil once from an existing clean BlochMetricBasis.

    The source owns all unprojected physics. This adapter owns only the
    projected matrices and independent operator/engine caches. It never
    modifies the source or invalidates the original oracle.
    """

    def __init__(self, source, config):
        cfg = resolve_projection(config)
        if cfg["kind"] == "none":
            raise ValueError("use the source basis for projection none")
        if source.route == "anticommutator":
            raise ValueError("spectral projection supports the exact lowdin/generalised metric routes")
        self.model, self.supercell = source.model, source.supercell
        self.velocity_mode, self.route, self.seed = source.velocity_mode, source.route, source.seed
        self.k_points, self.num_cells = source.k_points, source.num_cells
        self.num_wann, self.dimension = source.num_wann, source.dimension
        self.sk_array = source.sk_array
        self.inverse_array, self.inverse_sqrt_array = source.inverse_array, source.inverse_sqrt_array
        self.source, self.config = source, cfg
        self._blocks, self._op_cache, self.local_matrices = {}, {}, {}
        self._eigen_basis = None
        self.projection = (downfold_orthogonal(source.hk_array, source.sk_array, cfg["window"])
                           if cfg["kind"] == "orthogonal" else
                           spectral_projection(source.hk_array, source.sk_array, cfg["window"], kind=cfg["kind"]))
        self.hk_array = self.projection.hamiltonian
        if cfg["kind"] == "orthogonal":
            self.num_wann = self.hk_array.shape[-1]
            self.dimension = self.num_cells*self.num_wann
            self.sk_array = np.broadcast_to(np.eye(self.num_wann), self.hk_array.shape).copy()
            self.inverse_array = self.sk_array.copy()
            self.inverse_sqrt_array = self.sk_array.copy()
            self.centres = source.model.tau[self.projection.trial_columns]
        else:
            self.centres = source.model.tau
        if cfg["radius"] is not None:
            local_h = self._local("H", self.hk_array)
            self.hk_array = local_h.blocks(self.k_points)
        self._h_op = _core.BlockDiagonalOperator(np.ascontiguousarray(self.hk_array))
        self._s_op = _core.BlockDiagonalOperator(np.ascontiguousarray(self.sk_array))
        self.metric = SimpleNamespace(is_orthogonal=False, overlap=self._s_op,
                                      inverse_overlap=_core.BlockDiagonalOperator(self.inverse_array),
                                      inverse_sqrt_overlap=_core.BlockDiagonalOperator(self.inverse_sqrt_array))
        Y = self.inverse_sqrt_array
        lowdin_h = physical_blocks(Y @ self.hk_array @ Y, name="projected Lowdin H")
        self.eigenvalues, u = np.linalg.eigh(lowdin_h)
        self.eigenvectors = np.ascontiguousarray(Y @ u)
        self.physical_eigenvalues, self.physical_eigenvectors = self.eigenvalues, self.eigenvectors
        applied_h = lowdin_h if self.route == "lowdin" else self.inverse_array @ self.hk_array
        self._applied_h_op = _core.BlockDiagonalOperator(np.ascontiguousarray(applied_h))
        self._lowdin_seed = _core.BlockDiagonalOperator(np.ascontiguousarray(u))
        self._coefficient_seed = _core.BlockDiagonalOperator(self.eigenvectors)
        self.metric_residuals = dict(source.metric_residuals)

    def _local(self, name, blocks):
        if name not in self.local_matrices:
            self.local_matrices[name] = LocalMatrix(blocks, self.model.uc, self.centres,
                                                    self.supercell, self.k_points,
                                                    radius=self.config["radius"])
        return self.local_matrices[name]

    def operator_blocks(self, name):
        key = name.upper()
        if key in ("H", "HAM"):
            return self.hk_array
        if key == "S":
            return self.sk_array
        if key not in self._blocks:
            original = (self.source.sk_array if key in ("1", "I", "IDENTITY", "WINDOW")
                        else self.source.operator_blocks(key))
            projected = self.projection.project_operator(original)
            self._blocks[key] = (self._local(key, projected).blocks(self.k_points) if self.config["radius"] is not None
                                  else projected)
        return self._blocks[key]

    def operator(self, name):
        # Base class special-cases identity as S. Here explicit identity must
        # be projected too; otherwise a user word silently restores discarded modes.
        key = name.upper()
        if key not in self._op_cache:
            self._op_cache[key] = _core.BlockDiagonalOperator(self.operator_blocks(key))
        return self._op_cache[key]

    def eigenbasis(self):
        if self._eigen_basis is None:
            self._eigen_basis = ProjectedEigenBasis(self)
        return self._eigen_basis

    def projection_receipt(self):
        """Report actual (possibly truncated) spectral bounds and local errors."""
        if self.config["kind"] == "orthogonal":
            counts = self.projection.kept.sum(axis=1)
            local_h = self._local("H", self.projection.hamiltonian)
            def reference(k):
                return eigvalsh(self.source.model.Hk(k), self.source.model.Sk(k))
            result = dict(kind="orthogonal", window=self.config["window"],
                          kept_bands_per_k=counts.tolist(),
                          minimum_trial_singular_value=self.projection.minimum_singular_value,
                          band_validation=validate_bands(local_h, reference, self.config["window"],
                                                         self.model.uc, self.supercell),
                          original_bounds=[float(self.source.eigenvalues.min()), float(self.source.eigenvalues.max())],
                          metric="reduced identity; optional orthogonal tool")
        else:
            result = self.projection.receipt()
            result["metric"] = "original S; same cached inverse and inverse square root"
        result["spectral_bounds"] = [float(self.eigenvalues.min()), float(self.eigenvalues.max())]
        result["local_matrices"] = {name: matrix.receipt() for name, matrix in self.local_matrices.items()}
        result["radius"] = self.config["radius"]
        result["route"] = self.route
        return result

    def local_matrix(self, name):
        """Lazily build local physical hoppings; identity stores S P.

        The corresponding local coefficient projector is S^-1 times the
        selected-identity matrix. After radius truncation it need not be
        idempotent: use the reported matrix/word errors, not exact identities.
        """
        key = name.upper()
        blocks = self.hk_array if key == "H" else self.operator_blocks(key)
        return self._local(key, blocks)

    def supercell_operators(self, names):
        """Local CSR export including the unchanged, untruncated overlap.

        A disorder term added later can couple the retained and discarded
        subspaces; clean projection does not automatically commute with disorder.
        """
        result = {"H": self.local_matrix("H").supercell_csr()}
        overlap = LocalMatrix(self.sk_array, self.model.uc, self.centres,
                               self.supercell, self.k_points)
        result["S"] = overlap.supercell_csr()
        for name in names:
            result[name] = self.local_matrix(name).supercell_csr()
        return result


class ProjectedRealSpaceBasis:
    """Apply local hard/shift matrices as CSR with the original overlap.

    LocalMatrix supplies H and every physical operator in the same cell
    gauge, including the requested hopping-radius cutoff. RealSpaceMetricBasis
    supplies S, its metric appliers and the B-mapped phase stream.
    No velocity is differentiated from the projected Hamiltonian.
    """

    name = "real"

    def __init__(self, projected, *, applier="exact", order=256,
                 coeff_tol=1e-14, max_dim=4096, tol=1e-11):
        from ..disorder_bases import RealSpaceMetricBasis
        if projected.config["kind"] == "orthogonal":
            raise ValueError("local nonorthogonal metric adapter requires hard or shift projection")
        self.projected = projected
        self._base = RealSpaceMetricBasis(
            projected.model, projected.supercell, velocity_mode=projected.velocity_mode,
            route=projected.route, seed=projected.seed, applier=applier, order=order,
            coeff_tol=coeff_tol, max_dim=max_dim, tol=tol)
        for name in ("model", "supercell", "num_cells", "num_wann", "dimension",
                     "route", "seed", "metric", "metric_residuals", "transform"):
            setattr(self, name, getattr(self._base, name))
        self._operators, self._owners = {}, {}
        self._h_op = self.operator("H")
        self._applied_h_op = (
            _core.LowdinOperator(self.metric.inverse_sqrt_overlap, self._h_op)
            if self.route == "lowdin" else
            _core.MetricInverseProduct(self.metric.inverse_overlap, self._h_op))

    def operator(self, name):
        key = name.upper()
        if key == "S":
            return self.overlap()
        if key == "HAM":
            key = "H"
        if key not in self._operators:
            from ..models.disorder import checked_csr
            csr = checked_csr(self.projected.local_matrix(key).supercell_csr(),
                              name="local projected "+key)
            matrix = _core.SparseMatrix.from_csr(
                self.dimension, np.asarray(csr.indptr, np.int64),
                np.asarray(csr.indices, np.int64), np.ascontiguousarray(csr.data))
            self._owners[key] = (csr, matrix)
            self._operators[key] = _core.SparseMatrixOperator(matrix)
        return self._operators[key]

    def hamiltonian(self):
        return self._h_op

    def overlap(self):
        return self._base.overlap()

    def rescaled_hamiltonian(self, window):
        return _core.RescaledOperator(self._applied_h_op, window)

    def window(self, pad=0.10):
        return self.projected.window(pad)

    def states(self, trace_mode="random", num_states=None, seed=None):
        return self._base.states(trace_mode, num_states, seed)

    def bind_eigen_trace(self, plan, operators, trace_mode="random"):
        return self._base.bind_eigen_trace(plan, operators, trace_mode)
