"""Basis adapters (RealBasis, BlochBasis, EigenBasis) for the Chebyshev engine.

All three bases expose the same interface over a :class:`~lsquant.models.wannier.WannierModel`
and supercell ``(N_1, N_2, N_3)``:

- ``.hamiltonian() -> LinearOperator``
- ``.operator(name) -> LinearOperator``
- ``.states(trace_mode, num_states=None, seed=12345) -> StateGenerator``
- ``.window(pad=0.10) -> SpectralWindow``

The C++ Chebyshev engine (:class:`~lsquant._core.MomentEngine`) receives only
:class:`~lsquant._core.LinearOperator` and :class:`~lsquant._core.StateGenerator`
objects and is completely agnostic of which basis adapter produced them.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import shlex
from typing import Sequence, Protocol
from types import SimpleNamespace
import numpy as np
import numpy.typing as npt

from . import _core
from .moments import window_from_gershgorin
from .spectral_functions import GaussianDelta

class TightBindingModel(Protocol):
    """Structural model interface; no file format is required by an adapter."""
    num_wann: int
    seed_name: str
    model_dir: str | None
    uc: np.ndarray
    tau: np.ndarray

    def Hk(self, k_points): ...
    def dHk(self, k_points, direction): ...
    def k_mesh(self, supercell): ...
    def bloch_phases(self, k_points): ...
    def hamiltonian_csr(self, supercell): ...
    def velocity_csr(self, supercell, direction): ...


_W2S_ALL_OPERATORS: tuple[str, ...] = (
    "VX",
    "VY",
    "VZ",
    "SX",
    "SY",
    "SZ",
    "VXSX",
    "VXSY",
    "VXSZ",
    "VYSX",
    "VYSY",
    "VYSZ",
    "VZSX",
    "VZSY",
    "VZSZ",
)


def _parse_w2s_command_line(
    cmd: str, out_path: str
) -> tuple[
    str | None,
    str | None,
    tuple[int, int, int] | None,
    set[str] | None,
    str | None,
]:
    """Parse ``invocation.command_line`` from a ``wannier2sparse`` ``.out`` receipt."""
    try:
        tokens = shlex.split(cmd)
    except ValueError as exc:
        raise ValueError(
            f"RealBasis: malformed command_line in receipt '{out_path}': {exc}"
        ) from exc
    if not tokens:
        raise ValueError(f"RealBasis: empty command_line in receipt '{out_path}'")

    one_arg_opts = {
        "-o",
        "--output-dir",
        "-p",
        "--project",
        "--seed",
        "--mode",
        "-x",
        "--input-file",
        "--create",
        "-inp",
        "--inp",
        "--win",
        "--qe-bands",
        "--log-level",
        "--log-file",
        "--velocity-mode",
        "--tb",
        "--r-dat",
    }
    zero_arg_opts = {
        "--write",
        "--verbose",
        "--quiet",
        "--no-log-file",
        "--bounds",
        "--exact-spin",
        "--covariant-velocity",
        "--orbital-l",
        "--orbital-L",
    }
    known_checks = {
        "all",
        "hermiticity",
        "sum_rules",
        "algebra",
        "aliasing",
        "bounds",
    }

    config_path: str | None = None
    explicit_seed: str | None = None
    want_all = False
    extra_ops: set[str] = set()
    positional: list[str] = []

    i = 1
    n = len(tokens)
    while i < n:
        tok = tokens[i]
        if tok in one_arg_opts:
            if i + 1 >= n:
                raise ValueError(
                    f"RealBasis: missing value after '{tok}' in command_line of '{out_path}'"
                )
            val = tokens[i + 1]
            if tok == "--seed":
                explicit_seed = val
            elif tok in ("-x", "--input-file"):
                config_path = val
            i += 2
        elif tok == "--op-file":
            if i + 2 >= n:
                raise ValueError(
                    f"RealBasis: missing arguments after '--op-file' in command_line of '{out_path}'"
                )
            extra_ops.add(tokens[i + 1].strip().upper())
            i += 3
        elif tok == "--spin-current":
            if i + 2 >= n:
                raise ValueError(
                    f"RealBasis: missing arguments after '--spin-current' in command_line of '{out_path}'"
                )
            extra_ops.add(
                f"J{tokens[i + 1].strip().upper()}S{tokens[i + 2].strip().upper()}"
            )
            i += 3
        elif tok == "--check":
            if i + 1 < n and tokens[i + 1] in known_checks:
                i += 2
            else:
                i += 1
        elif tok in zero_arg_opts:
            if tok == "--exact-spin":
                extra_ops.update(("SXEXACT", "SYEXACT", "SZEXACT"))
            elif tok in ("--orbital-l", "--orbital-L"):
                extra_ops.update(("LX", "LY", "LZ"))
            i += 1
        elif tok == "all":
            want_all = True
            i += 1
        elif len(tok) > 1 and tok.startswith("-") and not tok[1].isdigit():
            raise ValueError(
                f"RealBasis: unknown option '{tok}' in command_line of '{out_path}'"
            )
        else:
            positional.append(tok)
            i += 1

    if config_path is not None:
        cmd_label = positional[0] if positional else None
        return cmd_label, explicit_seed, None, None, config_path

    if len(positional) < 4:
        raise ValueError(
            f"RealBasis: command_line in receipt '{out_path}' is missing positional LABEL N1 N2 N3"
        )
    cmd_label = positional[0]
    try:
        n1, n2, n3 = int(positional[1]), int(positional[2]), int(positional[3])
    except ValueError as exc:
        raise ValueError(
            f"RealBasis: non-integer supercell dimensions in command_line of '{out_path}'"
        ) from exc
    if n1 < 1 or n2 < 1 or n3 < 1:
        raise ValueError(
            f"RealBasis: invalid supercell dimensions {(n1, n2, n3)} in command_line of '{out_path}'"
        )
    cmd_sc = (n1, n2, n3)
    pos_ops = (
        set(_W2S_ALL_OPERATORS)
        if want_all
        else {op.strip().upper() for op in positional[4:]}
    )
    cmd_ops = {"HAM"} | pos_ops | extra_ops
    return cmd_label, explicit_seed, cmd_sc, cmd_ops, None


class W2sCsrError(ValueError, FileNotFoundError):
    """Raised when a ``wannier2sparse`` CSR directory is missing or disagrees with the requested supercell/seed."""


def format_w2s_command(
    model_dir: str | None,
    seed_name: str,
    supercell: Sequence[int],
    csr_dir: str,
    operators: Sequence[str] = ("VX",),
) -> str:
    """Return the exact ``wannier2sparse`` command that produces ``csr_dir``."""
    n1, n2, n3 = int(supercell[0]), int(supercell[1]), int(supercell[2])
    ops = [str(op) for op in operators if str(op).strip().upper() not in ("", "1", "HAM")]
    ops_str = " ".join(ops) if ops else "VX"
    proj_dir = os.path.abspath(model_dir) if model_dir else "<model_dir>"
    out_dir = os.path.abspath(csr_dir)
    return f"wannier2sparse {seed_name} {n1} {n2} {n3} {ops_str} -p {proj_dir} -o {out_dir}"


def validate_w2s_csr_dir(
    csr_dir: str,
    seed_name: str,
    supercell: Sequence[int],
    *,
    model_dir: str | None = None,
    operators: Sequence[str] = ("VX",),
) -> set[str]:
    """Validate ``csr_dir`` for ``<seed>.HAM.CSR`` and ``<seed>.out`` before any load.

    Refuses with a message naming ``csr_dir``, the mismatch, and the exact
    ``wannier2sparse`` command that produces it.
    """
    sc_tup = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
    w2s_cmd = format_w2s_command(
        model_dir, seed_name, sc_tup, csr_dir, operators=operators
    )
    ham_path = os.path.join(csr_dir, f"{seed_name}.HAM.CSR")
    if not os.path.isfile(ham_path):
        raise W2sCsrError(
            f"RealBasis: csr_dir '{csr_dir}' lacks '{seed_name}.HAM.CSR' for supercell {sc_tup}; "
            f"generate it with: {w2s_cmd}"
        )

    out_path = os.path.join(csr_dir, f"{seed_name}.out")
    if not os.path.isfile(out_path):
        raise W2sCsrError(
            f"RealBasis: csr_dir '{csr_dir}' has '{seed_name}.HAM.CSR' but lacks receipt '{out_path}'; "
            f"generate it with: {w2s_cmd}"
        )
    try:
        with open(out_path, encoding="utf-8") as fh:
            meta = json.load(fh)
    except Exception as exc:
        raise W2sCsrError(
            f"RealBasis: failed to read wannier2sparse receipt '{out_path}' in csr_dir '{csr_dir}': {exc}; "
            f"generate it with: {w2s_cmd}"
        ) from exc
    if not isinstance(meta, dict):
        raise W2sCsrError(
            f"RealBasis: wannier2sparse receipt '{out_path}' in csr_dir '{csr_dir}' must be a JSON object; "
            f"generate it with: {w2s_cmd}"
        )

    if "status" in meta and meta["status"] != "success":
        raise W2sCsrError(
            f"RealBasis: wannier2sparse receipt '{out_path}' in csr_dir '{csr_dir}' has status {meta['status']!r} (expected 'success'); "
            f"generate it with: {w2s_cmd}"
        )

    inv = meta.get("invocation")
    if inv is not None and not isinstance(inv, dict):
        raise W2sCsrError(
            f"RealBasis: wannier2sparse receipt '{out_path}' in csr_dir '{csr_dir}' field 'invocation' must be an object; "
            f"generate it with: {w2s_cmd}"
        )
    inv_dict = inv if isinstance(inv, dict) else {}
    if "mode" in inv_dict and inv_dict["mode"] != "sparse":
        raise W2sCsrError(
            f"RealBasis: wannier2sparse receipt '{out_path}' in csr_dir '{csr_dir}' has mode {inv_dict['mode']!r} (expected 'sparse'); "
            f"generate it with: {w2s_cmd}"
        )

    cmd_label: str | None = None
    cmd_seed: str | None = None
    cmd_sc: tuple[int, int, int] | None = None
    cmd_ops: set[str] | None = None
    if "command_line" in inv_dict and inv_dict["command_line"] is not None:
        if (
            not isinstance(inv_dict["command_line"], str)
            or not inv_dict["command_line"].strip()
        ):
            raise W2sCsrError(
                f"RealBasis: wannier2sparse receipt '{out_path}' in csr_dir '{csr_dir}' has invalid 'command_line'; "
                f"generate it with: {w2s_cmd}"
            )
        cmd_label, cmd_seed, cmd_sc, cmd_ops, _ = _parse_w2s_command_line(
            inv_dict["command_line"], out_path
        )

    labels_found: list[str] = []
    if "label" in meta:
        if not isinstance(meta["label"], str) or not meta["label"].strip():
            raise W2sCsrError(
                f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has invalid 'label'; "
                f"generate it with: {w2s_cmd}"
            )
        labels_found.append(meta["label"].strip())
    if "label" in inv_dict:
        if not isinstance(inv_dict["label"], str) or not inv_dict["label"].strip():
            raise W2sCsrError(
                f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has invalid invocation 'label'; "
                f"generate it with: {w2s_cmd}"
            )
        labels_found.append(inv_dict["label"].strip())
    if cmd_label is not None:
        labels_found.append(cmd_label)
    if not labels_found:
        raise W2sCsrError(
            f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' does not specify a label; "
            f"generate it with: {w2s_cmd}"
        )
    for lbl in labels_found:
        if lbl != seed_name:
            raise W2sCsrError(
                f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' label {lbl!r} does not match model seed {seed_name!r}; "
                f"generate it with: {w2s_cmd}"
            )
    if cmd_seed is not None and cmd_seed != seed_name:
        raise W2sCsrError(
            f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' --seed {cmd_seed!r} does not match model seed {seed_name!r}; "
            f"generate it with: {w2s_cmd}"
        )

    supercells_found: list[tuple[int, int, int]] = []
    for sc_raw in (meta.get("supercell"), inv_dict.get("supercell")):
        if sc_raw is not None:
            if not isinstance(sc_raw, (list, tuple)) or len(sc_raw) != 3:
                raise W2sCsrError(
                    f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has invalid 'supercell' {sc_raw!r}; "
                    f"generate it with: {w2s_cmd}"
                )
            try:
                sc_item = (int(sc_raw[0]), int(sc_raw[1]), int(sc_raw[2]))
            except (TypeError, ValueError) as exc:
                raise W2sCsrError(
                    f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has non-integer 'supercell' {sc_raw!r}; "
                    f"generate it with: {w2s_cmd}"
                ) from exc
            supercells_found.append(sc_item)
    if cmd_sc is not None:
        supercells_found.append(cmd_sc)
    if not supercells_found:
        raise W2sCsrError(
            f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' does not specify a supercell; "
            f"generate it with: {w2s_cmd}"
        )
    for sc_item in supercells_found:
        if sc_item != sc_tup:
            raise W2sCsrError(
                f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' supercell {sc_item} does not match requested supercell {sc_tup}; "
                f"generate it with: {w2s_cmd}"
            )

    op_sets_found: list[set[str]] = []
    if "outputs" in meta:
        outputs_raw = meta["outputs"]
        if not isinstance(outputs_raw, list) or not outputs_raw:
            raise W2sCsrError(
                f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' field 'outputs' must be a non-empty list; "
                f"generate it with: {w2s_cmd}"
            )
        out_ops: set[str] = set()
        for entry in outputs_raw:
            if isinstance(entry, dict):
                op_nm = entry.get("name")
                if not isinstance(op_nm, str) or not op_nm.strip():
                    raise W2sCsrError(
                        f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has an output entry without a valid 'name'; "
                        f"generate it with: {w2s_cmd}"
                    )
                op_key = op_nm.strip().upper()
                op_path = entry.get("path")
                if isinstance(op_path, str) and op_path.strip():
                    expected_base = f"{seed_name}.{op_key}.CSR"
                    if os.path.basename(op_path.strip()) != expected_base:
                        raise W2sCsrError(
                            f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' output '{op_key}' path '{op_path}' does not match '{expected_base}'; "
                            f"generate it with: {w2s_cmd}"
                        )
                out_ops.add(op_key)
            elif isinstance(entry, str) and entry.strip():
                out_ops.add(entry.strip().upper())
            else:
                raise W2sCsrError(
                    f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has invalid entry in 'outputs'; "
                    f"generate it with: {w2s_cmd}"
                )
        op_sets_found.append(out_ops)

    if "operators" in meta:
        ops_raw = meta["operators"]
        if not isinstance(ops_raw, list) or not ops_raw:
            raise W2sCsrError(
                f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' field 'operators' must be a non-empty list; "
                f"generate it with: {w2s_cmd}"
            )
        op_list_set: set[str] = set()
        for item in ops_raw:
            if not isinstance(item, str) or not item.strip():
                raise W2sCsrError(
                    f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has invalid operator name in 'operators'; "
                    f"generate it with: {w2s_cmd}"
                )
            op_list_set.add(item.strip().upper())
        op_sets_found.append(op_list_set)

    if cmd_ops is not None:
        op_sets_found.append(cmd_ops)

    if not op_sets_found:
        raise W2sCsrError(
            f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' does not specify an operator list; "
            f"generate it with: {w2s_cmd}"
        )

    canonical_ops = op_sets_found[0]
    if "HAM" not in canonical_ops:
        raise W2sCsrError(
            f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' operator list {sorted(canonical_ops)} does not include 'HAM'; "
            f"generate it with: {w2s_cmd}"
        )
    for other_ops in op_sets_found[1:]:
        if other_ops != canonical_ops:
            raise W2sCsrError(
                f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' has inconsistent operator lists: {sorted(canonical_ops)} vs {sorted(other_ops)}; "
                f"generate it with: {w2s_cmd}"
            )

    prefix = f"{seed_name}."
    suffix = ".CSR"
    disk_ops: set[str] = set()
    for fname in os.listdir(csr_dir):
        if fname.startswith(prefix) and fname.endswith(suffix):
            op_part = fname[len(prefix) : -len(suffix)].strip().upper()
            if op_part:
                disk_ops.add(op_part)
    if disk_ops != canonical_ops:
        raise W2sCsrError(
            f"RealBasis: receipt '{out_path}' in csr_dir '{csr_dir}' operator list {sorted(canonical_ops)} does not match CSR files on disk {sorted(disk_ops)}; "
            f"generate it with: {w2s_cmd}"
        )
    return canonical_ops


def make_anderson_disorder(
    dimension: int, disorder_w: float, seed: int = 12345
) -> npt.NDArray[np.complex128]:
    """Generate deterministic real-space box disorder in ``[-W/2, +W/2]``.

    Parameters
    ----------
    dimension : int
        Total Hilbert-space dimension ``N_k * W``.
    disorder_w : float
        Box disorder width ``W >= 0`` in eV.
    seed : int, optional
        Deterministic PRNG seed (default ``12345``).

    Returns
    -------
    numpy.ndarray
        1D ``complex128`` array of length ``dimension`` with real entries drawn
        uniformly from ``[-W/2, +W/2]``.
    """
    if disorder_w < 0.0:
        raise ValueError("make_anderson_disorder: disorder_w must be non-negative")
    if disorder_w == 0.0:
        return np.zeros(dimension, dtype=np.complex128)
    rng = np.random.default_rng(int(seed))
    vals = rng.uniform(-0.5 * float(disorder_w), 0.5 * float(disorder_w), size=dimension)
    return np.ascontiguousarray(vals, dtype=np.complex128)


class RealBasis:
    """Real-space supercell basis adapter.

    Reads ``<seed>.HAM.CSR`` and ``<seed>.VX.CSR`` (and ``VY``, ``VZ``) written
    by ``wannier2sparse`` when they exist in ``csr_dir`` or ``model.model_dir``
    for the requested supercell ``N_1 * N_2 * N_3 * W``, and builds them
    identically in the ``wannier2sparse`` cell numbering
    ``iR = n_1 + N_1 * n_2 + N_1 * N_2 * n_3`` otherwise. When ``disorder_w > 0``
    (or an explicit ``disorder`` vector is given), ``H`` is the matrix-free
    composition ``OperatorSum(SparseMatrixOperator(H_csr), DiagonalOperator(v_real))``.
    """

    def __init__(
        self,
        model: TightBindingModel,
        supercell: Sequence[int],
        *,
        disorder_w: float = 0.0,
        seed: int = 12345,
        disorder_seed: int | None = None,
        disorder: npt.NDArray[np.complex128] | npt.NDArray[np.float64] | None = None,
        csr_dir: str | None = None,
        use_w2s_csr: bool = True,
    ) -> None:
        if len(supercell) != 3:
            raise ValueError("RealBasis: supercell must have 3 dimensions (N1, N2, N3)")
        self.model = model
        self.supercell = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
        self.num_cells = self.supercell[0] * self.supercell[1] * self.supercell[2]
        self.num_wann = model.num_wann
        self.dimension = self.num_cells * self.num_wann
        self.disorder_w = float(disorder_w)
        self.seed = int(seed)
        self.disorder_seed = self.seed if disorder_seed is None else int(disorder_seed)
        self._csr_dir = csr_dir
        self._use_w2s_csr = bool(use_w2s_csr)

        self._op_sources: dict[str, str] = {}
        self._w2s_active_dir: str | None = None
        self._w2s_receipt_ops: set[str] | None = None

        if self._use_w2s_csr:
            seed_name = self.model.seed_name
            if self._csr_dir is not None:
                self._w2s_active_dir = self._csr_dir
            elif self.model.model_dir is not None:
                ham_candidate = os.path.join(
                    self.model.model_dir, f"{seed_name}.HAM.CSR"
                )
                if os.path.isfile(ham_candidate):
                    self._w2s_active_dir = self.model.model_dir

        if self._w2s_active_dir is not None:
            self._w2s_receipt_ops = self._validate_w2s_receipt(self._w2s_active_dir)
            self._h_csr, ham_src = self._load_w2s_csr("HAM")
            self._op_sources["HAM"] = ham_src
        else:
            self._h_csr = model.hamiltonian_csr(self.supercell)
            self._op_sources["HAM"] = "built by lsquant"
        self._h_op_clean = _core.SparseMatrixOperator(self._h_csr)

        if disorder is not None:
            v_arr = np.ascontiguousarray(disorder, dtype=np.complex128).ravel()
            if v_arr.size != self.dimension:
                raise ValueError("RealBasis: disorder array length must equal dimension")
            self._disorder_vec: npt.NDArray[np.complex128] | None = v_arr
            self._v_op: _core.DiagonalOperator | None = _core.DiagonalOperator(v_arr)
            self._h_op: _core.LinearOperator = _core.OperatorSum(
                self._h_op_clean, self._v_op
            )
        elif self.disorder_w > 0.0:
            v_arr = make_anderson_disorder(
                self.dimension, self.disorder_w, self.disorder_seed
            )
            self._disorder_vec = v_arr
            self._v_op = _core.DiagonalOperator(v_arr)
            self._h_op = _core.OperatorSum(self._h_op_clean, self._v_op)
        else:
            self._disorder_vec = None
            self._v_op = None
            self._h_op = self._h_op_clean

        self._vel_csr_cache: dict[str, _core.SparseMatrix] = {}
        self._op_cache: dict[str, _core.LinearOperator] = {}

    def _validate_w2s_receipt(self, d: str) -> set[str]:
        """Strictly validate ``<seed>.HAM.CSR`` and ``<seed>.out`` in ``d``."""
        return validate_w2s_csr_dir(
            d,
            self.model.seed_name,
            self.supercell,
            model_dir=self.model.model_dir,
        )

    def _load_w2s_csr(self, op_name: str) -> tuple[_core.SparseMatrix, str]:
        """Load ``<seed>.<op_name>.CSR`` from the validated ``wannier2sparse`` directory and return ``(mat, source_str)``."""
        assert self._w2s_active_dir is not None
        assert self._w2s_receipt_ops is not None
        key = op_name.strip().upper()
        seed = self.model.seed_name
        out_path = os.path.join(self._w2s_active_dir, f"{seed}.out")
        if key not in self._w2s_receipt_ops:
            raise ValueError(
                f"RealBasis: operator '{key}' is not in wannier2sparse receipt '{out_path}' (operators: {sorted(self._w2s_receipt_ops)})"
            )
        csr_path = os.path.join(self._w2s_active_dir, f"{seed}.{key}.CSR")
        if not os.path.isfile(csr_path):
            raise FileNotFoundError(
                f"RealBasis: CSR file '{csr_path}' listed in receipt '{out_path}' is missing"
            )
        h = hashlib.sha256()
        with open(csr_path, "rb") as fh:
            while chunk := fh.read(1024 * 1024):
                h.update(chunk)
        sha256_hex = h.hexdigest()
        mat = _core.ReadCsrText(csr_path)
        if int(mat.dimension) != self.dimension:
            raise ValueError(
                f"RealBasis: CSR '{csr_path}' dimension {mat.dimension} does not match expected supercell dimension {self.dimension}"
            )
        return mat, f"wannier2sparse {csr_path} {sha256_hex}"

    @property
    def hamiltonian_csr(self) -> _core.SparseMatrix:
        """Return the clean real-space supercell Hamiltonian :class:`~lsquant._core.SparseMatrix`."""
        return self._h_csr

    @property
    def disorder_values(self) -> npt.NDArray[np.complex128] | None:
        """Return the real-space diagonal disorder array if present."""
        return self._disorder_vec

    @property
    def operator_sources(self) -> dict[str, str]:
        """Return a copy of the provenance dictionary for all materialised operators."""
        return dict(self._op_sources)

    def operator_source(self, name: str) -> str:
        """Return the provenance string (``wannier2sparse <path> <sha256>`` or ``built by lsquant``) for ``name``."""
        key = name.strip().upper()
        if key in ("H", "HAM"):
            return self._op_sources["HAM"]
        if key in ("1", "I", "ID", "IDENTITY"):
            return "built by lsquant"
        if key in ("VX", "VY", "VZ"):
            if key not in self._op_sources:
                self.velocity_csr(key)
            return self._op_sources[key]
        raise ValueError(f"RealBasis.operator_source: unknown operator '{name}'")

    def hamiltonian(self) -> _core.LinearOperator:
        """Return the Hamiltonian :class:`~lsquant._core.LinearOperator` in real space."""
        return self._h_op

    def velocity_csr(self, name: str) -> _core.SparseMatrix:
        """Return the real-space supercell velocity :class:`~lsquant._core.SparseMatrix`."""
        key = name.strip().upper()
        if key not in ("VX", "VY", "VZ"):
            raise ValueError(f"RealBasis.velocity_csr: unsupported velocity '{name}'")
        if key not in self._vel_csr_cache:
            if self._w2s_active_dir is not None:
                mat, src = self._load_w2s_csr(key)
                self._vel_csr_cache[key] = mat
                self._op_sources[key] = src
            else:
                self._vel_csr_cache[key] = self.model.velocity_csr(self.supercell, key)
                self._op_sources[key] = "built by lsquant"
        return self._vel_csr_cache[key]

    def operator(self, name: str) -> _core.LinearOperator:
        """Return the named static operator as a :class:`~lsquant._core.LinearOperator`."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY"):
            if "1" not in self._op_cache:
                self._op_cache["1"] = _core.IdentityOperator(self.dimension)
            return self._op_cache["1"]
        if key in ("H", "HAM"):
            return self._h_op
        if key in ("VX", "VY", "VZ"):
            if key not in self._op_cache:
                mat = self.velocity_csr(key)
                self._op_cache[key] = _core.SparseMatrixOperator(mat)
            return self._op_cache[key]
        raise ValueError(f"RealBasis.operator: unknown operator '{name}'")

    def states(
        self,
        trace_mode: str = "exact",
        num_states: int | None = None,
        seed: int | None = None,
    ) -> _core.StateGenerator:
        """Return a :class:`~lsquant._core.StateGenerator` for ``trace_mode``."""
        mode = trace_mode.strip().lower()
        rng_seed = self.seed if seed is None else int(seed)
        if mode == "exact":
            n_st = self.dimension if num_states is None else int(num_states)
            if n_st <= 0 or n_st > self.dimension:
                raise ValueError(
                    "RealBasis.states: num_states must lie in 1..dimension for exact trace"
                )
            return _core.LocalStates(self.dimension, list(range(n_st)))
        if mode == "random":
            n_st = 1 if num_states is None else int(num_states)
            if n_st <= 0:
                raise ValueError("RealBasis.states: num_states must be positive")
            return _core.RandomPhaseStates(self.dimension, n_st, rng_seed)
        raise ValueError(
            f"RealBasis.states: unsupported trace_mode '{trace_mode}' (expected 'exact' or 'random')"
        )

    def window(self, pad: float = 0.10) -> _core.SpectralWindow:
        """Return the Gershgorin :class:`~lsquant._core.SpectralWindow` for this Hamiltonian.

        Kept on Gershgorin as the real-space default; on Wannier models with
        long-range hoppings, Gershgorin bounds are several times wider than the
        true spectrum and waste Chebyshev moments, so :class:`BlochBasis` and
        :class:`EigenBasis` use the block eigenvalue range instead.
        """
        win = window_from_gershgorin(self._h_csr, pad=float(pad))
        if self.disorder_w > 0.0:
            extra = 0.5 * self.disorder_w * (1.0 + float(pad))
            return _core.SpectralWindow(win.lower - extra, win.upper + extra)
        return win


class BlochBasis:
    """Bloch k-space basis adapter with optional FFT-sandwiched real-space disorder.

    Builds ``H_k = BlockDiagonalOperator(H(k))`` and ``v_a(k) = BlockDiagonalOperator(dH_a(k))``
    over the ``N_k = N_1 * N_2 * N_3`` supercell k-points, and owns the batched
    FFTW :class:`~lsquant._core.BlochTransform` ``B`` (`real -> Bloch`). When
    ``disorder_w > 0`` (or ``disorder`` is supplied), ``H`` is the composite
    ``OperatorSum(H_k, SandwichOperator(B, DiagonalOperator(v_real)))``.
    """

    def __init__(
        self,
        model: TightBindingModel,
        supercell: Sequence[int],
        *,
        disorder_w: float = 0.0,
        seed: int = 12345,
        disorder_seed: int | None = None,
        disorder: npt.NDArray[np.complex128] | npt.NDArray[np.float64] | None = None,
    ) -> None:
        if len(supercell) != 3:
            raise ValueError("BlochBasis: supercell must have 3 dimensions (N1, N2, N3)")
        self.model = model
        self.supercell = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
        self.num_cells = self.supercell[0] * self.supercell[1] * self.supercell[2]
        self.num_wann = model.num_wann
        self.dimension = self.num_cells * self.num_wann
        self.disorder_w = float(disorder_w)
        self.seed = int(seed)
        self.disorder_seed = self.seed if disorder_seed is None else int(disorder_seed)

        self.k_points = model.k_mesh(self.supercell)
        self.phases = model.bloch_phases(self.k_points)
        self.transform = _core.BlochTransform(
            self.supercell[0],
            self.supercell[1],
            self.supercell[2],
            self.num_wann,
            self.phases.ravel(),
            _core.BlochDirection.RealToBloch,
        )

        self.hk_array = model.Hk(self.k_points)
        self._hk_op = _core.BlockDiagonalOperator(self.hk_array)

        if disorder is not None:
            v_arr = np.ascontiguousarray(disorder, dtype=np.complex128).ravel()
            if v_arr.size != self.dimension:
                raise ValueError("BlochBasis: disorder array length must equal dimension")
            self._disorder_vec: npt.NDArray[np.complex128] | None = v_arr
            self._v_op: _core.DiagonalOperator | None = _core.DiagonalOperator(v_arr)
            self._sandwich_v: _core.SandwichOperator | None = _core.SandwichOperator(
                self.transform, self._v_op
            )
            self._h_op: _core.LinearOperator = _core.OperatorSum(
                self._hk_op, self._sandwich_v
            )
        elif self.disorder_w > 0.0:
            v_arr = make_anderson_disorder(
                self.dimension, self.disorder_w, self.disorder_seed
            )
            self._disorder_vec = v_arr
            self._v_op = _core.DiagonalOperator(v_arr)
            self._sandwich_v = _core.SandwichOperator(self.transform, self._v_op)
            self._h_op = _core.OperatorSum(self._hk_op, self._sandwich_v)
        else:
            self._disorder_vec = None
            self._v_op = None
            self._sandwich_v = None
            self._h_op = self._hk_op

        self._op_cache: dict[str, _core.LinearOperator] = {}
        self._eig_bounds: tuple[float, float] | None = None

    @property
    def disorder_values(self) -> npt.NDArray[np.complex128] | None:
        """Return the real-space diagonal disorder array if present."""
        return self._disorder_vec

    @property
    def operator_sources(self) -> dict[str, str]:
        """Return the provenance dictionary for operators on this Bloch basis."""
        out = {"HAM": "built by lsquant"}
        for key in self._op_cache:
            if key in ("VX", "VY", "VZ"):
                out[key] = "built by lsquant"
        return out

    def operator_source(self, name: str) -> str:
        """Return the provenance string for operator ``name``."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY", "H", "HAM", "VX", "VY", "VZ"):
            return "built by lsquant"
        raise ValueError(f"BlochBasis.operator_source: unknown operator '{name}'")

    def hamiltonian(self) -> _core.LinearOperator:
        """Return the Hamiltonian :class:`~lsquant._core.LinearOperator` in Bloch space."""
        return self._h_op

    def operator(self, name: str) -> _core.LinearOperator:
        """Return the named static operator in Bloch space."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY"):
            if "1" not in self._op_cache:
                self._op_cache["1"] = _core.IdentityOperator(self.dimension)
            return self._op_cache["1"]
        if key in ("H", "HAM"):
            return self._h_op
        if key in ("VX", "VY", "VZ"):
            if key not in self._op_cache:
                vel_blocks = self.model.dHk(self.k_points, key)
                self._op_cache[key] = _core.BlockDiagonalOperator(vel_blocks)
            return self._op_cache[key]
        raise ValueError(f"BlochBasis.operator: unknown operator '{name}'")

    def states(
        self,
        trace_mode: str = "exact",
        num_states: int | None = None,
        seed: int | None = None,
    ) -> _core.StateGenerator:
        """Return a :class:`~lsquant._core.StateGenerator` for ``trace_mode``.

        For ``trace_mode == "exact"`` at ``disorder_w == 0``, returns
        :class:`~lsquant._core.BlockLocalStates` ``(N_k, W)``, which evaluates
        the exact full-mesh trace of any block-diagonal word in ``W`` states.
        For ``trace_mode == "random"``, applies ``B`` to
        :class:`~lsquant._core.RandomPhaseStates` via
        :class:`~lsquant._core.TransformedStates` so that real-space and
        Bloch-space stochastic traces use the exact same physical state vectors.
        """
        mode = trace_mode.strip().lower()
        rng_seed = self.seed if seed is None else int(seed)
        if mode == "exact":
            if self._disorder_vec is None and (
                num_states is None or int(num_states) == self.num_wann
            ):
                return _core.BlockLocalStates(self.num_cells, self.num_wann)
            n_st = self.dimension if num_states is None else int(num_states)
            if n_st <= 0 or n_st > self.dimension:
                raise ValueError(
                    "BlochBasis.states: num_states must lie in 1..dimension for full exact trace"
                )
            return _core.LocalStates(self.dimension, list(range(n_st)))
        if mode == "random":
            n_st = 1 if num_states is None else int(num_states)
            if n_st <= 0:
                raise ValueError("BlochBasis.states: num_states must be positive")
            inner = _core.RandomPhaseStates(self.dimension, n_st, rng_seed)
            return _core.TransformedStates(self.transform, inner)
        raise ValueError(
            f"BlochBasis.states: unsupported trace_mode '{trace_mode}' (expected 'exact' or 'random')"
        )

    def window(self, pad: float = 0.10) -> _core.SpectralWindow:
        """Return the :class:`~lsquant._core.SpectralWindow` enclosing the block eigenvalues.

        Diagonalises the ``N_k`` blocks ``H(k)`` once and caches ``(emin, emax)``,
        returning ``[emin - pad * width, emax + pad * width]`` with
        ``width = emax - emin`` identical to :meth:`EigenBasis.window`.
        On Wannier models with long-range hoppings, Gershgorin bounds are
        several times wider than the true spectrum and waste Chebyshev moments,
        so ``BlochBasis`` uses the block eigenvalues rather than Gershgorin.
        """
        if pad < 0.0:
            raise ValueError("BlochBasis.window: pad must be non-negative")
        if self._eig_bounds is None:
            eigvals, _ = np.linalg.eigh(self.hk_array)
            self._eig_bounds = (float(np.min(eigvals)), float(np.max(eigvals)))
        emin, emax = self._eig_bounds
        if self.disorder_w > 0.0:
            emin -= 0.5 * self.disorder_w
            emax += 0.5 * self.disorder_w
        width = emax - emin
        if width <= 0.0:
            emin -= 1.0
            emax += 1.0
            width = 2.0
        margin = float(pad) * width
        return _core.SpectralWindow(emin - margin, emax + margin)


class EigenBasis:
    """Per-k exact eigenbasis adapter (clean systems, ``W = 0``).

    Diagonalises each ``W x W`` block ``H(k) = U(k) diag(eps(k)) U(k)^+`` via
    :func:`numpy.linalg.eigh`, builds ``H = DiagonalOperator(eps)`` and rotated
    operators ``O_eig(k) = U(k)^+ O(k) U(k)`` as
    :class:`~lsquant._core.BlockDiagonalOperator`, and uses
    :class:`~lsquant._core.BlockLocalStates` ``(N_k, W)`` so the standard
    :class:`~lsquant._core.MomentEngine` evaluates exact full-mesh moments in
    ``W`` state evaluations.
    """

    def __init__(
        self,
        model: TightBindingModel,
        supercell: Sequence[int],
        *,
        disorder_w: float = 0.0,
        seed: int = 12345,
    ) -> None:
        if float(disorder_w) != 0.0:
            raise ValueError(
                "EigenBasis: non-zero disorder_w is not supported in the per-k eigenbasis"
            )
        if len(supercell) != 3:
            raise ValueError("EigenBasis: supercell must have 3 dimensions (N1, N2, N3)")
        self.model = model
        self.supercell = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
        self.num_cells = self.supercell[0] * self.supercell[1] * self.supercell[2]
        self.num_wann = model.num_wann
        self.dimension = self.num_cells * self.num_wann
        self.disorder_w = 0.0
        self.seed = int(seed)

        self.k_points = model.k_mesh(self.supercell)
        self.hk_array = model.Hk(self.k_points)

        # Batched eigh over all N_k blocks of shape (N_k, W, W)
        eigvals, eigvecs = np.linalg.eigh(self.hk_array)
        self.eigenvalues: npt.NDArray[np.float64] = np.ascontiguousarray(
            eigvals, dtype=np.float64
        )
        self.eigenvectors: npt.NDArray[np.complex128] = np.ascontiguousarray(
            eigvecs, dtype=np.complex128
        )
        self.eigenvectors_dagger: npt.NDArray[np.complex128] = np.ascontiguousarray(
            np.conj(np.swapaxes(self.eigenvectors, -1, -2)), dtype=np.complex128
        )

        self._h_op = _core.DiagonalOperator(
            np.ascontiguousarray(self.eigenvalues.ravel(), dtype=np.complex128)
        )
        self.u_op = _core.BlockDiagonalOperator(self.eigenvectors)
        self.u_dagger_op = _core.BlockDiagonalOperator(self.eigenvectors_dagger)
        self._op_cache: dict[str, _core.LinearOperator] = {}

    @property
    def operator_sources(self) -> dict[str, str]:
        """Return the provenance dictionary for operators on this eigenbasis."""
        out = {"HAM": "built by lsquant"}
        for key in self._op_cache:
            if key in ("VX", "VY", "VZ"):
                out[key] = "built by lsquant"
        return out

    def operator_source(self, name: str) -> str:
        """Return the provenance string for operator ``name``."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY", "H", "HAM", "VX", "VY", "VZ"):
            return "built by lsquant"
        raise ValueError(f"EigenBasis.operator_source: unknown operator '{name}'")

    def hamiltonian(self) -> _core.LinearOperator:
        """Return ``DiagonalOperator(eps)``."""
        return self._h_op

    def rotate_blocks(
        self, bloch_blocks: npt.NDArray[np.complex128]
    ) -> npt.NDArray[np.complex128]:
        """Rotate ``(N_k, W, W)`` Bloch blocks into the eigenbasis: ``U(k)^+ O(k) U(k)``."""
        rotated = self.eigenvectors_dagger @ bloch_blocks @ self.eigenvectors
        return np.ascontiguousarray(rotated, dtype=np.complex128)

    def operator(self, name: str) -> _core.LinearOperator:
        """Return the named static operator rotated into the eigenbasis."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY"):
            if "1" not in self._op_cache:
                self._op_cache["1"] = _core.IdentityOperator(self.dimension)
            return self._op_cache["1"]
        if key in ("H", "HAM"):
            return self._h_op
        if key in ("VX", "VY", "VZ"):
            if key not in self._op_cache:
                vel_bloch = self.model.dHk(self.k_points, key)
                vel_eig = self.rotate_blocks(vel_bloch)
                self._op_cache[key] = _core.BlockDiagonalOperator(vel_eig)
            return self._op_cache[key]
        raise ValueError(f"EigenBasis.operator: unknown operator '{name}'")

    def states(
        self,
        trace_mode: str = "exact",
        num_states: int | None = None,
        seed: int | None = None,
    ) -> _core.StateGenerator:
        """Return a :class:`~lsquant._core.StateGenerator` for ``trace_mode``.

        Defaults to :class:`~lsquant._core.BlockLocalStates` ``(N_k, W)``, which
        evaluates the exact eigenbasis trace in ``W`` states.
        """
        mode = trace_mode.strip().lower()
        rng_seed = self.seed if seed is None else int(seed)
        if mode == "exact":
            if num_states is None or int(num_states) == self.num_wann:
                return _core.BlockLocalStates(self.num_cells, self.num_wann)
            n_st = int(num_states)
            if n_st <= 0 or n_st > self.dimension:
                raise ValueError(
                    "EigenBasis.states: num_states must lie in 1..dimension"
                )
            return _core.LocalStates(self.dimension, list(range(n_st)))
        if mode == "random":
            n_st = 1 if num_states is None else int(num_states)
            if n_st <= 0:
                raise ValueError("EigenBasis.states: num_states must be positive")
            phases = self.model.bloch_phases(self.k_points)
            b_op = _core.BlochTransform(
                self.supercell[0],
                self.supercell[1],
                self.supercell[2],
                self.num_wann,
                phases.ravel(),
                _core.BlochDirection.RealToBloch,
            )
            u_dag_b = _core.OperatorProduct(self.u_dagger_op, b_op)
            inner = _core.RandomPhaseStates(self.dimension, n_st, rng_seed)
            return _core.TransformedStates(u_dag_b, inner)
        raise ValueError(
            f"EigenBasis.states: unsupported trace_mode '{trace_mode}' (expected 'exact' or 'random')"
        )

    def window(self, pad: float = 0.10) -> _core.SpectralWindow:
        """Return the :class:`~lsquant._core.SpectralWindow` enclosing the block eigenvalues.

        Computes ``[emin - pad * width, emax + pad * width]`` from the minimum and
        maximum eigenvalues ``(emin, emax)`` of all ``N_k`` blocks ``H(k)`` with
        ``width = emax - emin``. On Wannier models with long-range hoppings,
        Gershgorin bounds are several times wider than the true spectrum and
        waste Chebyshev moments, so ``EigenBasis`` and :class:`BlochBasis`
        take their window from the block eigenvalues rather than from Gershgorin.
        """
        if pad < 0.0:
            raise ValueError("EigenBasis.window: pad must be non-negative")
        emin = float(np.min(self.eigenvalues))
        emax = float(np.max(self.eigenvalues))
        width = emax - emin
        if width <= 0.0:
            emin -= 1.0
            emax += 1.0
            width = 2.0
        margin = float(pad) * width
        return _core.SpectralWindow(emin - margin, emax + margin)

    def exact_dos(
        self,
        energies: Sequence[float] | npt.NDArray[np.float64],
        eta: float | Sequence[float] | npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate the exact per-cell Gaussian density of states on the block eigenvalues.

        Computes ``(1 / N_k) sum_{k, n} g(E - eps_n(k); eta)`` where the
        Gaussian ``g(x; eta) = exp(-x^2 / (2 * eta^2)) / (sqrt(2 * pi) * eta)``
        is evaluated via the SymPy
        :data:`~lsquant.spectral_functions.GaussianDelta`
        :class:`~lsquant.spectral_functions.SpectralFunction`.

        Parameters
        ----------
        energies : Sequence[float] or numpy.ndarray
            1D energy grid in eV.
        eta : float or numpy.ndarray
            Gaussian broadening width in eV (``eta > 0``).

        Returns
        -------
        numpy.ndarray
            1D ``float64`` array of length ``len(energies)`` giving the density
            of states per unit cell in ``1 / eV``.
        """
        e_arr = np.asarray(energies, dtype=np.float64).ravel()
        eta_arr = np.asarray(eta, dtype=np.float64)
        if np.any(eta_arr <= 0.0):
            raise ValueError("EigenBasis.exact_dos: eta must be positive")
        if eta_arr.ndim == 1:
            eta_arr = eta_arr[:, None, None]
        diff = e_arr[:, None, None] - self.eigenvalues[None, :, :]  # (N_E, N_k, W)
        g_vals = np.asarray(GaussianDelta.evaluate(diff, eta_arr), dtype=np.float64)
        dos_per_cell = np.sum(g_vals, axis=(1, 2)) / float(self.num_cells)
        return np.ascontiguousarray(dos_per_cell, dtype=np.float64)

    def exact_sigma_xx(
        self,
        energies: Sequence[float] | npt.NDArray[np.float64],
        eta: float | Sequence[float] | npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate the exact Gaussian Kubo-Greenwood conductivity ``sigma_xx(E)``.

        Computes ``pi * sum_{k, n, m} |v_{x, nm}(k)|^2 g(E - eps_n(k); eta) g(E - eps_m(k); eta)``
        with the same normalisation as :func:`~lsquant.reconstruct.kubo_greenwood`
        so that the two curves are directly comparable (multiply by
        ``2 * pi / (A_cell * N_k)`` to convert to ``hr_exactdiag.py``'s
        ``e^2 / h`` per-area convention).

        Parameters
        ----------
        energies : Sequence[float] or numpy.ndarray
            1D energy grid in eV.
        eta : float or numpy.ndarray
            Gaussian broadening width in eV (``eta > 0``).

        Returns
        -------
        numpy.ndarray
            1D ``float64`` array of length ``len(energies)``.
        """
        e_arr = np.asarray(energies, dtype=np.float64).ravel()
        eta_arr = np.asarray(eta, dtype=np.float64)
        if np.any(eta_arr <= 0.0):
            raise ValueError("EigenBasis.exact_sigma_xx: eta must be positive")
        if eta_arr.ndim == 1:
            eta_arr = eta_arr[:, None, None]
        vx_eig = self.rotate_blocks(self.model.dHk(self.k_points, "VX"))
        vx2 = np.abs(vx_eig) ** 2  # (N_k, W, W)
        diff = e_arr[:, None, None] - self.eigenvalues[None, :, :]  # (N_E, N_k, W)
        g_vals = np.asarray(GaussianDelta.evaluate(diff, eta_arr), dtype=np.float64)
        sigma_ext = math.pi * np.einsum(
            "ekn,knm,ekm->e", g_vals, vx2, g_vals, optimize=True
        )
        return np.ascontiguousarray(sigma_ext, dtype=np.float64)


def make_basis(
    basis_name: str,
    model: TightBindingModel,
    supercell: Sequence[int],
    *,
    disorder_w: float = 0.0,
    seed: int = 12345,
    disorder_seed: int | None = None,
    disorder: npt.NDArray[np.complex128] | npt.NDArray[np.float64] | None = None,
    csr_dir: str | None = None,
) -> RealBasis | BlochBasis | EigenBasis:
    """Construct one of :class:`RealBasis`, :class:`BlochBasis`, or :class:`EigenBasis` by name."""
    key = basis_name.strip().lower()
    if key == "real":
        return RealBasis(
            model,
            supercell,
            disorder_w=disorder_w,
            seed=seed,
            disorder_seed=disorder_seed,
            disorder=disorder,
            csr_dir=csr_dir,
        )
    if key == "bloch":
        return BlochBasis(
            model,
            supercell,
            disorder_w=disorder_w,
            seed=seed,
            disorder_seed=disorder_seed,
            disorder=disorder,
        )
    if key == "eigen":
        return EigenBasis(model, supercell, disorder_w=disorder_w, seed=seed)
    raise ValueError(
        f"make_basis: unknown basis '{basis_name}' (expected 'real', 'bloch', or 'eigen')"
    )


class BlochMetricBasis:
    """Clean nonorthogonal Bloch blocks with a single cached exact metric.

    ``X=S^-1`` and ``Y=S^-1/2`` are built once, O(N_k W^3), and exposed as
    ordinary C++ BlockDiagonalOperators. H and every observable are physical
    matrices; the word compiler inserts X or Y for the selected route. No
    iterative solver, sparse supercell inverse, or disorder is used here.
    The model protocol is Hk, Sk, derivatives, Ak, k_mesh, tau and uc; it has
    no dependency on a Wannier reader.
    """
    name = "bloch"

    def __init__(self, model, supercell, *, velocity_mode="berry_connection",
                 route="lowdin", seed=12345):
        from .models.nonorthogonal import adjoint, exact_block_metric
        if route not in ("lowdin", "generalised", "anticommutator"):
            raise ValueError("route must be lowdin, generalised or anticommutator")
        if velocity_mode not in ("bare", "berry_connection", "covariant"):
            raise ValueError("unknown velocity mode")
        if velocity_mode == "covariant" and model.dipoles is None:
            raise ValueError("covariant velocity requires a dipole matrix")
        self.model, self.supercell = model, tuple(map(int, supercell))
        self.velocity_mode, self.route, self.seed = velocity_mode, route, int(seed)
        self.k_points = model.k_mesh(self.supercell)
        self.num_cells, self.num_wann = len(self.k_points), model.num_wann
        self.dimension = self.num_cells*self.num_wann
        self.hk_array, self.sk_array = model.Hk(self.k_points), model.Sk(self.k_points)
        self.inverse_array, self.inverse_sqrt_array = exact_block_metric(self.sk_array)
        self._h_op = _core.BlockDiagonalOperator(self.hk_array)
        self._s_op = _core.BlockDiagonalOperator(self.sk_array)
        self.metric = SimpleNamespace(is_orthogonal=False, overlap=self._s_op,
                                      inverse_overlap=_core.BlockDiagonalOperator(self.inverse_array),
                                      inverse_sqrt_overlap=_core.BlockDiagonalOperator(self.inverse_sqrt_array))
        Y = self.inverse_sqrt_array
        lowdin_h = Y @ self.hk_array @ Y
        self.eigenvalues, u = np.linalg.eigh(lowdin_h)
        self.eigenvectors = np.ascontiguousarray(Y @ u)  # C^+ S C=I
        self.physical_eigenvalues = self.eigenvalues
        self.physical_eigenvectors = self.eigenvectors
        # Clean blocks permit the metric composition to be cached once.
        # Each recursion step then enters one C++ block matvec, rather than
        # allocating full-size intermediates for two/three metric matvecs.
        applied_h = lowdin_h if route == "lowdin" else self.inverse_array @ self.hk_array
        if route == "anticommutator":
            # metric_forms.hamiltonian_forms: this is a Hermitian approximation,
            # NOT a similarity transform of the generalised Hamiltonian.
            symmetric_h = 0.5*(self.inverse_array @ self.hk_array + self.hk_array @ self.inverse_array)
            self._symmetric_h_op = _core.BlockDiagonalOperator(symmetric_h)
            self.eigenvalues, u = np.linalg.eigh(symmetric_h)
            self.eigenvectors = np.ascontiguousarray(u)
            applied_h = symmetric_h
        self._applied_h_op = _core.BlockDiagonalOperator(np.ascontiguousarray(applied_h))
        self._lowdin_seed = _core.BlockDiagonalOperator(u)
        self._coefficient_seed = _core.BlockDiagonalOperator(self.eigenvectors)
        identity = np.eye(self.num_wann)
        self.metric_residuals = {
            "inverse": float(np.max(np.abs(self.sk_array @ self.inverse_array-identity))),
            "inverse_sqrt": float(np.max(np.abs(Y @ self.sk_array @ Y-identity))),
            "inverse_hermiticity": float(np.max(np.abs(self.inverse_array-adjoint(self.inverse_array)))),
            "generalised_eigenvectors": float(np.max(np.abs(
                adjoint(self.physical_eigenvectors) @ self.sk_array @ self.physical_eigenvectors-identity))),
        }
        self._blocks, self._op_cache = {}, {}
        self._eigen_basis = None

    def hamiltonian(self):
        """Return the physical H(k) matrix, BEFORE the metric policy acts."""
        return self._h_op

    def overlap(self):
        """Return S(k) as an ordinary block-diagonal operator."""
        return self._s_op

    def rescaled_hamiltonian(self, window):
        """Return ``a Y H Y+b I`` or ``a X H+b I = X(aH+bS)``.

        Symbolic source: metric_forms.hamiltonian_forms. Clean applied blocks
        are cached once; all per-moment block loops run in C++.
        The shift multiplies the identity in coefficient representation. A
        shift by bX would change every generalised eigenvalue incorrectly.
        """
        return _core.RescaledOperator(self._applied_h_op, window)

    def operator_blocks(self, name):
        """Return physical matrices before inserting any metric factors.

        Symbolic sources: metric_forms.velocity_forms, spin_forms and
        spin_current_forms. Hermiticity is checked in this physical form;
        the generalised applied form is only S-Hermitian.
        """
        from .models.nonorthogonal import (bloch_velocity, spin_operators,
                                           applied_spin_operators, spin_current,
                                           checked_hermitian)
        key = name.upper()
        if key in ("H", "HAM"):
            return self.hk_array
        if key == "S":
            return self.sk_array
        if key not in self._blocks:
            if key in ("VX", "VY", "VZ"):
                axis = "XYZ".index(key[1])
                dH, dS = self.model.derivatives(self.k_points, axis, cell=self.velocity_mode == "bare")
                A = self.model.Ak(self.k_points, axis) if self.velocity_mode == "covariant" else None
                self._blocks[key] = bloch_velocity(self.hk_array, self.sk_array, dH, dS, A,
                                                   self.velocity_mode, inverse_overlap=self.inverse_array)
            elif key in ("SX", "SY", "SZ"):
                if self.model.physical_spin is not None:
                    self._blocks.update(self.model.physical_spin)
                else:
                    self._blocks.update(spin_operators(self.sk_array, "interleaved"))
            elif len(key) == 4 and key[0] == "J" and key[1] in "XYZ" and key[2] == "S" and key[3] in "XYZ":
                V = self.operator_blocks("V"+key[1])
                if self.model.physical_spin is None:
                    q = applied_spin_operators(self.num_wann)["S"+key[3]]
                    self._blocks[key] = spin_current(V, q)
                else:
                    Sigma = self.operator_blocks("S"+key[3])
                    self._blocks[key] = 0.5*(V @ self.inverse_array @ Sigma +
                                              Sigma @ self.inverse_array @ V)
            else:
                raise ValueError(f"unknown BlochMetricBasis operator {name!r}")
            self._blocks[key] = checked_hermitian(self._blocks[key], name=f"physical {key}")
        return self._blocks[key]

    def operator(self, name):
        """Return a physical observable, or O_sym on the anticommutator route.

        The latter is already an applied Hermitian matrix. Its compiler uses
        the orthogonal policy, so no second inverse is inserted.
        """
        key = name.upper()
        if key in ("1", "I", "IDENTITY"):
            # The matrix of physical identity in a nonorthogonal basis is S.
            return (_core.IdentityOperator(self.dimension) if self.route == "anticommutator"
                    else self._s_op)
        if key not in self._op_cache:
            blocks = self.operator_blocks(key)
            if self.route == "anticommutator":
                blocks = self.anticommutator_blocks(blocks)
            self._op_cache[key] = _core.BlockDiagonalOperator(blocks)
        return self._op_cache[key]

    def anticommutator_blocks(self, blocks):
        """metric_forms.operator_forms: ``O_sym=(X O+O X)/2``.

        Only the single-operator trace identity is universal. Products and
        functions of H_sym are generally approximations; no plateau equality
        with the physical oracle is asserted on this route.
        """
        return 0.5*(self.inverse_array @ blocks + blocks @ self.inverse_array)

    def states(self, trace_mode="random", num_states=None, seed=None):
        """Generate orbital Bloch phases, or an exact block-local trace.

        Symbolic source: metric_forms.trace_equivalence. Random phases live
        in the orthonormal orbital (Lowdin) representation; the adapter maps
        them to coefficient vectors on the generalised route. Exact seeds
        span all bands at once in every independent k block, using W states.
        """
        if trace_mode == "random":
            return _core.RandomPhaseStates(self.dimension, 1 if num_states is None else int(num_states),
                                           self.seed if seed is None else int(seed))
        if trace_mode == "exact":
            if num_states is not None and int(num_states) != self.num_wann:
                raise ValueError("block exact trace requires W states")
            return _core.BlockLocalStates(self.num_cells, self.num_wann)
        raise ValueError("trace must be random or exact")

    def bind_eigen_trace(self, plan, operators, trace_mode="exact"):
        """Prepare an exact trace or an orbital-phase production estimator.

        Symbolic source: metric_forms.trace_equivalence. Exact seeds use C
        on the ket and SC on the dual bra, with C^+SC=I. Orbital phases use
        Y and SY instead, where Y=S^-1/2; Lowdin needs no seed transform.
        A rank-1 plan gets an explicit unexpanded bra.
        """
        operators = list(operators)
        if trace_mode == "exact":
            ket = self._coefficient_seed if self.route == "generalised" else self._lowdin_seed
        elif trace_mode == "random":
            ket = (self.metric.inverse_sqrt_overlap if self.route == "generalised"
                   else _core.IdentityOperator(self.dimension))
        else:
            raise ValueError("trace must be random or exact")
        # MirrorToBra starts from the SAME original r, not the prepared ket.
        # Prepare K r explicitly on both sides, then insert S on the bra.
        # With the compiler's X, X*A, T_m, S this gives equations (2a-c)
        # in the metric dual of K r: left_m = S*T_m*X*A*K*r.
        bra = self._s_op if self.route == "generalised" else _core.IdentityOperator(self.dimension)
        ki, bi = len(operators), len(operators)+1
        operators.extend([ket, bra])
        code, mirrored = [_core.ApplyOperator(ki)], False
        for instruction in plan.code:
            if isinstance(instruction, _core.MirrorToBra):
                code.extend([instruction, _core.ApplyOperator(ki), _core.ApplyOperator(bi)])
                mirrored = True
            elif isinstance(instruction, _core.Trace) and not mirrored:
                code.extend([_core.MirrorToBra(), _core.ApplyOperator(ki),
                             _core.ApplyOperator(bi), instruction])
            else:
                code.append(instruction)
        return _core.Plan(code, plan.word), operators

    def window(self, pad=0.10):
        """Return [emin-pad*width,emax+pad*width] from generalised eigenvalues."""
        if pad < 0:
            raise ValueError("window padding must be nonnegative")
        lo, hi = float(self.eigenvalues.min()), float(self.eigenvalues.max())
        width = max(hi-lo, 1.0)
        return _core.SpectralWindow(lo-pad*width, hi+pad*width)

    def eigenbasis(self):
        """Return the cached exact adapter ``diag(E), C^+ O C`` on the SAME blocks."""
        if self._eigen_basis is None:
            self._eigen_basis = EigenMetricBasis(self)
        return self._eigen_basis


class EigenMetricBasis:
    """Generalised eigenbasis oracle adapter, reusing the cached C^+SC=I solve."""
    name = "eigen"

    def __init__(self, bloch):
        self.bloch = bloch
        for name in ("dimension", "num_cells", "num_wann", "eigenvalues", "eigenvectors"):
            setattr(self, name, getattr(bloch, name))
        self._h_op = _core.DiagonalOperator(np.asarray(self.eigenvalues.ravel(), complex))
        self._op_cache = {}
        self.hermiticity_residuals = {}

    def hamiltonian(self):
        """Return diag(E), the operator in its exact eigenbasis."""
        return self._h_op

    def operator(self, name):
        """Return the Hermitian eigen representation of a physical operator.

        Symbolic source: metric_forms.operator_forms. Validate before the
        metric transform and remove only accepted transform round-off.
        """
        from .models.nonorthogonal import adjoint, checked_hermitian
        key = name.upper()
        if key in ("1", "I", "IDENTITY"):
            return _core.IdentityOperator(self.dimension)
        if key not in self._op_cache:
            C = self.eigenvectors
            blocks = adjoint(C) @ (
                    self.bloch.anticommutator_blocks(self.bloch.operator_blocks(key))
                    if self.bloch.route == "anticommutator" else self.bloch.operator_blocks(key)) @ C
            # The physical input was checked by operator_blocks. Only accepted
            # transform round-off may be folded; report the norm actually removed.
            anti_norm = float(np.max(0.5*np.linalg.norm(blocks-adjoint(blocks), axis=(-2, -1))))
            scale = max(float(np.max(np.linalg.norm(blocks, axis=(-2, -1)))), np.finfo(float).tiny)
            folded = checked_hermitian(blocks, name=f"eigen {key}")
            self.hermiticity_residuals[key] = dict(removed_anti_hermitian_norm=anti_norm,
                                                 operator_scale=scale, relative_norm=anti_norm/scale)
            print(f"eigen {key}: Hermitian fold removed anti-Hermitian norm={anti_norm:.17g} "
                  f"operator_scale={scale:.17g} relative_norm={anti_norm/scale:.17g}", flush=True)
            self._op_cache[key] = _core.BlockDiagonalOperator(folded)
        return self._op_cache[key]

    def states(self, trace_mode="random", num_states=None, seed=None):
        """Expose exact block trace, or transform orbital Bloch phases once.

        Symbolic source: metric_forms.trace_equivalence. Exact uses W
        block-local states with diagonal H. Random uses U^+ r, where r has
        orbital Bloch phases and U is the cached Lowdin eigenvector matrix.
        The transform and all k loops execute inside C++ TransformedStates.
        """
        if trace_mode == "exact":
            if num_states is not None and int(num_states) != self.num_wann:
                raise ValueError("block exact trace requires W states")
            return _core.BlockLocalStates(self.num_cells, self.num_wann)
        if trace_mode == "random":
            orbital = self.bloch.states("random", num_states, seed)
            return _core.TransformedStates(self.bloch._lowdin_seed.adjoint(), orbital)
        raise ValueError("trace must be random or exact")

    def window(self, pad=0.10):
        """Use the same generalised spectral window as the Bloch calculation."""
        return self.bloch.window(pad)
