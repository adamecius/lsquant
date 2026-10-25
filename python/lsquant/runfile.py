"""Run-file schema, strict validation, execution, reconstruction, and inspection."""

from __future__ import annotations

from dataclasses import dataclass, field
import glob
import json
import os
import shutil
from typing import Any, Mapping, Sequence

import numpy as np

from . import _core
from .projection.spectral import resolve_projection
from .bases import (
    BlochBasis,
    EigenBasis,
    RealBasis,
    make_basis,
    validate_w2s_csr_dir,
)
from .models.wannier import WannierModel, load_wannier_model
from .moments import run_plan, window_from_gershgorin
from .reconstruct import (
    dos_curve,
    shc,
    kubo_bastin,
    kubo_bastin_sea,
    kubo_bastin_surf,
    kubo_fft_grid,
    kubo_greenwood,
    time_correlation_curve,
)
from .words import compile as compile_word, parse as parse_word
from .words.catalogue import SHAPE_ALLOWED_FORMULAS


ALLOWED_TOP_KEYS: frozenset[str] = frozenset(
    {
        "model",
        "basis",
        "metric",
        "projection",
        "word",
        "words",
        "expansion",
        "time",
        "reconstruction",
        "output",
        "receipt",
        "reference",
        "build_profile",
        "runtime",
    }
)

ALLOWED_MODEL_KEYS: frozenset[str] = frozenset(
    {
        "kind",
        "dir",
        "label",
        "seed",
        "operators",
        "csr_dir",
        "name",
        "tutorial",
        "params",
    }
)

GENERATOR_CATALOGUE: dict[str, str] = {
    "make_chain": "examples/01_dos_1d_chain/make_chain.py",
    "make_disordered_chain": "examples/02_localization_1d_chain/make_disordered_chain.py",
    "make_haldane": "examples/03_haldane_hall/make_haldane.py",
    "make_chain_mag": "examples/04_spin_precession_magnetic_chain/lsqspin.py",
    "make_square": "examples/05_charge_transport_square/make_square.py",
    "make_ey": "examples/06_spin_relaxation_square/make_ey.py",
    "make_rashba": "examples/06_spin_relaxation_square/make_rashba.py",
    "gen_models": "python/lsquant/models/wannier.py",
}

_GENERATOR_TO_TUTORIAL: dict[str, int] = {
    "make_chain": 1,
    "make_disordered_chain": 2,
    "make_haldane": 3,
    "make_chain_mag": 4,
    "make_square": 5,
    "make_ey": 6,
    "make_rashba": 6,
}

ALLOWED_GENERATOR_PARAMS: dict[str, frozenset[str]] = {
    "make_chain": frozenset({"sizes", "hopping_t"}),
    "make_disordered_chain": frozenset(
        {
            "sizes",
            "W",
            "realisations",
            "readme_size",
            "readme_disorder_seed",
        }
    ),
    "make_haldane": frozenset(),
    "make_chain_mag": frozenset(),
    "make_square": frozenset({"L", "W", "disorder_seed", "dos_moments"}),
    "make_ey": frozenset(
        {
            "ey_l",
            "ey_ds",
            "dp_l",
            "dp_lambda",
            "dp_ws",
            "dp_moments",
            "dp_steps",
            "dp_tmax",
            "ratelaw_l",
            "ratelaw_ds",
            "ratelaw_nr",
            "crosscheck_l",
            "crosscheck_d",
            "crosscheck_nr",
            "crosscheck_moments",
            "crosscheck_steps",
            "crosscheck_tmax",
            "disorder_seed",
        }
    ),
    "make_rashba": frozenset(
        {
            "ey_l",
            "ey_ds",
            "dp_l",
            "dp_lambda",
            "dp_ws",
            "dp_moments",
            "dp_steps",
            "dp_tmax",
            "ratelaw_l",
            "ratelaw_ds",
            "ratelaw_nr",
            "crosscheck_l",
            "crosscheck_d",
            "crosscheck_nr",
            "crosscheck_moments",
            "crosscheck_steps",
            "crosscheck_tmax",
            "disorder_seed",
        }
    ),
    "gen_models": frozenset(),
}


def resolve_generator_path(
    name: str, repo_root: str | None = None
) -> str:
    """Resolve a generator catalogue name to its module path under the repository root."""
    if name not in GENERATOR_CATALOGUE:
        raise ValueError(
            f"Unknown generator name '{name}'; allowed generators are {sorted(GENERATOR_CATALOGUE)}"
        )
    root = (
        os.path.abspath(repo_root)
        if repo_root is not None
        else os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..")
        )
    )
    return os.path.abspath(os.path.join(root, GENERATOR_CATALOGUE[name]))

ALLOWED_BASIS_KEYS: frozenset[str] = frozenset(
    {"name", "supercell", "disorder"}
)

ALLOWED_DISORDER_KEYS: frozenset[str] = frozenset({"W", "seed"})

ALLOWED_METRIC_KEYS: frozenset[str] = frozenset(
    {
        "policy",
        "route",
        "applier",
        "order",
        "coeff_tol",
        "tol",
        "max_iter",
        "max_dim",
    }
)

ALLOWED_EXPANSION_KEYS: frozenset[str] = frozenset(
    {
        "moments",
        "window",
        "trace",
        "states",
        "seed",
        "batch",
        "state_label",
    }
)

ALLOWED_TIME_KEYS: frozenset[str] = frozenset(
    {"steps", "tmax", "dt", "bessel"}
)

ALLOWED_RECONSTRUCTION_KEYS: frozenset[str] = frozenset(
    {
        "kernel",
        "broadening_meV",
        "fermi_energies",
        "formulas",
        "num_points",
        "fft_grid",
        "legacy_text_round_trip",
        "per_site_dos",
        "targets",
        "num_cells",
        "cell_area",
        "energies",
        "energy_grid",
        "bastin_backend",
    }
)

ALLOWED_RECON_TARGET_KEYS: frozenset[str] = frozenset(
    {
        "formulas",
        "broadening_meV",
        "fermi_energies",
        "subdir",
        "also_top_dir",
        "legacy_text_round_trip",
    }
)

ALLOWED_OUTPUT_KEYS: frozenset[str] = frozenset(
    {
        "dir",
        "legacy_names",
        "legacy_layout",
        "write_chebmom",
        "write_basis_json",
        "capture_tutorial",
    }
)

def default_window_for_basis(basis_name: str) -> str:
    """Return the default spectral window mode for ``basis_name``.

    Defaults to ``"eigen"`` on ``"bloch"`` and ``"eigen"`` bases (Stage 5,
    review 1) and to ``"gershgorin"`` only on ``"real"``.
    """
    if basis_name in ("bloch", "eigen"):
        return "eigen"
    return "gershgorin"


DEFAULT_RUNFILE_DEFAULTS: dict[str, Any] = {
    "model": {
        "kind": "csr",
        "operators": [],
        "csr_dir": None,
        "name": None,
        "tutorial": None,
        "params": {},
    },
    "basis": {
        "name": "real",
        "supercell": [1, 1, 1],
        "disorder": {"W": 0.0, "seed": 7},
    },
    "metric": {
        "policy": "orthogonal",
        "route": "lowdin",
        "applier": "exact",
        "order": 60,
        "coeff_tol": 1e-12,
        "tol": 1e-13,
        "max_iter": 200,
        "max_dim": 4096,
    },
    "expansion": {
        "moments": 64,
        "window": {"real": "gershgorin", "bloch": "eigen", "eigen": "eigen"},
        "trace": "exact",
        "states": None,
        "seed": 12345,
        "batch": None,
        "state_label": None,
    },
    "time": None,
    "reconstruction": {
        "kernel": "jackson",
        "broadening_meV": None,
        "fermi_energies": [],
        "formulas": None,
        "num_points": None,
        "fft_grid": False,
        "legacy_text_round_trip": False,
        "per_site_dos": False,
        "targets": [],
    },
    "output": {
        "dir": ".",
        "legacy_names": True,
        "legacy_layout": True,
        "write_chebmom": False,
        "write_basis_json": False,
        "capture_tutorial": False,
    },
}


def defaults_for_basis(basis_name: str = "real") -> dict[str, Any]:
    """Return a copy of ``DEFAULT_RUNFILE_DEFAULTS`` with ``expansion.window`` resolved for ``basis_name``."""
    out = {
        k: (dict(v) if isinstance(v, Mapping) else v)
        for k, v in DEFAULT_RUNFILE_DEFAULTS.items()
    }
    out["expansion"] = dict(
        DEFAULT_RUNFILE_DEFAULTS["expansion"],
        window=default_window_for_basis(basis_name),
    )
    return out


def _get_or(mapping: Mapping[str, Any], key: str, default: Any = None) -> Any:
    """Return ``mapping[key]`` when present and not ``None``, otherwise ``default``."""
    val = mapping.get(key)
    return default if val is None else val


def _check_keys(
    section_name: str, data: Mapping[str, Any], allowed: frozenset[str]
) -> None:
    """Refuse any unknown key in ``data`` with a descriptive ``ValueError``."""
    for key in data.keys():
        if key not in allowed:
            raise ValueError(
                f"Unknown key '{key}' in run file ({section_name}); "
                f"allowed keys are {sorted(allowed)}"
            )


def validate_runfile_dict(raw: Mapping[str, Any]) -> None:
    """Strictly validate a raw run-file dictionary before any file or model is loaded."""
    if not isinstance(raw, Mapping):
        raise ValueError("Run file must be a JSON object")

    _check_keys("top-level", raw, ALLOWED_TOP_KEYS)
    projection = resolve_projection(raw.get("projection"))
    # Item 3: a resolved receipt cannot silently select a different binary profile.
    # Tester: replay matching receipts and verify mismatched profiles fail before work.
    if "build_profile" in raw and raw["build_profile"] != _core.build_profile:
        raise ValueError(
            f"Run file build_profile={raw['build_profile']!r} differs from "
            f"loaded engine profile={_core.build_profile!r}"
        )

    for req in ("model", "word", "expansion"):
        if req not in raw:
            raise ValueError(f"Missing required top-level key '{req}' in run file")

    model = raw["model"]
    if not isinstance(model, Mapping):
        raise ValueError("Run file 'model' must be a JSON object")
    _check_keys("model", model, ALLOWED_MODEL_KEYS)
    kind = str(_get_or(model, "kind", "csr"))
    reference = raw.get("reference")
    if reference is not None and not isinstance(reference, bool):
        raise ValueError("reference must be a boolean when supplied")
    if reference is not None and kind != "siesta":
        raise ValueError("reference is supported only for SIESTA stage6b runs")
    if projection["kind"] != "none" and kind != "siesta":
        raise ValueError("stage6d run-file projection requires a SIESTA model")
    if kind not in ("csr", "wannier", "siesta", "generator"):
        raise ValueError(
            f"Unknown model.kind '{kind}'; expected 'csr', 'wannier', 'siesta', or 'generator'"
        )
    if kind == "csr":
        if not model.get("dir") or not model.get("label"):
            raise ValueError(
                "Run file with model.kind='csr' requires 'dir' and 'label'"
            )
    elif kind == "wannier":
        if not model.get("dir") or not model.get("seed"):
            raise ValueError(
                "Run file with model.kind='wannier' requires 'dir' and 'seed'"
            )
    elif kind == "siesta":
        if not model.get("dir") or not model.get("label"):
            raise ValueError("SIESTA model requires dir and label")
        if model.get("params") is not None and not isinstance(model["params"], Mapping):
            raise ValueError("SIESTA model.params must be a JSON object")
        _check_keys("model.params (siesta)", model.get("params") or {}, frozenset({
            "velocity_mode", "oracle", "oracle_degeneracy_tolerance", "bloch_memory_mb",
            "quadrature_order", "quadrature_max_order", "quadrature_tolerance",
            "nonorthogonal_disorder"}))
        if "nonorthogonal_disorder" in (model.get("params") or {}) and not isinstance(
                model["params"]["nonorthogonal_disorder"], bool):
            raise ValueError("model.params.nonorthogonal_disorder must be a boolean")
        degeneracy_tolerance = float((model.get("params") or {}).get("oracle_degeneracy_tolerance", 1e-6))
        if not np.isfinite(degeneracy_tolerance) or degeneracy_tolerance <= 0:
            raise ValueError("oracle_degeneracy_tolerance must be positive and finite (eV)")
        bloch_memory_mb = float((model.get("params") or {}).get("bloch_memory_mb", 256.0))
        if not np.isfinite(bloch_memory_mb) or bloch_memory_mb <= 0:
            raise ValueError("bloch_memory_mb must be positive and finite (MiB)")
        if (model.get("params") or {}).get("velocity_mode", "berry_connection") not in (
                "bare", "berry_connection", "covariant"):
            raise ValueError("unknown SIESTA velocity_mode")
    elif kind == "generator":
        gen_name = model.get("name")
        if not gen_name or not isinstance(gen_name, str):
            raise ValueError(
                "Run file with model.kind='generator' requires 'name' from "
                f"{sorted(GENERATOR_CATALOGUE)}"
            )
        if gen_name not in GENERATOR_CATALOGUE:
            raise ValueError(
                f"Unknown generator name '{gen_name}' in run file (model.name); "
                f"allowed generators are {sorted(GENERATOR_CATALOGUE)}"
            )
        raw_params = model.get("params")
        if raw_params is not None:
            if not isinstance(raw_params, Mapping):
                raise ValueError("model.params must be a JSON object")
            _check_keys(
                f"model.params ({gen_name})",
                raw_params,
                ALLOWED_GENERATOR_PARAMS.get(gen_name, frozenset()),
            )
        if "basis" in raw and isinstance(raw["basis"], Mapping):
            gen_bname = str(_get_or(raw["basis"], "name", "real"))
            if gen_bname != "real":
                raise ValueError(
                    f"model.kind='generator' only supports basis.name='real' (got '{gen_bname}')"
                )
            if "disorder" in raw["basis"] and raw["basis"]["disorder"] is not None:
                raise ValueError(
                    "basis.disorder is not supported for model.kind='generator' "
                    "(pass disorder parameters in model.params)"
                )
        if "metric" in raw and isinstance(raw["metric"], Mapping):
            gen_pol = str(_get_or(raw["metric"], "policy", "orthogonal"))
            if gen_pol != "orthogonal":
                raise ValueError(
                    f"model.kind='generator' only supports metric.policy='orthogonal' (got '{gen_pol}')"
                )

    if "basis" in raw and raw["basis"] is not None:
        basis = raw["basis"]
        if not isinstance(basis, Mapping):
            raise ValueError("Run file 'basis' must be a JSON object")
        _check_keys("basis", basis, ALLOWED_BASIS_KEYS)
        bname = str(_get_or(basis, "name", "real"))
        if bname not in ("real", "bloch", "eigen"):
            raise ValueError(
                f"Unknown basis.name '{bname}'; expected 'real', 'bloch', or 'eigen'"
            )
        if "supercell" in basis and basis["supercell"] is not None:
            sc = basis["supercell"]
            if (
                not isinstance(sc, (list, tuple))
                or len(sc) != 3
                or not all(isinstance(x, int) and x > 0 for x in sc)
            ):
                raise ValueError(
                    "basis.supercell must be a list of 3 positive integers [N1, N2, N3]"
                )
        if "disorder" in basis and basis["disorder"] is not None:
            dis = basis["disorder"]
            if not isinstance(dis, Mapping):
                raise ValueError("basis.disorder must be a JSON object")
            _check_keys("basis.disorder", dis, ALLOWED_DISORDER_KEYS)
            if float(_get_or(dis, "W", 0.0)) < 0.0:
                raise ValueError("basis.disorder.W must be non-negative")

    if "metric" in raw and raw["metric"] is not None:
        metric = raw["metric"]
        if not isinstance(metric, Mapping):
            raise ValueError("Run file 'metric' must be a JSON object")
        _check_keys("metric", metric, ALLOWED_METRIC_KEYS)
        pol = str(_get_or(metric, "policy", "orthogonal"))
        if pol not in ("orthogonal", "nonorthogonal"):
            raise ValueError(
                f"Unknown metric.policy '{pol}'; expected 'orthogonal' or 'nonorthogonal'"
            )
        rte = str(_get_or(metric, "route", "lowdin"))
        if rte not in ("lowdin", "generalised", "anticommutator"):
            raise ValueError(
                f"Unknown metric.route '{rte}'; expected 'lowdin', 'generalised' or 'anticommutator'"
            )
        if rte == "anticommutator" and kind != "siesta":
            raise ValueError("anticommutator metric is an explicitly approximate SIESTA route")
        app = str(_get_or(metric, "applier", "exact"))
        if app not in ("exact", "chebyshev", "ischo"):
            raise ValueError(
                f"Unknown metric.applier '{app}'; expected 'exact', 'chebyshev', or 'ischo'"
            )
        if pol == "nonorthogonal" and rte == "lowdin" and app == "ischo":
            raise ValueError(
                "metric.applier='ischo' provides S^-1 only and cannot be used with route='lowdin'"
            )

    exp = raw["expansion"]
    if not isinstance(exp, Mapping):
        raise ValueError("Run file 'expansion' must be a JSON object")
    _check_keys("expansion", exp, ALLOWED_EXPANSION_KEYS)
    if _get_or(exp, "moments") is None or int(exp["moments"]) <= 0:
        raise ValueError("expansion.moments must be a positive integer")
    tr = str(_get_or(exp, "trace", "exact"))
    if tr not in ("exact", "random"):
        raise ValueError(
            f"Unknown expansion.trace '{tr}'; expected 'exact' or 'random'"
        )
    if "window" in exp and exp["window"] is not None:
        win = exp["window"]
        if isinstance(win, str):
            if win not in ("eigen", "gershgorin"):
                raise ValueError(
                    f"Unknown expansion.window '{win}'; expected 'eigen', 'gershgorin', or [lo, hi]"
                )
            if kind in ("csr", "generator") and win == "eigen":
                raise ValueError(
                    f"expansion.window='eigen' is not supported for model.kind='{kind}' "
                    "(no k blocks; use gershgorin or an explicit window)"
                )
        elif (
            isinstance(win, (list, tuple))
            and len(win) == 2
            and all(isinstance(v, (int, float)) for v in win)
        ):
            if not (float(win[1]) > float(win[0])):
                raise ValueError(
                    "expansion.window [lo, hi] requires hi > lo"
                )
        else:
            raise ValueError(
                "expansion.window must be 'eigen', 'gershgorin', or [lo, hi]"
            )

    if "time" in raw and raw["time"] is not None:
        tm = raw["time"]
        if not isinstance(tm, Mapping):
            raise ValueError("Run file 'time' must be a JSON object")
        _check_keys("time", tm, ALLOWED_TIME_KEYS)
        if _get_or(tm, "steps") is None or int(tm["steps"]) <= 1:
            raise ValueError("time.steps must be an integer > 1")
        if _get_or(tm, "tmax") is None and _get_or(tm, "dt") is None:
            raise ValueError("time requires 'tmax' or 'dt'")

    if "reconstruction" in raw and raw["reconstruction"] is not None:
        rec = raw["reconstruction"]
        if not isinstance(rec, Mapping):
            raise ValueError("Run file 'reconstruction' must be a JSON object")
        _check_keys("reconstruction", rec, ALLOWED_RECONSTRUCTION_KEYS)
        if rec.get("bastin_backend", "tables") not in ("tables", "folded"):
            raise ValueError("reconstruction.bastin_backend must be tables or folded")
        if rec.get("bastin_backend") == "folded" and (
                raw.get("model", {}).get("kind") != "siesta" or rec.get("fft_grid", False)):
            raise ValueError("folded run-file reconstruction requires SIESTA and the uniform grid")
        if rec.get("energies") is not None and rec.get("energy_grid") is not None:
            raise ValueError("choose reconstruction.energies or energy_grid, not both")
        if rec.get("energies") is not None:
            energies = np.asarray(rec["energies"], dtype=float)
            if energies.ndim != 1 or energies.size == 0 or not np.all(np.isfinite(energies)):
                raise ValueError("reconstruction.energies must be a nonempty finite list")
            if np.any(np.diff(energies) <= 0):
                raise ValueError("reconstruction.energies must be strictly increasing")
        if rec.get("energy_grid") is not None:
            grid = np.asarray(rec["energy_grid"], dtype=float)
            if (grid.shape != (3,) or not np.all(np.isfinite(grid)) or grid[0] >= grid[1]
                    or grid[2] < 2 or grid[2] != int(grid[2])):
                raise ValueError("reconstruction.energy_grid must be [E_min, E_max, n], n >= 2")
        for tgt in rec.get("targets") or ():
            if not isinstance(tgt, Mapping):
                raise ValueError(
                    "Each entry in reconstruction.targets must be a JSON object"
                )
            _check_keys("reconstruction.targets", tgt, ALLOWED_RECON_TARGET_KEYS)

    if "output" in raw and raw["output"] is not None:
        out = raw["output"]
        if not isinstance(out, Mapping):
            raise ValueError("Run file 'output' must be a JSON object")
        _check_keys("output", out, ALLOWED_OUTPUT_KEYS)

    # Validate word(s) syntax, catalogue membership, basis consistency, and time requirements
    word_str = raw["word"]
    if not isinstance(word_str, str) or not word_str.strip():
        raise ValueError("Run file 'word' must be a non-empty string")
    words_list = [word_str]
    if "words" in raw and raw["words"] is not None:
        if not isinstance(raw["words"], (list, tuple)) or not raw["words"]:
            raise ValueError("Run file 'words' must be a non-empty list of strings")
        words_list = [str(w) for w in raw["words"]]

    basis_name = str(_get_or(raw.get("basis") or {}, "name", "real"))
    model_ops = model.get("operators")
    basis_spec: dict[str, Any] = {"name": basis_name}
    if model_ops is not None:
        basis_spec["operators"] = list(model_ops)
    if "metric" in raw and raw["metric"] is not None:
        basis_spec["metric"] = dict(raw["metric"])

    raw_rec_obj = (
        dict(raw["reconstruction"])
        if isinstance(raw.get("reconstruction"), Mapping)
        else None
    )
    req_formulas_raw = (
        raw_rec_obj.get("formulas") if raw_rec_obj is not None else None
    )
    if req_formulas_raw is not None and not isinstance(
        req_formulas_raw, (list, tuple)
    ):
        raise ValueError("reconstruction.formulas must be a list of strings")

    union_allowed_formulas: set[str] = set()
    for w in words_list:
        parsed_w = parse_word(w)
        comp_base = compile_word(
            parsed_w,
            basis=basis_spec,
            params={
                "expansion": dict(exp),
                "time": dict(raw["time"]) if raw.get("time") else None,
                "metric": dict(raw["metric"]) if raw.get("metric") else None,
            },
        )
        w_allowed = SHAPE_ALLOWED_FORMULAS.get(
            comp_base.match.entry.shape_id, frozenset()
        )
        union_allowed_formulas.update(w_allowed)
        if raw_rec_obj is not None:
            rec_for_w = dict(raw_rec_obj)
            if req_formulas_raw is not None:
                rec_for_w["formulas"] = [
                    str(f)
                    for f in req_formulas_raw
                    if str(f).lower() in w_allowed
                ]
            compile_word(
                parsed_w,
                basis=basis_spec,
                params={
                    "expansion": dict(exp),
                    "time": dict(raw["time"]) if raw.get("time") else None,
                    "metric": dict(raw["metric"]) if raw.get("metric") else None,
                    "reconstruction": rec_for_w,
                },
            )

    if req_formulas_raw is not None:
        for f in req_formulas_raw:
            if str(f).lower() not in union_allowed_formulas:
                raise ValueError(
                    f"Reconstruction formula '{f}' in reconstruction.formulas "
                    f"does not match any word in {words_list} "
                    f"(allowed formulas for given word(s): {sorted(union_allowed_formulas)})"
                )


@dataclass
class RunFile:
    """Validated and defaults-resolved Stage-7 run configuration."""

    model: dict[str, Any]
    basis: dict[str, Any]
    metric: dict[str, Any]
    word: str
    words: list[str]
    expansion: dict[str, Any]
    time: dict[str, Any] | None
    reconstruction: dict[str, Any]
    output: dict[str, Any]
    source_path: str | None = field(default=None, repr=False)
    reference: bool | None = None
    projection: dict[str, Any] = field(default_factory=lambda: {"kind": "none", "window": None, "radius": None})

    @classmethod
    def from_dict(
        cls,
        raw: Mapping[str, Any],
        *,
        source_path: str | None = None,
        repo_root: str | None = None,
    ) -> RunFile:
        """Validate ``raw`` and populate defaults to construct a :class:`RunFile`."""
        validate_runfile_dict(raw)

        root = (
            os.path.abspath(repo_root)
            if repo_root is not None
            else os.path.abspath(
                os.path.join(os.path.dirname(__file__), "..", "..")
            )
        )

        def _resolve_path(p: str | None) -> str | None:
            if p is None:
                return None
            p = os.path.expanduser(p)
            if os.path.isabs(p):
                return p
            cand = os.path.abspath(os.path.join(root, p))
            if os.path.exists(cand):
                return cand
            if source_path is not None:
                rel_src = os.path.abspath(
                    os.path.join(os.path.dirname(source_path), p)
                )
                if os.path.exists(rel_src):
                    return rel_src
            return cand

        raw_model = dict(raw["model"])
        kind = str(_get_or(raw_model, "kind", "csr"))
        gen_name = raw_model.get("name")
        if gen_name is not None:
            _ = resolve_generator_path(str(gen_name), root)
        model_cfg: dict[str, Any] = {
            "kind": kind,
            "dir": _resolve_path(raw_model.get("dir")),
            "label": raw_model.get("label"),
            "seed": raw_model.get("seed"),
            "operators": (
                list(raw_model["operators"])
                if raw_model.get("operators") is not None
                else None
            ),
            "csr_dir": _resolve_path(raw_model.get("csr_dir")),
            "name": str(gen_name) if gen_name is not None else None,
            "tutorial": raw_model.get("tutorial"),
            "params": dict(raw_model.get("params") or {}),
        }

        raw_basis = dict(raw.get("basis") or {})
        raw_dis = dict(raw_basis.get("disorder") or {})
        basis_cfg: dict[str, Any] = {
            "name": str(_get_or(raw_basis, "name", "real")),
            "supercell": [
                int(x) for x in _get_or(raw_basis, "supercell", [1, 1, 1])
            ],
            "disorder": {
                "W": float(_get_or(raw_dis, "W", 0.0)),
                "seed": int(_get_or(raw_dis, "seed", 7)),
            },
        }

        raw_metric = dict(raw.get("metric") or {})
        metric_cfg: dict[str, Any] = {
            "policy": str(_get_or(raw_metric, "policy", "orthogonal")),
            "route": str(_get_or(raw_metric, "route", "lowdin")),
            "applier": str(_get_or(raw_metric, "applier", "exact")),
            "order": int(_get_or(raw_metric, "order", 60)),
            "coeff_tol": float(_get_or(raw_metric, "coeff_tol", 1e-12)),
            "tol": float(_get_or(raw_metric, "tol", 1e-13)),
            "max_iter": int(_get_or(raw_metric, "max_iter", 200)),
            "max_dim": int(_get_or(raw_metric, "max_dim", 4096)),
        }

        raw_exp = dict(raw["expansion"])
        default_win = default_window_for_basis(basis_cfg["name"])
        win_val = raw_exp.get("window")
        if win_val is None:
            win_val = default_win
        elif isinstance(win_val, (list, tuple)):
            win_val = [float(win_val[0]), float(win_val[1])]
        moments_val = int(raw_exp["moments"])
        expansion_cfg: dict[str, Any] = {
            "moments": moments_val,
            "window": win_val,
            "trace": str(_get_or(raw_exp, "trace", "exact")),
            "states": (
                int(raw_exp["states"])
                if raw_exp.get("states") is not None
                else None
            ),
            "seed": int(_get_or(raw_exp, "seed", 12345)),
            "batch": (
                int(raw_exp["batch"])
                if raw_exp.get("batch") is not None
                else moments_val
            ),
            "state_label": raw_exp.get("state_label"),
        }

        time_cfg: dict[str, Any] | None = None
        if raw.get("time") is not None:
            rt = dict(raw["time"])
            time_cfg = {
                "steps": int(rt["steps"]),
                "tmax": float(rt["tmax"]) if rt.get("tmax") is not None else None,
                "dt": float(rt["dt"]) if rt.get("dt") is not None else None,
                "bessel": str(_get_or(rt, "bessel", "scipy")),
            }

        raw_rec = dict(raw.get("reconstruction") or {})
        recon_cfg: dict[str, Any] = {
            "kernel": str(_get_or(raw_rec, "kernel", "jackson")).lower(),
            "broadening_meV": (
                float(raw_rec["broadening_meV"])
                if raw_rec.get("broadening_meV") is not None
                else None
            ),
            "fermi_energies": [
                float(e) for e in (raw_rec.get("fermi_energies") or [])
            ],
            "formulas": (
                [str(f) for f in raw_rec["formulas"]]
                if raw_rec.get("formulas") is not None
                else None
            ),
            "num_points": (
                int(raw_rec["num_points"])
                if raw_rec.get("num_points") is not None
                else None
            ),
            "fft_grid": bool(_get_or(raw_rec, "fft_grid", False)),
            **({"bastin_backend": raw_rec["bastin_backend"]} if "bastin_backend" in raw_rec else {}),
            "legacy_text_round_trip": bool(
                _get_or(raw_rec, "legacy_text_round_trip", False)
            ),
            "per_site_dos": bool(_get_or(raw_rec, "per_site_dos", False)),
            "targets": [dict(t) for t in (raw_rec.get("targets") or [])],
            **({key: raw_rec[key] for key in ("num_cells", "cell_area") if key in raw_rec}),
            **({key: raw_rec[key] for key in ("energies", "energy_grid") if key in raw_rec}),
        }

        raw_out = dict(raw.get("output") or {})
        default_write_cheb = (
            kind == "wannier"
            or metric_cfg["policy"] == "nonorthogonal"
            or not bool(_get_or(raw_out, "legacy_layout", True))
        )
        output_cfg: dict[str, Any] = {
            "dir": str(_get_or(raw_out, "dir", ".")),
            "legacy_names": bool(_get_or(raw_out, "legacy_names", True)),
            "legacy_layout": bool(_get_or(raw_out, "legacy_layout", True)),
            "write_chebmom": bool(
                _get_or(raw_out, "write_chebmom", default_write_cheb)
            ),
            "write_basis_json": bool(
                _get_or(raw_out, "write_basis_json", kind == "wannier")
            ),
            "capture_tutorial": bool(_get_or(raw_out, "capture_tutorial", False)),
        }

        word_canonical = parse_word(str(raw["word"])).canonical()
        if raw.get("words") is not None:
            words_canonical = [
                parse_word(str(w)).canonical() for w in raw["words"]
            ]
        else:
            words_canonical = [word_canonical]

        return cls(
            model=model_cfg,
            basis=basis_cfg,
            metric=metric_cfg,
            word=word_canonical,
            words=words_canonical,
            expansion=expansion_cfg,
            time=time_cfg,
            reconstruction=recon_cfg,
            output=output_cfg,
            source_path=source_path,
            reference=raw.get("reference"),
            projection=resolve_projection(raw.get("projection")),
        )

    @classmethod
    def from_file(
        cls, path: str, *, repo_root: str | None = None
    ) -> RunFile:
        """Load and validate a :class:`RunFile` from a JSON file on disk."""
        with open(path, encoding="utf-8") as fh:
            raw = json.load(fh)
        return cls.from_dict(
            raw, source_path=os.path.abspath(path), repo_root=repo_root
        )

    def to_resolved_dict(self) -> dict[str, Any]:
        """Return the fully resolved run-file dictionary written to ``run.resolved.json``."""
        out: dict[str, Any] = {
            "build_profile": _core.build_profile,
            # Item 6: record effective launch choices, never use receipt values
            # to override the caller's environment. Tester archives with timings.
            "runtime": {key: os.environ.get(key) for key in (
                "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "LSQUANT_SPMV_MIN_NNZ",
                "LSQUANT_KERNEL_BATCH", "LSQUANT_FFT_BACKEND")},
            "model": dict(self.model),
            "basis": dict(self.basis),
            "metric": dict(self.metric),
            "word": self.word,
            "expansion": dict(self.expansion),
            "time": dict(self.time) if self.time is not None else None,
            "reconstruction": dict(self.reconstruction),
            "output": dict(self.output),
        }
        if self.reference is not None:
            out["reference"] = self.reference
        if len(self.words) > 1 or self.words != [self.word]:
            out["words"] = list(self.words)
        if self.projection["kind"] != "none":
            out["projection"] = dict(self.projection)
        return out


def _compute_padded_overlap_bounds(
    s_op: _core.LinearOperator,
) -> _core.SpectralBounds:
    """Compute Lanczos bounds on overlap ``S`` with the 5% safety pad of ``stage6_nonorth``."""
    ritz = _core.LanczosBounds(s_op, 80, 12345)
    span = max(1e-6, float(ritz.upper - ritz.lower))
    lower = max(0.5 * float(ritz.lower), float(ritz.lower) - 0.05 * span)
    upper = float(ritz.upper) + 0.05 * span
    return _core.SpectralBounds(lower, upper)


def _discover_csr_operators(model_dir: str, label: str) -> list[str]:
    """List available static operator names ``<label>.<OP>.CSR`` in ``model_dir``."""
    prefix = os.path.join(model_dir, f"{label}.")
    suffix = ".CSR"
    ops: list[str] = []
    for path in sorted(glob.glob(f"{prefix}*{suffix}")):
        base = os.path.basename(path)
        op_name = base[len(label) + 1 : -len(suffix)]
        if op_name and op_name not in ("HAM", "S"):
            ops.append(op_name)
    return ops


def _reconstruct_single_formula(
    formula: str,
    mu: _core.ChebyshevMoments,
    broadening_meV: float,
    *,
    fermi_energies: Sequence[float] = (),
    num_points: int | None = None,
    num_cells: int | None = None,
    cell_area: float | None = None,
) -> list[tuple[str, str, np.ndarray, np.ndarray]]:
    """Run one named reconstruction formula on ``mu``.

    Returns a list of ``(legacy_prefix, legacy_suffix, x_grid, values)`` tuples.
    """
    f = formula.lower()
    if f == "shc":
        if num_cells is None or cell_area is None:
            raise ValueError("shc reconstruction requires num_cells and cell_area in the resolved run file")
        energies, values = shc(mu, broadening_meV, num_cells, cell_area)
        return [("SpinHall_", "JACKSON.dat", energies, values)]
    if f in ("dos", "dos_curve"):
        energies, values = dos_curve(mu, broadening_meV)
        return [("mean", "JACKSON.dat", energies, values)]
    if f in ("greenwood", "kubo_greenwood"):
        energies, values = kubo_greenwood(mu, broadening_meV)
        return [("KuboGreenWood_", "JACKSON.dat", energies, values)]
    if f in ("bastin", "kubo_bastin"):
        energies, values = kubo_bastin(mu, broadening_meV)
        return [("KuboBastin_", "JACKSON.dat", energies, values)]
    if f in ("sea", "kubo_bastin_sea"):
        energies, values = kubo_bastin_sea(mu, broadening_meV)
        return [("KuboBastinSea_", "JACKSON.dat", energies, values)]
    if f in ("surf", "kubo_bastin_surf"):
        energies, values = kubo_bastin_surf(mu, broadening_meV)
        return [("KuboBastinSurf_", "JACKSON.dat", energies, values)]
    if f in ("greenwood_fft", "kubo_greenwood_fft"):
        energies, values = kubo_fft_grid(mu, broadening_meV, formula="greenwood")
        return [("KuboGreenWood_", "JACKSON_FFTgrid.dat", energies, values)]
    if f in ("bastin_fft", "kubo_bastin_fft"):
        energies, values = kubo_fft_grid(mu, broadening_meV, formula="bastin")
        return [("KuboBastin_", "JACKSON_FFTgrid.dat", energies, values)]
    if f in ("time_correlation", "time_correlation_curve"):
        out: list[tuple[str, str, np.ndarray, np.ndarray]] = []
        for ef in fermi_energies:
            times, values = time_correlation_curve(
                mu, broadening_meV, float(ef)
            )
            out.append(
                ("mean", f"EF{float(ef):.6f}JACKSON.dat", times, values)
            )
        return out
    if f in ("row_m0", "row0"):
        return []
    raise ValueError(f"Unknown reconstruction formula '{formula}'")


def _execute_generator_runfile(
    rf: RunFile, out_dir: str, repo_root: str | None = None
) -> list[str]:
    """Execute a tutorial generator run file while writing ``run.resolved.json``."""
    gen_name = str(rf.model.get("name") or "")
    gen_path = resolve_generator_path(gen_name, repo_root)
    if not os.path.isfile(gen_path):
        raise FileNotFoundError(
            f"Generator module '{gen_name}' not found at '{gen_path}'"
        )
    tut = int(
        rf.model.get("tutorial")
        or _GENERATOR_TO_TUTORIAL.get(gen_name, 0)
    )
    params = dict(rf.model.get("params") or {})
    if tut == 1:
        from .apps.tutorial01 import run_tutorial01

        return run_tutorial01(
            out_dir,
            sizes=params.get("sizes", (256, 1024, 4096)),
            hopping_t=float(_get_or(params, "hopping_t", -1.0)),
            num_moments=int(rf.expansion["moments"]),
            window=rf.expansion.get("window"),
            word=rf.word,
            broadening_meV=float(
                _get_or(rf.reconstruction, "broadening_meV", 20.0)
            ),
        )
    if tut == 2:
        from .apps.tutorial02 import run_tutorial02

        assert rf.time is not None
        return run_tutorial02(
            out_dir,
            sizes=params.get("sizes", (256, 1024, 4096)),
            w_disorder=float(_get_or(params, "W", 1.0)),
            realisations=int(_get_or(params, "realisations", 16)),
            readme_size=int(_get_or(params, "readme_size", 512)),
            readme_disorder_seed=int(
                _get_or(params, "readme_disorder_seed", 1)
            ),
            kpm_seed=int(_get_or(rf.expansion, "seed", 12345)),
            num_moments=int(rf.expansion["moments"]),
            num_states=int(_get_or(rf.expansion, "states", 1)),
            window=rf.expansion.get("window"),
            ntime=int(rf.time["steps"]),
            tmax=float(_get_or(rf.time, "tmax", 600.0)),
            bessel=str(_get_or(rf.time, "bessel", "legacy")),
            word=rf.word,
            broadening_meV=float(
                _get_or(rf.reconstruction, "broadening_meV", 10.0)
            ),
            energies=(
                rf.reconstruction["fermi_energies"]
                if rf.reconstruction.get("fermi_energies")
                else (0.0, 1.0, -1.0, 1.9)
            ),
            legacy_text_round_trip=bool(
                _get_or(rf.reconstruction, "legacy_text_round_trip", True)
            ),
        )
    if tut == 3:
        from .apps.tutorial03 import run_tutorial03

        paths, _ = run_tutorial03(
            out_dir,
            operators_dir=rf.model.get("dir"),
            num_moments=int(rf.expansion["moments"]),
            broadening_meV=float(
                _get_or(rf.reconstruction, "broadening_meV", 10.0)
            ),
        )
        return paths
    if tut == 4:
        from .apps.tutorial04 import run_tutorial04

        paths, _ = run_tutorial04(
            out_dir, operators_dir=rf.model.get("dir")
        )
        return paths
    if tut == 5:
        from .apps.tutorial05 import run_tutorial05

        assert rf.time is not None
        raw_w = _get_or(params, "W", (0.5, 1.0, 2.0, 4.0))
        w_list = (
            [float(x) for x in raw_w]
            if isinstance(raw_w, (list, tuple))
            else [float(raw_w)]
        )
        return run_tutorial05(
            out_dir,
            l_side=int(_get_or(params, "L", 200)),
            w_list=w_list,
            disorder_seed=int(_get_or(params, "disorder_seed", 1)),
            dos_moments=int(_get_or(params, "dos_moments", 1024)),
            kpm_seed=int(_get_or(rf.expansion, "seed", 1)),
            num_moments=int(rf.expansion["moments"]),
            num_states=int(_get_or(rf.expansion, "states", 1)),
            window=rf.expansion.get("window"),
            ntime=int(rf.time["steps"]),
            tmax=float(_get_or(rf.time, "tmax", 150.0)),
            bessel=str(_get_or(rf.time, "bessel", "legacy")),
            words=rf.words,
            broadening_meV=float(
                _get_or(rf.reconstruction, "broadening_meV", 80.0)
            ),
            fermi_energies=(
                rf.reconstruction["fermi_energies"]
                if rf.reconstruction.get("fermi_energies")
                else (2.0,)
            ),
            legacy_text_round_trip=bool(
                _get_or(rf.reconstruction, "legacy_text_round_trip", True)
            ),
        )
    if tut == 6:
        from .apps.tutorial06 import run_tutorial06

        assert rf.time is not None
        ey_steps_val = int(rf.time["steps"])
        ey_tmax_raw = rf.time.get("tmax")
        if ey_tmax_raw is None and rf.time.get("dt") is not None:
            ey_tmax_val = float(rf.time["dt"]) * (ey_steps_val - 1)
        else:
            ey_tmax_val = float(_get_or(rf.time, "tmax", 260.0))
        return run_tutorial06(
            out_dir,
            ey_l=int(_get_or(params, "ey_l", 20)),
            ey_ds=_get_or(params, "ey_ds", (0.05, 0.07, 0.10)),
            ey_moments=int(rf.expansion["moments"]),
            ey_steps=ey_steps_val,
            ey_tmax=ey_tmax_val,
            dp_l=int(_get_or(params, "dp_l", 20)),
            dp_lambda=float(_get_or(params, "dp_lambda", 0.2)),
            dp_ws=_get_or(params, "dp_ws", (0.0, 6.0)),
            dp_moments=int(_get_or(params, "dp_moments", 160)),
            dp_steps=int(_get_or(params, "dp_steps", 120)),
            dp_tmax=float(_get_or(params, "dp_tmax", 45.0)),
            ratelaw_l=int(_get_or(params, "ratelaw_l", 50)),
            ratelaw_ds=_get_or(
                params, "ratelaw_ds", (0.04, 0.06, 0.08, 0.10, 0.14)
            ),
            ratelaw_nr=int(_get_or(params, "ratelaw_nr", 20)),
            ratelaw_moments=int(rf.expansion["moments"]),
            ratelaw_steps=ey_steps_val,
            ratelaw_tmax=ey_tmax_val,
            crosscheck_l=int(_get_or(params, "crosscheck_l", 8)),
            crosscheck_d=float(_get_or(params, "crosscheck_d", 0.10)),
            crosscheck_nr=int(_get_or(params, "crosscheck_nr", 20)),
            crosscheck_moments=int(_get_or(params, "crosscheck_moments", 96)),
            crosscheck_steps=int(_get_or(params, "crosscheck_steps", 60)),
            crosscheck_tmax=float(_get_or(params, "crosscheck_tmax", 220.0)),
            disorder_seed=int(_get_or(params, "disorder_seed", 1)),
            window=rf.expansion.get("window"),
            trace=str(_get_or(rf.expansion, "trace", "exact")),
            num_states=rf.expansion.get("states"),
            kpm_seed=int(_get_or(rf.expansion, "seed", 12345)),
            kernel=str(_get_or(rf.reconstruction, "kernel", "jackson")),
            broadening_meV=float(
                _get_or(rf.reconstruction, "broadening_meV", 140.0)
            ),
            fermi_energies=(
                rf.reconstruction["fermi_energies"]
                if rf.reconstruction.get("fermi_energies")
                else (2.0,)
            ),
            formulas=rf.reconstruction.get("formulas"),
            bessel=str(_get_or(rf.time, "bessel", "legacy")),
            legacy_text_round_trip=bool(
                _get_or(rf.reconstruction, "legacy_text_round_trip", True)
            ),
            capture_tutorial=bool(
                _get_or(rf.output, "capture_tutorial", True)
            ),
            words=rf.words,
        )
    raise ValueError(f"Unsupported generator '{gen_name}' (tutorial={tut})")


def execute_runfile(
    runfile: RunFile | Mapping[str, Any] | str,
    *,
    out_dir_override: str | None = None,
    repo_root: str | None = None,
) -> dict[str, Any]:
    """Execute a :class:`RunFile` and write moments, curves, and ``run.resolved.json``."""
    if isinstance(runfile, str):
        rf = RunFile.from_file(runfile, repo_root=repo_root)
    elif isinstance(runfile, Mapping):
        rf = RunFile.from_dict(runfile, repo_root=repo_root)
    else:
        rf = runfile

    out_dir = os.path.abspath(
        out_dir_override if out_dir_override is not None else rf.output["dir"]
    )
    kind = str(rf.model["kind"])
    if kind == "siesta":
        if rf.projection["kind"] != "none":
            from .projection.runfile import execute_projection_runfile
            return execute_projection_runfile(rf, out_dir)
        if (rf.model["params"].get("nonorthogonal_disorder", False) or
                rf.basis["name"] == "real" or rf.basis["disorder"]["W"] != 0):
            from .apps.stage6c_disorder import execute_disorder_runfile
            return execute_disorder_runfile(rf, out_dir)
        from .apps.stage6b_siesta import execute_siesta_runfile
        return execute_siesta_runfile(rf, out_dir)
    if kind == "wannier":
        basis_name = str(rf.basis["name"])
        csr_dir_val = rf.model.get("csr_dir")
        if basis_name == "real" and csr_dir_val is not None:
            validate_w2s_csr_dir(
                str(csr_dir_val),
                str(rf.model["seed"]),
                tuple(int(x) for x in rf.basis["supercell"]),
                model_dir=str(rf.model["dir"]),
                operators=("VX",),
            )
    elif kind == "csr":
        model_dir = str(rf.model["dir"])
        label = str(rf.model["label"])
        ham_path = os.path.join(model_dir, f"{label}.HAM.CSR")
        if not os.path.isfile(ham_path):
            raise FileNotFoundError(
                f"Hamiltonian CSR file not found: {ham_path}"
            )
        if str(rf.metric["policy"]) != "orthogonal":
            s_path = os.path.join(model_dir, f"{label}.S.CSR")
            if not os.path.isfile(s_path):
                raise FileNotFoundError(
                    f"Overlap matrix file required for nonorthogonal metric not found: {s_path}"
                )

    os.makedirs(out_dir, exist_ok=True)

    resolved_dict = rf.to_resolved_dict()
    resolved_dict["output"] = dict(resolved_dict["output"], dir=out_dir)
    resolved_path = os.path.join(out_dir, "run.resolved.json")
    with open(resolved_path, "w", encoding="utf-8") as fh:
        json.dump(resolved_dict, fh, indent=2)
        fh.write("\n")

    if kind == "generator":
        written_files = _execute_generator_runfile(
            rf, out_dir, repo_root=repo_root
        )
        return {
            "output_dir": out_dir,
            "resolved_runfile": resolved_path,
            "moment_files": [
                p
                for p in written_files
                if p.endswith((".chebmom", ".chebmom1D", ".chebmom2D", ".chebmomTD"))
            ],
            "curve_files": [p for p in written_files if p.endswith(".dat")],
            "written_files": written_files,
        }

    # Build list of words to run from rf.words (no implicit companion words).
    words_to_run: list[tuple[str, list[str] | None]] = []
    req_formulas = rf.reconstruction.get("formulas")
    for w_str in rf.words:
        parsed_w = parse_word(w_str)
        comp_preview = compile_word(
            parsed_w,
            basis={"name": rf.basis["name"]},
            params={
                "expansion": rf.expansion,
                "time": rf.time,
                "metric": rf.metric,
            },
        )
        allowed_f = SHAPE_ALLOWED_FORMULAS.get(
            comp_preview.match.entry.shape_id, frozenset()
        )
        w_formulas = (
            [f for f in req_formulas if f.lower() in allowed_f]
            if req_formulas is not None
            else None
        )
        words_to_run.append((w_str, w_formulas))

    moment_files: list[str] = []
    curve_files: list[str] = []
    written_files: list[str] = []
    curves_json_payload: dict[str, Any] = {}

    if kind == "wannier":
        model_dir = str(rf.model["dir"])
        seed_name = str(rf.model["seed"])
        supercell = tuple(int(x) for x in rf.basis["supercell"])
        disorder_w = float(rf.basis["disorder"]["W"])
        disorder_seed = int(rf.basis["disorder"]["seed"])
        rng_seed = int(rf.expansion["seed"])
        basis_name = str(rf.basis["name"])

        wmodel = load_wannier_model(model_dir, seed=seed_name)
        basis_obj = make_basis(
            basis_name,
            wmodel,
            supercell,
            disorder_w=disorder_w,
            seed=rng_seed,
            disorder_seed=disorder_seed,
            csr_dir=rf.model.get("csr_dir"),
        )
        win_spec = rf.expansion["window"]
        if win_spec == "eigen":
            if isinstance(basis_obj, (BlochBasis, EigenBasis)):
                window = basis_obj.window()
            else:
                clean_win = EigenBasis(wmodel, supercell).window(pad=0.0)
                emin = float(clean_win.lower) - 0.5 * disorder_w
                emax = float(clean_win.upper) + 0.5 * disorder_w
                width = emax - emin
                if width <= 0.0:
                    emin -= 1.0
                    emax += 1.0
                    width = 2.0
                margin = 0.10 * width
                window = _core.SpectralWindow(emin - margin, emax + margin)
            window_source = "eigen"
        elif win_spec == "gershgorin":
            if isinstance(basis_obj, RealBasis):
                window = basis_obj.window()
            else:
                window = RealBasis(
                    wmodel,
                    supercell,
                    disorder_w=disorder_w,
                    seed=rng_seed,
                    disorder_seed=disorder_seed,
                ).window()
            window_source = "gershgorin"
        else:
            window = _core.SpectralWindow(float(win_spec[0]), float(win_spec[1]))
            window_source = "explicit"

        trace_mode = str(rf.expansion["trace"])
        state_gen = basis_obj.states(
            trace_mode=trace_mode,
            num_states=rf.expansion.get("states"),
            seed=rng_seed,
        )
        num_states = int(state_gen.number_of_states)
        num_moments = int(rf.expansion["moments"])
        batch_size = int(rf.expansion["batch"] or num_moments)

        concrete_ops: dict[str, _core.LinearOperator] = {}
        for op_name in ("VX", "VY", "VZ"):
            try:
                concrete_ops[op_name] = basis_obj.operator(op_name)
            except Exception:  # noqa: BLE001
                pass

        for w_str, w_formulas in words_to_run:
            compiled = compile_word(
                w_str,
                basis={
                    "name": basis_name,
                    "operators": concrete_ops,
                    "metric": rf.metric,
                },
                params={
                    "expansion": dict(
                        rf.expansion,
                        states=num_states,
                        window=window,
                    ),
                    "time": rf.time,
                    "metric": rf.metric,
                    "reconstruction": dict(
                        rf.reconstruction, formulas=w_formulas
                    ),
                },
            )
            symmetrise = (
                compiled.match.entry.shape_id == "kubo_greenwood"
                and len(compiled.operator_names) == 1
            )
            mu = run_plan(
                compiled.plan,
                basis_obj.hamiltonian(),
                window,
                state_gen,
                operators=compiled.operators,
                batch_size=batch_size,
                symmetrise=symmetrise,
                dt=compiled.dt,
            )
            if trace_mode == "exact":
                mu.number_of_states = basis_obj.dimension

            label_tag = f"{wmodel.seed_name}_{basis_name}_"
            prefix = compiled.reconstruction.legacy_moment_prefix
            ext = compiled.reconstruction.legacy_moment_ext
            if mu.rank == 2:
                stem = (
                    f"{prefix}{label_tag}KPM_M{num_moments}x{num_moments}"
                    f"_state{mu.state_label}"
                )
            else:
                stem = f"{prefix}{label_tag}KPM_M{num_moments}_state{mu.state_label}"

            if rf.output["legacy_layout"]:
                leg_path = os.path.join(out_dir, f"{stem}{ext}")
                if ext == ".chebmom1D":
                    _core.WriteLegacyChebmom1D(leg_path, mu)
                elif ext == ".chebmom2D":
                    _core.WriteLegacyChebmom2D(leg_path, mu)
                else:
                    _core.WriteLegacyChebmomTD(leg_path, mu)
                moment_files.append(leg_path)
                written_files.append(leg_path)

            if rf.output["write_chebmom"] or not rf.output["legacy_layout"]:
                cheb_path = os.path.join(out_dir, f"{stem}.chebmom")
                _core.WriteChebmom(cheb_path, mu)
                moment_files.append(cheb_path)
                written_files.append(cheb_path)

            broad = compiled.reconstruction.broadening_meV
            if broad is not None:
                for form in compiled.reconstruction.formulas:
                    for c_pref, c_suf, x_arr, y_arr in _reconstruct_single_formula(
                        form,
                        mu,
                        broad,
                        fermi_energies=compiled.reconstruction.fermi_energies,
                        num_points=compiled.reconstruction.num_points,
                    ):
                        if rf.output["legacy_names"]:
                            c_path = os.path.join(
                                out_dir, f"{c_pref}{stem}{c_suf}"
                            )
                            _core.WriteCurveText(c_path, x_arr, y_arr)
                            curve_files.append(c_path)
                            written_files.append(c_path)
                        else:
                            curves_json_payload[f"{stem}:{form}"] = {
                                "grid": [float(v) for v in x_arr],
                                "values": [float(v) for v in y_arr],
                            }

            if rf.output["write_basis_json"]:
                used_ops = ["HAM"] + list(compiled.operator_names)
                word_tag = (
                    "dos"
                    if compiled.match.entry.shape_id.startswith("dos")
                    else (
                        "kuboxx"
                        if len(compiled.operator_names) == 1
                        else "kuboxy"
                    )
                )
                basis_meta = {
                    "basis": basis_name,
                    "seed": wmodel.seed_name,
                    "supercell": list(basis_obj.supercell),
                    "mesh": list(basis_obj.supercell),
                    "dimension": int(basis_obj.dimension),
                    "trace": trace_mode,
                    "states": int(mu.number_of_states),
                    "state_seed": int(rng_seed),
                    "window_source": window_source,
                    "window": [float(window.lower), float(window.upper)],
                    "disorder_w": float(disorder_w),
                    "disorder_seed": (
                        int(disorder_seed) if disorder_w > 0.0 else None
                    ),
                    "moments": int(num_moments),
                    "broadening_meV": float(broad) if broad is not None else None,
                    "word": word_tag,
                    "operators": {
                        op_n: {"source": basis_obj.operator_source(op_n)}
                        for op_n in used_ops
                    },
                }
                bjson_path = os.path.join(out_dir, f"{stem}.basis.json")
                with open(bjson_path, "w", encoding="utf-8") as fh:
                    json.dump(basis_meta, fh, indent=2)
                    fh.write("\n")
                written_files.append(bjson_path)

    else:
        # kind == "csr"
        model_dir = str(rf.model["dir"])
        label = str(rf.model["label"])
        ham_path = os.path.join(model_dir, f"{label}.HAM.CSR")
        if not os.path.isfile(ham_path):
            raise FileNotFoundError(f"Hamiltonian CSR file not found: {ham_path}")
        ham = _core.ReadCsrText(ham_path)
        dim = int(ham.dimension)

        win_spec = rf.expansion["window"]
        if win_spec == "eigen":
            raise ValueError(
                "expansion.window='eigen' is not supported for model.kind='csr' "
                "(no k blocks; use gershgorin or an explicit window)"
            )
        elif win_spec == "gershgorin":
            window = window_from_gershgorin(ham)
        elif isinstance(win_spec, (list, tuple)):
            window = _core.SpectralWindow(float(win_spec[0]), float(win_spec[1]))
        else:
            window = window_from_gershgorin(ham)

        trace_mode = str(rf.expansion["trace"])
        rng_seed = int(rf.expansion["seed"])
        req_states = rf.expansion.get("states")
        num_states = (
            int(req_states)
            if req_states is not None
            else (dim if trace_mode == "exact" else 1)
        )
        if trace_mode == "exact":
            if num_states > dim:
                raise ValueError(
                    f"expansion.states ({num_states}) exceeds Hamiltonian dimension ({dim})"
                )
            state_gen: _core.StateGenerator = _core.LocalStates(
                dim, list(range(num_states))
            )
            state_label = str(rf.expansion.get("state_label") or "exact")
        else:
            state_gen = _core.RandomPhaseStates(dim, num_states, rng_seed)
            state_label = str(
                rf.expansion.get("state_label") or state_gen.label
            )

        # Metric policy setup (orthogonal or nonorthogonal)
        policy_name = str(rf.metric["policy"])
        route_name = str(rf.metric["route"])
        applier_name = str(rf.metric["applier"])
        metric_obj: _core.MetricPolicy
        s_mat: _core.SparseMatrix | None = None
        if policy_name == "orthogonal":
            metric_obj = _core.OrthogonalMetric(dim)
            rescaled_h: _core.LinearOperator = _core.RescaledOperator(
                ham, window
            )
        else:
            s_path = os.path.join(model_dir, f"{label}.S.CSR")
            if not os.path.isfile(s_path):
                raise FileNotFoundError(
                    f"Overlap matrix file required for nonorthogonal metric not found: {s_path}"
                )
            s_mat = _core.ReadCsrText(s_path)
            s_op = _core.SparseMatrixOperator(s_mat)
            order = int(rf.metric["order"])
            coeff_tol = float(rf.metric["coeff_tol"])
            tol = float(rf.metric["tol"])
            max_iter = int(rf.metric["max_iter"])
            max_dim = int(rf.metric["max_dim"])

            if route_name == "lowdin":
                inv_app: _core.LinearOperator = _core.IdentityOperator(dim)
                if applier_name == "exact":
                    inv_sqrt_app: _core.LinearOperator = (
                        _core.ExactInverseSqrtOverlap(s_mat, max_dim)
                    )
                else:
                    bounds_s = _compute_padded_overlap_bounds(s_op)
                    inv_sqrt_app = _core.ChebyshevInverseSqrtOverlap(
                        s_mat, bounds_s, order, coeff_tol
                    )
                metric_obj = _core.NonOrthogonalMetric(
                    s_mat, inv_app, inv_sqrt_app
                )
                ham_op = _core.SparseMatrixOperator(ham)
                lowdin_h = _core.LowdinOperator(
                    metric_obj.inverse_sqrt_overlap, ham_op
                )
                rescaled_h = _core.RescaledOperator(lowdin_h, window)
            else:
                if applier_name == "exact":
                    inv_app = _core.ExactInverseOverlap(s_mat, max_dim)
                elif applier_name == "ischo":
                    inv_app = _core.IschoInverseOverlap(s_mat, tol, max_iter)
                else:
                    bounds_s = _compute_padded_overlap_bounds(s_op)
                    inv_app = _core.ChebyshevInverseOverlap(
                        s_mat, bounds_s, order, coeff_tol
                    )
                inv_sqrt_app = _core.IdentityOperator(dim)
                metric_obj = _core.NonOrthogonalMetric(
                    s_mat, inv_app, inv_sqrt_app
                )
                rescaled_h = _core.GeneralisedRescaledOperator(
                    ham, s_mat, metric_obj.inverse_overlap, window
                )

        discovered_ops = (
            list(rf.model["operators"])
            if rf.model.get("operators") is not None
            else _discover_csr_operators(model_dir, label)
        )
        raw_csr_mats: dict[str, _core.SparseMatrix] = {}
        for op_n in discovered_ops:
            op_p = os.path.join(model_dir, f"{label}.{op_n}.CSR")
            if os.path.isfile(op_p):
                raw_csr_mats[op_n] = _core.ReadCsrText(op_p)

        num_moments = int(rf.expansion["moments"])
        batch_size = int(rf.expansion["batch"] or num_moments)
        captured_idx = 0

        for w_str, w_formulas in words_to_run:
            parsed_w = parse_word(w_str)
            # For MSD, the velocity operator passed to the engine is RescaledOperator(V, window)
            concrete_ops = {}
            for op_n, mat in raw_csr_mats.items():
                if any(
                    getattr(s, "velocity_operator", None) == op_n
                    for s in parsed_w.slots
                ):
                    concrete_ops[op_n] = _core.RescaledOperator(mat, window)
                else:
                    concrete_ops[op_n] = _core.SparseMatrixOperator(mat)
            concrete_ops["1"] = _core.IdentityOperator(dim)

            metric_with_obj = dict(rf.metric, metric_object=metric_obj)
            compiled = compile_word(
                parsed_w,
                basis={
                    "name": rf.basis["name"],
                    "operators": concrete_ops,
                    "metric": metric_with_obj,
                },
                params={
                    "expansion": dict(
                        rf.expansion,
                        states=num_states,
                        window=window,
                    ),
                    "time": rf.time,
                    "metric": metric_with_obj,
                    "reconstruction": dict(
                        rf.reconstruction, formulas=w_formulas
                    ),
                },
            )

            engine = _core.MomentEngine(
                rescaled_h, compiled.operators, state_gen, batch_size
            )
            mu = engine.run(compiled.plan)
            mu.window = window
            mu.system_size = dim
            mu.number_of_states = num_states
            mu.state_label = state_label
            mu.word = compiled.plan.word
            if compiled.dt is not None:
                mu.time_step = float(compiled.dt)
            if mu.rank == 1 or mu.number_of_time_steps > 0:
                mu.kpm_weighted = True
                mu.symmetrised = False
            elif mu.rank == 2 and not mu.kpm_weighted:
                _core.ApplyKpmWeights(mu, num_states, True)

            prefix = compiled.reconstruction.legacy_moment_prefix
            ext = compiled.reconstruction.legacy_moment_ext
            if mu.number_of_time_steps > 0:
                if compiled.match.entry.shape_id in ("msd", "time_correlation"):
                    assert rf.time is not None
                    tmax_val = rf.time.get("tmax")
                    if tmax_val is None and compiled.dt is not None:
                        tmax_val = compiled.dt * (int(rf.time["steps"]) - 1)
                    tmax_str = f"{float(tmax_val):g}"
                    stem = (
                        f"{prefix}{label}KPM_M{num_moments}"
                        f"_state{tmax_str}{state_label}"
                    )
                else:
                    stem = (
                        f"{prefix}{label}KPM_M{num_moments}_state{state_label}"
                    )
            elif mu.rank == 2:
                stem = (
                    f"{prefix}{label}KPM_M{num_moments}x{num_moments}"
                    f"_state{state_label}"
                )
            else:
                stem = f"{prefix}{label}KPM_M{num_moments}_state{state_label}"

            leg_path: str | None = None
            if (
                rf.output["legacy_layout"]
                or compiled.reconstruction.legacy_text_round_trip
            ):
                leg_path = os.path.join(out_dir, f"{stem}{ext}")
                if ext == ".chebmom1D":
                    _core.WriteLegacyChebmom1D(leg_path, mu)
                elif ext == ".chebmom2D":
                    _core.WriteLegacyChebmom2D(leg_path, mu)
                else:
                    _core.WriteLegacyChebmomTD(leg_path, mu)
                moment_files.append(leg_path)
                written_files.append(leg_path)

                if (
                    rf.output["capture_tutorial"]
                    and ext == ".chebmomTD"
                    and compiled.match.entry.shape_id == "projected_evolution"
                ):
                    cap_dir = os.path.join(out_dir, "captured")
                    os.makedirs(cap_dir, exist_ok=True)
                    cap_p = os.path.join(
                        cap_dir,
                        f"{captured_idx:04d}_{os.path.basename(leg_path)}",
                    )
                    shutil.copyfile(leg_path, cap_p)
                    written_files.append(cap_p)
                    captured_idx += 1

            if rf.output["write_chebmom"] or not rf.output["legacy_layout"]:
                cheb_path = os.path.join(out_dir, f"{stem}.chebmom")
                _core.WriteChebmom(cheb_path, mu)
                moment_files.append(cheb_path)
                written_files.append(cheb_path)

            mu_for_recon = (
                _core.ReadChebmom(leg_path)
                if (
                    compiled.reconstruction.legacy_text_round_trip
                    and leg_path is not None
                )
                else mu
            )

            broad = compiled.reconstruction.broadening_meV
            if broad is not None:
                for form in compiled.reconstruction.formulas:
                    for c_pref, c_suf, x_arr, y_arr in _reconstruct_single_formula(
                        form,
                        mu_for_recon,
                        broad,
                        fermi_energies=compiled.reconstruction.fermi_energies,
                        num_points=compiled.reconstruction.num_points,
                    ):
                        if rf.output["legacy_names"]:
                            c_path = os.path.join(
                                out_dir, f"{c_pref}{stem}{c_suf}"
                            )
                            _core.WriteCurveText(c_path, x_arr, y_arr)
                            curve_files.append(c_path)
                            written_files.append(c_path)
                            if rf.output["capture_tutorial"]:
                                cap_dir = os.path.join(out_dir, "captured")
                                os.makedirs(cap_dir, exist_ok=True)
                                cap_p = os.path.join(
                                    cap_dir,
                                    f"{captured_idx:04d}_{os.path.basename(c_path)}",
                                )
                                shutil.copyfile(c_path, cap_p)
                                written_files.append(cap_p)
                                captured_idx += 1
                        else:
                            curves_json_payload[f"{stem}:{form}"] = {
                                "grid": [float(v) for v in x_arr],
                                "values": [float(v) for v in y_arr],
                            }

            if rf.reconstruction.get("per_site_dos") and mu.rank == 1:
                # Reconstruct per-site DOS on 30 * M points with full M moments (broadening = 1/M)
                full_m_meV = 1.0e-6
                e_ps, v_ps = dos_curve(mu_for_recon, full_m_meV)
                ps_path = os.path.join(out_dir, "DOS_reconstructed.dat")
                tmp_ps_path = os.path.join(out_dir, ".tmp_dos_standalone.dat")
                _core.WriteCurveText(tmp_ps_path, e_ps, v_ps)
                with (
                    open(tmp_ps_path, encoding="utf-8") as fin,
                    open(ps_path, "w", encoding="utf-8") as fout,
                ):
                    for line in fin:
                        parts = line.split()
                        if len(parts) >= 2 and all(
                            tok.lower() not in ("nan", "inf", "-inf", "+inf")
                            for tok in parts[:2]
                        ):
                            c1 = float(parts[0])
                            c2 = float(parts[1]) / float(dim)
                            fout.write(f"{c1:.10g} {c2:.10g}\n")
                os.remove(tmp_ps_path)
                curve_files.append(ps_path)
                written_files.append(ps_path)

            for tgt in rf.reconstruction.get("targets") or ():
                t_broad = float(tgt["broadening_meV"])
                t_formulas = [str(f) for f in tgt["formulas"]]
                t_fermi = [
                    float(e) for e in (tgt.get("fermi_energies") or ())
                ]
                t_rt = bool(
                    _get_or(
                        tgt,
                        "legacy_text_round_trip",
                        compiled.reconstruction.legacy_text_round_trip,
                    )
                )
                mu_tgt = (
                    _core.ReadChebmom(leg_path)
                    if (t_rt and leg_path is not None)
                    else mu
                )
                sub = tgt.get("subdir")
                dest_d = (
                    os.path.join(out_dir, str(sub)) if sub else out_dir
                )
                os.makedirs(dest_d, exist_ok=True)
                also_top = bool(_get_or(tgt, "also_top_dir", False))
                for form in t_formulas:
                    for c_pref, c_suf, x_arr, y_arr in _reconstruct_single_formula(
                        form,
                        mu_tgt,
                        t_broad,
                        fermi_energies=t_fermi,
                    ):
                        fname = f"{c_pref}{stem}{c_suf}"
                        c_path = os.path.join(dest_d, fname)
                        _core.WriteCurveText(c_path, x_arr, y_arr)
                        curve_files.append(c_path)
                        written_files.append(c_path)
                        if also_top and dest_d != out_dir:
                            top_p = os.path.join(out_dir, fname)
                            if not os.path.exists(top_p):
                                shutil.copyfile(c_path, top_p)
                                written_files.append(top_p)

    if not rf.output["legacy_names"] and curves_json_payload:
        cjson_path = os.path.join(out_dir, "curves.json")
        with open(cjson_path, "w", encoding="utf-8") as fh:
            json.dump(curves_json_payload, fh, indent=2)
            fh.write("\n")
        curve_files.append(cjson_path)
        written_files.append(cjson_path)

    return {
        "output_dir": out_dir,
        "resolved_runfile": resolved_path,
        "moment_files": moment_files,
        "curve_files": curve_files,
        "written_files": written_files,
    }


def reconstruct_target(
    target_path: str,
    *,
    formula: str | None = None,
    broadening_meV: float | None = None,
    fermi_energies: Sequence[float] | None = None,
    out_dir: str | None = None,
) -> dict[str, Any]:
    """Reconstruct curve(s) from a moment file or an output directory containing moment files."""
    abs_target = os.path.abspath(target_path)
    if os.path.isdir(abs_target):
        dest_dir = os.path.abspath(out_dir) if out_dir else abs_target
        mom_candidates = sorted(
            glob.glob(os.path.join(abs_target, "*.chebmom1D"))
            + glob.glob(os.path.join(abs_target, "*.chebmom2D"))
            + glob.glob(os.path.join(abs_target, "*.chebmomTD"))
        )
        if not mom_candidates:
            mom_candidates = sorted(
                glob.glob(os.path.join(abs_target, "*.chebmom"))
            )
        if not mom_candidates:
            raise FileNotFoundError(
                f"No moment files found in directory '{abs_target}'"
            )
        resolved_json = os.path.join(abs_target, "run.resolved.json")
        rf_rec: dict[str, Any] = {}
        if os.path.isfile(resolved_json):
            with open(resolved_json, encoding="utf-8") as fh:
                rf_rec = (json.load(fh).get("reconstruction") or {})
        eff_broad = (
            float(broadening_meV)
            if broadening_meV is not None
            else float(rf_rec.get("broadening_meV") or 10.0)
        )
        eff_fermi = (
            list(fermi_energies)
            if fermi_energies is not None
            else list(rf_rec.get("fermi_energies") or [0.0])
        )
        written: list[str] = []
        for mp in mom_candidates:
            res = reconstruct_target(
                mp,
                formula=formula,
                broadening_meV=eff_broad,
                fermi_energies=eff_fermi,
                out_dir=dest_dir,
            )
            written.extend(res["curve_files"])
        return {"output_dir": dest_dir, "curve_files": written}

    if not os.path.isfile(abs_target):
        raise FileNotFoundError(f"Moment file not found: '{abs_target}'")

    dest_dir = (
        os.path.abspath(out_dir)
        if out_dir is not None
        else os.path.dirname(abs_target)
    )
    os.makedirs(dest_dir, exist_ok=True)
    mu = _core.ReadChebmom(abs_target)
    base = os.path.basename(abs_target)
    stem, _ = os.path.splitext(base)

    eff_broad = float(broadening_meV) if broadening_meV is not None else 10.0
    eff_fermi = (
        [float(x) for x in fermi_energies]
        if fermi_energies is not None
        else [0.0]
    )
    if formula is None:
        if mu.number_of_time_steps > 0:
            formula = "time_correlation"
        elif mu.rank == 2:
            formula = "bastin"
        else:
            formula = "dos"

    shc_metadata = {}
    resolved_path = os.path.join(os.path.dirname(abs_target), "run.resolved.json")
    if formula == "shc" and os.path.isfile(resolved_path):
        with open(resolved_path, encoding="utf-8") as fh:
            reconstruction = json.load(fh).get("reconstruction", {})
        shc_metadata = {key: reconstruction.get(key) for key in ("num_cells", "cell_area")}
    written_curves: list[str] = []
    for c_pref, c_suf, x_arr, y_arr in _reconstruct_single_formula(
        formula, mu, eff_broad, fermi_energies=eff_fermi, **shc_metadata
    ):
        c_path = os.path.join(dest_dir, f"{c_pref}{stem}{c_suf}")
        if formula == "shc":
            np.savetxt(c_path, np.column_stack((x_arr, y_arr)), fmt="%.17g",
                       header="energy_eV sigma_z_xy [e^2/h, sigma_z Pauli]")
        else:
            _core.WriteCurveText(c_path, x_arr, y_arr)
        written_curves.append(c_path)

    return {
        "moment_file": abs_target,
        "output_dir": dest_dir,
        "formula": formula,
        "broadening_meV": eff_broad,
        "fermi_energies": eff_fermi,
        "curve_files": written_curves,
    }


def inspect_target(path: str) -> dict[str, Any]:
    """Inspect a moment file, a run-file JSON, or a run output directory."""
    abs_path = os.path.abspath(path)
    if not os.path.exists(abs_path):
        raise FileNotFoundError(f"Path not found: '{path}'")

    if os.path.isdir(abs_path):
        resolved_path = os.path.join(abs_path, "run.resolved.json")
        resolved_data = None
        if os.path.isfile(resolved_path):
            with open(resolved_path, encoding="utf-8") as fh:
                resolved_data = json.load(fh)
        mom_files = sorted(
            glob.glob(os.path.join(abs_path, "*.chebmom*"))
        )
        moments_info = [
            inspect_target(mp) for mp in mom_files if os.path.isfile(mp)
        ]
        dir_basis = str(
            _get_or((resolved_data or {}).get("basis") or {}, "name", "real")
        )
        return {
            "kind": "directory",
            "path": abs_path,
            "resolved_runfile": resolved_data,
            "defaults": defaults_for_basis(dir_basis),
            "moment_files": moments_info,
        }

    if abs_path.endswith(".json"):
        with open(abs_path, encoding="utf-8") as fh:
            raw = json.load(fh)
        rf = RunFile.from_dict(raw, source_path=abs_path)
        return {
            "kind": "runfile",
            "path": abs_path,
            "resolved_runfile": rf.to_resolved_dict(),
            "defaults": defaults_for_basis(rf.basis["name"]),
        }

    mu = _core.ReadChebmom(abs_path)
    half_width = 0.5 * float(mu.window.upper - mu.window.lower)
    band_center = 0.5 * float(mu.window.upper + mu.window.lower)
    return {
        "kind": "moments",
        "path": abs_path,
        "rank": int(mu.rank),
        "dims": [int(d) for d in mu.dims],
        "system_size": int(mu.system_size),
        "number_of_states": int(mu.number_of_states),
        "state_label": str(mu.state_label),
        "word": str(mu.word),
        "window": [float(mu.window.lower), float(mu.window.upper)],
        "band_width": float(mu.window.upper - mu.window.lower),
        "half_width": half_width,
        "band_center": band_center,
        "kpm_weighted": bool(mu.kpm_weighted),
        "symmetrised": bool(mu.symmetrised),
        "number_of_time_steps": int(mu.number_of_time_steps),
        "time_step": float(mu.time_step),
    }
