"""Compiler from parsed observable :class:`~lsquant.words.grammar.Word` to bytecode :class:`~lsquant._core.Plan`."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterator, Mapping, Sequence

from .. import _core
from ..moments import (
    dos_plan,
    kubo_plan,
    msd_plan,
    projected_evolution_plan,
    time_correlation_plan,
    time_evolution_series,
    time_step_from_tmax,
)
from .catalogue import (
    CATALOGUE,
    CatalogueMatch,
    SHAPE_ALLOWED_FORMULAS,
    match_catalogue,
    validate_word_bases,
)
from .grammar import OperatorSlot, Word, parse


DEFAULT_KNOWN_OPERATORS: frozenset[str] = frozenset(
    {
        "1",
        "VX",
        "VY",
        "VZ",
        "SX",
        "SY",
        "SZ",
        "JXSX",
        "JXSY",
        "JXSZ",
        "JYSX",
        "JYSY",
        "JYSZ",
        "JZSX",
        "JZSY",
        "JZSZ",
        "VXSX",
        "VXSY",
        "VXSZ",
        "VYSX",
        "VYSY",
        "VYSZ",
        "VZSX",
        "VZSY",
        "VZSZ",
        "P",
        "PX",
        "PY",
        "PZ",
        "A",
        "B",
        "O",
    }
)


@dataclass(frozen=True)
class ReconstructionSpec:
    """Specification of post-moment reconstruction emitted by :func:`compile`."""

    shape_id: str
    default_reconstruction: str
    formulas: tuple[str, ...]
    kernel: str = "jackson"
    broadening_meV: float | None = None
    fermi_energies: tuple[float, ...] = ()
    num_points: int | None = None
    use_fft_grid: bool = False
    legacy_moment_prefix: str = ""
    legacy_moment_ext: str = ".chebmom1D"
    legacy_text_round_trip: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "shape_id": self.shape_id,
            "default_reconstruction": self.default_reconstruction,
            "formulas": list(self.formulas),
            "kernel": self.kernel,
            "broadening_meV": self.broadening_meV,
            "fermi_energies": list(self.fermi_energies),
            "num_points": self.num_points,
            "use_fft_grid": self.use_fft_grid,
            "legacy_moment_prefix": self.legacy_moment_prefix,
            "legacy_moment_ext": self.legacy_moment_ext,
            "legacy_text_round_trip": self.legacy_text_round_trip,
        }


@dataclass
class CompiledWord:
    """Result of compiling a :class:`Word` for a given basis and parameter set.

    Supports both attribute access and 3-tuple unpacking
    ``plan, operators, reconstruction = compile(word, basis, params)``.
    """

    word: Word
    match: CatalogueMatch
    plan: _core.Plan
    operator_names: tuple[str, ...]
    operators: list[_core.LinearOperator]
    reconstruction: ReconstructionSpec
    dt: float | None = None

    def __iter__(self) -> Iterator[Any]:
        yield self.plan
        yield self.operators if self.operators else list(self.operator_names)
        yield self.reconstruction

    def to_dict(self) -> dict[str, object]:
        return {
            "word": self.word.to_dict(),
            "catalogue": self.match.to_dict(),
            "plan": {
                "word": self.plan.word,
                "instructions": format_plan_instructions(
                    self.plan, self.operator_names
                ),
            },
            "operators": list(self.operator_names),
            "reconstruction": self.reconstruction.to_dict(),
            "dt": self.dt,
        }


def format_plan_instructions(
    plan: _core.Plan, operator_names: Sequence[str] = ()
) -> list[dict[str, object]]:
    """Serialize bytecode instructions in ``plan.code`` to JSON-friendly dicts."""
    items: list[dict[str, object]] = []
    for instr in plan.code:
        if isinstance(instr, _core.Expand):
            items.append(
                {"op": "Expand", "num_moments": int(instr.num_moments)}
            )
        elif isinstance(instr, _core.ApplyOperator):
            idx = int(instr.operator_index)
            name = (
                operator_names[idx]
                if 0 <= idx < len(operator_names)
                else f"op[{idx}]"
            )
            items.append(
                {
                    "op": "ApplyOperator",
                    "operator_index": idx,
                    "operator_name": name,
                }
            )
        elif isinstance(instr, _core.MirrorToBra):
            items.append({"op": "MirrorToBra"})
        elif isinstance(instr, _core.ApplySeries):
            items.append(
                {
                    "op": "ApplySeries",
                    "num_coeffs": len(instr.coeffs),
                    "scale": [
                        float(instr.scale.real),
                        float(instr.scale.imag),
                    ],
                }
            )
        elif isinstance(instr, _core.ApplyDerivedSeries):
            idx = int(instr.operator_index)
            name = (
                operator_names[idx]
                if 0 <= idx < len(operator_names)
                else f"op[{idx}]"
            )
            items.append(
                {
                    "op": "ApplyDerivedSeries",
                    "operator_index": idx,
                    "operator_name": name,
                    "num_coeffs": len(instr.coeffs),
                    "scale": [
                        float(instr.scale.real),
                        float(instr.scale.imag),
                    ],
                }
            )
        elif isinstance(instr, _core.TimeAxis):
            steps_desc: list[dict[str, object]] = []
            for ev in instr.evolve:
                if isinstance(ev, _core.EvolveRegister):
                    steps_desc.append(
                        {
                            "step": "EvolveRegister",
                            "target": str(ev.target).split(".")[-1],
                            "num_coeffs": len(ev.coeffs),
                        }
                    )
                elif isinstance(ev, _core.EvolveDerived):
                    idx = int(ev.operator_index)
                    name = (
                        operator_names[idx]
                        if 0 <= idx < len(operator_names)
                        else f"op[{idx}]"
                    )
                    steps_desc.append(
                        {
                            "step": "EvolveDerived",
                            "target": str(ev.target).split(".")[-1],
                            "source": str(ev.source).split(".")[-1],
                            "operator_index": idx,
                            "operator_name": name,
                            "num_coeffs": len(ev.coeffs),
                        }
                    )
                elif isinstance(ev, _core.CopyRegister):
                    steps_desc.append(
                        {
                            "step": "CopyRegister",
                            "from": str(ev.from_reg).split(".")[-1],
                            "to": str(ev.to_reg).split(".")[-1],
                        }
                    )
            items.append(
                {
                    "op": "TimeAxis",
                    "num_steps": int(instr.num_steps),
                    "evolve": steps_desc,
                }
            )
        elif isinstance(instr, _core.Trace):
            items.append(
                {
                    "op": "Trace",
                    "num_weights": len(instr.weights),
                }
            )
        else:
            items.append({"op": type(instr).__name__})
    return items


def _build_generalised_kubo_plan(
    num_moments: int,
    *,
    sinv_left_index: int,
    sinv_right_index: int,
    sinv_index: int,
    s_index: int,
    word: str,
) -> _core.Plan:
    """Build the non-orthogonal generalised-route Kubo bytecode plan."""
    code = [
        _core.Expand(num_moments),
        _core.ApplyOperator(sinv_right_index),
        _core.MirrorToBra(),
        _core.ApplyOperator(sinv_index),
        _core.ApplyOperator(sinv_left_index),
        _core.Expand(num_moments),
        _core.ApplyOperator(s_index),
        _core.Trace([complex(1.0, 0.0)]),
    ]
    return _core.Plan(code, word)


def _get_or(mapping: Mapping[str, Any], key: str, default: Any = None) -> Any:
    """Return ``mapping[key]`` when present and not ``None``, otherwise ``default``."""
    val = mapping.get(key)
    return default if val is None else val


def _extract_basis_info(
    basis: Any,
) -> tuple[
    str,
    set[str] | None,
    dict[str, _core.LinearOperator],
    dict[str, str],
    dict[str, Any],
]:
    """Extract basis name, allowed operator names, concrete operators, operator bases, and metric info."""
    if basis is None:
        return ("real", None, {}, {}, {"policy": "orthogonal", "route": "lowdin"})
    if isinstance(basis, str):
        return (basis, None, {}, {}, {"policy": "orthogonal", "route": "lowdin"})
    if isinstance(basis, Mapping):
        bname = str(_get_or(basis, "name", "real"))
        op_map: dict[str, _core.LinearOperator] = {}
        allowed: set[str] | None = None
        raw_ops = basis.get("operators")
        if isinstance(raw_ops, Mapping):
            allowed = {"1"} | {str(k) for k in raw_ops.keys()}
            for k, v in raw_ops.items():
                if isinstance(v, _core.LinearOperator):
                    op_map[str(k)] = v
        elif isinstance(raw_ops, (list, tuple, set, frozenset)):
            allowed = {"1"} | {str(x) for x in raw_ops}
        op_bases = {
            str(k): str(v)
            for k, v in (basis.get("operator_bases") or {}).items()
        }
        metric = dict(
            basis.get("metric") or {"policy": "orthogonal", "route": "lowdin"}
        )
        return (bname, allowed, op_map, op_bases, metric)

    # Support custom BasisContext-like objects
    bname = str(getattr(basis, "name", "real"))
    raw_ops = getattr(basis, "operators", None)
    op_map = {}
    allowed = None
    if isinstance(raw_ops, Mapping):
        allowed = {"1"} | {str(k) for k in raw_ops.keys()}
        for k, v in raw_ops.items():
            if isinstance(v, _core.LinearOperator):
                op_map[str(k)] = v
    op_bases = dict(getattr(basis, "operator_bases", {}) or {})
    metric = dict(
        getattr(basis, "metric", None)
        or {"policy": "orthogonal", "route": "lowdin"}
    )
    return (bname, allowed, op_map, op_bases, metric)


def compile(
    word: str | Word,
    basis: Any = None,
    params: Mapping[str, Any] | None = None,
    *,
    _allow_default_time: bool = False,
) -> CompiledWord:
    """Compile an observable word into a :class:`CompiledWord` ``(Plan, operators, reconstruction)``.

    Parameters
    ----------
    word : str or Word
        Observable word string or parsed :class:`Word` AST.
    basis : str, Mapping, or None, optional
        Basis specification (`"real"`, `"bloch"`, `"eigen"`, or dict with
        `"name"`, `"operators"`, `"operator_bases"`, `"metric"`).
    params : Mapping or None, optional
        Compilation parameters (`moments`, `states`, `window`, `time`,
        `metric`, `reconstruction`).

    Returns
    -------
    CompiledWord
        Compiled plan, ordered operator list, and reconstruction spec.
    """
    parsed = parse(word) if isinstance(word, str) else word
    basis_name, allowed_ops, concrete_ops, op_bases, basis_metric = (
        _extract_basis_info(basis)
    )

    if basis_name not in ("real", "bloch", "eigen"):
        raise ValueError(
            f"Unknown basis '{basis_name}'; expected 'real', 'bloch', or 'eigen'"
        )

    validate_word_bases(
        parsed, calculation_basis=basis_name, operator_bases=op_bases
    )
    match = match_catalogue(parsed)

    cfg = dict(params or {})
    exp_cfg = (
        dict(cfg["expansion"])
        if isinstance(cfg.get("expansion"), Mapping)
        else cfg
    )
    time_cfg = (
        dict(cfg["time"]) if isinstance(cfg.get("time"), Mapping) else None
    )
    if time_cfg is None and ("steps" in cfg or "tmax" in cfg or "dt" in cfg):
        time_cfg = {
            k: cfg[k] for k in ("steps", "tmax", "dt", "bessel") if k in cfg
        }
    metric_cfg = (
        dict(cfg["metric"])
        if isinstance(cfg.get("metric"), Mapping)
        else basis_metric
    )
    recon_cfg = (
        dict(cfg["reconstruction"])
        if isinstance(cfg.get("reconstruction"), Mapping)
        else {}
    )

    # Validate static operators against known/registered operators
    known_ops = (
        allowed_ops if allowed_ops is not None else set(DEFAULT_KNOWN_OPERATORS)
    )
    projector_aliases: dict[str, str] = {}
    if "P" not in known_ops and "PX" in known_ops:
        projector_aliases["P"] = "PX"

    for op_slot in match.operators:
        key = op_slot.operator_key
        resolved_key = projector_aliases.get(key, key)
        if resolved_key not in known_ops:
            raise ValueError(
                f"Unknown operator '{key}' in word '{parsed.canonical()}' "
                f"(registered operators: {sorted(known_ops)})"
            )
        if resolved_key in concrete_ops:
            op_obj = concrete_ops[resolved_key]
            if not op_obj.is_hermitian:
                raise ValueError(
                    f"Operator '{resolved_key}' in word '{parsed.canonical()}' is not Hermitian"
                )

    num_moments = int(_get_or(exp_cfg, "moments", _get_or(cfg, "moments", 64)))
    num_states = int(_get_or(exp_cfg, "states", _get_or(cfg, "states", 1)))
    if num_moments <= 0:
        raise ValueError("compile: moments must be positive")
    if num_states <= 0:
        raise ValueError("compile: states must be positive")

    win_raw = _get_or(exp_cfg, "window", cfg.get("window"))
    if isinstance(win_raw, _core.SpectralWindow):
        window = win_raw
    elif (
        isinstance(win_raw, (list, tuple))
        and len(win_raw) == 2
        and isinstance(win_raw[0], (int, float))
    ):
        window = _core.SpectralWindow(float(win_raw[0]), float(win_raw[1]))
    else:
        window = _core.SpectralWindow(-1.0, 1.0)

    # Check time parameters for time-dependent words
    dt: float | None = None
    num_steps: int = 0
    use_legacy_bessel: bool = False
    if parsed.is_time_dependent:
        if time_cfg is None or not time_cfg:
            if not _allow_default_time:
                raise ValueError(
                    f"Word '{parsed.canonical()}' contains time evolution U(t) "
                    "without 'time' parameters (steps and tmax required)"
                )
            num_steps = 2
            dt = 1.0
        else:
            if _get_or(time_cfg, "steps") is None:
                raise ValueError(
                    f"Word '{parsed.canonical()}' requires 'steps' in 'time' parameters"
                )
            num_steps = int(time_cfg["steps"])
            if "dt" in time_cfg and time_cfg["dt"] is not None:
                dt = float(time_cfg["dt"])
            elif "tmax" in time_cfg and time_cfg["tmax"] is not None:
                dt = time_step_from_tmax(num_steps, float(time_cfg["tmax"]))
            else:
                raise ValueError(
                    f"Word '{parsed.canonical()}' requires 'tmax' or 'dt' in 'time' parameters"
                )
            bessel_mode = str(_get_or(time_cfg, "bessel", "scipy"))
            if bessel_mode not in ("scipy", "legacy"):
                raise ValueError(
                    f"Unknown time.bessel mode '{bessel_mode}'; expected 'scipy' or 'legacy'"
                )
            use_legacy_bessel = bessel_mode == "legacy"

    policy = str(_get_or(metric_cfg, "policy", "orthogonal"))
    route = str(_get_or(metric_cfg, "route", "lowdin"))
    metric_policy_obj = metric_cfg.get("metric_object")

    shape_id = match.entry.shape_id
    op_names: list[str] = []
    op_objs: list[_core.LinearOperator] = []

    def _resolve_op_key(slot: OperatorSlot) -> str:
        k = slot.operator_key
        return projector_aliases.get(k, k)

    if shape_id == "dos_identity":
        plan = dos_plan(num_moments, num_states, match.engine_word)

    elif shape_id == "dos_operator":
        if not match.operators:
            plan = dos_plan(num_moments, num_states, match.engine_word)
        else:
            key_a = _resolve_op_key(match.operators[0])
            op_names = [key_a]
            if key_a in concrete_ops:
                raw_a = concrete_ops[key_a]
                if (
                    policy == "nonorthogonal"
                    and metric_policy_obj is not None
                    and not metric_policy_obj.is_orthogonal
                ):
                    if route == "lowdin":
                        op_objs.append(
                            _core.LowdinOperator(
                                metric_policy_obj.inverse_sqrt_overlap, raw_a
                            )
                        )
                    else:
                        op_objs.append(
                            _core.MetricInverseProduct(
                                metric_policy_obj.inverse_overlap, raw_a
                            )
                        )
                else:
                    op_objs.append(raw_a)
            plan = dos_plan(
                num_moments,
                num_states,
                match.engine_word,
                operator_index=0,
            )

    elif shape_id in ("kubo_greenwood", "kubo_bastin"):
        key_a = _resolve_op_key(match.operators[0])
        key_b = _resolve_op_key(match.operators[1])

        if policy == "nonorthogonal" and route == "generalised":
            if key_a == key_b:
                op_names = [f"S^-1*{key_a}", "S^-1", "S"]
                sinv_left_idx = 0
                sinv_right_idx = 0
                sinv_idx = 1
                s_idx = 2
                if (
                    key_a in concrete_ops
                    and metric_policy_obj is not None
                    and not metric_policy_obj.is_orthogonal
                ):
                    sinv_a = _core.MetricInverseProduct(
                        metric_policy_obj.inverse_overlap, concrete_ops[key_a]
                    )
                    op_objs = [
                        sinv_a,
                        metric_policy_obj.inverse_overlap,
                        metric_policy_obj.overlap,
                    ]
            else:
                op_names = [f"S^-1*{key_a}", f"S^-1*{key_b}", "S^-1", "S"]
                sinv_left_idx = 0
                sinv_right_idx = 1
                sinv_idx = 2
                s_idx = 3
                if (
                    key_a in concrete_ops
                    and key_b in concrete_ops
                    and metric_policy_obj is not None
                    and not metric_policy_obj.is_orthogonal
                ):
                    sinv_a = _core.MetricInverseProduct(
                        metric_policy_obj.inverse_overlap, concrete_ops[key_a]
                    )
                    sinv_b = _core.MetricInverseProduct(
                        metric_policy_obj.inverse_overlap, concrete_ops[key_b]
                    )
                    op_objs = [
                        sinv_a,
                        sinv_b,
                        metric_policy_obj.inverse_overlap,
                        metric_policy_obj.overlap,
                    ]
            plan = _build_generalised_kubo_plan(
                num_moments,
                sinv_left_index=sinv_left_idx,
                sinv_right_index=sinv_right_idx,
                sinv_index=sinv_idx,
                s_index=s_idx,
                word=match.engine_word,
            )
        else:
            if key_a == key_b:
                op_names = [key_a]
                left_idx = 0
                right_idx = 0
            else:
                op_names = [key_a, key_b]
                left_idx = 0
                right_idx = 1

            if all(k in concrete_ops for k in op_names):
                for k in op_names:
                    raw_op = concrete_ops[k]
                    if (
                        policy == "nonorthogonal"
                        and metric_policy_obj is not None
                        and not metric_policy_obj.is_orthogonal
                    ):
                        if route == "lowdin":
                            op_objs.append(
                                _core.LowdinOperator(
                                    metric_policy_obj.inverse_sqrt_overlap,
                                    raw_op,
                                )
                            )
                        else:
                            op_objs.append(
                                _core.MetricInverseProduct(
                                    metric_policy_obj.inverse_overlap, raw_op
                                )
                            )
                    else:
                        op_objs.append(raw_op)

            plan = kubo_plan(
                num_moments,
                left_index=left_idx,
                right_index=right_idx,
                word=match.engine_word,
            )

    elif shape_id == "msd":
        vel_key = _resolve_op_key(match.operators[0])
        op_names = [vel_key]
        if vel_key in concrete_ops:
            op_objs = [concrete_ops[vel_key]]
        assert dt is not None
        coeffs, phase = time_evolution_series(
            window, dt, use_legacy_bessel=use_legacy_bessel
        )
        plan = msd_plan(
            num_moments,
            num_steps,
            num_states,
            coeffs,
            velocity_index=0,
            phase=phase,
            word=match.engine_word,
        )

    elif shape_id == "projected_evolution":
        obs_key = _resolve_op_key(match.operators[0])
        proj_key = _resolve_op_key(match.operators[1])
        op_names = [obs_key, proj_key]
        if obs_key in concrete_ops and proj_key in concrete_ops:
            op_objs = [concrete_ops[obs_key], concrete_ops[proj_key]]
        assert dt is not None
        coeffs, phase = time_evolution_series(
            window, dt, use_legacy_bessel=use_legacy_bessel
        )
        legacy_prefix = f"TimeProj{obs_key}-{proj_key}"
        match = CatalogueMatch(
            entry=match.entry,
            operators=match.operators,
            engine_word=f"{proj_key} U(t)^+ {obs_key} U(t) delta {proj_key}",
            legacy_moment_prefix=legacy_prefix,
            legacy_moment_ext=match.legacy_moment_ext,
            default_formulas=match.default_formulas,
        )
        plan = projected_evolution_plan(
            num_moments,
            num_steps,
            num_states,
            coeffs,
            operator_index=0,
            projector_index=1,
            phase=phase,
            word=match.engine_word,
        )

    elif shape_id == "time_correlation":
        key_a = _resolve_op_key(match.operators[0])
        key_b = _resolve_op_key(match.operators[1])
        if key_a == key_b:
            op_names = [key_a]
            left_idx = 0
            right_idx = 0
        else:
            op_names = [key_a, key_b]
            left_idx = 0
            right_idx = 1
        if all(k in concrete_ops for k in op_names):
            op_objs = [concrete_ops[k] for k in op_names]
        assert dt is not None
        coeffs, phase = time_evolution_series(
            window, dt, use_legacy_bessel=use_legacy_bessel
        )
        plan = time_correlation_plan(
            num_moments,
            num_steps,
            num_states,
            coeffs,
            left_index=left_idx,
            right_index=right_idx,
            phase=phase,
            word=match.engine_word,
        )
    else:
        raise ValueError(f"Unhandled catalogue shape '{shape_id}'")

    raw_formulas = recon_cfg.get("formulas")
    if raw_formulas is not None:
        formulas = tuple(str(f) for f in raw_formulas)
        allowed_f = SHAPE_ALLOWED_FORMULAS.get(shape_id, frozenset())
        for f in formulas:
            if f.lower() not in allowed_f:
                raise ValueError(
                    f"Reconstruction formula '{f}' does not match word "
                    f"'{parsed.canonical()}' (shape '{shape_id}'; allowed formulas: "
                    f"{sorted(allowed_f)})"
                )
    else:
        formulas = match.default_formulas

    kernel = str(_get_or(recon_cfg, "kernel", "jackson")).lower()
    broadening = (
        float(recon_cfg["broadening_meV"])
        if recon_cfg.get("broadening_meV") is not None
        else None
    )
    fermi_list = tuple(
        float(x) for x in (recon_cfg.get("fermi_energies") or ())
    )
    num_points = (
        int(recon_cfg["num_points"])
        if recon_cfg.get("num_points") is not None
        else None
    )
    use_fft_grid = bool(_get_or(recon_cfg, "fft_grid", False))
    legacy_rt = bool(_get_or(recon_cfg, "legacy_text_round_trip", False))

    recon = ReconstructionSpec(
        shape_id=shape_id,
        default_reconstruction=match.entry.default_reconstruction,
        formulas=formulas,
        kernel=kernel,
        broadening_meV=broadening,
        fermi_energies=fermi_list,
        num_points=num_points,
        use_fft_grid=use_fft_grid,
        legacy_moment_prefix=match.legacy_moment_prefix,
        legacy_moment_ext=match.legacy_moment_ext,
        legacy_text_round_trip=legacy_rt,
    )

    return CompiledWord(
        word=parsed,
        match=match,
        plan=plan,
        operator_names=tuple(op_names),
        operators=op_objs,
        reconstruction=recon,
        dt=dt,
    )


def describe(
    word: str | Word,
    basis: Any = None,
    params: Mapping[str, Any] | None = None,
) -> dict[str, object]:
    """Return a structured description of ``word`` (parse AST, plan instructions, operators, reconstruction)."""
    compiled = compile(
        word,
        basis=basis,
        params=params,
        _allow_default_time=(params is None),
    )
    return compiled.to_dict()


def list_catalogue() -> list[dict[str, object]]:
    """Return the list of all recognised catalogue shapes as dictionaries."""
    return [entry.to_dict() for entry in CATALOGUE]
