"""Catalogue of recognised observable word shapes, legacy names, and output conventions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from .grammar import (
    CommutatorSlot,
    OperatorSlot,
    SpectralSlot,
    TRANSFORM_PAIRS,
    Word,
)


@dataclass(frozen=True)
class CatalogueEntry:
    """Metadata for one recognised word shape in the LSQUANT word catalogue."""

    shape_id: str
    canonical_example: str
    pattern_description: str
    plan_builder: str
    default_reconstruction: str
    legacy_moment_prefix: str
    legacy_curve_prefix: str
    rank: int
    time_dependent: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "shape_id": self.shape_id,
            "canonical_example": self.canonical_example,
            "pattern_description": self.pattern_description,
            "plan_builder": self.plan_builder,
            "default_reconstruction": self.default_reconstruction,
            "legacy_moment_prefix": self.legacy_moment_prefix,
            "legacy_curve_prefix": self.legacy_curve_prefix,
            "rank": self.rank,
            "time_dependent": self.time_dependent,
        }


CATALOGUE: tuple[CatalogueEntry, ...] = (
    CatalogueEntry(
        shape_id="dos_identity",
        canonical_example="Tr[delta(E-H)]",
        pattern_description="Tr[delta(E-H)]",
        plan_builder="dos_plan",
        default_reconstruction="dos_curve",
        legacy_moment_prefix="SpectralOp1",
        legacy_curve_prefix="mean",
        rank=1,
        time_dependent=False,
    ),
    CatalogueEntry(
        shape_id="dos_operator",
        canonical_example="Tr[A delta(E-H)]",
        pattern_description="Tr[A delta(E-H)]",
        plan_builder="dos_plan",
        default_reconstruction="dos_curve",
        legacy_moment_prefix="SpectralOp<A>",
        legacy_curve_prefix="mean",
        rank=1,
        time_dependent=False,
    ),
    CatalogueEntry(
        shape_id="kubo_greenwood",
        canonical_example="Tr[VX delta(E-H) VX delta(E-H)]",
        pattern_description="Tr[A delta(E-H) B delta(E-H)]",
        plan_builder="kubo_plan",
        default_reconstruction="kubo_greenwood",
        legacy_moment_prefix="NonEqOp<A>-<B>",
        legacy_curve_prefix="KuboGreenWood_",
        rank=2,
        time_dependent=False,
    ),
    CatalogueEntry(
        shape_id="kubo_bastin",
        canonical_example="Tr[VX delta(E-H) VY dG-(E)/dE]",
        pattern_description="Tr[A delta(E-H) B dG-(E)/dE] (or sea[...])",
        plan_builder="kubo_plan",
        default_reconstruction="kubo_bastin",
        legacy_moment_prefix="NonEqOp<A>-<B>",
        legacy_curve_prefix="KuboBastin_",
        rank=2,
        time_dependent=False,
    ),
    CatalogueEntry(
        shape_id="msd",
        canonical_example="Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]",
        pattern_description="Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]",
        plan_builder="msd_plan",
        default_reconstruction="time_correlation_curve",
        legacy_moment_prefix="Correlation<V>",
        legacy_curve_prefix="mean",
        rank=1,
        time_dependent=True,
    ),
    CatalogueEntry(
        shape_id="projected_evolution",
        canonical_example="Tr[P U(t)^+ SX U(t) delta(E-H) P]",
        pattern_description="Tr[P U(t)^+ A U(t) delta(E-H) P]",
        plan_builder="projected_evolution_plan",
        default_reconstruction="row_m0",
        legacy_moment_prefix="TimeProj<A>-<P>",
        legacy_curve_prefix="TimeProj",
        rank=1,
        time_dependent=True,
    ),
    CatalogueEntry(
        shape_id="time_correlation",
        canonical_example="Tr[A delta(E-H) U(t)^+ B U(t)]",
        pattern_description="Tr[A delta(E-H) U(t)^+ B U(t)]",
        plan_builder="time_correlation_plan",
        default_reconstruction="time_correlation_curve",
        legacy_moment_prefix="TimeCorr<A>-<B>",
        legacy_curve_prefix="mean",
        rank=1,
        time_dependent=True,
    ),
)

CATALOGUE_BY_ID: dict[str, CatalogueEntry] = {
    entry.shape_id: entry for entry in CATALOGUE
}

SHAPE_ALLOWED_FORMULAS: dict[str, frozenset[str]] = {
    "dos_identity": frozenset({"dos", "dos_curve"}),
    "dos_operator": frozenset({"dos", "dos_curve"}),
    "kubo_greenwood": frozenset(
        {"greenwood", "kubo_greenwood", "greenwood_fft", "kubo_greenwood_fft"}
    ),
    "kubo_bastin": frozenset(
        {
            "shc",
            "bastin",
            "kubo_bastin",
            "sea",
            "kubo_bastin_sea",
            "surf",
            "kubo_bastin_surf",
            "bastin_fft",
            "kubo_bastin_fft",
        }
    ),
    "msd": frozenset({"time_correlation", "time_correlation_curve"}),
    "projected_evolution": frozenset(
        {"row_m0", "row0", "time_correlation", "time_correlation_curve"}
    ),
    "time_correlation": frozenset(
        {"time_correlation", "time_correlation_curve"}
    ),
}


@dataclass(frozen=True)
class CatalogueMatch:
    """Resolved match between a parsed :class:`Word` and a :class:`CatalogueEntry`."""

    entry: CatalogueEntry
    operators: tuple[OperatorSlot, ...]
    engine_word: str
    legacy_moment_prefix: str
    legacy_moment_ext: str
    default_formulas: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "shape_id": self.entry.shape_id,
            "canonical_example": self.entry.canonical_example,
            "plan_builder": self.entry.plan_builder,
            "default_reconstruction": self.entry.default_reconstruction,
            "default_formulas": list(self.default_formulas),
            "engine_word": self.engine_word,
            "legacy_moment_prefix": self.legacy_moment_prefix,
            "legacy_moment_ext": self.legacy_moment_ext,
            "operators": [op.to_dict() for op in self.operators],
        }


def validate_word_bases(
    word: Word,
    calculation_basis: str | None = None,
    operator_bases: Mapping[str, str] | None = None,
) -> None:
    """Check that all operator slots in ``word`` are in a consistent basis or explicitly transformed.

    Raises ``ValueError`` with ``'mixed basis without explicit transform'`` if:
    - an unmatched transform operator (``B^+``, ``U^+``) appears in the word;
    - an operator has an explicit ``@<basis>`` tag or registered basis that
      disagrees with the active basis without an enclosing ``B^+ ... B`` or
      ``U^+ ... U`` transform sandwich;
    - two operators in the word carry conflicting ``@<basis>`` tags without
      transforms bringing them into a common basis.
    """
    inferred_basis = calculation_basis
    for slot in word.slots:
        if not isinstance(slot, OperatorSlot):
            continue
        if slot.dagger and slot.name in ("B", "U") and slot.transform is None:
            raise ValueError(
                f"Mixed basis without explicit transform: unmatched transform operator "
                f"'{slot.canonical()}' in '{word.canonical()}'"
            )

        op_basis = slot.basis
        if op_basis is None and operator_bases is not None:
            op_basis = operator_bases.get(slot.operator_key) or operator_bases.get(
                slot.name
            )

        if slot.transform is not None:
            inner_basis, outer_basis = TRANSFORM_PAIRS[slot.transform]
            if op_basis is not None and op_basis != inner_basis:
                raise ValueError(
                    f"Mixed basis without explicit transform: operator '{slot.name}' "
                    f"is in '{op_basis}' basis, but transform {slot.transform[0]} ... "
                    f"{slot.transform[1]} expects inner basis '{inner_basis}'"
                )
            effective_basis = outer_basis
        else:
            effective_basis = op_basis

        if effective_basis is not None:
            if inferred_basis is None:
                inferred_basis = effective_basis
            elif effective_basis != inferred_basis:
                raise ValueError(
                    f"Mixed basis without explicit transform: operator '{slot.canonical()}' "
                    f"is in '{effective_basis}' basis while calculation/active basis is "
                    f"'{inferred_basis}' (insert explicit B^+ ... B or U^+ ... U transform)"
                )


def nearest_catalogue_entry(word: Word) -> CatalogueEntry:
    """Return the nearest :class:`CatalogueEntry` to ``word`` for diagnostic messages."""
    has_comm = any(isinstance(s, CommutatorSlot) for s in word.slots)
    if has_comm:
        return CATALOGUE_BY_ID["msd"]

    has_u = any(
        isinstance(s, SpectralSlot) and s.kind in ("U", "U^+") for s in word.slots
    )
    if has_u:
        ops = [s for s in word.slots if isinstance(s, OperatorSlot)]
        if (
            len(ops) >= 2
            and isinstance(word.slots[0], OperatorSlot)
            and isinstance(word.slots[-1], OperatorSlot)
            and word.slots[0].operator_key == word.slots[-1].operator_key
        ):
            return CATALOGUE_BY_ID["projected_evolution"]
        return CATALOGUE_BY_ID["time_correlation"]

    has_green = any(
        isinstance(s, SpectralSlot) and s.kind in ("G+", "G-", "dG-/dE", "dG+/dE")
        for s in word.slots
    )
    if has_green or word.sea:
        return CATALOGUE_BY_ID["kubo_bastin"]

    if word.rank >= 2:
        return CATALOGUE_BY_ID["kubo_greenwood"]

    ops = [s for s in word.slots if isinstance(s, OperatorSlot)]
    if ops:
        return CATALOGUE_BY_ID["dos_operator"]
    return CATALOGUE_BY_ID["dos_identity"]


def match_catalogue(word: Word) -> CatalogueMatch:
    """Match ``word`` against the word catalogue or raise a descriptive ``ValueError``."""
    validate_word_bases(word)

    if word.rank >= 3:
        nearest = nearest_catalogue_entry(word)
        raise ValueError(
            f"Unsupported word '{word.canonical()}' of rank {word.rank} "
            f"({word.rank} expanded spectral slots): rank >= 3 moment tensors are out of scope; "
            f"nearest recognised shape is '{nearest.canonical_example}' ({nearest.shape_id})"
        )

    slots = word.slots

    # 1. Tr[delta(E-H)]
    if (
        not word.sea
        and len(slots) == 1
        and isinstance(slots[0], SpectralSlot)
        and slots[0].kind == "delta"
    ):
        entry = CATALOGUE_BY_ID["dos_identity"]
        return CatalogueMatch(
            entry=entry,
            operators=(),
            engine_word="1 delta",
            legacy_moment_prefix="SpectralOp1",
            legacy_moment_ext=".chebmom1D",
            default_formulas=("dos",),
        )

    # 2. Tr[A delta(E-H)]
    if (
        not word.sea
        and len(slots) == 2
        and isinstance(slots[0], OperatorSlot)
        and isinstance(slots[1], SpectralSlot)
        and slots[1].kind == "delta"
    ):
        op_a = slots[0]
        if op_a.dagger:
            raise ValueError(
                f"Operator '{op_a.canonical()}' in '{word.canonical()}' must be Hermitian (no '^+')"
            )
        entry = CATALOGUE_BY_ID["dos_operator"]
        key_a = op_a.operator_key
        ops_tuple = () if op_a.is_identity else (op_a,)
        return CatalogueMatch(
            entry=entry,
            operators=ops_tuple,
            engine_word=f"{key_a} delta",
            legacy_moment_prefix=f"SpectralOp{key_a}",
            legacy_moment_ext=".chebmom1D",
            default_formulas=("dos",),
        )

    # 3 & 4. Rank-2 Kubo words: Tr[A delta(E-H) B delta(E-H)] or Tr[A delta(E-H) B dG-(E)/dE]
    if (
        len(slots) == 4
        and isinstance(slots[0], OperatorSlot)
        and isinstance(slots[1], SpectralSlot)
        and slots[1].kind == "delta"
        and isinstance(slots[2], OperatorSlot)
        and isinstance(slots[3], SpectralSlot)
    ):
        op_a = slots[0]
        op_b = slots[2]
        spec_2 = slots[3]
        if op_a.dagger or op_b.dagger:
            raise ValueError(
                f"Kubo operators '{op_a.canonical()}', '{op_b.canonical()}' must be Hermitian"
            )
        key_a = op_a.operator_key
        key_b = op_b.operator_key
        if spec_2.kind == "delta" and not word.sea:
            entry = CATALOGUE_BY_ID["kubo_greenwood"]
            return CatalogueMatch(
                entry=entry,
                operators=(op_a, op_b),
                engine_word=f"{key_a} delta {key_b} delta",
                legacy_moment_prefix=f"NonEqOp{key_a}-{key_b}",
                legacy_moment_ext=".chebmom2D",
                default_formulas=("greenwood",),
            )
        if spec_2.kind == "dG-/dE":
            entry = CATALOGUE_BY_ID["kubo_bastin"]
            default_f = ("sea",) if word.sea else ("bastin",)
            return CatalogueMatch(
                entry=entry,
                operators=(op_a, op_b),
                engine_word=f"{key_a} delta {key_b} delta",
                legacy_moment_prefix=f"NonEqOp{key_a}-{key_b}",
                legacy_moment_ext=".chebmom2D",
                default_formulas=default_f,
            )

    # 5. Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]
    if (
        not word.sea
        and len(slots) == 3
        and isinstance(slots[0], CommutatorSlot)
        and slots[0].dagger
        and isinstance(slots[1], SpectralSlot)
        and slots[1].kind == "delta"
        and isinstance(slots[2], CommutatorSlot)
        and not slots[2].dagger
        and slots[0].axis == slots[2].axis
    ):
        comm = slots[0]
        vel_name = comm.velocity_operator
        vel_op = OperatorSlot(name=vel_name)
        entry = CATALOGUE_BY_ID["msd"]
        return CatalogueMatch(
            entry=entry,
            operators=(vel_op,),
            engine_word=f"[{comm.axis}, U(t)]^+ delta [{comm.axis}, U(t)]",
            legacy_moment_prefix=f"Correlation{vel_name}",
            legacy_moment_ext=".chebmomTD",
            default_formulas=("time_correlation",),
        )

    # 6. Tr[P U(t)^+ SX U(t) delta(E-H) P]
    if (
        not word.sea
        and len(slots) == 6
        and isinstance(slots[0], OperatorSlot)
        and isinstance(slots[1], SpectralSlot)
        and slots[1].kind == "U^+"
        and isinstance(slots[2], OperatorSlot)
        and isinstance(slots[3], SpectralSlot)
        and slots[3].kind == "U"
        and isinstance(slots[4], SpectralSlot)
        and slots[4].kind == "delta"
        and isinstance(slots[5], OperatorSlot)
        and slots[0].canonical() == slots[5].canonical()
    ):
        proj = slots[0]
        obs = slots[2]
        if proj.dagger or obs.dagger:
            raise ValueError(
                f"Projected evolution operators '{proj.canonical()}', '{obs.canonical()}' must be Hermitian"
            )
        proj_key = proj.operator_key
        obs_key = obs.operator_key
        entry = CATALOGUE_BY_ID["projected_evolution"]
        return CatalogueMatch(
            entry=entry,
            operators=(obs, proj),
            engine_word=f"{proj_key} U(t)^+ {obs_key} U(t) delta {proj_key}",
            legacy_moment_prefix=f"TimeProj{obs_key}-{proj_key}",
            legacy_moment_ext=".chebmomTD",
            default_formulas=("row_m0",),
        )

    # 7. Tr[A delta(E-H) U(t)^+ B U(t)]
    if (
        not word.sea
        and len(slots) == 5
        and isinstance(slots[0], OperatorSlot)
        and isinstance(slots[1], SpectralSlot)
        and slots[1].kind == "delta"
        and isinstance(slots[2], SpectralSlot)
        and slots[2].kind == "U^+"
        and isinstance(slots[3], OperatorSlot)
        and isinstance(slots[4], SpectralSlot)
        and slots[4].kind == "U"
    ):
        op_a = slots[0]
        op_b = slots[3]
        if op_a.dagger or op_b.dagger:
            raise ValueError(
                f"Time-correlation operators '{op_a.canonical()}', '{op_b.canonical()}' must be Hermitian"
            )
        key_a = op_a.operator_key
        key_b = op_b.operator_key
        entry = CATALOGUE_BY_ID["time_correlation"]
        return CatalogueMatch(
            entry=entry,
            operators=(op_a, op_b),
            engine_word=f"{key_a} delta U(t)^+ {key_b} U(t)",
            legacy_moment_prefix=f"TimeCorr{key_a}-{key_b}",
            legacy_moment_ext=".chebmomTD",
            default_formulas=("time_correlation",),
        )

    nearest = nearest_catalogue_entry(word)
    raise ValueError(
        f"Unrecognised word shape '{word.canonical()}'; "
        f"nearest recognised shape is '{nearest.canonical_example}' ({nearest.shape_id})"
    )
