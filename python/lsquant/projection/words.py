"""Projection-aware preparation on top of the unchanged stage-6b compiler."""
from dataclasses import replace

from .. import _core
from ..apps.stage6b_siesta import prepare_moments, compute_moments
from ..words import parse as parse_word
from ..words.grammar import OperatorSlot


def prepare_projection_moments(basis, word, **options):
    """Preserve an EXPLICIT projected identity slot through compilation.

    The general word catalogue removes literal 1 before compiling. Once
    selected it represents P, not the full identity, so give it the existing
    named-projector spelling P:WINDOW. A DOS word without any static slot
    still counts every zero/edge complement state, as the brief requires.
    """
    parsed = parse_word(word)
    selected = replace(parsed, slots=tuple(
        replace(slot, name="P:WINDOW") if isinstance(slot, OperatorSlot) and slot.is_identity else slot
        for slot in parsed.slots))
    prepared = prepare_moments(basis, selected.canonical(), **options)
    engine, plan, compiled, generator, window, operators, h = prepared
    plan = _core.Plan(plan.code, parsed.canonical())
    return engine, plan, compiled, generator, window, operators, h


def compute_projection_moments(basis, word, *, _prepared=None, **options):
    """Execute a physical projected word, retaining its original receipt text."""
    prepared = prepare_projection_moments(basis, word, **options) if _prepared is None else _prepared
    mu, compiled = compute_moments(basis, word, _prepared=prepared)
    mu.word = parse_word(word).canonical()
    return mu, compiled
