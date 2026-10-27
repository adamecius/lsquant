"""Observable word language: grammar, catalogue, and compiler to bytecode plans."""

from __future__ import annotations

from .catalogue import (
    CATALOGUE,
    CATALOGUE_BY_ID,
    CatalogueEntry,
    CatalogueMatch,
    match_catalogue,
    nearest_catalogue_entry,
    validate_word_bases,
)
from .compiler import (
    CompiledWord,
    ReconstructionSpec,
    compile,
    describe,
    format_plan_instructions,
    list_catalogue,
)
from .grammar import (
    CommutatorSlot,
    OperatorSlot,
    SpectralSlot,
    Token,
    Word,
    parse,
    tokenize,
)

__all__ = [
    "CATALOGUE",
    "CATALOGUE_BY_ID",
    "CatalogueEntry",
    "CatalogueMatch",
    "CommutatorSlot",
    "CompiledWord",
    "OperatorSlot",
    "ReconstructionSpec",
    "SpectralSlot",
    "Token",
    "Word",
    "compile",
    "describe",
    "format_plan_instructions",
    "list_catalogue",
    "match_catalogue",
    "nearest_catalogue_entry",
    "parse",
    "tokenize",
    "validate_word_bases",
]
