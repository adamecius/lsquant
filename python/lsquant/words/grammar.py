"""Tokenizer, AST, and hand-written recursive-descent parser for observable words.

Grammar (EBNF)::

    word          ::= "Tr[" body "]" | "sea[" "Tr[" product "]" "]" ;
    body          ::= "sea[" product "]" | product ;
    product       ::= factor ( ["*"] factor )* ;
    factor        ::= spectral | commutator | operator ;
    spectral      ::= "delta(E-H)"
                    | "G+(E)" | "G-(E)"
                    | "dG-(E)/dE" | "dG+(E)/dE"
                    | "U(t)^+" | "U(t)" ;
    commutator    ::= "[" axis "," "U(t)" "]^+"
                    | "[" axis "," "U(t)" "]" ;
    axis          ::= "X" | "Y" | "Z" | "VX" | "VY" | "VZ" ;
    operator      ::= ("1" | "P:" ident | ident) ["^+"] ["@" basis_tag] ;
    ident         ::= LETTER (LETTER | DIGIT | "_")* ;
    basis_tag     ::= "real" | "bloch" | "eigen" ;
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union


VALID_BASES: frozenset[str] = frozenset({"real", "bloch", "eigen"})
VALID_AXES: frozenset[str] = frozenset({"X", "Y", "Z", "VX", "VY", "VZ"})
TRANSFORM_PAIRS: dict[tuple[str, str], tuple[str, str]] = {
    ("B^+", "B"): ("bloch", "real"),
    ("B", "B^+"): ("real", "bloch"),
    ("U^+", "U"): ("eigen", "bloch"),
    ("U", "U^+"): ("bloch", "eigen"),
}


@dataclass(frozen=True)
class Token:
    """Lexical token emitted by :func:`tokenize`."""

    kind: str
    value: str
    pos: int


@dataclass(frozen=True)
class OperatorSlot:
    """Static operator slot in a parsed :class:`Word`."""

    name: str
    dagger: bool = False
    basis: str | None = None
    transform: tuple[str, str] | None = None

    @property
    def is_identity(self) -> bool:
        return self.name == "1" and not self.dagger and self.transform is None

    @property
    def is_projector(self) -> bool:
        return (
            self.name == "P"
            or self.name.startswith("P:")
            or self.name in ("PX", "PY", "PZ")
        )

    @property
    def operator_key(self) -> str:
        """Name used to look up the operator in a model or basis."""
        if self.name.startswith("P:"):
            return self.name[2:]
        return self.name

    def canonical(self) -> str:
        core = self.name
        if self.dagger:
            core += "^+"
        if self.basis is not None:
            core += f"@{self.basis}"
        if self.transform is not None:
            return f"{self.transform[0]} {core} {self.transform[1]}"
        return core

    def to_dict(self) -> dict[str, object]:
        out: dict[str, object] = {
            "type": "operator",
            "name": self.name,
            "operator_key": self.operator_key,
        }
        if self.dagger:
            out["dagger"] = True
        if self.basis is not None:
            out["basis"] = self.basis
        if self.transform is not None:
            out["transform"] = list(self.transform)
        return out


@dataclass(frozen=True)
class SpectralSlot:
    """Spectral function slot ``f(H)`` in a parsed :class:`Word`."""

    kind: str  # "delta", "G+", "G-", "dG-/dE", "dG+/dE", "U", "U^+"

    @property
    def mode(self) -> str:
        """Evaluation mode: ``'expand'`` for energy slots, ``'apply'`` for ``U(t)``."""
        if self.kind in ("U", "U^+"):
            return "apply"
        return "expand"

    def canonical(self) -> str:
        mapping = {
            "delta": "delta(E-H)",
            "G+": "G+(E)",
            "G-": "G-(E)",
            "dG-/dE": "dG-(E)/dE",
            "dG+/dE": "dG+(E)/dE",
            "U": "U(t)",
            "U^+": "U(t)^+",
        }
        return mapping[self.kind]

    def to_dict(self) -> dict[str, object]:
        return {
            "type": "spectral",
            "kind": self.kind,
            "token": self.canonical(),
            "mode": self.mode,
        }


@dataclass(frozen=True)
class CommutatorSlot:
    """Position-commutator time-evolution slot ``[X,U(t)]`` or ``[X,U(t)]^+``."""

    axis: str
    dagger: bool = False

    @property
    def mode(self) -> str:
        return "apply"

    @property
    def velocity_operator(self) -> str:
        if self.axis in ("X", "VX"):
            return "VX"
        if self.axis in ("Y", "VY"):
            return "VY"
        if self.axis in ("Z", "VZ"):
            return "VZ"
        raise ValueError(f"Unsupported commutator axis '{self.axis}'")

    def canonical(self) -> str:
        suffix = "^+" if self.dagger else ""
        return f"[{self.axis},U(t)]{suffix}"

    def to_dict(self) -> dict[str, object]:
        return {
            "type": "commutator",
            "axis": self.axis,
            "dagger": self.dagger,
            "velocity_operator": self.velocity_operator,
            "token": self.canonical(),
            "mode": self.mode,
        }


Slot = Union[OperatorSlot, SpectralSlot, CommutatorSlot]


@dataclass(frozen=True)
class Word:
    """Parsed observable word ``Tr[A_0 f_1(H) ... f_p(H) A_p]``."""

    slots: tuple[Slot, ...]
    sea: bool = False

    @property
    def rank(self) -> int:
        """Number of expanded spectral slots in the word."""
        return sum(
            1
            for s in self.slots
            if isinstance(s, SpectralSlot) and s.mode == "expand"
        )

    @property
    def is_time_dependent(self) -> bool:
        """Whether the word contains ``U(t)``, ``U(t)^+``, or ``[X,U(t)]``."""
        for s in self.slots:
            if isinstance(s, CommutatorSlot):
                return True
            if isinstance(s, SpectralSlot) and s.kind in ("U", "U^+"):
                return True
        return False

    def canonical(self) -> str:
        """Return the canonical string representation of the word."""
        inner = " ".join(s.canonical() for s in self.slots)
        if self.sea:
            return f"Tr[sea[{inner}]]"
        return f"Tr[{inner}]"

    def __str__(self) -> str:
        return self.canonical()

    def to_dict(self) -> dict[str, object]:
        return {
            "canonical": self.canonical(),
            "sea": self.sea,
            "rank": self.rank,
            "time_dependent": self.is_time_dependent,
            "slots": [s.to_dict() for s in self.slots],
        }


def _match_spectral_at(text: str, pos: int) -> tuple[str, int] | None:
    """Match a spectral function token starting at ``pos``."""
    rest = text[pos:]
    # Check delta(E-H) with optional spaces around '-'
    if rest.startswith("delta("):
        close = rest.find(")")
        if close != -1:
            inside = "".join(rest[6:close].split())
            if inside == "E-H":
                return ("delta", pos + close + 1)
            raise ValueError(
                f"Malformed delta spectral token '{rest[:close + 1]}' at position {pos}; "
                "expected 'delta(E-H)'"
            )
    for literal, kind in (
        ("dG-(E)/dE", "dG-/dE"),
        ("dG+(E)/dE", "dG+/dE"),
        ("dG-(E-H)/dE", "dG-/dE"),
        ("dG+(E-H)/dE", "dG+/dE"),
        ("G+(E-H)", "G+"),
        ("G-(E-H)", "G-"),
        ("G+(E)", "G+"),
        ("G-(E)", "G-"),
        ("U(t)^+", "U^+"),
        ("U(t)", "U"),
    ):
        if rest.startswith(literal):
            return (kind, pos + len(literal))
    return None


def _match_commutator_at(text: str, pos: int) -> tuple[tuple[str, bool], int] | None:
    """Match ``[<axis>,U(t)]`` or ``[<axis>,U(t)]^+`` starting at ``pos``."""
    if text[pos] != "[":
        return None
    close = text.find("]", pos)
    if close == -1:
        raise ValueError(f"Unclosed '[' in commutator at position {pos}")
    inside = "".join(text[pos + 1 : close].split())
    if "," not in inside:
        raise ValueError(
            f"Malformed commutator '[{text[pos + 1 : close]}]' at position {pos}; "
            "expected '[X,U(t)]' or '[X,U(t)]^+'"
        )
    axis, rhs = inside.split(",", 1)
    if axis not in VALID_AXES or rhs != "U(t)":
        raise ValueError(
            f"Malformed commutator '[{text[pos + 1 : close]}]' at position {pos}; "
            "expected '[X,U(t)]' with axis in {X, Y, Z, VX, VY, VZ}"
        )
    end = close + 1
    dagger = False
    if text[end:].startswith("^+"):
        dagger = True
        end += 2
    return ((axis, dagger), end)


def tokenize(text: str) -> list[Token]:
    """Tokenize an observable word string into a sequence of :class:`Token` objects."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Word text must be a non-empty string")

    tokens: list[Token] = []
    n = len(text)
    pos = 0

    while pos < n:
        ch = text[pos]
        if ch.isspace():
            pos += 1
            continue
        if ch == "*":
            tokens.append(Token("STAR", "*", pos))
            pos += 1
            continue
        if text[pos:].startswith("Tr["):
            tokens.append(Token("TR_OPEN", "Tr[", pos))
            pos += 3
            continue
        if text[pos:].startswith("sea["):
            tokens.append(Token("SEA_OPEN", "sea[", pos))
            pos += 4
            continue
        if ch == "]":
            tokens.append(Token("RBRACKET", "]", pos))
            pos += 1
            continue

        spec = _match_spectral_at(text, pos)
        if spec is not None:
            kind, new_pos = spec
            tokens.append(Token("SPECTRAL", kind, pos))
            pos = new_pos
            continue

        comm = _match_commutator_at(text, pos)
        if comm is not None:
            (axis, dagger), new_pos = comm
            val = f"{axis}:{'+' if dagger else '-'}"
            tokens.append(Token("COMMUTATOR", val, pos))
            pos = new_pos
            continue

        # Operator token: "1" or identifier, optional ":ident" for P:<name>,
        # optional "^+", optional "@<basis>".
        if ch == "1" and (
            pos + 1 == n or text[pos + 1].isspace() or text[pos + 1] in "*]@^"
        ):
            op_name = "1"
            end = pos + 1
        elif ch.isalpha() or ch == "_":
            end = pos + 1
            while end < n and (text[end].isalnum() or text[end] == "_"):
                end += 1
            op_name = text[pos:end]
            if op_name == "P" and end < n and text[end] == ":":
                end += 1
                sub_start = end
                while end < n and (text[end].isalnum() or text[end] == "_"):
                    end += 1
                if end == sub_start:
                    raise ValueError(
                        f"Expected projector name after 'P:' at position {pos}"
                    )
                op_name = f"P:{text[sub_start:end]}"
        else:
            raise ValueError(
                f"Unexpected character '{ch}' at position {pos} in word '{text}'"
            )

        if end < n and text[end] == "(":
            close = text.find(")", end)
            snippet = text[pos : (close + 1 if close != -1 else n)]
            raise ValueError(
                f"Unknown spectral function or malformed token '{snippet}' at position {pos}"
            )

        dagger = False
        if text[end:].startswith("^+"):
            dagger = True
            end += 2

        basis_tag: str | None = None
        if end < n and text[end] == "@":
            end += 1
            b_start = end
            while end < n and (text[end].isalnum() or text[end] == "_"):
                end += 1
            basis_tag = text[b_start:end]
            if not basis_tag or basis_tag not in VALID_BASES:
                raise ValueError(
                    f"Unknown basis tag '@{basis_tag}' at position {pos}; "
                    f"expected one of {sorted(VALID_BASES)}"
                )

        encoded = f"{op_name}|{'1' if dagger else '0'}|{basis_tag or ''}"
        tokens.append(Token("OPERATOR", encoded, pos))
        pos = end

    return tokens


class _Parser:
    """Hand-written recursive-descent parser over :class:`Token` stream."""

    def __init__(self, tokens: list[Token], source: str) -> None:
        self._tokens = tokens
        self._source = source
        self._idx = 0

    def _peek(self) -> Token | None:
        if self._idx < len(self._tokens):
            return self._tokens[self._idx]
        return None

    def _consume(self, expected_kind: str | None = None) -> Token:
        tok = self._peek()
        if tok is None:
            raise ValueError(
                f"Unexpected end of word in '{self._source}'"
                + (f" (expected {expected_kind})" if expected_kind else "")
            )
        if expected_kind is not None and tok.kind != expected_kind:
            raise ValueError(
                f"Expected token {expected_kind} but found '{tok.value}' ({tok.kind}) "
                f"at position {tok.pos} in '{self._source}'"
            )
        self._idx += 1
        return tok

    def parse_word(self) -> Word:
        tok = self._peek()
        if tok is None:
            raise ValueError("Empty word")

        sea = False
        outer_sea = False
        if tok.kind == "SEA_OPEN":
            self._consume("SEA_OPEN")
            sea = True
            outer_sea = True
            self._consume("TR_OPEN")
        elif tok.kind == "TR_OPEN":
            self._consume("TR_OPEN")
            if self._peek() is not None and self._peek().kind == "SEA_OPEN":
                self._consume("SEA_OPEN")
                sea = True
        else:
            raise ValueError(
                f"Word must start with 'Tr[' (or 'sea[Tr[...]]'), found '{tok.value}' "
                f"at position {tok.pos} in '{self._source}'"
            )

        raw_slots = self._parse_product()
        if not raw_slots:
            raise ValueError(f"Empty trace body in '{self._source}'")

        if sea and not outer_sea:
            self._consume("RBRACKET")
        self._consume("RBRACKET")
        if outer_sea:
            self._consume("RBRACKET")

        if self._peek() is not None:
            extra = self._peek()
            raise ValueError(
                f"Unexpected trailing token '{extra.value}' at position {extra.pos} "
                f"in '{self._source}'"
            )

        slots = _collapse_transform_sandwiches(raw_slots)
        return Word(slots=tuple(slots), sea=sea)

    def _parse_product(self) -> list[Slot]:
        slots: list[Slot] = []
        while True:
            tok = self._peek()
            if tok is None or tok.kind == "RBRACKET":
                break
            if tok.kind == "STAR":
                if not slots:
                    raise ValueError(
                        f"Leading '*' at position {tok.pos} in '{self._source}'"
                    )
                self._consume("STAR")
                nxt = self._peek()
                if nxt is None or nxt.kind in ("RBRACKET", "STAR"):
                    raise ValueError(
                        f"Trailing or consecutive '*' at position {tok.pos} in '{self._source}'"
                    )
            slots.append(self._parse_factor())
        return slots

    def _parse_factor(self) -> Slot:
        tok = self._peek()
        if tok is None:
            raise ValueError(f"Unexpected end of word in '{self._source}'")
        if tok.kind == "SPECTRAL":
            self._consume("SPECTRAL")
            return SpectralSlot(kind=tok.value)
        if tok.kind == "COMMUTATOR":
            self._consume("COMMUTATOR")
            axis, dag_flag = tok.value.split(":", 1)
            return CommutatorSlot(axis=axis, dagger=(dag_flag == "+"))
        if tok.kind == "OPERATOR":
            self._consume("OPERATOR")
            op_name, dag_str, basis_str = tok.value.split("|", 2)
            return OperatorSlot(
                name=op_name,
                dagger=(dag_str == "1"),
                basis=(basis_str if basis_str else None),
            )
        raise ValueError(
            f"Unexpected token '{tok.value}' ({tok.kind}) at position {tok.pos} "
            f"in '{self._source}'"
        )


def _collapse_transform_sandwiches(slots: list[Slot]) -> list[Slot]:
    """Collapse explicit basis transform sandwiches ``B^+ O B`` / ``U^+ O U``."""
    out: list[Slot] = []
    i = 0
    n = len(slots)
    while i < n:
        if i + 2 < n:
            left, mid, right = slots[i], slots[i + 1], slots[i + 2]
            if (
                isinstance(left, OperatorSlot)
                and isinstance(mid, OperatorSlot)
                and isinstance(right, OperatorSlot)
                and left.basis is None
                and right.basis is None
                and left.transform is None
                and mid.transform is None
                and right.transform is None
            ):
                left_tok = left.canonical()
                right_tok = right.canonical()
                if (left_tok, right_tok) in TRANSFORM_PAIRS:
                    out.append(
                        OperatorSlot(
                            name=mid.name,
                            dagger=mid.dagger,
                            basis=mid.basis,
                            transform=(left_tok, right_tok),
                        )
                    )
                    i += 3
                    continue
        out.append(slots[i])
        i += 1
    return out


def parse(text: str) -> Word:
    """Parse an observable word string into a :class:`Word` AST."""
    tokens = tokenize(text)
    parser = _Parser(tokens, text)
    return parser.parse_word()
