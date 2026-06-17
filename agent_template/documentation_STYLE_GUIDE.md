# How to document LSQUANT code

This guide explains how to turn `DOCSTRING_TEMPLATE.md` into finished in-code
documentation, and it is opinionated about voice on purpose, in the same spirit
as the tutorial guide. "Documentation" here means the comments that live inside
the source: docstrings in Python, Doxygen blocks in C and C++. It is not the
rendered website. Read a well documented module in the codebase, then read this.

## The voice

Write the documentation as if you were explaining to a competent colleague who
will call your code tomorrow and has never seen it, in the spirit of Sagan,
Feynman, or Tyson. That is not decoration. It is a discipline with consequences
for every block you write.

Lead with the problem, not the mechanism. The first paragraph after the summary
says what computational or physical problem this unit solves and why it exists,
before a single parameter is named. "Given a sparse Hamiltonian, estimate the
spectral density by a stochastic trace over Chebyshev moments" earns attention.
"This function takes `H`, `M`, and `n_states` and returns an array" does not.

Carry one responsibility. A docstring documents one unit doing one thing. If you
cannot state that thing in a single summary line, the unit is doing too much,
and the fix is in the code, not in longer prose.

Earn the block; never restate the signature. The reader can already see
`def dos(H, M)`. The documentation earns its place by saying what `H` and `M`
mean physically, in what units, with what shape and sign convention, and why the
function exists. Types and names that introspection already exposes are not
content.

Respect the reader's expertise, and only theirs. Assume a programmer fluent in
the language. Do not explain a list comprehension, a pointer, or what a trace
is. Do explain anything that steps outside that vocabulary: the physics, the
provenance of the algorithm, the units, the broadening, the storage order.
Define a domain term the first time it earns its place.

Trust the code to be legible. Compact documentation reads as confidence. Do not
narrate the implementation line by line; the source is right there. Document the
contract and the traps, not the control flow.

Name the trap. The most useful sentence is often the one that stops a competent
caller from the obvious misuse: the array that must be C-contiguous, the energy
that is in Rydberg and not eV, the routine that mutates its input in place, the
estimator whose variance grows when the moment count is too low. Find that
sentence for your unit and give it a warning block.

## What the generator does, and when

There are two layers, and they are not the same task.

Layer one is the in-code documentation: the docstrings and Doxygen blocks this
guide is about. Writing them is part of writing the code. They are always
present.

Layer two is the rendered output: the HTML or PDF that Sphinx (with numpydoc)
or Doxygen builds from layer one. This is a separate, explicit action.

The rule: always write layer one. Never produce layer two on your own
initiative. Do not create `conf.py`, a `Doxyfile`, a `docs/` build tree, or run
`sphinx-build` or `doxygen` unless a task asks for the rendered documentation in
words. When the task does not say otherwise, the deliverable is source whose
in-code documentation is complete, not a generated site.

Default depth, when nothing else is specified: every public unit gets an
extensive, problem-oriented description. The extended summary is where that
description lives. Terse one-liners are for genuinely obvious units only, and
only when the task asks for brevity.

## Hard rules

These are not stylistic preferences. Apply them without exception.

- Python follows PEP 8 for layout and PEP 257 for docstring conventions. The
  docstring body uses the numpydoc section format (Sphinx reads it through the
  numpydoc or napoleon extension). C and C++ use Doxygen blocks in Javadoc
  style: `/** ... */` with `@command`. Pick one comment convention per language
  and never mix it within a project.
- The summary line is one line, in the imperative mood, ends with a period, fits
  on the line, and does not restate the signature. "Estimate the density of
  states.", not "Estimates the DOS by taking H, M and returning an array."
- All mathematics is LaTeX. In Python docstrings use `` :math:`...` `` inline and
  a `.. math::` directive for displayed equations. In Doxygen use `\f$ ... \f$`
  inline and `\f[ ... \f]` displayed. Never write equations as ASCII in prose.
- No em dashes anywhere. Use a comma, a colon, parentheses, or a full stop.
- Document the contract, not the implementation. Parameters, returns, raised
  errors, units, array shapes, ownership, and side effects belong in the
  documentation. Implementation detail the caller does not need goes in `Notes`
  or `@details`, never in the summary.
- Every public symbol is documented: the module or file, every public function,
  class, method, and any module-level constant whose meaning is not obvious.
  Private helpers (leading underscore in Python, `static` or internal in C) get
  documentation only where the reason for their existence is non-obvious.
- Keep documentation lines within the project line length. numpydoc targets 75
  characters for terminal readability; if the project sets nothing, use PEP 8
  (79 for code, 72 for flowing text in docstrings and comments). State the limit
  once and hold it.
- Minimal formatting. The numpydoc section underlines (`-----`) and the Doxygen
  command markers are the only structure you need. No decorative banners, no
  ASCII separators, no bold scattered through prose.

## Pitfalls, TODOs, and warnings

Each kind of hazard has one home, so a reader scanning for danger finds it in
the same place every time, and so the generator can collect the machine-readable
ones into their own lists.

- A pitfall the caller must know to use the unit correctly goes in a `Warnings`
  section (numpydoc free text, distinct from `Warns`, which lists the
  `warnings.warn` categories the unit may raise) or a `@warning` block in
  Doxygen. Examples: not thread safe, mutates its argument, requires a
  contiguous array, energies in Rydberg.
- An unfinished piece of work is an inline `# TODO(owner): ...` in Python or a
  `@todo` in Doxygen. Tie each to an issue and an owner; a bare `# TODO` with no
  owner is noise.
- A known defect is `@bug` in Doxygen; in Python state it in `Warnings` and tag
  the offending line `# FIXME(owner): ...`. Both should reference a filed issue.
- A useful remark that is not a hazard goes in `Notes` or `@note`.
- On the C side, a precondition, postcondition, or invariant goes in `@pre`,
  `@post`, or `@invariant`. In Python state it in the extended summary, or in
  `Raises` if violating it raises.

One rule of phrasing for warnings: state the consequence first, then the
condition. "Returns garbage if `H` is not Hermitian" beats "If `H` is not
Hermitian, then ...". The reader scanning for what can break should learn the
cost in the first clause.

## Structure of a documentation block

The shape is fixed because it works, and numpydoc prescribes the section set and
order, so do not invent sections.

1. Summary line. The single responsibility, imperative, one line.
2. Extended summary. The problem-oriented description: the problem solved, the
   physics or algorithm, its provenance, and the conventions (units, sign,
   shape, storage order). This is the part that is extensive by default.
3. The contract. `Parameters` / `@param`, `Returns` / `@return` (plus `@retval`
   for status codes), `Yields`, and the error conditions in `Raises`, with type
   and shape stated.
4. Hazards and notes. `Warnings` / `@warning`, `Notes` / `@note`, and on the C
   side `@pre`, `@post`, `@invariant`.
5. Cross references. `See Also` / `@see` to the related units a reader will want
   next.
6. References. `References` / `@cite` for the one equation or paper the unit
   implements. Cite the source the way the tutorial guide cites a method: the
   single reference the code rests on, not a bibliography.
7. Examples. `Examples` / `@code`, written as doctests where they can run.

The module or file header carries its own block. A Python module docstring says
what the module is for in terms of the objects it provides. Every C and C++ `.c`
and `.h` file opens with a `@file` and a `@brief`; without them Doxygen will not
emit documentation for the file's globals.

## Naming, so callers can trust the contract

Keep the parameter names in the documentation identical to the signature, letter
for letter. Keep unit names, symbol names, and shape conventions consistent with
the papers the code implements and consistent across modules, so a reader who
learned them in one function reads the next without relearning. When you document
a new unit in an existing module, match the section set the module already uses
rather than inventing a new one.

## Before you commit

A finished documentation pass passes all of these:

- Every public symbol has a block whose summary is one imperative line, ending
  in a period, that does not restate the signature.
- The extended summary states the problem the code solves, with units, sign, and
  shape conventions, not a walk through the implementation.
- Every parameter, return value, and raised error is documented, and the names
  match the signature exactly.
- Every pitfall is in a `Warnings` or `@warning` block, consequence first; every
  TODO has an owner and an issue; every known defect is `@bug` or `FIXME`.
- All mathematics is LaTeX (`` :math:`` ``/`.. math::` or `\f$`/`\f[`), no ASCII
  math, and there are no em dashes.
- Python passes `pydocstyle` or `ruff` for PEP 257 and `pycodestyle` or `ruff`
  for PEP 8, and the numpydoc sections are valid and in the prescribed order.
  C and C++ produce no Doxygen warnings.
- The rendered documentation was not built, because no task asked for it.
- Nothing repeats, and nothing documents what introspection already shows.
