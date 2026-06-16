# documentation_incode — in-code comment template

Scope: **in-code documentation (comments/docstrings)**. These comments are the
*source* an automated generator (Doxygen for C++ — the project default; Sphinx /
NumPy-style for Python) can later turn into rendered docs.

## When the generator runs
- **Default: never.** Do not run a doc generator or emit generated/rendered
  documentation at any stage **unless explicitly requested.**
- **When explicitly requested:** read this file and `documentation_pr.md`, then
  follow them.

## What to write by default (not requested)
For the code you wrote, write **an extensive, problem-oriented description** as
in-code comments — describe the *problem* the code solves, the approach taken,
assumptions, units/conventions, and complexity. Do **not** restate what the code
literally does line by line. Always include, where applicable: **Pitfall**,
**To-do**, **Warning**.

## C++ skeleton (Doxygen — matches the existing `@brief` style)
```cpp
/**
 * @brief One-line summary.
 *
 * Problem-oriented description: the problem this solves, the chosen approach,
 * assumptions, units/conventions, complexity. Not a line-by-line restatement.
 *
 * @param  x   meaning + units/range
 * @return     meaning + units
 *
 * @warning Hard constraint; violating it is UB / a correctness break.
 * @par Pitfall
 *      The subtle gotcha that has bitten (or will bite) a caller here.
 * @todo Concrete follow-up, ideally tied to a plan phase.
 */
```
> `@par Pitfall` works as-is. To get a first-class `@pitfall` tag, add to the
> Doxyfile: `ALIASES += pitfall="@par Pitfall:\n"`.

## Python skeleton (NumPy-style docstring)
```python
def f(x):
    """One-line summary.

    Problem-oriented description: the problem solved, the approach, assumptions,
    units/conventions, complexity. Not a restatement of the code.

    Parameters
    ----------
    x : <type>
        Meaning + units/range.

    Returns
    -------
    <type>
        Meaning + units.

    Warnings
    --------
    Hard constraint; violating it breaks correctness.

    Notes
    -----
    Pitfall: the subtle gotcha a caller will hit here.

    .. todo:: Concrete follow-up.
    """
```
