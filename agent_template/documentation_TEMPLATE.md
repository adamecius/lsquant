<!--
  LSQUANT in-code documentation template.
  Copy the skeleton for your language into your source file and fill the
  [[ ... ]] slots. The <!-- ... --> comments are guidance for the author; they
  never reach the reader, so delete them or leave them in this file.

  Read DOCSTRING_STYLE_GUIDE.md once before writing. The short version:
    - lead with the problem, earn the block, never restate the signature
    - one responsibility per unit; extensive problem-oriented extended summary
    - Python: PEP 8 + PEP 257, numpydoc sections; C/C++: Doxygen, Javadoc style
    - all math in LaTeX, no em dashes
    - pitfalls in Warnings/@warning (consequence first), TODOs owned, bugs filed
    - write the in-code docs; never build the rendered site unless asked
-->

# Documentation template

This file holds one fillable skeleton per language. Copy the relevant one, fill
the `[[ ... ]]` slots, and delete the slots that do not apply.

## Python (PEP 257 + numpydoc)

### Module

<!--
  The module docstring is the first statement in the file. It says what the
  module is for in terms of the objects it provides, not how they are built.
  Triple double quotes; summary line on its own; closing quotes on their own
  line.
-->

```python
"""[[One imperative line naming what this module provides.]]

[[Extended summary: the problem this module addresses and the physical or
computational objects it exposes. Name the conventions that hold module wide
(units, sign, storage order) so each function does not repeat them.]]

Notes
-----
[[Optional: provenance, the paper or method the module implements, in LaTeX
where there is math, e.g. the moments :math:`\\mu_m = \\mathrm{Tr}\\,T_m(H)`.]]
"""
```

### Function

<!--
  Summary in the imperative, ending with a period, not a restatement of the
  signature. Then the problem-oriented extended summary. numpydoc section order
  is fixed: Parameters, Returns (or Yields), Raises, Warns, Warnings, See Also,
  Notes, References, Examples. Use only the sections you need, in that order.
-->

```python
def [[name]]([[args]]):
    r"""[[Imperative one-line summary of the single responsibility.]]

    [[Extended summary: the problem this solves and why it exists. State the
    physics or algorithm and the conventions. Put math in LaTeX, e.g.

    .. math:: \rho(E) = \sum_n \delta(E - E_n).
    ]]

    Parameters
    ----------
    [[name]] : [[type, with shape if an array, e.g. ndarray of shape (n, n)]]
        [[What it means physically and its units. Not just its type.]]
    [[name]] : [[type]], optional
        [[Meaning and the default's effect.]]

    Returns
    -------
    [[name]] : [[type and shape]]
        [[What the value represents physically and its units.]]

    Raises
    ------
    [[ExceptionType]]
        [[The condition that triggers it, consequence first.]]

    Warnings
    --------
    [[The pitfall a caller must know, consequence first: e.g. "Returns
    garbage if `H` is not Hermitian; the routine does not check."]]

    Notes
    -----
    [[Optional: derivation, complexity, or a subtlety that would otherwise be
    misread. The sentence that stops the obvious wrong conclusion lives here.]]

    References
    ----------
    [[The one method paper the unit implements, e.g. A. Weisse et al.,
    Rev. Mod. Phys. 78, 275 (2006).]]

    Examples
    --------
    >>> [[a runnable doctest that shows the intended call and result]]
    """
    # TODO([[owner]]): [[unfinished work, tied to issue #NNN]]  # remove if none
    [[...]]
```

### Class

<!--
  Document the class for the role it plays. Document each public method with its
  own function-style docstring. Attributes the caller reads or sets go in an
  Attributes section.
-->

```python
class [[Name]]:
    """[[Imperative one-line summary of the role this type plays.]]

    [[Extended summary: the problem this type models and the invariant it
    maintains across its methods.]]

    Parameters
    ----------
    [[name]] : [[type]]
        [[Constructor argument: meaning and units.]]

    Attributes
    ----------
    [[name]] : [[type]]
        [[A public attribute the caller may read; meaning and units.]]

    Warnings
    --------
    [[A type-wide hazard, e.g. "Instances are not thread safe; do not share one
    across workers."]]
    """
```

## C and C++ (Doxygen, Javadoc style)

### File header

<!--
  Every .c, .h, .cpp, and .hpp opens with this. Without @file and @brief,
  Doxygen will not document the file's global functions, typedefs, or enums.
  Math in LaTeX with \f$ ... \f$ inline and \f[ ... \f] displayed.
-->

```c
/**
 * @file [[name.h]]
 * @brief [[One line naming what this file provides.]]
 *
 * @details [[The problem this file addresses and the objects it exposes,
 *           with the conventions (units, sign, storage order) that hold
 *           across the file. Math in LaTeX, e.g. \f$ H = H^\dagger \f$.]]
 *
 * @author  [[name]]
 */
```

### Function

<!--
  @brief is the imperative one-line summary; @details is the problem-oriented
  description. Use @param[in], @param[out], @param[in,out] to mark direction.
  @retval documents specific status codes; @return documents the return value.
  Hazards and contract: @warning (consequence first), @pre, @post, @note.
-->

```c
/**
 * @brief [[Imperative one-line summary of the single responsibility.]]
 *
 * @details [[The problem this solves and why it exists, with the algorithm and
 *           its conventions. Displayed math in LaTeX:
 *           \f[ \mu_m = \mathrm{Tr}\, T_m(H). \f] ]]
 *
 * @param[in]     [[name]] [[Meaning and units, with array length if a pointer.]]
 * @param[out]    [[name]] [[What the routine writes here and its units.]]
 * @param[in,out] [[name]] [[State on entry and on exit.]]
 *
 * @return [[What the return value represents.]]
 * @retval [[CODE]] [[The condition that yields this status code.]]
 *
 * @pre  [[A precondition the caller must satisfy, e.g. n > 0.]]
 * @post [[A guarantee that holds on success.]]
 *
 * @warning [[The pitfall, consequence first: e.g. "Writes out of bounds if
 *          `n` exceeds the allocation; the routine does not check."]]
 * @note    [[A useful remark that is not a hazard.]]
 *
 * @todo [[Unfinished work, tied to a filed issue.]]   // remove if none
 * @bug  [[A known defect, tied to a filed issue.]]    // remove if none
 *
 * @see  [[a related function the reader will want next]]
 */
[[return_type]] [[name]]([[parameters]]);
```

### Struct or class

<!--
  Document the type for the role it plays. Document each member with a trailing
  /**< ... */ where its meaning is not obvious from the name.
-->

```c
/**
 * @brief [[One line naming the role this type plays.]]
 *
 * @details [[The problem this type models and the invariant its fields keep.]]
 *
 * @warning [[A type-wide hazard, if any.]]
 */
typedef struct {
    [[type]] [[name]];   /**< [[meaning and units]] */
    [[type]] [[name]];   /**< [[meaning and units]] */
} [[Name]];
```
