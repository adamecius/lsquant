<!--
  LinQT tutorial template.
  Copy this file to examples/NN_short_name/README.md and fill the [[ ... ]] slots.
  The <!-- ... --> comments are guidance for the author. They are HTML comments,
  so they stay invisible when the page renders on GitHub. Delete them or leave
  them; readers never see them.

  Read STYLE_GUIDE.md once before writing. The short version:
    - lead with the physics, earn every command with a physical reason
    - one controlling idea per tutorial
    - all math in LaTeX, all runnable lines in code blocks
    - no em dashes, no setup or PATH or pip paragraphs
    - embed every figure; its alt text states the takeaway
-->

# Tutorial [[N]]: [[the physical question, in plain words]]

<!--
  THE HOOK. Two short paragraphs, no headings.
  Paragraph 1: pose the phenomenon as a small puzzle a curious person would feel.
  Paragraph 2: name what we are after and the single idea this tutorial teaches.
  Do not mention the tool, the build, or the file formats here. Only the physics.
-->

[[Open with the phenomenon. Make the reader curious before they meet a command.]]

[[State plainly what we want to compute and the one lesson that will survive into
later tutorials. End with a sentence of the form "the lesson here is ...".]]

## The physics

<!--
  The minimum theory needed to read the steps, no more. Define any term that
  steps outside a working physicist's vocabulary. All equations in LaTeX.
  Put the ONE conceptual result the tutorial turns on in its own sentence, and
  say in words why it matters before any code appears.
-->

[[Define the system and its Hamiltonian.]]

$$ [[ H = \dots ]] $$

[[Give the closed form or known limit we will compare against, in LaTeX.]]

$$ [[ \text{closed form, e.g. } \rho(E) = \dots ]] $$

[[State the key conceptual result of the tutorial and what it means physically.]]

![[[alt text that states the takeaway of this opening figure]]](fig_[[name]].png)

## Step 1: build [[the system]]

<!--
  Light. One command per size or case. One sentence on what it produces, framed
  as physical objects (a Hamiltonian, its band edges, a state set), never as a
  walk through file formats. The reader does not open these files.
-->

```bash
[[python make_*.py <args>]]
```

[[One sentence: what physical objects this writes, and that the next step reads them.]]

## Step 2: [[compute the moments]]

<!--
  The expensive kernel call. State which operator and why, then the command.
  Show the one output the reader will refer to. Do not explain threading,
  C++, or "this is the heavy part".
-->

[[One sentence naming the operator and the quantity its spectral function gives.]]

```bash
[[inline_compute-kpm-spectralOp <label> <op> <M> <state>]]
```

[[One sentence: this evaluates <the trace or estimator>, and writes one file:]]

```
[[output file name pattern]]
```

## Step [[k]]: [[a physical statement, not a verb]]

<!--
  Each later step is titled with a physical claim the step demonstrates, e.g.
  "the moments are intensive", "resolution follows the moment count". Run a
  command, embed the figure, then read the figure as physics.
-->

```bash
[[python -c "import lsqplot; lsqplot.figure_*( ... )"   or   python lsqplot.py]]
```

![[[alt text stating what this figure shows]]](fig_[[name]].png)

[[Read the figure: what the curves do, and the physical reason. If there is a
subtlety the reader would otherwise get wrong, name it plainly here.]]

<!-- Repeat the "Step k" block as needed. Keep the count small. -->

## What to take away

<!-- Three to five bullets. Each is a physics statement, not a software note. -->

- [[the central result, stated as physics]]
- [[the scaling or invariance lesson]]
- [[the resolution or accuracy lesson]]
- [[the provenance or operator lesson, if relevant]]

[[One sentence linking forward: what the next tutorial reuses from this one.]]

## References and links

- LinQT source and documentation: https://github.com/adamecius/lsquant
- Methodology: Z. Fan, J. H. García, A. W. Cummings et al., *Linear Scaling
  Quantum Transport Methodologies*, arXiv:1811.07387.
- Installation: see the main README of the repository.

<!--
  Optional, only if the tutorial leans on a specific method or result:
  ## Further reading
  - [[author, title, venue, year]]
  e.g. A. Weisse et al., The kernel polynomial method, Rev. Mod. Phys. 78, 275 (2006).
-->
