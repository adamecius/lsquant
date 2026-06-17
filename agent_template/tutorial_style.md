# How to write a LSQUANT tutorial

This guide explains how to turn `TEMPLATE.md` into a finished tutorial, and it
is opinionated about voice on purpose. Tutorial 1 (`examples/01_dos_1d_chain`)
is the worked reference; read it, then read this.

## The voice

Write as if you were explaining the physics to a smart friend who has never met
the code, in the spirit of Sagan, Feynman, or Tyson. That is not decoration. It
is a discipline with consequences for every paragraph.

Lead with the phenomenon, not the tool. The reader should feel a small puzzle
before they ever see a command. "A finite chain has $N$ levels, yet we draw a
smooth curve; where did the curve come from?" earns attention. "This tutorial
computes the DOS with `inline_compute-kpm-spectralOp`" does not.

Carry one idea. A tutorial is not a tour of features. It is the unfolding of a
single physical lesson, and every step exists to advance that one lesson. If a
step does not move the idea forward, cut it. Tutorial 1 has exactly one spine:
the moment count is the resolution dial, and finiteness hides in one revival.

Earn every command. No command appears without a physical reason for running it.
The reader should always know what question the next line answers. The code is
the answer to a physics question, never the subject.

Respect the reader's expertise, and only theirs. Assume a working physicist who
programs. Do not explain Chebyshev polynomials or what a trace is. Do explain
anything that steps outside that vocabulary, and define a term the first time it
earns its place. Plain language everywhere that is not physics.

Trust the phenomenon to be interesting. You do not need to tell the reader that
a result is "striking" more than once, and you never need filler. Compact prose
reads as confidence. Repetition reads as padding.

Name the trap. The most useful sentence in a tutorial is often the one that
stops a smart reader from drawing the obvious wrong conclusion. In tutorial 1
that sentence is: while only $\mu_0$ is nonzero, the broadening does nothing,
because there is nothing above $m=0$ to damp. Find that sentence in your topic
and give it room.

## Hard rules

These are not stylistic preferences. Apply them without exception.

- All mathematics is LaTeX: `$...$` inline, `$$...$$` for displayed equations.
  Never write equations as plain ASCII in prose or in code fences.
- Everything the reader runs goes in a code block. Outputs the reader will refer
  to (a file name, a header) also go in a plain code block.
- No em dashes anywhere. Use a comma, a colon, parentheses, or a full stop.
- No setup paragraphs. Do not write about `PATH`, `LSQUANT_BIN`, `pip install`,
  build directories, threading, or which language the kernel is in. Installation
  lives in the main README, and the tutorial points there once, in References.
- Every figure is embedded with a relative path, and its alt text states the
  takeaway, not the file name. `![All sizes overlap at the m=0 moment](fig_moments.png)`,
  not `![figure 2](fig_moments.png)`.
- Minimal formatting. No bold scattered through prose, no nested bullets. Headings
  and the occasional short list, nothing more.

## Structure

Follow `TEMPLATE.md`. The shape is fixed because it works:

1. Title as a physical question in plain words.
2. A two paragraph hook: the puzzle, then what we are after and the one lesson.
3. `## The physics`: the minimum theory, the closed form to compare against, the
   one conceptual result stated in words, and an opening figure that shows the
   problem.
4. `## Step 1: build ...`: light, the system as physical objects.
5. `## Step 2: ...`: the kernel call that produces the moments.
6. `## Step k: <a physical claim>`: later steps are titled with the claim they
   demonstrate, run a command, embed a figure, read the figure as physics.
7. `## What to take away`: three to five bullets, each a physics statement.
8. One sentence forward to the next tutorial.
9. `## References and links`: the fixed footer below.

Step titles are claims, not verbs. "The moments are intensive" tells the reader
what they are about to learn. "Plot the moments" tells them nothing.

## Figures

Figures are the argument, so make them carry it.

Ship a small script with each tutorial that produces every embedded PNG, either
by calling the compiled kernels through `lsqplot.py` or, where the physics has a
closed form, by computing it directly. The two agree to machine precision when
the moments come from an exact trace, so a closed-form figure is faithful, not a
mock-up. Commit the PNGs next to the README so the page renders offline.

Name figures for their content: `fig_levels.png`, `fig_moments.png`,
`fig_dos.png`. Keep one idea per figure or per panel. A two panel figure should
make a single comparison the reader can state in one sentence. Put the physical
parameters the reader needs (the size, the moment count, the broadening) into
the legend, so the figure is readable on its own.

## Naming, so commands stay copy-pasteable

Keep labels and file patterns consistent across tutorials, because readers paste
commands and expect the output names to match. Tutorial 1 uses `chain1d_N<N>`
for the system label, operator label `1` for the identity, and the moment file
pattern `SpectralOp<op><label>KPM_M<M>_state<state>.chebmom1D`. When you add a
system, extend the convention rather than inventing a new one.

## References

Every tutorial ends with the same footer, nothing more by default:

```markdown
## References and links

- LSQUANT source and documentation: https://github.com/adamecius/lsquant
- Methodology: Z. Fan, J. H. García, A. W. Cummings et al., *Linear Scaling
  Quantum Transport Methodologies*, arXiv:1811.07387.
- Installation: see the main README of the repository.
```

If, and only if, the tutorial leans on a specific method or result, add a short
`## Further reading` with that one citation. The KPM kernel paper (Weisse et al.,
Rev. Mod. Phys. 78, 275, 2006) belongs there when a tutorial uses kernel details.
Do not build a bibliography. Two or three links is the ceiling.

## Before you commit

A finished tutorial passes all of these:

- A reader could state the one lesson in a single sentence after reading.
- Every command has a physical reason stated before it.
- Every equation is LaTeX, every runnable line is in a code block.
- Every figure is embedded, with alt text that states its takeaway, and the PNG
  is committed.
- There are no em dashes and no setup, PATH, pip, or build paragraphs.
- The references footer is the standard three lines.
- Nothing repeats, and nothing is there that the lesson does not need.
