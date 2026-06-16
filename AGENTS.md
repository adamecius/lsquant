# AGENTS.md — agent operating contract (project-agnostic)

This file defines how an agent **behaves**, not what any one project contains.
Drop the same file into every repo; per-project facts live in that repo's own
context doc (see §1), never here. Auto-read by Claude Code; symlink
`CLAUDE.md -> AGENTS.md` and point other agents (Kimi, etc.) at it explicitly.
**Read it fully before acting.**

---

## §0 — Session preamble (run first, every session)

You operate **in place** inside an existing checkout — do not clone elsewhere.
Verify the working tree is the right one and not stale, then read context (§1).

```bash
git rev-parse --is-inside-work-tree            # must print: true
git remote get-url origin                      # must match the remote in the project context doc
git rev-parse --abbrev-ref HEAD                # branch
git rev-parse HEAD                             # HEAD — record it in the PR body
git fetch --quiet origin
git rev-list --left-right --count HEAD...@{u}  # right number = commits behind upstream
```

If the remote is not the expected one, or you are behind upstream, **stop and
report** before doing anything. (Past failure: an agent worked in a stale clone
of a *renamed/old* remote and declared the plan "fiction.")

---

## §1 — Read project context before assuming anything

Do not reconstruct state from the diff. Read the project's context docs first;
they are the source of truth for all project facts (backends, placeholders,
conventions, banned files, phases). On conflict: **this file governs behavior,
the project docs govern facts.** A missing doc → say so, don't invent state.

Expect **one living doc per role**, not a pile of overlapping handoffs:

| Role | Holds | When |
|------|-------|------|
| Orient | build/run + repo map | once, for bearings |
| **State** | current phase · invariants · open items; dated + `verified @ <HEAD>` | **first, every session** |
| Plan | phased plan + STOP conditions (committed, not local-only) | before starting a phase |
| Reference | domain conventions / closed forms / oracle | before any domain change |
| Tests | living gate: green vs pending, target numbers | before & after changes |
| Plotting | figure specs: widths · fonts · palette · caption format | before producing any figure |

**Anti-rot:** the State doc must carry `verified against HEAD <hash> on <date>`.
If that hash is behind the §0 HEAD, treat the State doc as **possibly stale** —
reconcile or flag before trusting its invariants. Point-in-time discovery notes
and superseded docs go in an archive dir, never the read-first set.

---

## §2 — Language, communication, and documentation

- Write **everything in English**: PR bodies, commit messages, code comments,
  docstrings, docs. (The maintainer may converse in another language; agent
  output is English regardless.)
- **Templates (match them; never invent a structure when one exists):**
  - Communications / tutorials → `agents_templates/tutorial_*`
  - In-code comments / PR text → `agents_templates/documentation_*`
  - If a referenced template is absent, flag it rather than guessing.
- **Documentation = in-code comments by default.** They are the source an
  automated generator (Doxygen for C++, Sphinx/NumPy-style for Python) can later
  turn into rendered docs.
  - **Do not run a doc generator or emit generated/rendered docs at any stage
    unless explicitly requested.** When explicitly requested, read both
    `documentation_*` templates first and follow them.
  - When not requested, the in-code comments you write must be an **extensive,
    problem-oriented description** of the code: the problem solved, the approach,
    assumptions, units/conventions, complexity — not a line-by-line restatement.
  - Always include, where applicable: **Pitfall**, **To-do**, **Warning**.

---

## §3 — Change discipline

- **Coverage before change.** If the project has a regression / golden / oracle
  suite, never refactor a path that lacks coverage — add coverage first. A
  refactor step without a safety net is unsafe.
- **Green before commit.** Run the full suite (with any fixed seed the project
  specifies) before every commit; state the pass count in the PR body. Keep
  exact-match goldens byte-identical unless the change is explicitly intended to
  move them, and then justify it.
- **No generated or bulk artifacts.** No build caches, no large binaries, no
  data blobs above the project's size threshold. If you find committed cruft of
  this kind, propose removing and gitignoring it — don't add more.
- **Respect placeholders.** Intentional throwing stubs / not-implemented guards
  stay as they are unless the task is to implement them. Don't silently revive
  or re-stub them.

---

## §4 — PR and merge-commit contract

The durable record lives in the **squash-merge commit body**, not only the
GitHub PR UI (PR text can be rate-limited or unreachable; `git log` is forever).
Every PR body **and** the squash-merge message must contain all four:

1. **Why these files changed** — per file or logical group, the *reason*, not the
   *what* (the diff shows what).
2. **Tests** — which targets ran, the seed, the pass count, and confirmation that
   exact-match goldens stayed byte-identical. Name any new gate added.
3. **Reason for acceptance** — why this is safe to merge now (behavior-preserving?
   covered? gated?), tied to a project invariant (§1) or a plan phase.
4. **To-do / follow-ups** — what this PR deliberately leaves open and what it
   unblocks. Note any placeholder touched.

Record the HEAD hash from §0 and the relevant commit hashes.

---

## §5 — Reference material (mild RAG, never committed)

Reference books/papers live in `agent_references/` — **gitignored, local only.**
Only the manifest `agent_references.md` is committed; it carries title + topics +
a "use-for" tag per entry so you can grep it to pick what to open. Never
`git add` a reference file; never paste its content into the repo. When a change
needs a literature justification, grep the manifest by topic, open only the
matching local file, and cite it in the PR body by title.

`.gitignore`:  `agent_references/`

Manifest entry format (`agent_references.md`):
```
## <Author, Venue Year> — <Title>
- file:    agent_references/<filename>   (local only, not committed)
- topics:  <semicolon-separated topics>
- use-for: <when an agent should reach for this>
```

---

## §6 — Definition of done

§0 verification recorded · full suite green at the project's fixed seed with
goldens byte-identical · PR body + squash-merge message satisfy all four points
of §4 · no banned artifacts committed (§3, §5) · output in English (§2) ·
project context doc / plan updated if the change advances a phase.

## §7 — Figures and plots

Behavior here; the **specifications** (column widths, fonts, palette, encoding,
caption format) live in the project's `plotting_style.md`, a Reference-role doc
(§1). On conflict: this section governs *when/how you act*, the style doc governs
*the numbers*.

- **Trigger.** Any figure you generate or edit, for any surface — article,
  GitHub Pages, slides, README. No exceptions for "quick" or "throwaway" plots.
- **Read the style doc first.** Before producing a figure, read the project's
  `plotting_style.md` and use its presets. Missing doc → **say so, don't invent
  defaults** (§1). Do not reconstruct style from existing figures.
- **Figure source is code, not raster.** Commit the generating script; rendered
  figures are regenerable cache. Do not commit rasters above the project size
  threshold or any bulk image blob (§3, §5). Export final-size artifacts only
  where the project explicitly collects them.
- **No title; the caption carries everything.** Never set a figure/panel title.
  Every figure ships **with** a self-contained caption (begins `FIG. N.`, defines
  each colour/linestyle/marker, and lists the parameters needed to reproduce the
  panel). **A figure delivered without its caption is not done** — emit the
  caption in the `.tex`, the PR body, or a sidecar per the style doc.
- **Redundant encoding is mandatory.** Traces separable by colour **and**
  linestyle **and** marker; legible in grayscale and under colour-vision
  deficiency. Colour-only encoding, `jet`/rainbow, and sole red–green contrast
  are rejected.
- **Self-check gate.** Run the style doc's acceptance checklist and report
  PASS/FAIL in the PR body next to the test gate (§3, §4). Any failed box →
  not shipped. Treat a non-compliant figure like a red suite: fix or escalate.
- **No post-export rescaling.** Export at the target width and include at natural
  size; rescaling on `\includegraphics` breaks the lettering-height rule.