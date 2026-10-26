# Paper Writing Pipeline — High-Performance Configuration

**Author:** José Hugo García Aguilar (ICN2)  
**Date:** 2026-10-26  
**Status:** Optimized for maximum quality and minimal context contamination

## Overview

This pipeline uses specialized Grok Bots for scientific paper writing. It follows the best-performing multi-agent patterns from PaperOrchestra, Story2Proposal, and PaperOrchestrator (2026).

**Core principle:** Specialization + shared state + human gates. Never let one agent do everything.

## Agents

| Bot | Role | Responsibilities | Key Constraints |
|-----|------|------------------|-----------------|
| **Paper Architect** | Planner | Generate structured outline, claims, figure plan, citation needs. Create paper folder. | Does **not** write full prose. |
| **Article Writer Clean** (Drafter) | Writer | Write full draft from outline + technical notes. Match author's style. | Clean context. Follow outline closely. |
| **Paper Refiner** | Reviewer | 1–2 focused passes: coherence, clarity, style, technical precision, consistency. | Only keep real improvements. |
| **Literature Scout** | Support | Find and verify canonical references. Provide BibTeX. | High-impact works only. |

## Shared State

All materials live in:

```
/root/Agents/Article Writer/Draft/<PaperShortName>/
```

Expected files:
- `outline.md`
- `draft_v1.md`
- `draft_v2.md` (after refinement)
- `references.bib`
- `review_notes.md` (optional)
- `figures/` (optional)

## Recommended Flow

1. **User → Architect**: Provide idea, results, key claims.
2. **Architect** creates folder + `outline.md`.
3. **Human gate** (short): Review and adjust outline.
4. **Drafter** writes `draft_v1.md`.
5. **Refiner** produces improved version.
6. **Human gate** (final review).

Optional: Call **Literature Scout** at any point for citations.

## Performance Rules

- Each bot must stay focused on its role.
- Do not regurgitate previous conversation context.
- Prefer precise scientific English.
- Quantify claims when possible.
- Keep sections modular.
- Always save outputs in the paper’s folder.

## How to Use

1. Start with the Architect for a new paper.
2. After outline approval, send outline + notes to the Drafter.
3. Send the draft to the Refiner.
4. Iterate only if needed.

This configuration prioritizes quality, consistency, and reduced confirmation bias.
