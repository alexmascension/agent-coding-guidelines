# Maintainability resources and adaptation notes

Reviewed on 2026-11-16 for `CODING_GUIDELINES.md` 4.2.0. These are engineering recommendations and
project-specific design choices, not evidence of a measured improvement in this repository's usability.
The additions use original wording; no external source text or code is copied.

## Review of the existing guidelines

The baseline already prioritizes readability over brevity/performance (§1), cohesive functions and descriptive
names (§4), contracts and rationale (§9), explicit side effects (§20), and fresh-kernel notebooks (§22).
`AGENTIC_GUIDELINES.md` already incorporates Co-Scientist toolkit and superpowers practices for scoped work,
review and verification; their upstream provenance remains in `SOURCE_REPOS.toml`. Another agent workflow
would duplicate those rules. The gap addressed here is the human reader's navigation and next edit.

## Sources consulted

| Resource | Material read | Application |
|---|---|---|
| [Aleksey Kladov: ARCHITECTURE.md (2021)](https://matklad.github.io/2021/02/06/ARCHITECTURE.md.html) | Entire short article | A coarse code map, important boundaries and named modules; `CODE-PROJ-010` to `012`. |
| [Google Engineering Practices: What to look for in a code review](https://google.github.io/eng-practices/review/reviewer/looking-for.html) | Design, complexity, tests, naming, comments and documentation sections | Evaluate how easily a reader can understand and modify code; avoid speculative abstractions; `CODE-PRIN-006` to `008`, `CODE-VAL-005`. |
| [Rule et al. (2019): Ten simple rules for writing and sharing computational analyses in Jupyter Notebooks](https://doi.org/10.1371/journal.pcbi.1007007) | Full-text HTML, especially rules 1–5 and 7 | Audience, narrative, meaningful cells, modularization and execution; `CODE-NB-006` to `009`. |
| [Wilson et al. (2017): Good enough practices in scientific computing](https://doi.org/10.1371/journal.pcbi.1005510) | Full-text HTML, software, collaboration and project organization sections | Small comprehensible steps, explicit requirements and a runnable small example; `CODE-DOC-017`, reading-guide walkthrough. |

## Deliberate adaptations

- Kladov discusses larger projects and recommends infrequent maintenance and symbol names instead of links.
  Here, a separate map is conditional; stable file links and symbols are checked when the relevant code changes.
  This follows the repository's existing documentation synchronization rule and the owner's navigation needs.
- The template separates calls/dependencies from data movement. The common-change table, external dependency
  boundary table and handoff reading route are local design choices, not requirements claimed by the sources.
- The notebook paper permits more exploratory styles. These additions retain this repository's stricter
  package-logic and stripped-output rules. Inspect execution before stripping outputs, and report unavailable
  data or unexecuted expensive steps rather than weakening validation claims.
- No fixed function-length, cell-count or complexity score is introduced. A coherent reading path matters more
  than splitting code into many small pieces. These defaults need judgement and feedback from the maintainer.

## Access and evidence limits

All four sources above were accessible as HTML; no source needed for these changes remains blocked.
The initial web fetch of the target GitHub repository failed, but GitHub connector access and a local clone
succeeded. No paywalled book was used or represented as having been read. The changes have not been evaluated
in a controlled study or a fresh agent session; structural checks cannot establish human comprehension.
