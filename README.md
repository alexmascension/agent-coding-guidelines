# agent-coding-guidelines

Project-agnostic guidelines for AI coding agents (and the people who review their work) in bioinformatics and
computational biology, with tested reference implementations of the non-obvious procedures.

The guides are written for agents but are meant to be readable by a person with intermediate statistics. Each rule
states what to do and why, its cost and when to deviate; claims cite the literature (`OMICS_REFERENCES.md`,
`STATISTICS_REFERENCES.md`).

## Contents

| File | Scope | Status |
|---|---|---|
| `guidelines/CODING_GUIDELINES.md` | Code, data, provenance, reproducibility, testing, project documentation (`AGENTS.md`, `PROJECT_CONTEXT.md`) | stable |
| `guidelines/AGENTIC_GUIDELINES.md` | Agent working discipline: scope, evidence, verification, claims, decisions, unattended batches, delegation, confidential material, handoffs, iteration | new |
| `guidelines/MANUSCRIPT_GUIDELINES.md` | Manuscripts built on the analyses: traceability of figures, tables and numbers, claims, citations, pre-submission review, revisions, writing | new |
| `guidelines/STATISTICS_GUIDELINES.md` | Design, test choice, multiple testing, resampling, reporting | stable |
| `guidelines/RNASEQ_BULK_GUIDELINES.md` | Bulk RNA-seq from reads to differential expression and enrichment | stable |
| `guidelines/SCRNASEQ_GUIDELINES.md` | Single-cell and single-nucleus RNA-seq (10x, STARsolo, scverse) | stable |
| `guidelines/SPATIAL_GUIDELINES.md` | Visium, Xenium, MERFISH | stable |
| `guidelines/META_ANALYSIS_GUIDELINES.md` | Combining studies, general and omics | stable |
| `guidelines/SEX_INFERENCE_GUIDELINES.md` | Sex-chromosome complement from expression; reference routing | stable |
| `guidelines/SAMPLE_IDENTITY_GUIDELINES.md` | Duplicates, swaps, mixtures, genotype fingerprints | stable |
| `guidelines/DECONVOLUTION_GUIDELINES.md` | Cell-type composition of bulk samples | stable |
| `guidelines/HPC_GUIDELINES.md` | Shared clusters with Slurm: storage roles, staging and promotion, containers, jobs, agents on clusters | new |
| `guidelines/hpc_sites/` | Site profiles (facts of one cluster, dated): `DIPC_HYPERION.md` | — |
| `guidelines/FIGURE_GUIDELINES.md` | Figures | planned |
| `references/` | `OMICS_REFERENCES.md` / `.bib`, `STATISTICS_REFERENCES.md` / `.bib`: references with keys used in the guides | — |
| `recipes/` | Reference implementations and tests; coverage in `recipes/README.md` | 0.4.0 |
| `CHANGELOG.md` | Changes per guide (semantic versioning per file) | — |
| `RULE_IDS.tsv` | Registry of rule identifiers (active and retired) | — |

Scope: human and mouse; bulk RNA-seq, sc/snRNA-seq, spatial transcriptomics. Python first (scverse); R through
scripts with a file interface when no validated Python implementation exists. Not covered: CITE-seq, multiome,
TCR/BCR.

## Use in a project

1. Copy `guidelines/CODING_GUIDELINES.md` and the domain guides the project needs, **verbatim**, into `guidelines/`
   at the project root; for clusters, the site profile into `guidelines/hpc_sites/`; the reference files from
   `references/` into `guidelines/references/` (`CODING_GUIDELINES.md` §0, `CODE-USE-001`):

   ```bash
   UPSTREAM=/path/to/agent-coding-guidelines   # a clone at a pinned commit
   mkdir -p guidelines/references
   cp "$UPSTREAM"/guidelines/{CODING,AGENTIC,STATISTICS,SCRNASEQ}_GUIDELINES.md guidelines/
   cp "$UPSTREAM"/references/* guidelines/references/
   mkdir -p guidelines/hpc_sites && cp "$UPSTREAM"/guidelines/hpc_sites/DIPC_HYPERION.md guidelines/hpc_sites/  # clusters only
   ```
2. Create `AGENTS.md` (entry point for agents; `CLAUDE.md` may contain `@AGENTS.md`) and `PROJECT_CONTEXT.md` with
   the structure of `CODING_GUIDELINES.md` §2. Adaptations go to *Overrides* in `PROJECT_CONTEXT.md`, never into the
   copied guides; *Overrides* starts with the list of adopted guides and their versions, and refers to rules by
   their stable identifier (e.g. `SC-QC-003`; registry in `RULE_IDS.tsv`).
3. For recipes, install this repository at a pinned commit (see `recipes/README.md`) or vendor the modules as
   described in `CODING_GUIDELINES.md` §26 (recipes).
4. To update, compare the project's copies with upstream using `CHANGELOG.md`; MAJOR versions may renumber or remove
   rules that overrides point to.

## Development

```bash
pip install -e ".[scrnaseq,bulk,dev]" pre-commit
pre-commit install              # once per clone: the checks below run on every commit
pre-commit run --all-files      # line endings, YAML/TOML, large files, ruff, shellcheck, rule IDs, guides, versions
pytest -m "not r"               # Python recipe tests
pytest -m r                     # R-backed tests: need RECIPES_RSCRIPT (or Rscript on PATH) with the R packages
```

The pre-commit hooks run `tools/rule_ids.py assign` and `check` (identifiers and `RULE_IDS.tsv`),
`tools/check_guides.py` (citation keys, cross-references, recipe paths, versions) and `tools/check_versions.py`
(a changed guide bumps its version and has a `CHANGELOG.md` entry, `CODE-USE-007`).

CI (`.github/workflows/ci.yml`) runs on pushes to `main`, pull requests, weekly and on demand:

| Job | Runs when | What |
|---|---|---|
| Guides and lint | always | the pre-commit hooks on all files; version and changelog check against the previous commit |
| Recipes (Python 3.11, 3.12) | a file under `recipes/` or `tools/`, `pyproject.toml` or the CI changed | `pytest -m "not r"` (recipe and tool tests) |
| Install from the built wheel | as above | install in a clean environment; import and package data (R scripts) |
| Recipes backed by R | an `*.R` file changed, or a file that depends on R (`tools/ci_scope.py`) | R with scran, scDblFinder, tximport, DESeq2; `pytest -m r`, where a skip is a failure |

`tools/ci_scope.py` decides the last two: a file depends on R when it calls R or names an R script, when a test
marked `r` uses it, or when such a file imports it; it checks the files before and after the change, so adding or
removing a call to R also triggers the R job. Scheduled and manual runs run every job.

## License

MIT (`LICENSE`), for the guidelines, references and code.

`guidelines/AGENTIC_GUIDELINES.md` adapts text from the Co-Scientist toolkit by the Shocklab, University of Cape
Town (<https://github.com/shocklab/co-scientist-toolkit>), licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/): its working rules (`guidance/CLAUDE-starter.md`,
`guidance/paper-CLAUDE.md`, `guidance/handoff-CURRENT.md`) and the skills `fable-mode`, `pre-register`,
`convergence-check`, `decision-page` and `claim-fact-checker`. The text was rewritten as rules, generalised and
extended; §0 of the guide lists the sources.

`guidelines/MANUSCRIPT_GUIDELINES.md` adapts text from the same toolkit (CC BY 4.0): `guidance/paper-CLAUDE.md`
(itself adapted from Dominik Lukeš's agents-for-reproducibility course), the *Writing* and *Research papers*
sections of `guidance/CLAUDE-starter.md`, and the skills `paper-review` (review templates adapted from Dominik
Lukeš's paper-reviewer skill, MIT), `claim-fact-checker` (adapted from Dominik Lukeš's agent skills, MIT),
`reference-verifier` and `convergence-check`. §0 of the guide lists the sources.

The rules added to `guidelines/CODING_GUIDELINES.md` in 3.7.0 (`CODE-LANG-004` shell traps, `CODE-DOC-005`,
`CODE-NUM-013`, `CODE-PROV-006`, `CODE-TEST-009`, `CODE-TEST-010`, `CODE-VAL-004`) adapt the *Shell commands*,
*Verifying work* and *Writing* sections of `guidance/CLAUDE-starter.md` and the results archive of
`guidance/paper-CLAUDE.md` from the same toolkit (CC BY 4.0).

`STAT-DES-010` to `STAT-DES-012` in `guidelines/STATISTICS_GUIDELINES.md` 2.1.0 adapt the forking-paths guard, the
append-only deviations and the pilot decision rule of the toolkit's `pre-register` skill (CC BY 4.0).
