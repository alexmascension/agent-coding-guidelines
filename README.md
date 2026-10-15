# agent-coding-guidelines

Project-agnostic guidelines for AI coding agents (and the people who review their work) in bioinformatics and
computational biology, with tested reference implementations of the non-obvious procedures.

The guides are written for agents but are meant to be readable by a person with intermediate statistics. Each rule
states what to do and why, its cost and when to deviate; claims cite the literature (`OMICS_REFERENCES.md`,
`STATISTICS_REFERENCES.md`).

## Contents

| File | Scope | Status |
|---|---|---|
| `CODING_GUIDELINES.md` | Code, data, provenance, reproducibility, testing, project documentation (`AGENTS.md`, `PROJECT_CONTEXT.md`) | stable |
| `STATISTICS_GUIDELINES.md` | Design, test choice, multiple testing, resampling, reporting | stable |
| `RNASEQ_BULK_GUIDELINES.md` | Bulk RNA-seq from reads to differential expression and enrichment | stable |
| `SCRNASEQ_GUIDELINES.md` | Single-cell and single-nucleus RNA-seq (10x, STARsolo, scverse) | stable |
| `SPATIAL_GUIDELINES.md` | Visium, Xenium, MERFISH | stable |
| `META_ANALYSIS_GUIDELINES.md` | Combining studies, general and omics | stable |
| `SEX_INFERENCE_GUIDELINES.md` | Sex-chromosome complement from expression; reference routing | stable |
| `SAMPLE_IDENTITY_GUIDELINES.md` | Duplicates, swaps, mixtures, genotype fingerprints | stable |
| `DECONVOLUTION_GUIDELINES.md` | Cell-type composition of bulk samples | stable |
| `FIGURE_GUIDELINES.md` | Figures | planned |
| `OMICS_REFERENCES.md` / `.bib`, `STATISTICS_REFERENCES.md` / `.bib` | References with keys used in the guides | — |
| `recipes/` | Reference implementations and tests; coverage in `recipes/README.md` | 0.2.0 |
| `CHANGELOG.md` | Changes per guide (semantic versioning per file) | — |
| `RULE_IDS.tsv` | Registry of rule identifiers (active and retired) | — |

Scope: human and mouse; bulk RNA-seq, sc/snRNA-seq, spatial transcriptomics. Python first (scverse); R through
scripts with a file interface when no validated Python implementation exists. Not covered: CITE-seq, multiome,
TCR/BCR.

## Use in a project

1. Copy `CODING_GUIDELINES.md` and the domain guides the project needs, **verbatim**, into the project root, with the
   reference files they cite.
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
pip install -e ".[scrnaseq,bulk,dev]"
pytest                          # R tests run if RECIPES_RSCRIPT points to Rscript with the R packages
ruff check recipes tools && ruff format --check recipes tools
python tools/check_guides.py    # citation keys, cross-references, recipe paths, versions
python tools/rule_ids.py check  # every rule has a unique, registered identifier
python tools/rule_ids.py assign # after adding rules: assign identifiers, update RULE_IDS.tsv
```

CI (`.github/workflows/ci.yml`) runs these checks, the tests with and without R, and an installation test from the
built wheel.

## License

MIT (`LICENSE`), for the guidelines, references and code.
