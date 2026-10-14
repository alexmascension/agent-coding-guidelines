# Changelog

All notable changes to the guideline files in this repository are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and each guideline file is
versioned with [Semantic Versioning](https://semver.org/) as defined in `CODING_GUIDELINES.md`
§0 (versioning).

## 2026-10-14 (second release of the day)

### SEX_INFERENCE_GUIDELINES.md, SAMPLE_IDENTITY_GUIDELINES.md, DECONVOLUTION_GUIDELINES.md [1.0.0]

#### Added

- Three new guidelines, based on a literature review summarized in `REVIEW_DECISIONS.md`.

### SCRNASEQ_GUIDELINES.md [2.0.0]

#### Added

- Rationale blocks (*Why*, *Trade-off*, *Deviate when*) for every non-obvious rule.
- §1.1 named layers instead of `adata.raw`; §1.3 handling of Ensembl IDs and symbols in marker dictionaries.
- §2 single-donor designs and several samples per donor.
- §3 gene-level quantification with spliced/unspliced/ambiguous layers.
- §4.1 nuclear fraction from USA layers; §4.3 two doublet detectors (Scrublet and DoubletDetection) with consensus
  and union columns, flagged first and removed after a first clustering.
- §6 procedure to choose the number of PCs; §8 stability-based choice of clustering resolution.
- §9 annotation with marker dictionaries (CellAssign, decoupler); §11 differential abundance in detail.
- §12 cell-cell communication as its own section; §13 trajectories, CellRank 2 kernels, PAGA and velocity caveats.

#### Changed

- §4.1 CellBender is the default for cell calling and ambient RNA (SoupX as alternative).
- §4.2 mitochondrial genes identified by chromosome, not by prefix.
- §5 scran pooled size factors (or median-based `normalize_total`) instead of a fixed target sum.
- §7 Harmony as default, scVI/scANVI as second opinion.
- Sections renumbered (MAJOR).

### RNASEQ_BULK_GUIDELINES.md [2.0.0]

#### Added

- Rationale blocks for every non-obvious rule.
- §4 quantification in detail: `star_salmon`, Salmon `--seqBias --gcBias --numGibbsSamples 20`, tximport modes
  (`lengthScaledTPM` for `pydeseq2`).
- §5 sample-specific gene length bias; §6 VST as default transformation.
- §7 known and hidden factors: qSVA, RUV, SVA and PCA with criteria and number of factors.
- §9 transcript-level analyses with inferential replicates; §10 enrichment biases, camera and decoupler footprints.

#### Removed

- Sex check (moved to `SEX_INFERENCE_GUIDELINES.md`) and deconvolution (moved to `DECONVOLUTION_GUIDELINES.md`).

### CODING_GUIDELINES.md [2.2.0], STATISTICS_GUIDELINES.md [1.1.0], SPATIAL_GUIDELINES.md [1.1.0], META_ANALYSIS_GUIDELINES.md [1.1.0]

#### Added

- Expanded rationale for rules whose justification was implicit or terse; no rule changed and no section renumbered.
- CODING §9: rule documents give each non-obvious rule its rationale; §26 lists the new guidelines.
- SPATIAL §7 links to `DECONVOLUTION_GUIDELINES.md`.

### OMICS_REFERENCES.md

#### Added

- 112 new Crossref-verified references (249 in total).

#### Changed

- Keys renamed to resolve collisions: `Love 2016` → `Love 2016a`, `Zhao 2021` → `Zhao 2021a`.

## RNASEQ_BULK_GUIDELINES.md, SCRNASEQ_GUIDELINES.md, SPATIAL_GUIDELINES.md, META_ANALYSIS_GUIDELINES.md

### [1.0.0] - 2026-10-14

#### Added

- First versions of the four omics guidelines (bulk RNA-seq; single-cell and single-nucleus RNA-seq; spatial
  transcriptomics for Visium, Xenium and MERFISH; general and omics meta-analysis).
- Companion `OMICS_REFERENCES.md` (137 Crossref-verified references with notes from their abstracts) and
  `OMICS_REFERENCES.bib`.

## STATISTICS_GUIDELINES.md (patch)

### [1.0.1] - 2026-10-14

#### Changed

- References to `RNASEQ_GUIDELINES.md` replaced by the new guideline names (§0, §3.2, §13).

## CODING_GUIDELINES.md (patch)

### [2.1.1] - 2026-10-14

#### Changed

- §2 (documentation set) and §26 (related guidelines) list the new guideline and reference files.

## STATISTICS_GUIDELINES.md

### [1.0.0] - 2026-10-13

#### Added

- First version: design, exploration, model choice, hierarchical data, batch effects, circular
  analysis and leakage, multiple testing, effect sizes, enrichment, resampling, optional
  Bayesian analysis, reporting checklist, Python implementations and glossary.
- Companion `STATISTICS_REFERENCES.md` (88 Crossref-verified references and practitioner
  sources) and `STATISTICS_REFERENCES.bib`.

## CODING_GUIDELINES.md

### [2.1.0] - 2026-10-13

#### Added

- `AGENTS.md` as a required entry-point file, with template; `CLAUDE.md` bridge for agents
  that do not read `AGENTS.md` (§2, documentation set).
- Optional `LICENSE`, `CITATION.cff` and `CHANGELOG.md` for published repositories (§2).
- *Overrides* section of `PROJECT_CONTEXT.md` records the adopted guideline versions (§2).
- Semantic versioning rules for guideline files (§0, versioning).
- Bash rules: 100-line limit before rewriting in Python, `$(…)`, `[[ … ]]`, `main` function,
  mandatory `shellcheck` (§3, languages and tooling).
- R style: tidyverse style guide, `styler`, `lintr` (§3).
- Function-length signal (~40 lines) and comprehension complexity limit (§4, principles).
- Short names allowed for notation from a cited reference paper (§6, naming).
- Self-contained comments and documents: state the reason, not the development history; no
  commented-out code; `TODO` format; minimal and current documentation (§9, documentation).
- Markdown rules: single H1, unique ATX headings, fenced code with language, `$...$` and
  `$$...$$` for equations, descriptive cross-references (§9).
- No `assert` for input validation, no mutable default arguments, `with` for resources,
  avoid power features, precise error messages with `{value=}` (§10, robustness).
- `is None` checks instead of truthiness (§11, numerical code).
- Runnable-from-fresh-clone requirement; scripts use `main()` with `__name__` guard (§14,
  configuration and paths).
- Avoid mutable global state (§20, side effects).
- Commit message format (50/72), one logical change per commit, user approval before
  committing and pushing (§23, version control).
- `shellcheck` in pre-commit and validation (§17, §24).

#### Changed

- f-strings are used for all interpolation, including `logging` calls; ruff `G004` is
  disabled (§17, style).
- Project root resolved from `pyproject.toml`, with `PROJECT_ROOT` as optional override;
  external data is linked into `data/` (§14).
- "Preserve existing behaviour" moved from §4 (principles) to §5 (changes to existing code).

#### Removed

- Untracked local configuration file for machine-specific settings (§14, §22, §23).

### [2.0.0] - 2026-10-13

#### Added

- First project-agnostic version: project-specific content moved to `PROJECT_CONTEXT.md`, with
  a required nine-section structure.
