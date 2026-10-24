# Changelog

All notable changes to the guideline files in this repository are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and each guideline file is
versioned with [Semantic Versioning](https://semver.org/) as defined in `CODING_GUIDELINES.md`
§0 (versioning).

## 2026-10-24 (continuous integration)

### Repository

#### Added

- CI (`.github/workflows/ci.yml`): pre-commit hooks and version check on every push and pull request; Python recipe
  tests (3.11, 3.12), installation from the built wheel and R-backed tests, each only when the change needs it
  (`tools/ci_scope.py`); every job weekly and on demand.
- `.pre-commit-config.yaml`: line endings, YAML/TOML, large files (`CODE-PROV-002`), private keys, ruff, shellcheck,
  rule identifiers, guide checks and `tools/check_versions.py`.
- `tools/check_versions.py`: a changed guide MUST bump its version and have a `### <file> [<version>]` entry here
  (`CODE-USE-007`). `tools/ci_scope.py`: selects the test suites a change needs.
- `.gitattributes`: LF line endings on every platform.
- Pytest marker `r` for R-backed tests; with `RECIPES_REQUIRE_R=1` a skipped R-backed test fails.

#### Fixed

- `tools/rule_ids.py`, `tools/check_guides.py`: read and write UTF-8 with LF line endings on every platform (on
  Windows they wrote CRLF and depended on the locale encoding).
- Recipes 0.3.1: `pyarrow` is a dependency (`polars.DataFrame.to_pandas()` needs it; two bulk tests failed in a clean
  environment).
- `README.md`, *Development*: describes the CI and pre-commit that now exist.

## 2026-10-23 (proposals from the snRNA-seq reprocessing project on DIPC Hyperion)

### CODING_GUIDELINES.md [3.2.0]

#### Added

- `CODE-DL-009`: waits inside download jobs are bounded and counted in the time request; the transfer time is sized
  from throughput measured on the nodes that run the jobs, plus the back-off of `CODE-DL-005`, and grows with the
  attempt number.
- `CODE-DL-010`: order of transfer methods per repository (ENA, SRA-only runs, Synapse, EGA, other HTTPS sources) and
  how each is verified.

#### Changed

- `CODE-DL-001`: now an explicit MUST that covers every code path yielding a file (resume, retry, multi-connection
  transfer, file already present); the pipeline compares checksums itself and records observed and expected values
  and the method; a file without a published checksum fails unless a container integrity check (SRA) is run or the
  acceptance is recorded with its reason; URL and checksum tables are matched one to one.
- `CODE-DL-006`: the connection budget counts files in parallel times connections per file, per archive and source
  address (3–4 by default, up to about 8 with measured low per-connection throughput); *Why* added.

### HPC_GUIDELINES.md [1.3.0]

#### Changed

- `HPC-FLOW-006`: CRAM version and profile are set explicitly and chosen by the project's priority: CRAM 3.0 for
  compatibility with readers without CRAM 3.1 support (htsjdk before 4.3.0, GATK before 4.7.0.0); CRAM 3.1 `normal`
  as the default once every reader is known to read 3.1; CRAM 3.1 `archive` to minimise volume, after measuring
  region extraction. A converted file is read once with the oldest reader in use; profile and writing tool are
  recorded in provenance. *Why* notes that the default output version depends on the writer (samtools/htslib 1.22,
  htsjdk 5.0.0) and cites [Bonfield 2022].

### OMICS_REFERENCES.md

- [Bonfield 2022] (CRAM 3.1), new subsection *File formats*; `HPC_GUIDELINES.md` listed among the guides it serves.

### hpc_sites/DIPC_HYPERION.md

- §8: download throughput from compute nodes (ENA with `aria2c` and `wget`, SRA `prefetch`), one size-correct file
  with a wrong MD5, and download jobs timing out while waiting inside the job.

## 2026-10-18 (proposals from the snRNA-seq project on DIPC Hyperion)

### HPC_GUIDELINES.md [1.2.0]

#### Added

- `HPC-FLOW-006`: alignments that later questions may need are archived as CRAM 3.0 against the exact mapping
  reference (record count checked, reference and removed tags in provenance) instead of being deleted; region BAMs
  are re-extracted from the CRAMs when the list of regions changes.

#### Changed

- `HPC-BATCH-005`: "before BAMs are deleted" now covers conversion to CRAM (`HPC-FLOW-006`).

### SCRNASEQ_GUIDELINES.md [3.2.0]

#### Added

- `SC-QUANT-005`: archived droplet alignments keep `CB`/`UB` and the gene and mapping tags; the raw barcode tags
  (`sS`, `sQ`, `CR`, `UR`) MAY be dropped when the FASTQs stay retrievable (pilot: 3.63 → 1.99 GB per library).

### RNASEQ_BULK_GUIDELINES.md [3.2.0]

#### Added

- `BULK-QUANT-005`: archived genome alignments keep the tags downstream steps read (`NH`, `HI`, `AS`, `nM`, `XS`,
  `RG`; `RX`/`MI` with UMIs) and MAY drop raw UMI tags (`OX`, `QX`); the transcriptome-projected BAM is not archived.

### SAMPLE_IDENTITY_GUIDELINES.md [2.2.2]

#### Fixed

- `ID-GENO-003`: IBS0 and IBS2 are numbers of sites (somalier `ibs0`, `ibs2`, out of `n`), not proportions; divide by
  `n` to compare pairs.

### SAMPLE_IDENTITY_GUIDELINES.md [2.2.1]

#### Changed

- `ID-GENO-001` *Trade-off*: fingerprints can also be extracted later from archived CRAMs (`HPC-FLOW-006`).

## 2026-10-17 (proposals from a bulk RNA-seq project on DIPC Hyperion)

### CODING_GUIDELINES.md [3.1.0]

#### Added

- §13, new subsection *Downloads from public archives* (`CODE-DL-001`–`008`): checksum verification; whole-unit
  retries, never client retries into a pipe; record-by-record validation of streamed reads; redirects treated as
  failures for archives that do not redirect; back-off over minutes; few concurrent connections per archive; a
  validated, tracked run table; testing downloads from the nodes that run them.
- `CODE-LANG-004`: `set -E` and an `ERR` trap reporting line, command and exit code; no discarded STDERR on commands
  whose failure aborts the script.
- `CODE-PROV-002`: *Deviate when* for small derived tables tracked in git, with their regeneration script.

#### Changed

- `CODE-USE-001`, §2 (documentation set and `AGENTS.md` template): copied guidelines live in `guidelines/` at the
  project root (harmonizes the unversioned edit of commit `3585824` with the template and `README.md`);
  `AGENTS.md`, `PROJECT_CONTEXT.md` and `README.md` stay in the root.
- `CODE-USE-008`: prefix list includes `HPC`.

### HPC_GUIDELINES.md [1.1.0]

#### Added

- §2, subsection *Datasets larger than the work area* (`HPC-BATCH-001`–`006`): batches sized to a space budget,
  footprint calibrated on a pilot, stage markers, validation before deletion (MUST), extraction of what needs an
  intermediate before deleting it, `--dependency=singleton`.
- `HPC-WF-005` parameters in a tracked `-params-file`; `HPC-WF-006` prebuilt images with Apptainer
  (`-profile singularity` for nf-core).
- `HPC-FLOW-001`: *Deviate when* for large read-only references read from `ARCHIVE`, after measurement.

#### Changed

- `HPC-ENV-003`: `NXF_HOME` among the caches; caches written by jobs must be writable from compute nodes.
- `HPC-JOB-001`: `set -Eeuo pipefail` and `ERR` trap. `recipes/hpc/job.slurm` updated; `recipes/hpc/sites/hyperion.env`
  falls back to `WORK` when `/lscratch` is absent and sets `NXF_HOME`.
- `hpc_sites/DIPC_HYPERION.md`: §8, observations on the cluster (Nextflow version and Java, `NXF_HOME`, nf-core
  profile, `/data` read times, quota commands, network address of compute nodes, `/lscratch`, group partitions,
  `/dipc` → `/home`).

### SEX_INFERENCE_GUIDELINES.md [2.1.0]

- `SEX-MARK-001`: *TSIX* no longer an optional marker (barely detected in adult tissue; antisense to *XIST*).
- `SEX-MARK-002`, `SEX-MARK-003`: *TXLNGY* among brain MSY genes; *TXLNG/TXLNGY* among homolog pairs.
- `SEX-COUNT-003`: expected informative reads with gene-body counting (hundreds to thousands); screening reads
  MAY be streamed, with validation and the caveat that the first reads are not a random sample.
- `SEX-COUNT-005` (new): marker counting MUST use the library's strandedness.
- `SEX-LOY-001`: comparison across a donor's samples.

### RNASEQ_BULK_GUIDELINES.md [3.1.0]

- `BULK-RQC-004` (new): per-sample strandedness from a screening alignment (STAR `GeneCounts`), with the caveat that
  its fraction is lower than Salmon's and tissue-dependent, so Salmon cut-offs are not transferable.

### SAMPLE_IDENTITY_GUIDELINES.md [2.2.0]

- `ID-DUP-004` (new): re-deposited libraries, processed once, with a few pairs kept as positive controls.
- `ID-DUP-005` (new): metadata discordance within duplicate pairs is flagged and resolved in the decisions log.

### Not adopted

- Reference ranges of sex markers from a single cohort (one depth, tissue and protocol) as guide content: thresholds
  stay data-driven (`SEX-THR-001`); the values can serve as an example in the future `recipes/sex/` documentation.
- "MUST NOT follow redirects" for every download: adopted only for archives that serve files without redirection,
  because several services (signed cloud-storage URLs, data repositories) redirect by design.
- `SAMPLE_IDENTITY`: re-deposits as duplicates "by default" without confirmation; adopted as nomination plus
  positive controls.

## 2026-10-16

### HPC_GUIDELINES.md [1.0.0]

#### Added

- New guide for shared clusters with Slurm (40 rules, prefix `HPC`): storage roles (`CODE`, `WORK`, `ARCHIVE`,
  `NODE`) with the same path `<root>/$USER/<repo>` on every filesystem, role variables from a tracked site file,
  repository links; staging and promotion with manifests; no-backup assumption; Apptainer images from tracked
  definitions; job conventions (`--export=NONE`, `--no-requeue`, measured resources, arrays with throttle,
  idempotent tasks); workflow managers; small files and node-local scratch; provenance; rules for agents on
  login nodes and for mass submission or deletion.
- Site profile `hpc_sites/DIPC_HYPERION.md` (checked against the DIPC documentation on 2026-10-16): filesystems and
  quotas, general partitions and QoS (further ones can be granted per group), GPU constraints, modules.
- Recipes `recipes/hpc/` (`setup_storage.sh`, `stage.sh`, `promote.sh`, `load_site.sh`, `job.slurm`,
  `sites/local.env`, `sites/hyperion.env`) with tests in `recipes/tests/test_hpc_scripts.py`.

### CODING_GUIDELINES.md [3.0.1]

#### Changed

- §26 (related guidelines) lists `HPC_GUIDELINES.md`.

### Tools

- `tools/rule_ids.py`: prefix `HPC` and its topics. `tools/check_guides.py`: also checks cross-references and
  recipe paths in `hpc_sites/*.md`, and recipe paths ending in `.sh`, `.slurm`, `.env`.

### SCRNASEQ_GUIDELINES.md [3.1.0]

#### Changed

- `SC-QUANT-001`: the main matrix for nuclei (and for cells when introns are counted) is `GeneFull_Ex50pAS`
  instead of `GeneFull` (`--soloFeatures Gene GeneFull_Ex50pAS Velocyto`, STAR ≥ 2.7.10a). `GeneFull` treats any
  overlap with another gene's locus as ambiguity, so genes nested in introns of other genes on the same strand lose
  their counts; `GeneFull_Ex50pAS` assigns reads that are more than 50% exonic in one gene to that gene. `GeneFull`
  remains an optional comparison layer. MINOR: the rule keeps its scope (full-gene counting for nuclei) and
  refines the counting option.
## 2026-10-17 (labelling of single-cell pools)

Prompted by reprocessing public pooled snRNA-seq (Masrori et al. 2025, Nat Neurosci,
doi:10.1038/s41593-025-02075-1), where donors were assigned to freemuxlet clusters by matching to
per-patient bulk RNA-seq, and per-cell labels were only available in a SCope loom.

### SAMPLE_IDENTITY_GUIDELINES.md [2.1.0]

#### Added

- §6: order of evidence for labelling clusters (author labels, genotype references of the same donors from
  other data, sex in XX+XY pools, unnamed individuals); per-cluster fingerprints MUST be stored before BAM
  deletion; acceptance conditions for cluster-to-donor matches; unlabelled clusters of mixed-condition pools
  MUST NOT enter condition comparisons; genotype-free demultiplexers interchangeable on common transcribed variants.

### META_ANALYSIS_GUIDELINES.md [2.1.0]

#### Added

- §10: when reprocessing public pools, re-run demultiplexing, compare with author assignments, look for labels in
  processed objects and browser exports, and record dependence on undeposited reference data.

## 2026-10-15 (external review integrated; rule identifiers)

Changes after an external review of commit `97ae0ce`. Each code finding has a regression test in
`recipes/tests/test_review_regressions.py`.

### All guidelines [MAJOR]

- **Stable rule identifiers** (`PREFIX-TOPIC-NNN`, e.g. `SC-QC-004`) on every rule, with the registry
  `RULE_IDS.tsv` (563 active rules) and `tools/rule_ids.py` (`assign`, `check`). Overrides now refer to rules by
  identifier (`CODE-USE-008`, `CODE-PROJ-008`); existing overrides that point to section numbers should be rewritten
  with identifiers. Versioning redefined (`CODE-USE-007`): renumbering sections is no longer MAJOR.
- Key words: MUST reserved for checkable validity conditions; a rule without a key word has the force of SHOULD
  (header of every guide). MUST rules whose *Deviate when* described cases outside their scope now say *Scope*
  (`SEX-WHY-001`, `STAT-DES-004`, `STAT-MOD-007`); `CODE-LANG-004` (shellcheck) is SHOULD; `ID-POOL-001` covers
  hashtag pools.
- `LICENSE` (MIT) for the whole repository.

### SEX_INFERENCE_GUIDELINES.md [2.0.0]

#### Changed (MAJOR: a rule reversed)

- §6 "within a donor, all samples MUST agree; discordance indicates a swap" replaced by: discordance MUST trigger
  review and MUST NOT by itself trigger reassignment; other explanations (mosaic loss of Y, pathology, quality,
  transplantation) listed. Expression calls are screening results; a complement is confirmed only with DNA.
- §2 *XIST* limits (XXY, some tumours of male patients); the homolog ratio cancels shared scale factors, not
  differential effects of tissue or quality.
- §3, §5 depth, minimum reads and ratio thresholds are starting values that MUST be calibrated.
- §4 single routing rule (SHOULD when sex-chromosome or PAR genes are analysed or sex is a variable; uncertain calls to
  the XY reference), shared with RNASEQ_BULK §2.

### DECONVOLUTION_GUIDELINES.md [2.0.0]

#### Changed (MAJOR: a rule reversed)

- §3 "linear scale, not log" replaced by "use the input scale each method documents" (hspe takes log-scale input);
  normalization only as each method requires (Bisque and BayesPrism take raw counts).
- §1 the estimand depends on the method (RNA fraction, cell proportion, fraction of reads); method contract table.
- §4 agreement between methods is necessary, not sufficient (shared biases); §5 validation adds bias per cell type,
  abundance range and an incomplete reference; §8 aligned.

### STATISTICS_GUIDELINES.md [2.0.0]

- §10 exchangeability is a property of the design; permuting within strata does not preserve associations with
  other covariates; two-sided permutation p-values with the SciPy convention; Monte Carlo uncertainty of permutation
  frequencies.
- §7 Bonferroni is valid under any dependence (union bound); FDR is an expectation; p-value histograms: a peak near 1
  can come from conservative or discrete tests.
- §3 Student vs Welch as models; Mann-Whitney as the score test of the proportional-odds model; rank tests stated as
  $P(X>Y)$ hypotheses; ANCOVA limited to randomized before/after comparisons; Wald/LRT/score disagreement needs a
  method for the cause (Firth, parametric bootstrap, exact), not a switch to the LRT; beta regression support (0, 1).
- §1 low power and false positives qualified by significance level and proportion of true effects.

### SCRNASEQ_GUIDELINES.md [3.0.0]

- §1.1 table of the layer each step reads (`counts` for QC and doublets; ambient-corrected counts downstream).
- §2 donor minimums are planning minimums; aggregate to the donor only when it contributes to one level of every
  compared factor; §10 aggregation unit defined by the design.
- §4.2 scaled MAD defined (5 scaled = 7.41 unscaled); non-finite metrics and MAD = 0 handled explicitly; panels.
- §6 PCA recomputed instead of truncating the proposal; "no reproducible partition" as an outcome; check runs on the
  embedding used for clustering; removed "extra PCs add little noise"; stability is not biological validity.
- §7 the Harmony benchmark described with its scope.

### RNASEQ_BULK_GUIDELINES.md [3.0.0]

- §8 permutation check renamed global-null diagnostic, with exchangeability and Monte Carlo uncertainty; it does not
  establish FDR control under alternatives.
- §7 residual PCs described as exploratory, not qSVA/RUV/SVA; §2 routing deferred to SEX_INFERENCE §4; §1 replicate
  minimum is a planning default.

### META_ANALYSIS_GUIDELINES.md [2.0.0]

- §2 independence between studies (shared donors, reused controls, re-deposited data) checked before combining.
- §3 unshrunken log fold changes and their SE; harmonized estimands.
- §5 prediction interval only for $k \ge 3$, imprecise with few studies.
- §9 harmonic mean p-value: asymptotically exact procedure, approximate FWER control.

### SPATIAL_GUIDELINES.md [2.0.0]

- §1 standard Visium vs Visium HD geometry; enrichment universe = eligible panel genes. §4 scRNA-seq defaults that
  do not apply to panels (mitochondrial and nuclear filters, 2,000 HVGs).

### SAMPLE_IDENTITY_GUIDELINES.md [2.0.0]

- §7 sex mismatch with matching genotype: review; correct the annotation only with independent evidence.

### CODING_GUIDELINES.md [3.0.0]

- Header: MUST reserved for checkable validity conditions; SHOULD for defaults; scope of *Deviate when*.
- §0 precedence between domain guidelines; adopted guides and versions listed in *Overrides*.
- §7 aggregates may combine disjoint groups of decomposable statistics; §14 installed library code locates its data
  with `importlib.resources`; §15 bit-identical results not required, random draws identical and floating-point
  results within a stated tolerance; §21 and §26 vendoring of recipes allowed with dependencies, tests, license and
  commit; §26 recipes MAY exist, coverage in `recipes/README.md`; FIGURE_GUIDELINES marked as planned.

### STATISTICS_REFERENCES.md

- Explicit bold keys (`[Author Year]`) for all 88 references, as in OMICS_REFERENCES.md.

### recipes/ [0.2.0]

#### Fixed

- `pseudobulk`, `abundance`: covariates keep their dtype (numeric covariates are slopes); missing values rejected;
  report lists absent and below-threshold combinations; empty results keep the schema; abundance accepts a sample
  manifest (samples without cells flagged, not analysed).
- `stats.permutation`: non-finite statistics rejected; SciPy two-sided convention; Clopper-Pearson bound.
- `doublets`: co-expression measured per cell (a cluster of mixed singlets is no longer removed); enrichment against
  capture-specific call rates; DoubletDetection ambiguous calls kept as `ambiguous`; detector parameters reported.
- `qc`: MAD = 0 and non-finite values reported (`not_evaluated_mad_zero`, `qc_invalid`, `qc_fail`); optional
  absolute mitochondrial limit; panels without mitochondrial genes.
- `normalization`: size factors validated (finite, positive, one per cell); median factors no longer re-centred to
  the mean (now equal to `normalize_total(target_sum=None)`); R script located as package data; unique run
  directories.
- `pc_robustness`: proposal not truncated; status `recompute_pca`.
- `deseq2`: filter ignores unused categories; exact coefficient names; designs with interactions rejected for
  shrinkage; count validation (`allow_non_integer` only for lengthScaledTPM).
- `hidden_factors`: rank checks with tolerance; aligned factors; k range validated; documented as exploratory.
- `fdr_calibration`: global-null diagnostic with exchangeability caveats, adjusted p-values validated, Monte Carlo
  bound, distinct permutations counted.
- `marker_annotation`: contradictory signs rejected; ties and unscored cells labelled; index alignment checked.
- `gene_ids`: rule names validated; deterministic choice among duplicated IDs; warning for non-Ensembl IDs.
- `_clustering`: Jaccard from the contingency table; stability reports `n_present` and `not_assessed`.
- R scripts: `--overwrite` required to replace outputs, atomic writes, provenance file.

#### Added

- Packaging: setuptools build with R scripts as package data; extras `scrnaseq`, `bulk`, `dev`; line length 99 as in
  CODING §17.
- Test of pydeseq2 against R DESeq2 (`recipes/tests/r/deseq2_reference.R`); `recipes/README.md` coverage table.
- `tools/check_guides.py` (citation keys, cross-references, recipe paths, versions); CI workflow with Python, R and
  wheel-installation jobs; repository `README.md`.

## 2026-10-14 (third release of the day)

### SCRNASEQ_GUIDELINES.md [2.1.0]

#### Added

- §3 STARsolo as the default mapper and gene-level quantifier, with parameters to reproduce Cell Ranger (≥ 4)
  cell calling (`EmptyDrops_CR`), `Gene` for cells and `GeneFull` for nuclei, and `Velocyto` output for the
  spliced/unspliced/ambiguous layers. Isoform-level work: Salmon preferred over alevin, and only with full-length or
  long-read data.
- §6 detailed procedure to test the robustness of the number of PCs ($p \pm 5$, ARI and per-cluster Jaccard), with
  the reference implementation `recipes/scrnaseq/pc_robustness.py`.
- §8 number of subsampling iterations by dataset size (10 / 5 / 1) and reference implementation
  `recipes/scrnaseq/clustering_stability.py`.

#### Changed

- §1.3 the `var` index is the stable Ensembl ID without version; the versioned ID is kept in
  `var["gene_id_version"]`.
- §4.1 CellBender runs on the raw STARsolo matrix and its calls are compared with `EmptyDrops_CR`; the nuclear
  fraction is computed from the STARsolo `Velocyto` layers.
- §8 resolution grid 0.1–2.0; stability reading made explicit for the interval [0.5, 0.75).
- §5 pooled size factors in Python with `scranPY` (sfortma2/scranPY, pinned commit) instead of the Bioconductor
  `scranpy`, which has no pooled size factors; R scran kept as cross-check.
- §4.3 citation fixed: DoubletDetection is [Gayoso 2020]; [Wolock 2019] is Scrublet.
- Links to the recipes of §1.3, §4, §5, §9, §10 and §11.

### RNASEQ_BULK_GUIDELINES.md [2.1.0]

#### Added

- Links to the recipes of §4 (count import), §5 (filter, length bias), §7 (hidden factors and $k$ sensitivity) and §8
  (pydeseq2 wrapper, FDR calibration by permutation).
- §8 note: `pydeseq2` 0.5 has no likelihood-ratio test; use DESeq2 in R.

### CODING_GUIDELINES.md [2.3.0]

#### Added

- §26 *Recipes*: use the linked reference implementation instead of re-implementing; adapt through parameters;
  changes of logic are recorded as overrides.

### All other guidelines [PATCH]

STATISTICS_GUIDELINES.md [1.1.1],
SPATIAL_GUIDELINES.md [1.1.1], META_ANALYSIS_GUIDELINES.md [1.1.1], SEX_INFERENCE_GUIDELINES.md [1.0.1],
SAMPLE_IDENTITY_GUIDELINES.md [1.0.1], DECONVOLUTION_GUIDELINES.md [1.0.1].

#### Changed

- Acronyms expanded at their first occurrence in every guideline, including those listed above (e.g. RUV, RIN, qSVA, LOY, VST, TMM). No rule changed.

### recipes/ [0.1.0]

#### Added

- Reference implementations with tests on synthetic data (`pytest`, 42 tests; R-dependent tests are skipped
  without R):
  - `recipes/common/gene_ids.py`;
  - `recipes/scrnaseq/`: `pc_robustness.py`, `clustering_stability.py`, `qc.py`, `doublets.py`,
    `normalization.py` (+ `scran_size_factors.R`), `marker_annotation.py`, `pseudobulk.py`, `abundance.py`;
  - `recipes/bulk/`: `counts_io.py`, `tximport_lengthscaled.R`, `deseq2.py`, `length_bias.py`,
    `hidden_factors.py`, `fdr_calibration.py`;
  - `recipes/stats/permutation.py`.
- `pyproject.toml` with dependencies, pytest and ruff configuration.
- STATISTICS_GUIDELINES.md §10 links `recipes/stats/permutation.py` (PATCH).
- OMICS_REFERENCES.md: two software records, [Fortmann 2023] (scranPY) and [Gayoso 2020] (DoubletDetection).

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
