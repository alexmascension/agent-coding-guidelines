# Changelog

All notable changes to the guideline files in this repository are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and each guideline file is
versioned with [Semantic Versioning](https://semver.org/) as defined in `CODING_GUIDELINES.md`
§0 (versioning).


## 2026-11-05 (external audit of c981ac89; recipes 0.8.0; CODING 4.1.0, AGENTIC 1.2.0, MANUSCRIPT 1.2.0, STATISTICS 2.3.0, SCRNASEQ 5.3.0, RNASEQ_BULK 3.4.0, SEX_INFERENCE 2.2.0, SAMPLE_IDENTITY 3.1.0, DECONVOLUTION 2.0.1, SPATIAL 3.1.0, META_ANALYSIS 2.3.0, HPC 1.5.0)

An external audit of commit `c981ac89` raised 99 observations on the guides and 23 on the recipes and tooling. Each
item was checked against the rule text or the code before it was adopted (`AGT-REV-004`, `AGT-REV-006`); the items
not adopted are listed at the end. The reproduced recipe failures have regression tests that fail on `c981ac89`.
Reference notes state whether they rest on the full text or on the abstract.

### CODING_GUIDELINES.md [4.1.0]

#### Changed

- `CODE-PROV-001`: the copy of record of a raw input is read-only and MUST NOT be modified, overwritten or deleted; a
  replica (a copy staged on `WORK` or `NODE`, a cache) MAY be deleted once the copy of record has been verified to
  exist and to match its recorded checksum (`HPC-FLOW-001`, `HPC_GUIDELINES.md` §2). *Why* added.
- `CODE-DL-003`: a streamed head is a partial extract: validated record by record, recorded as a partial extract in the
  provenance file and never recorded as a checksum-verified file (`CODE-DL-001`).
- `CODE-PIPE-001`: *Why* no longer says that the outputs of a hand-built equivalent are not comparable; they match the
  pipeline's only if every tool, version, parameter and reference matches, which the project would have to show and a
  pinned release with tracked parameters fixes.
- `CODE-LANG-007`: the `conda-lock` lock file is the reproducible record; an export with `--no-builds` records package
  versions, not builds, and still lists platform-specific dependencies.

### AGENTIC_GUIDELINES.md [1.2.0]

#### Changed

- `AGT-DEL-003`: do not assume that a subagent sees the project's instructions; unless the harness documents what it
  inherits, the brief names the rules that apply.
- `AGT-USER-001`: scoped to questions about the work; a request phrased as a question ("can you fix X?") is a request
  from the user; a request relayed from a co-author or reviewer stays the user's to decide (`MS-RESP-001`).
- `AGT-USER-005`: the fix for the task at hand goes ahead; the generalised rule is confirmed before it is recorded
  beyond the task.
- `AGT-REV-004`: an unclear item blocks itself and the items that depend on it or touch the same code, data or text;
  clear independent items may proceed.
- `AGT-EVID-006`: a figure, or a table or equation the conversion may have garbled, is checked on the page image when a
  claim depends on it.
- `AGT-CLAIM-003`: the source of a regenerated number is the stored result file of the run (`CODE-PROV-003`), not a
  re-run of the analysis.
- `AGT-INTENT-005`: playback length scales with the part (about 200 to 300 words for a complex part, fewer for a simple
  one).
- `AGT-LOOP-006`: an attempt is a fix based on a hypothesis about the cause; automatic retries of a transient failure
  (`CODE-DL-002`, `CODE-DL-005`) are not attempts.
- `AGT-LOOP-002`: three or more framings of unchanged results are a signal to examine the evidence, the argument and
  the audience, not a diagnosis that the data are the bottleneck (consistent with `MS-RESP-004`).
- `AGT-INSTR-002`, `AGT-INSTR-003`: *Note* on the evidence (superpowers' wording tests: one model, one test harness,
  five samples per variant); the table is a set of heuristics to check against the failures observed.

### MANUSCRIPT_GUIDELINES.md [1.2.0]

#### Changed

- `MS-TRACE-001`: *Deviate when* for hand-drawn schematics and acquired images: keep the source file and record in the
  sidecar the source and, for an image, the acquisition and every processing step.
- `MS-RESP-004`: the reframing count is a prompt to examine the evidence, the argument and the audience, not a
  diagnosis that the bottleneck is the data or the question.

### STATISTICS_GUIDELINES.md [2.3.0]

#### Changed

- `STAT-MOD-006`: statistic, hypothesis and assumptions separated. Mann-Whitney $U$ estimates
  $\theta = P(X > Y) + \tfrac{1}{2}P(X = Y)$ but is a valid test of $\theta = 1/2$ only when the groups have the same
  distribution under the null (e.g. shift model); Brunner-Munzel (`scipy.stats.brunnermunzel`) when spread or shape may
  differ; Kruskal-Wallis and Friedman stated as tests of identical distributions; MUST NOT on medians kept. *Why* and
  *Note* (Brunner-Munzel is an approximation; group sizes) added [Fay 2010; Brunner 2000; Divine 2018]. §13 table:
  `scipy.stats.brunnermunzel`; the equal-medians wording of the `kruskal` docstring holds only under a shift model;
  `brunnermunzel` returns `nan` (`0.0` with `distribution="normal"`) for non-overlapping groups.
- `STAT-RES-003`: $(b+1)/(m+1)$ scoped to $m$ random permutations; with complete enumeration the exact p-value is $b/M$,
  the observed arrangement counted among the $M$. Both two-sided conventions (doubling the smaller one-sided p,
  comparing $|T|$) give valid tests under exchangeability; $|T|$ matches a two-sided alternative only when the null
  of $T$ is centred at 0 and roughly symmetric, so the recipe doubles; state the convention used. Replaces "comparing
  absolute values is valid only for statistics centred at 0 under the null". §13 table: `permutation_test` enumerates
  when `n_resamples` reaches the number of distinct permutations.
- `STAT-RES-004`: smallest attainable p-value $1/(m+1)$ one-sided and $2/(m+1)$ for the doubled two-sided p ($1/M$ and
  $2/M$ with complete enumeration).
- `STAT-RES-005`: units MUST be resampled; a two-stage bootstrap MAY additionally resample measurements within each
  resampled unit; resampling only the measurements is what is prohibited.
- `STAT-RES-008`: "so the estimates are upper bounds" replaced: recall and precision against the held-out reference
  measure agreement with an imperfect reference, not true power or precision (recall tends to be optimistic,
  precision pessimistic, neither a guaranteed bound).
- `STAT-DES-002`: observational studies: the unit is the unit of sampling and inference (the smallest entity sampled
  independently, e.g. a donor); repeated samples of one donor are dependent. Glossary entry *Experimental unit*
  aligned.
- `STAT-EXPL-005`: the prohibition covers choosing a transformation or deleting data after looking at the results (to
  pass a normality test or reach significance); a transformation chosen in advance for the estimand and data type is
  part of the model; *Why* added (`STAT-DES-010`).
- `STAT-MT-002`: *Trade-off* gives the FDR level as an average; one list can contain more.
- `STAT-MT-009`: name the column after the procedure; a Storey q-value and the BH-adjusted p of the same p-values are
  different numbers [Storey 2003].
- `STAT-EFF-005`: values near the point estimate, not the centre (intervals can be asymmetric).
- `STAT-USE-002`: installation follows `CODE-USE-001` (this file into `guidelines/`, the reference files into
  `guidelines/references/`) instead of "into the project root".
- `STAT-USE-003`, `STAT-REP-007`: `FIGURE_GUIDELINES.md` marked as planned.

### SCRNASEQ_GUIDELINES.md [5.3.0]

#### Changed

- `SC-DATA-004`: pseudobulk row extended to depth matching and empirical power (`counts_ambient_corrected` if it
  exists, else `counts`); recipes resolve the layer with one function, `layer=None` in normalization and pseudobulk
  recipes meaning the table's row; a non-integer corrected layer MUST NOT be rounded silently for count-based DE
  (explicit rounding recorded in `PROJECT_CONTEXT.md`, e.g. SoupX `adjustCounts(roundToInt = TRUE)` with its seed, or
  `layer="counts"`); *Recipe* line (`recipes/scrnaseq/layers.py`, `input_layer`).
- `SC-DE-001`, `SC-DES-004`: the prohibition is on tests that treat cells as independent replicates; mixed models
  with the donor (and the sample) as random effects remain allowed.
- `SC-DE-010`: the DE list is incomplete, enriched for large effects and highly expressed genes, and contains false
  discoveries at the controlled rate (was "a lower bound of the true changes"); genes for validation from the most
  significant, claims about single genes need family-wise control.
- `SC-DE-011`: `median` MUST NOT be the input of count-based DE tools (was "used for DE"); "No benchmark supports it"
  removed.
- `SC-DE-014`: zero-inflated mixtures (read counts without UMIs) distinguished from hurdle models, which separate
  detection from level (`SC-DE-012`) and, fitted to cells, need the donor as a random effect; Poisson or negative
  binomial remains the default for UMI counts; cites [Finak 2015] (MAST).
- `SC-QC-003`: MADs on logit(clip(f, ε, 1 − ε)), after checking that f lies in [0, 1] (was logit(f + ε)); ε = 10⁻⁴ in
  the recipe.
- `SC-QC-016`, §16: the first pass removes no cells (not only doublets): samplesheet values that disable each cell
  filter at commit `0048567`, the filters that cannot be disabled there, and a per-capture check that every barcode
  leaving cell calling reaches the first-pass object.
- `SC-QC-014`: *Deviate when* requires roughly equal donor proportions; same-donor doublet share under random pairing is
  the sum of squared donor fractions.
- `SC-PC-002`: when the full singular-value spectrum cannot be computed, the elbow alone, recording that Gavish–Donoho
  was not evaluated; `SC-PC-008` records the criteria evaluated.
- `SC-PC-004`: embeddings not ordered by variance (scVI latent space) vary the number of latent dimensions by
  refitting, not by truncation. `SC-PC-006`: a clustering that changes with one or two more PCs may depend on noise or
  on a weak signal (e.g. a rare population); check which.
- `SC-QUANT-002`: full-length isoforms are not quantified from 3' tag data; poly(A)-site usage (APA) can be studied
  with dedicated methods [Wu 2020; Bi 2026].
- `SC-QC-001`: *Why* and *Trade-off* cite a seven-tool decontamination benchmark with no single best tool [Nie 2026].
- `SC-DES-007` (*Trade-off*), `SC-DE-008`, `SC-DE-009`, `SC-DE-010` (*Note*): scope of [Crampton 2026†] (sex
  differences, three brain cell types, one cohort of 1,494 donors).
- `SC-DES-007` (*Trade-off*): recall and precision of the benchmark are read as `STATISTICS_GUIDELINES.md`
  `STAT-RES-008` describes.

### RNASEQ_BULK_GUIDELINES.md [3.4.0]

#### Changed

- `BULK-DES-001`: samples that share a donor (tissues, time points, aliquots, pools containing the same donor) are
  dependent; the donor is modelled (`STAT-DES-002`) with the designs of `SC-DES-005`.
- `BULK-HF-002`: *Note* separating PCA of expression from residual principal components with the design protected
  (`recipes/bulk/hidden_factors.py`), which are exploratory and not a substitute for qSVA, SVA or RUV.
- `BULK-HF-004`: MUST NOT kept; batch is modelled as a covariate on the original counts. *Why* from [Nygaard 2015†];
  [Zhang 2020†] now cited for what it proposes (adjusted counts for DE software, with its simulation results) instead
  of as support for the prohibition.

### SEX_INFERENCE_GUIDELINES.md [2.2.0]

#### Changed

- `SEX-MARK-003`: pairs with X + Y = 0 are left out of the median; the project states a minimum number of informative
  pairs, below which the sample is `LOW_SIGNAL`; "most sensitive signal" and "without a cohort-level threshold"
  replaced by statements consistent with `SEX-THR-002`, `SEX-THR-004` and `SEX-THR-006`.
- `SEX-COUNT-001`: introns counted for nuclei including 3' single-nucleus libraries; exonic counts for whole-cell
  poly(A) or 3' libraries; "about half of the reads are intronic" replaced by a hedged statement.
- `SEX-COUNT-002`: the log scale separates the modes and compresses the range ("stabilizes Poisson noise" removed).
- `SEX-COUNT-003`: calls from a streamed head only when it reaches the informative reads of `SEX-COUNT-004`; otherwise
  `LOW_SIGNAL`.
- `SEX-DEC-002`: RNA genotype evidence is partly independent of expression (shares reads, mapping and coverage), not
  orthogonal.
- `SEX-LOY-002`: genotype evidence suggests the explanation of `XIST+Y+`; a karyotype is confirmed only with DNA.
- `SEX-REP-001`: the table states that the call is an expression-based screening result; labels unchanged.

### SAMPLE_IDENTITY_GUIDELINES.md [3.1.0]

#### Changed

- `ID-DUP-001`: identical file hashes mean byte-identical files (compare decompressed content for recompressed files);
  identical count vectors flag a probable duplicate to confirm, empty or very sparse vectors excluded; flagged pairs
  MUST NOT be deleted automatically.
- `ID-GENO-002`: renaming contigs presupposes the same assembly; genome build and REF/ALT alleles of the sites are
  checked against the BAMs' reference.
- `ID-GENO-004`: two modes are expected only with same-individual pairs in the cohort; otherwise calibrate with known
  pairs, controls and expected relatedness.
- `ID-POOL-001`: the MUST applies to donors that differ genetically at the covered sites; *Scope* for pools of one
  inbred strain (hashtags or design, limit stated).
- `ID-POOL-004`: animals of one inbred strain cannot be distinguished by genotype; pooling different strains has been
  reported to work [Cardiello 2023].
- `ID-POOL-005`: conflicts between sources are recorded and checked for provenance and fingerprint quality, not
  resolved by the order; unresolved clusters stay unnamed.

### DECONVOLUTION_GUIDELINES.md [2.0.1]

#### Changed

- `DECON-EST-001`: comparisons between groups cancel only a bias common to the groups.
- `DECON-METH-002`: agreement between methods is neither necessary nor sufficient for correctness; a method that does
  not beat NNLS adds no accuracy on the validation criteria and is kept only for its other outputs.
- `DECON-DOWN-001`: *Note* that the estimation error of proportions is not propagated when they enter a model as
  fixed covariates.

### SPATIAL_GUIDELINES.md [3.1.0]

#### Added

- §1: table of which `SCRNASEQ_GUIDELINES.md` QC and normalization rules apply to Visium spots, imaging panels and
  whole-transcriptome data with segmented cells, with the spatial rule that replaces each; §0 points to it.

#### Changed

- `SP-DES-001`: sections from one block are nested in the donor; they can differ anatomically but add no independent
  donors (were "technical replicates").
- `SP-NORM-001`: comparing normalized and unnormalized results is a sensitivity analysis, judged with histology,
  technical controls and replicate sections.
- `SP-NBH-002`: donors are the units (title fixed); sections are aggregated to the donor or modelled as nested in it,
  and the aggregation or model is stated.

### META_ANALYSIS_GUIDELINES.md [2.3.0]

#### Changed

- `META-MOD-003`: DerSimonian-Laird can underestimate $\tau^2$.
- `META-PUB-001`: asymmetry tests only with at least 10 studies; with fewer, a funnel plot MAY be shown as a
  description and is not interpreted.
- `META-STRAT-001`: *Note* that an effect-size threshold used to select is prioritization (`STAT-MT-007`).
- `META-STRAT-002`: the two one-sided analyses form one family, corrected across both directions; a combined p-value
  shows an effect in at least one study, not replication.
- `META-SC-001`: *Why* corrected (pooling ignores study and donor; confounding when conditions are unbalanced across
  studies); *Deviate when* for a joint pseudobulk model when every study contains every compared condition.

### HPC_GUIDELINES.md [1.5.0]

#### Changed

- `HPC-FLOW-001`: staging is verified against checksums of the source (`HPC-FLOW-005`). *Trade-off*: a staged copy
  is a replica, which `CODE-PROV-001` lets be deleted, or overwritten with the source, once the copy of record in
  `ARCHIVE/data/raw` has been verified to exist and match its checksum; the copy of record is never modified or
  deleted.
- `HPC-FLOW-002`: the promotion manifest holds checksums computed on the run before the copy and checked on the copy.
- `HPC-FLOW-005`: verification compares content (`rsync --checksum`, or SHA-256 checksums of the source checked on
  the copy). *Why* added (by default rsync skips files whose size and modification time match), and *Trade-off*
  (checksums read every file on both sides).
- `HPC-FLOW-006`: an equal record count no longer suffices before the BAM is deleted. `samtools checksum -a` of the
  BAM and of the CRAM (samtools 1.22 or later, the same version for both; mate checksums changed in 1.24) MUST match
  in every column, except `+aux` and `combined` when tags were removed on purpose; the checksum reports go into
  provenance.
- `HPC-BATCH-001`: the inputs a batch deletes are replicas (`CODE-PROV-001`), whose copy of record stays in
  `ARCHIVE/data/raw` or is the public archive (`CODE-DL-001`).
- `HPC-JOB-001`: *Why* refers to the submission environment. *Note* added: with `--export=NONE`, only `SLURM_*` and
  SPANK variables come from the submitting shell, and Slurm then loads the user's login environment on the node, as
  with `--get-user-env` (`--export=NIL` does not). `srun` steps inherit `NONE` unless `SLURM_EXPORT_ENV=ALL` (or
  `srun --export=ALL`) is set. `--export=NONE` stays the default.
- `HPC-JOB-004`: *Why* corrected: an array is submitted, monitored and cancelled as one unit, but each task counts as
  a job toward per-user limits (`MaxJobs`, `MaxSubmitJobs`), and the throttle bounds the running tasks. This replaces
  "accounted as one job".
- `HPC-IO-001`: zarr is no longer listed as a single-file container. A zarr directory store writes one file per
  chunk, so zarr is used with sharding (zarr format 3, `shards=`) or in a zip store (`zarr.storage.ZipStore`), with
  chunk and shard shapes matched to access; *Why* added.
- §10 (implementation): the staging and promotion rows describe the checksum verification.

### STATISTICS_REFERENCES.md

#### Added

- [Fay 2010], [Brunner 2000]; metadata and BibTeX from Crossref (the `number={none}` field of the Crossref export
  dropped from [Fay 2010]).

#### Changed

- [Divine 2018]: note names the quantity tested (from the abstract) and the paper's remark on unequal variances and
  the Brunner-Munzel variant, instead of "it tests stochastic dominance".
- [Zhang 2020]: note rewritten from the full text (adjusted integer counts for DE software; power and false-positive
  rates in the authors' simulations, against batch as a covariate); it described the corrected data as being for
  visualisation, which the paper does not say.

### OMICS_REFERENCES.md

#### Added

- [Wu 2020] (scAPAtrap), [Bi 2026] (benchmark of poly(A)-site methods for 3' tag scRNA-seq), [Nie 2026] (benchmark of
  seven ambient-RNA tools), [Finak 2015] (MAST); DOIs verified against Crossref, BibTeX from Crossref (month of
  [Nie 2026] written as the macro `Sep`).

#### Changed

- [Sosina 2021]: truncated note completed from the abstract (DNA-methylation reference; RNA composition, not cellular
  composition).
- [Cardiello 2023]: note adds the inbred-strain result cited by `ID-POOL-004` (from the full text).

### recipes/ [0.8.0]

#### Added

- `scrnaseq/layers`: `input_layer(adata, step)`, the table of `SC-DATA-004` (steps `qc` and `doublets` read
  `counts`; `normalization` and `pseudobulk` read `counts_ambient_corrected` if it exists, else `counts`; an unknown
  step, or a resolved layer absent from `adata.layers`, raises; the choice is logged); `require_integer_counts`, which
  rejects a non-integer layer for count-based DE with an error naming the layer and the remedies (round explicitly and
  record it in `PROJECT_CONTEXT.md`, e.g. SoupX `adjustCounts(roundToInt = TRUE)`, or pass `layer="counts"`), never
  rounding silently; `check_output_layer` and `PROTECTED_LAYERS` (`counts`, `counts_ambient_corrected`, `spliced`,
  `unspliced`, `ambiguous`; `SC-DATA-003`).
- `scrnaseq/units`: `unit_ids` (missing identifiers, and distinct values that become the same string such as `1` and
  `"1"`, raise before any conversion), `unit_metadata` (moved from `scrnaseq/pseudobulk`, which still exports it) and
  `conflicting_values`.
- `scrnaseq/abundance`: `cell_count_table(..., subcomposition=False)` and output column `denominator`
  (`all_cell_types` or `listed_cell_types`).
- `scrnaseq/pc_robustness`: `propose_n_pcs` returns `criteria` (the criteria evaluated) and logs a warning when
  Gavish-Donoho is not evaluated (`SC-PC-002`).
- `scrnaseq/empirical_power`: run status `reference_not_testable` and summary column `n_reference_not_testable`.
- `tests/test_audit_regressions.py`.
- `common/build_gene_groups`: every response a build uses is saved byte for byte in a snapshot directory (default
  `build/gene_groups_sources/<date>`, gitignored), with `SOURCES.json` (URL, parameters, size, SHA-256, download
  date) and `SHA256SUMS`. `--from-snapshot DIR` rebuilds the lists offline after checking every checksum, and the
  result is identical to the original build. The lists record `sources_sha256`.
- `tests/test_build_gene_groups.py` (offline, with responses written by hand).

#### Changed

- Default input layer (behaviour change): `scrnaseq/normalization` (`median_size_factors`, `scranpy_size_factors`,
  `scran_size_factors`, `log_normalize`), `scrnaseq/pseudobulk.aggregate_counts`, `scrnaseq/depth_matching.match_depth`
  and `scrnaseq/empirical_power.empirical_power` take `layer=None`, resolved with `input_layer`: they now read
  `counts_ambient_corrected` when it exists (0.7.0 read `counts`). An explicit layer is honoured; the QC and doublet
  recipes keep `layer="counts"`.
- `stats/permutation.permute_within_strata`: strata are visited in order of first appearance instead of set order, so
  a seed gives the same permutations whatever `PYTHONHASHSEED` (`CODE-REPRO-002`). For the same seed the permutation
  stream differs from 0.7.0, also in `bulk/fdr_calibration` and `scrnaseq/pseudobulk.null_split_de`, which use it.
  `stats/heldout`: ties in the order of strata are broken by first appearance (splits unchanged otherwise).
- `stats/permutation`, `stats/heldout`, `scrnaseq/empirical_power` docstrings and `recipes/README.md`, worded as
  `STAT-RES-003` and `STAT-RES-008`: both two-sided conventions are valid under exchangeability; comparing |T| matches
  a two-sided alternative only when the null distribution of T is centred at 0 and roughly symmetric, so the recipe
  doubles the smaller one-sided p; recall and precision against a held-out reference measure agreement with an
  imperfect reference (recall tends to be optimistic, precision pessimistic), not upper bounds.
- `scrnaseq/doublets.run_scdblfinder`, called directly, applies `ambient_method_of_layer` as `detect_doublets` does
  (`SC-QC-012`).
- `scrnaseq/abundance` imports `scrnaseq/units` instead of `scrnaseq/pseudobulk`, so it no longer needs `pydeseq2`.
- `common/build_gene_groups`: existing lists and a non-empty snapshot are replaced only with `--overwrite`, so a
  rebuild into the shipped directory now needs the flag. Lists are written atomically (temporary file +
  `os.replace`) after every source has been read. `mouse_genes(human, sources=None, file=...)` takes the snapshot.
- `hpc/stage.sh`: `logs/staging.tsv` gains two columns (number of files, `sha256-verified`); `--help` documents that
  a staged replica differing from the source is overwritten with it, and that files present only in the replica
  are left in place.
- `hpc/job.slurm`: comment on `--export=NONE` and on `SLURM_EXPORT_ENV=ALL` before `srun` steps.

#### Fixed

- `bulk/deseq2.validated_counts`: the integer check is exact (`np.allclose` let `100000.4` through and returned
  `100000`); counts at or above 2**63 raise. `run_deseq2` raises a clear error when there are no genes.
- `scrnaseq/pseudobulk`: missing unit or sample identifiers, and identifiers that collide as strings, raise (they
  became a `"nan"` unit or were merged); a cell type with no gene passing the filter is reported in `not_testable`
  (pydeseq2 failed with an unrelated error).
- `scrnaseq/abundance.cell_count_table`: a `cell_types` list that omits observed types raises unless
  `subcomposition=True` (it changed the denominator silently: 1 A and 3 B gave A = 100%); manifest covariates that
  contradict `obs` raise (they overwrote them); the manifest index is validated (missing, repeated or colliding IDs).
- `scrnaseq/pc_robustness`: `elbow_index` raises on constant, non-finite, unordered or too short spectra (a constant
  spectrum returned 1 after a 0/0); `robustness_verdict` raises on clusters of interest absent from the reference, or
  on an empty list (an empty selection was read as robust).
- `scrnaseq/normalization.log_normalize` refuses `key_added` in `PROTECTED_LAYERS` or equal to the input layer
  (`key_added="counts"` overwrote the raw counts; `SC-DATA-003`).
- `scrnaseq/qc.mad_outlier_flags`: fraction metrics with finite values outside [0, 1] (e.g. percentages) raise before
  the logit clip (`SC-QC-003`).
- `scrnaseq/empirical_power`: a split whose reference has no unit with `min_cells` cells, or whose reference DE is not
  testable, is marked `reference_not_testable` (it crashed in `np.vstack`, or reported `ok` against an empty
  reference).
- `hpc/stage.sh`: a staged file with the size and modification time of the source but different content was kept,
  and the script exited 0. The script now copies with `rsync -a --checksum` and checks every file against SHA-256
  checksums computed on the source before the copy (exit 6 and no log line on mismatch). `--dry-run` creates,
  copies and logs nothing. Absolute paths, `.`/`..` components and destinations that resolve outside
  `WORK/data/raw` are refused.
- `hpc/promote.sh`: a run argument such as `../outside` copied and made read-only a directory outside the results
  roots. The run must now be one directory name, and source and destination must resolve directly under
  `WORK/results` and `ARCHIVE/results` before any copy or `chmod`. Checksums are computed on the source before the
  copy (the manifest was computed on the destination), and the copy is checked against them and against the file
  count. `PROMOTION.yaml` is appended to `MANIFEST.sha256`. Runs that already contain `MANIFEST.sha256` or
  `PROMOTION.yaml`, or that hold no files, are refused.
- `hpc/load_site.sh`, `hpc/sites/local.env`, `hpc/sites/hyperion.env`: the scripts failed when `USER` was unset,
  after copying in the case of `stage.sh`. The user is now resolved as `${USER:-$(id -un)}` before any side effect,
  with an early error when neither gives a name. `tests/test_hpc_scripts.py` runs the scripts without `USER` and
  adds regressions for each of these fixes.

### Tools

- `tools/check_sources.py`: a compare whose file list reached the API limit (`COMPARE_FILE_LIMIT`, 300 files) now
  needs review (`needs_review()` is true). `review_reasons()` lists why, and the report shows a "Review needed: ..."
  line. The unused `REVIEW_STATUSES` constant is removed. Tests in `tools/tests/test_check_sources.py`.

### Repository

#### Fixed

- `README.md`: recipes version 0.8.0 (the table said 0.6.0 while `pyproject.toml` said 0.7.0).
- `TODO.md`: follow-ups of this round (reference notes, gene-list snapshot, versioning text, checker coverage).

### Not adopted

- G-C04 (a named constant for every literal): `CODE-NAME-006` already applies only to controlled vocabulary.
- G-C06, G-W02 (moving CODING §§24–25, AGENTIC §§13–16 or writing for funders to other documents): organisation, not
  correctness; not in this round.
- G-B04 (moving the sex and reference rules out of `RNASEQ_BULK_GUIDELINES.md`): the bulk guide keeps the decision
  point and links to `SEX_INFERENCE_GUIDELINES.md`.
- G-SC13 (linear scale for τ): withdrawn by the audit; the cited benchmark evaluates the log scale.
- G-X06 (renaming the `XX`/`XY` labels): labels kept for compatibility; `SEX-REP-001` now states that they are
  expression-based screening results.
- G-I06 (one-to-one donor matching), G-D03 (moving the cross-platform rule to SPATIAL), G-M06 (prediction interval
  for $k \ge 3$): no error found.
- G-A09 (confidential review material, `AGT-PRIV-003`): the policy is kept.
- G-H06 (command and output for each observation of `hpc_sites/DIPC_HYPERION.md`): needs a session on the cluster;
  in `TODO.md`.
- C19 (`n_evaluated` in `scrnaseq/library_qc`): the any-flag column is correct for its contract.
- C22 (citation checks for `CODING_GUIDELINES.md` in `tools/check_guides.py`): in `TODO.md`.


## 2026-11-05 (recipes 0.7.0 for single-cell DE power and pseudoreplication; SCRNASEQ 5.2.1, STATISTICS 2.2.1)

Recipes for the rules added in SCRNASEQ 5.1.0 and STATISTICS 2.2.0 (`SC-DES-007`, `SC-DE-002`, `SC-DE-006` to
`SC-DE-010`, `STAT-RES-008`).

### SCRNASEQ_GUIDELINES.md [5.2.1]

#### Changed

- `SC-DES-007`, `SC-DE-002`, `SC-DE-003`, `SC-DE-006` to `SC-DE-010`: *Recipe* lines for the new functions and modules.

### STATISTICS_GUIDELINES.md [2.2.1]

#### Changed

- `STAT-RES-008`: *Recipe* line (`recipes/stats/heldout.py`).

### recipes/ [0.7.0]

#### Added

- `scrnaseq/pseudobulk`: `check_no_pseudo_replicates` and `aggregate_counts(..., sample_keys=...)`, which raise when
  the cells of one sample (of one donor, in pooled captures) fall in several aggregation units (`SC-DE-006`);
  `detectability_report`, status per gene of interest and cell type, with "no change of interest" only from an
  equivalence test (TOST, Holm-adjusted) against a margin fixed in advance (`SC-DE-009`); `null_split_de`, random
  splits of the control units or permuted labels per cell type through `bulk/fdr_calibration` (`SC-DE-007`);
  `min_cells_diagnostic`, principal components of pseudobulk samples against their number of cells, by bins
  (`SC-DE-002`); `power_inputs`, units, cells per unit and UMIs per cell per cell type and group (`SC-DE-008`).
- `scrnaseq/depth_matching`: `match_depth`, cells per unit subsampled and UMIs thinned binomially to a reference cell
  type, for comparisons between cell types (`SC-DE-008`).
- `scrnaseq/empirical_power`: recall and precision per effect-size bin over a grid of donors and cells per donor,
  against a held-out reference (`SC-DES-007`).
- `stats/heldout`: disjoint, stratified, seeded reference and validation splits (`STAT-RES-008`).
- `tests/test_scrnaseq_de_power.py`.

#### Changed

- `scrnaseq/pseudobulk`: `pseudobulk_de` adds `lfc_ci_low`, `lfc_ci_high` (Wald interval of the unshrunken log fold
  change), `rank` (by p-value within the cell type) and `significance_quartile` (`SC-DE-010`), and accepts `alpha` and
  `shrink`.

## 2026-11-04 (pseudobulk aggregation function, detection, distribution and variability; SCRNASEQ 5.2.0)

Which aggregation function to use for pseudobulk (`scanpy.get.aggregate`: `sum`, `mean`, `median`, `count_nonzero`,
`var`) and whether percentiles, per-donor models (ZINB) or Gini and skewness indices are alternatives. Full texts read:
[Crowell 2020], [Zimmerman 2021a], [Squair 2021], [Junttila 2022], [Murphy 2022], [Gilis 2025], [Bouland 2023],
[Hoffman 2026]; the other new references from abstracts.

### SCRNASEQ_GUIDELINES.md [5.2.0]

#### Added

- `SC-DE-011`: `mean`, `median` and `var` are not the default input for DE between conditions; mean of log-normalized
  values only with limma-trend, with log fold changes from summed counts; `median` not for DE.
- `SC-DE-012`: differential detection on binarized pseudobulk counts (`count_nonzero`), with cell numbers as trials and
  a normalization offset, combined with DE by stage-wise testing.
- `SC-DE-013`: changes in variability or distribution shape tested with methods that separate them from mean and depth
  (distinct, IDEAS, memento, BASiCS); no comparison of raw Gini indices, skewness, percentiles or fractions detected.
- `SC-DE-014`: UMI counts modelled without zero inflation; zero-inflated or hurdle models for read counts without UMIs.

#### Changed

- `SC-DE-002`: step 1 sums raw counts from the counts layer, never `X` after normalization (`sc.get.aggregate` uses `X`
  by default), and normalizes after aggregation; *Why* cites the benchmarks that compared sum and mean.
- `SC-DE-004`: precision-weighted pseudobulk (dreamlet) for unequal cell numbers or repeated measures.
- `SC-ANNOT-005`: specificity across cell types summarized with tau.
- `SC-REP-006`: the aggregation function and its layer are reported.
- §16 tools table: dreamlet; new row for differential detection, distribution and variability.

### STATISTICS_REFERENCES.md

#### Added

- [Murphy 2022], [Thurman 2021], [Gilis 2025], [Bouland 2021], [Bouland 2023], [Hoffman 2026], [Kim 2024],
  [Lun 2018] (preprint), [Svensson 2020], [Kim 2020], [Jiang 2022], [Tiberi 2023], [Zhang 2022], [Hafner 2025],
  [Eling 2018], [Jiang 2016]; metadata and BibTeX from Crossref, months as BibTeX macros.

#### Changed

- [Crowell 2020], [Junttila 2022]: notes state the result of the comparison of aggregation functions.

### OMICS_REFERENCES.md

#### Added

- [Kryuchkova-Mostacci 2017]; BibTeX from Crossref, year, month, volume and pages of the journal issue.

### Repository

#### Changed

- `references/SOURCE_REPOS.toml`: the deferred differential-expression chapter of the single-cell best-practices book
  notes that the aggregation function is settled (`SC-DE-002`, `SC-DE-011`).

## 2026-11-04 (single-cell DE: power and pseudoreplication; SCRNASEQ 5.1.0, STATISTICS 2.2.0, RNASEQ_BULK 3.3.0)

Rules from [Crampton 2026] (preprint: empirical power of pseudobulk DE in 1,494 human brain donors, with a held-out
reference set) and from a full reading of [Squair 2021] (why cell-level tests produce false discoveries).

### SCRNASEQ_GUIDELINES.md [5.1.0]

#### Added

- `SC-DES-007`: justify the number of donors against the smallest effect of interest and the expression of the genes of
  interest; simulation tools with parameters from large external data, not from a pilot; cells per donor of the cell
  type of interest before depth per cell; budget in UMIs per donor.
- `SC-DE-006`: no replicates made by partitioning the cells of one sample (pseudo-replicates).
- `SC-DE-007`: global-null diagnostic (random split of control donors, or permuted donor labels) with
  `recipes/bulk/fdr_calibration.py`.
- `SC-DE-008`: the number of DE genes per cell type is not a measure of response; report donors, cells per donor and
  UMIs per cell; compare at matched depth or by effect sizes.
- `SC-DE-009`: a gene absent from a DE list is not evidence of no change; report its expression and the interval of its
  log fold change.
- `SC-DE-010`: DE genes near the FDR cut-off are the least reproducible; ranks and adjusted p-values reported, genes for
  validation chosen from the most significant.

#### Changed

- `SC-DE-001`: names the marker functions that run cell-level tests (`scanpy.tl.rank_genes_groups`, Seurat
  `FindMarkers`); *Why* states the mechanism (variance between donors ignored, underestimated most for highly
  expressed genes).
- `SC-DE-002`: the minimum number of cells per pseudobulk sample is checked on the data (principal components against
  the number of cells), fixed before the results and applied to all groups.
- `SC-DE-004`: *Trade-off* and *Deviate when* for mixed models.
- `SC-DES-002`: points to `SC-DES-007`; cohorts of 3–5 donors are powered only for large effects.
- `SC-REP-006`: donors tested, cells per donor and UMIs per cell reported per cell type.
- Citation key [Zimmerman 2021†] renamed [Zimmerman 2021a†].

### STATISTICS_GUIDELINES.md [2.2.0]

#### Added

- `STAT-RES-008`: empirical power or replication estimated by subsampling defines the reference on disjoint units,
  stratified and repeated; the estimates are upper bounds.
- `STAT-ENR-006`: a coherent enrichment result does not validate the gene list it was computed from.

#### Changed

- `STAT-HIER-002`: citation key [Zimmerman 2021] renamed [Zimmerman 2021a].

### RNASEQ_BULK_GUIDELINES.md [3.3.0]

#### Changed

- `BULK-DE-003`: log fold changes compared across DE tools are shrunk with one common estimator.

### SPATIAL_GUIDELINES.md [3.0.1]

#### Changed

- Citation key [Zimmerman 2021†] renamed [Zimmerman 2021a†].

### STATISTICS_REFERENCES.md

#### Added

- [Crampton 2026] (preprint), [Schmid 2021] (scPower), [Zimmerman 2021b] (hierarchicell) and [Vieth 2017] (powsimR);
  metadata and BibTeX from Crossref, months as BibTeX macros, notes from the abstracts.

#### Changed

- [Squair 2021]: the note states the mechanism and the pseudo-replicate result.
- [Zimmerman 2021] renamed [Zimmerman 2021a] (`.bib` key `Zimmerman_2021a`), to tell it from [Zimmerman 2021b].


## 2026-11-04 (established pipelines: CODING 4.0.0, SCRNASEQ 5.0.0)

### CODING_GUIDELINES.md [4.0.0]

#### Added

- §12 *Established pipelines*: `CODE-PIPE-001` (run a processing stage with the established pipeline the domain guide
  names; what counts as established; unreleased pipelines pinned to a tested commit), `CODE-PIPE-002` (domain rules fix
  what is computed; five implementation levels, from parameters to an external implementation, for pipelines and
  recipes; levels 4 and 5 need confirmation), `CODE-PIPE-003` (a modified step mirrors the original; `ext.args` in a
  `withName` block replaces the pipeline's string), `CODE-PIPE-004` (`pipeline_deviations` YAML block with typed return
  conditions), `CODE-PIPE-005` (check the deviations at each release change).

#### Changed

- `CODE-PROJ-008`: deviations from a pipeline or recipe go in the `pipeline_deviations` block.
- `CODE-RCP-002`: adapting a recipe follows the levels of `CODE-PIPE-002`.

#### Removed

- `CODE-RCP-003` (retired): merged into `CODE-PIPE-002` (level 4) and `CODE-PIPE-004`.

### SCRNASEQ_GUIDELINES.md [5.0.0]

#### Added

- `SC-QC-016`: from matrix to clustered object, `nf-core/scdownstream` pinned to a commit of `dev` until its first
  release, in two passes with the review of the recipes between them.
- §16: *Baseline parameters for `nf-core/scdownstream`* (read at commit `0048567`) and its differences from this guide.
- `SC-QC-003`: minimum scaled MAD of 1.5 percentage points for fraction metrics; a scaled MAD of 0 in log counts or
  log genes is not evaluated.

#### Changed

- `SC-QUANT-001`: STARsolo through `nf-core/scrnaseq`; one gene feature plus `Velocyto` per run (`Gene` for cells,
  `GeneFull_Ex50pAS` for nuclei, through a pinned fork while the schema does not accept it); extra STAR options and
  `EmptyDrops_CR` in `ext.args`; `GeneFull` in a separate run.
- `SC-QC-006`: default ribosomal and haemoglobin lists of Ensembl IDs ship with the recipes; an empty group is not
  evaluated, never 0.
- §16 *Implementation*: quantification and matrix-to-object rows name the pipelines.
- [Curion 2024] is now [Curion 2024a].

#### Removed

- `SC-QC-003`: the absolute threshold as the way out of a MAD of 0 (replaced by the minimum MAD).

### SAMPLE_IDENTITY_GUIDELINES.md [3.0.0]

#### Changed

- `ID-POOL-001`: without donor genotypes, `nf-core/hadge` [Curion 2024b]; with genotypes, Demuxalot and Dropulation
  from the Demuxafy image combined with `Combine_Results.R` (level 5, until hadge includes them); hashtags through
  hadge; demuxlet runs inside hadge.
- §9 *Implementation*: demultiplexing row.

### SPATIAL_GUIDELINES.md [3.0.0]

#### Changed

- `SP-SEG-001`: segmentation through `nf-core/sopa` [Blampey 2024]; Mesmer through sopa's `custom_staining_based`
  (level 3).
- `SP-SEG-004`: SPLIT runs on the SpatialData output of sopa. `SP-SEG-005`: with nf-core/sopa, the release and its
  parameters file record the segmentation.
- §13 *Implementation*: segmentation row.

### META_ANALYSIS_GUIDELINES.md [2.2.0]

#### Changed

- `META-DATA-001`: propose the reprocessing pipeline to the user, checking at that moment the latest release, other
  nf-core pipelines and uniformly processed resources.

### HPC_GUIDELINES.md [1.4.0]

#### Changed

- `HPC-WF-003`: three tracked configuration files (site, parameters, project configuration of each pipeline); a
  patched pipeline runs from a pinned fork, never from an edited copy in `$NXF_HOME/assets`.

### recipes/ [0.6.0]

#### Added

- `common/gene_groups`: reader of the default gene lists (`recipes/common/data/gene_groups/`, human and mouse,
  ribosomal and haemoglobin; Ensembl ID, symbol, source, download date, Ensembl release, SHA-256 of the HGNC file) and
  of the symbols the guide excludes.
- `common/build_gene_groups`: builds the lists from the HGNC REST API (human gene groups) and, for mouse, from the MGI
  orthologues that HGNC records, mapped to Ensembl protein-coding genes with BioMart (Ensembl Compara lacks one-to-one
  orthologues for several core ribosomal proteins); fails when a group is empty or an excluded symbol remains.

#### Changed

- `scrnaseq/qc`: minimum scaled MAD for fraction metrics (`MIN_SCALED_MAD` = 0.015, `min_scaled_mad`); report
  statuses `mad_floor_applied` and `not_evaluated_no_values`; `absolute_threshold_only` removed.
- `scrnaseq/qc`: `flag_gene_groups(..., species=...)` uses the default lists and reports the genes of the list absent
  from the matrix; `qc_metrics` adds `ribo_fraction` and `hb_fraction`, missing when the group has no gene in the matrix.

### Repository

#### Added

- `references/OMICS_REFERENCES.md` and `.bib`: [Curion 2024b] (hadge) and [Blampey 2024] (sopa).
- `RULE_IDS.tsv`: `CODE-PIPE-001` to `CODE-PIPE-005`, `SC-QC-016`; `CODE-RCP-003` retired.
- `tools/rule_ids.py`: subtopic `PIPE` for *Established pipelines*.
- `TODO.md`: script to check the `pipeline_deviations` block; items for the next round.

## 2026-11-03 (SCRNASEQ 4.1.0: library QC and cell filters)

### SCRNASEQ_GUIDELINES.md [4.1.0]

#### Added

- `SC-QUANT-006`: check every library before its cells enter QC: FastQC and MultiQC; from STARsolo `Summary.csv`
  and `Barcodes.stats`, reads, valid barcodes, Q30, saturation, mapping, estimated cells, reads in cells, UMIs and
  genes per cell; libraries compared with each other and with the cells loaded.
- `SC-QUANT-007`: recipe `recipes/scrnaseq/library_qc.py`.
- `SC-DATA-011`: library-level metrics in `adata.uns["library_qc"]`, one row per `sample_id`, added after the
  captures are concatenated.
- `SC-QC-015`: check what the filters remove from each population: flagged fraction per cluster after the first
  clustering, inspection of clusters that lose much more than the others [Germain 2020].
- §0: upstream sources (single-cell best-practices book, panpipes [Curion 2024]).

#### Changed

- `SC-QC-003`: log total counts and log genes, both tails; mitochondrial fraction upper tail only, on the logit
  scale, flagged only above a floor stated by the project [Germain 2020]; fraction of counts in the 20 most expressed
  genes, both tails.
- `SC-QC-006`: ribosomal and haemoglobin fractions are diagnostics, not filters.
- `SC-QC-008`: the recipe description lists the logit mitochondrial rule with a floor, the top-20 fraction and the
  flagged fraction per cluster.
- §0 scope: pathway and transcription-factor activity listed; multimodal assays named. Pooled CRISPR screens, which
  the draft named with a reference to a section that does not exist, are moved to `TODO.md`.

### recipes/ [0.5.0]

#### Added

- `scrnaseq/library_qc`: reads STARsolo `Summary.csv` (feature name detected in the rows; `NoMulti` kept missing) and
  `Barcodes.stats` (STAR 2.7.11b names) in the nf-core/scrnaseq output layout; one row per library with the STAR
  version and source files; flags libraries more than 3 scaled MADs from the other libraries in the direction of a
  problem (not evaluated with fewer than 3 libraries or MAD = 0); stores the table in `uns["library_qc"]` with nullable
  booleans and categorical strings, which `write_h5ad` accepts in current and older anndata versions. Tests on
  synthetic outputs.
- `scrnaseq/qc`: `qc_metrics` adds `top20_fraction`; `mad_outlier_flags` flags it on both tails (`n_mads_top20`);
  new `flagged_fraction_by_cluster` (`SC-QC-015`), with clusters marked for inspection at twice the median flagged
  fraction and at least 5 points above it (project defaults).

#### Changed

- `scrnaseq/qc`: the mitochondrial MAD rule is computed on the logit scale (`MITO_LOGIT_EPS`) and, with the new
  `mito_floor`, flags a cell only above the floor; without a floor the report status is `evaluated_no_floor`. The
  report gains `scale` and `floor`; `median` and `upper` of the mitochondrial row are on the fraction scale.

### Repository

#### Added

- `references/OMICS_REFERENCES.md` and `.bib`: [Curion 2024] (panpipes) and [Germain 2020] (pipeComp), cited by
  SCRNASEQ 4.1.0.
- `RULE_IDS.tsv`: `SC-DATA-011`, `SC-QUANT-006`, `SC-QUANT-007`, `SC-QC-015`.
- `TODO.md`: planned work not yet in the guides.

## 2026-10-26 (single-cell sources)

### Repository

#### Added

- `references/SOURCE_REPOS.toml`: three single-cell sources, reviewed against `SCRNASEQ_GUIDELINES.md` and nothing
  adapted yet: the single-cell best-practices book (Theis lab, Apache-2.0, `dd6eb86`), panpipes (Dendrou lab,
  BSD-3-Clause, `81721cf`) and its inactive community edition (`4a4d914`, no licence stated). Chapters and files
  with proposed changes are `deferred` to a SCRNASEQ round or a multimodal round; the spatial and deconvolution
  chapters are `watch`, not yet compared with their guides.

## 2026-10-26 (SCRNASEQ citation fixes)

### SCRNASEQ_GUIDELINES.md [4.0.1]

#### Fixed

- `SC-QC-009`: DoubletDetection cites [Gayoso 2025], the key of the DoubletDetection entry since the BibTeX fixes of
  the same day (the rule cited the former key [Gayoso 2020]); Solo cites [Bernstein 2020] (Solo, *Cell Systems*
  2020), not [Bernstein 2025], which is in neither reference file.

## 2026-10-26 (agent behaviour from superpowers)

### AGENTIC_GUIDELINES.md [1.1.0]

#### Added

- §13 *Eliciting intent* (`AGT-INTENT-001` to `AGT-INTENT-010`): a private check of what is known about what the
  user wants and why; look before asking; a first message that restates the request and asks the one question that
  would most change the result; one open question per message; play-back in parts, with the agent's proposals and
  defaults marked as such; two or three approaches with a recommendation; sizing and a written design; a check of
  the design by a fresh reader playing the implementer; no building of a project before its design is approved,
  except after an opt-out or for an agreed throwaway spike.
- §14 *Plans* (`AGT-PLAN-001` to `AGT-PLAN-005`): what a plan records, task size, a specification's silence is not
  permission, self-review against the specification, execution from the plan text.
- §15 *Reviews* (`AGT-REV-001` to `AGT-REV-006`): fresh-context, read-only review from files; findings by severity
  with a verdict and a list of what was set aside; no pre-judging in the reviewer's brief; receiving feedback by
  checking it; no performative agreement; external suggestions checked before adoption.
- §16 *Writing instructions for agents* (`AGT-INSTR-001` to `AGT-INSTR-008`): instructions written for observed
  failures, in the form that fits the failure; no nuance or exemption clauses; wording tested against a
  no-guidance control; descriptions say when, not what; mechanical checks automated; persuasive framing only in the
  user's interest; proposals upstream state their motivating case and who drafted them.
- `AGT-EVID-007` (read the current text of an instruction), `AGT-VER-007` to `AGT-VER-009` (the check that proves a
  claim, completion item by item, checks on what is delivered), `AGT-CLAIM-008` (every observed failure reported),
  `AGT-BATCH-007` to `AGT-BATCH-009` (rulings recorded and reported, no check-ins within a batch, the ledger
  outranks memory), `AGT-DEL-007` to `AGT-DEL-011` (context constructed for subagents, report contract, workers do
  not dispatch agents, parallel work only on independent problems, models chosen by turns), `AGT-LOOP-006` (three
  failed attempts: question the approach), `AGT-REP-005`.

#### Changed

- `AGT-SCOPE-001`: points to §13 for open-ended work. `AGT-BATCH-002`: irreversible, security-sensitive and unnamed
  external actions, and a broken plan, also stop a batch; ambiguities of implementation do not.
- *What to record in `PROJECT_CONTEXT.md`* moves from §13 to §17; identifiers unchanged.
- Adapted in part from superpowers (Jesse Vincent and Prime Radiant, <https://github.com/obra/superpowers>, MIT
  licence, commit `2ce8aa1`): skills `brainstorming`, `writing-plans`, `executing-plans`, `requesting-code-review`,
  `receiving-code-review`, `verification-before-completion`, `subagent-driven-development`,
  `dispatching-parallel-agents`, `systematic-debugging` (one rule) and `writing-skills` with its references, and the
  design notes on brainstorming and on positive instructions. Test-driven development, test quality and debugging
  techniques are deferred to a later round of `CODING_GUIDELINES.md`.

### STATISTICS_REFERENCES.md

- [Meincke 2026] (persuasion principles and LLM compliance, PNAS), section L; metadata checked against Crossref,
  figures against the abstract (Europe PMC). Cited by `AGENTIC_GUIDELINES.md` §16.

### Repository

#### Changed

- `tools/rule_ids.py`: topics of `AGENTIC_GUIDELINES.md` §13 to §17.
- `README.md`: licence note with the attribution of the text adapted from superpowers.

## 2026-10-26 (upstream sources: registry and monthly check)

### Repository

#### Added

- `references/SOURCE_REPOS.toml`: registry of the upstream repositories the guides adapt (Co-Scientist toolkit at
  `d1d35f5`, superpowers at `2ce8aa1`), with the commit last reviewed and the files followed in each, by status
  (`integrated`, `deferred`, `watch`, `not-integrated`) and destination.
- `tools/check_sources.py`: compares each repository with its reviewed commit through the GitHub compare API and
  reports the followed files that changed; `mark-reviewed` records the commit a review covered. Tests in
  `tools/tests/test_check_sources.py` (offline, with compare responses written by hand).
- `.github/workflows/sources.yml`: runs the check on the first day of each month (and on demand) and opens, or
  comments on, an issue labelled `source-review` when a followed file with status `integrated`, `deferred` or
  `watch` changed.

#### Changed

- `README.md`: *Maintenance: upstream sources*; the copy example copies only the reference files
  (`references/*_REFERENCES.*`).

## 2026-10-26 (doublets and ambient RNA correction: [Meier 2026], [Neavin 2024])

### SCRNASEQ_GUIDELINES.md [4.0.0]
#### Changed (breaking)

- `SC-QC-009`: the two default detectors are Scrublet and scDblFinder (R, through a script); DoubletDetection or Solo
  MAY be added as a third detector in large captures. DoubletFinder is no longer listed as an alternative: on in silico
  pools with known doublets it had one of the lowest MCC values and flagged more singlets as doublets [Neavin 2024].
  New columns `doublet_n_evaluated` and `doublet_n_calls`; `doublet_consensus` requires every evaluated detector and at
  least two.

#### Added

- `SC-QC-012`: doublets are detected on counts without ambient correction; SoupX- or decontX-corrected counts MUST
  NOT be the input of the detectors or of the doublet review; with these methods the whole doublet procedure, first
  clustering included, runs on `counts` and the correction is applied afterwards; CellBender-corrected counts MAY be
  used as a recorded deviation [Meier 2026].
- `SC-QC-013`: the Scrublet threshold is checked per capture and set by hand when the histogram of simulated scores
  is not bimodal; a failed automatic threshold MUST NOT be read as a capture without doublets.
- `SC-QC-014`: in pooled captures, genotype-demultiplexing doublets are removed in addition to `SC-QC-010`, and
  transcriptional detection is still run (up to 20 donors per capture) [Neavin 2024].

#### Changed

- `SC-DATA-001`: `adata.uns["ambient_correction"]` records the method and version behind `counts_ambient_corrected`.
- `SC-DATA-004`: the `counts` row covers the co-expression review; *Why* cites [Meier 2026] and states that, with SoupX
 or decontX, the corrected layer is created after doublet removal.
- §4 introduction: order of the QC steps with CellBender and with SoupX or decontX.
- `SC-QC-001`: *Why* adds that CellBender is the only compared method that models doublets; *Deviate when* adds decontX
  for data without the raw matrix, and both SoupX and decontX run after doublet removal.
- `SC-QC-010`: step 1 on `counts` with SoupX or decontX; step 3, co-expression measured on `counts`; *Note* corrected:
  [Neavin 2024] compared both combination rules.
 `SC-QC-011`: the recipe runs scDblFinder through a script, handles Scrublet thresholds and checks the input layer.
- `SC-REP-003`: report the ambient-correction method, the layer read by doublet detection and the Scrublet thresholds.
- §16: doublet tools updated.



## 2026-10-26 (reference files: verification and BibTeX fixes)

Every entry of both `.bib` files was checked against Crossref, OpenAlex, Semantic Scholar and arXiv with the
`reference-verifier` skill of the Co-Scientist toolkit (Shocklab, University of Cape Town,
<https://github.com/shocklab/co-scientist-toolkit>), and Zenodo DOIs against DataCite: all 347 references exist and
match their titles and DOIs. Years differ from the databases' print dates in 39 entries only because the files
use the year of first publication; this is now stated in both headers.

### SCRNASEQ_GUIDELINES.md [3.2.1]

#### Changed

- `SC-QC-009` (doublets): cites [Gayoso 2025] (DoubletDetection 4.3.0.post1, the 4.x defaults the recipes use)
  instead of [Gayoso 2020].

### STATISTICS_REFERENCES.md

#### Fixed

- `.bib`: six Krzywinski & Altman columns shared two keys; they now use the keys of the Markdown list
  (`Krzywinski_2013a` to `Krzywinski_2014c`).
- `.bib`: month values `Sept`, `June` and `July`, which are not BibTeX month macros and stop some parsers, replaced
  by `sep`, `jun` and `jul`.
- Header: years are those of first publication, online or in print (Crossref's *issued* date).

### OMICS_REFERENCES.md

#### Fixed

- `.bib`: duplicated keys `Love_2016`, `Soneson_2016`, `Palla_2022` and `Zhao_2021` split into the `a`/`b` keys of the
  Markdown list; keys `2019` (two entries), `2024` and `2017` renamed `Lun_2019`, `Thompson_2019`, `CZI_2024` and `Tukiainen_2017`.
- `.bib`: four entries had an empty first author where Crossref lists a group author; the group is kept as first
  author for `CZI_2024` (cited as such in the Markdown list) and dropped from the author lists of `Lun_2019`,
  `Thompson_2019` and `Tukiainen_2017`.
- `.bib`: month values `Sept`, `June` and `July` replaced by `sep`, `jun` and `jul`.
- [Gayoso 2020] replaced by [Gayoso 2025]: the entry cited version 3.0 with five authors under the concept DOI of
  all releases, which resolves to the latest one; it now cites version 4.3.0.post1 (the 4.x defaults the recipes
  use) with its version DOI and the two authors of its DataCite record.
- [Fortmann 2023]: "version 0.1" is the package version; the entry now gives the commit the recipes pin (62292a0,
  five commits after release v1.1) and states that the DOI is the concept DOI.
- Header: software DOIs are checked against DataCite; years are those of first publication, online or in print.

### Repository

#### Added

- `tools/check_guides.py`: checks of the `.bib` files: unique keys, no empty author names, months as BibTeX macros,
  and the same DOIs as the Markdown list, one entry per DOI on each side. Run on the files before the fixes above,
  it reports every problem fixed here.
- `tools/tests/test_check_guides.py`: tests of these checks (`CODE-TEST-009`); `tools/tests` added to the pytest
  paths.

#### Changed

- `tools/ci_scope.py`: a change under `tools/` also runs the Python test job, which now includes the tool tests;
  `README.md` CI table updated.

## 2026-10-26 (pre-registration)

### STATISTICS_GUIDELINES.md [2.1.0]

#### Added

- `STAT-DES-009`: the analysis plan of a confirmatory claim (hypotheses with direction and rough size, primary
  outcome, model, inclusion and QC rules, multiple-testing family, decision rule) is fixed, dated and committed
  before the analyses that produce the outcome run; with existing data, state what had already been seen.
- `STAT-DES-010`: forking-paths audit of the choices that could be tuned on the outcome; choices made afterwards are
  justified independently or analysed under each alternative.
- `STAT-DES-011`: departures from the plan are appended with date and reason, never edited in; no hypothesis formed
  after the results presented as stated in advance.
- `STAT-DES-012`: pilots decide feasibility under a rule written beforehand; sample size is not based on the
  pilot's effect size.
- `STAT-REP-008`: report whether each analysis was planned or exploratory, the date of the plan and its departures.
- `STAT-DES-010` to `STAT-DES-012` adapt the `pre-register` skill of the Co-Scientist toolkit (Shocklab, University of
  Cape Town, <https://github.com/shocklab/co-scientist-toolkit>, CC BY 4.0).

### STATISTICS_REFERENCES.md

- [Nosek 2018], [Gelman & Loken 2014], [Kerr 1998] and [Albers 2018], section B; metadata checked against Crossref,
  content against the abstracts (Gelman & Loken 2014 against the authors' copy of the article).

### Repository

#### Changed

- `README.md`: licence note lists the `STATISTICS_GUIDELINES.md` rules adapted from the Co-Scientist toolkit.

## 2026-10-26 (verification, tests and shell rules)

### CODING_GUIDELINES.md [3.7.0]

#### Added

- `CODE-NUM-013`: round only when presenting; computations, bounds and stored tables use full precision, and rounded
  values are never read back as inputs.
- `CODE-PROV-006`: broken or superseded runs go to `results/_superseded/<run>/` with a `DO_NOT_USE.md`; loaders
  raise when given one, unless an explicit, logged argument allows it.
- `CODE-TEST-009`: every input validation has a test that triggers it; a pipeline's end-to-end test has a negative
  counterpart on a corrupted fixture.
- `CODE-TEST-010`: expected values are established independently of the code under test; a failing test's
  expectation is not changed to the observed value until it is established which is right; snapshot tests are
  marked as such.
- `CODE-VAL-004`: run the changed entry point and inspect its output where the next step reads it, not only the
  tests.

#### Changed

- `CODE-LANG-004`: quoted heredocs; waiting on background jobs by PID, not `pgrep -f`; `timeout` checked to exist
  (`timeout=120 command` is an environment assignment); background processes stopped on exit.
- `CODE-DOC-005`: labels coined during development are replaced by what they name.
- Adapted in part from `guidance/CLAUDE-starter.md` (*Shell commands*, *Verifying work*, *Writing*) and
  `guidance/paper-CLAUDE.md` of the Co-Scientist toolkit (Shocklab, University of Cape Town,
  <https://github.com/shocklab/co-scientist-toolkit>, CC BY 4.0).

### MANUSCRIPT_GUIDELINES.md [1.1.0]

#### Changed

- `MS-TRACE-003`: superseded runs go to `results/_superseded/<run>/` as `CODE-PROV-006` requires (previously
  `results/_ARCHIVE/`, a name easily confused with the `ARCHIVE` storage role of `HPC_GUIDELINES.md`, which holds
  promoted results).

### Repository

#### Changed

- `README.md`: licence note lists the `CODING_GUIDELINES.md` rules adapted from the Co-Scientist toolkit.

## 2026-10-26 (manuscript guidelines)

### MANUSCRIPT_GUIDELINES.md [1.0.0]

#### Added

- New guide (prefix `MS`) for manuscripts built on the project's analyses. Layout and inputs: manuscript in `paper/`
  in a text format, cited papers as read-only inputs with their converted text, one vocabulary for text, figures
  and code (§1); traceability: figures and tables from scripts with a sidecar naming run and commit, numbers in the
  text from a generated values file, superseded runs archived and refused by figure scripts, regeneration and tag
  before each submission, availability statements (§2); claims: a claims table tying load-bearing claims to
  results, scope matching the design, exploratory analyses labelled, limitations, reporting standards (§3);
  citations: no reference without its record, automated bibliography check, support of each load-bearing citation
  checked in the cited text, primary sources, reading limits of `AGENTIC_GUIDELINES.md` §10 (§4); pre-submission
  review along independent dimensions, findings with location and severity, venue criteria looked up, changes only
  with the author's approval and verified afterwards (§5); revisions and response letters (§6); writing (§7); what
  to record in `PROJECT_CONTEXT.md` (§8).
- Adapted in part from the Co-Scientist toolkit (Shocklab, University of Cape Town,
  <https://github.com/shocklab/co-scientist-toolkit>, CC BY 4.0): `guidance/paper-CLAUDE.md` (adapted there from
  Dominik Lukeš's agents-for-reproducibility course), `guidance/CLAUDE-starter.md` (*Writing*, *Research papers*)
  and the skills `paper-review` (templates adapted from Dominik Lukeš's paper-reviewer skill, MIT),
  `claim-fact-checker` (adapted from Dominik Lukeš's agent skills, MIT), `reference-verifier` and
  `convergence-check`.

### CODING_GUIDELINES.md [3.6.0]

#### Changed

- `CODE-USE-008`: prefix `MS` added to the list of guide prefixes.
- §2 (documentation set): the manuscript and its claims table in `paper/` added to the optional documents.
- §26 (related guidelines): lists `MANUSCRIPT_GUIDELINES.md`.

### STATISTICS_REFERENCES.md

- [Jergas 2015] (quotation accuracy) and [Walters 2023] (fabricated references generated by ChatGPT), section L;
  metadata checked against Crossref. Cited by `MANUSCRIPT_GUIDELINES.md` §4.

### Repository

#### Changed

- `tools/rule_ids.py`: prefix `MS` and its topics.
- `README.md`: contents table and licence note with the attribution of the adapted text.
## 2026-10-26 (doublets and ambient RNA correction: [Meier 2026], [Neavin 2024])

### SCRNASEQ_GUIDELINES.md [4.0.0]

#### Changed (breaking)

- `SC-QC-009`: the two default detectors are Scrublet and scDblFinder (R, through a script); DoubletDetection or Solo
  MAY be added as a third detector in large captures. DoubletFinder is no longer listed as an alternative: on in silico
  pools with known doublets it had one of the lowest MCC values and flagged more singlets as doublets [Neavin 2024].
  New columns `doublet_n_evaluated` and `doublet_n_calls`; `doublet_consensus` requires every evaluated detector and at
  least two.

#### Added

- `SC-QC-012`: doublets are detected on counts without ambient correction; SoupX- or decontX-corrected counts MUST
  NOT be the input of the detectors or of the doublet review; with these methods the whole doublet procedure, first
  clustering included, runs on `counts` and the correction is applied afterwards; CellBender-corrected counts MAY be
  used as a recorded deviation [Meier 2026].
- `SC-QC-013`: the Scrublet threshold is checked per capture and set by hand when the histogram of simulated scores
  is not bimodal; a failed automatic threshold MUST NOT be read as a capture without doublets.
- `SC-QC-014`: in pooled captures, genotype-demultiplexing doublets are removed in addition to `SC-QC-010`, and
  transcriptional detection is still run (up to 20 donors per capture) [Neavin 2024].

#### Changed

- `SC-DATA-001`: `adata.uns["ambient_correction"]` records the method and version behind `counts_ambient_corrected`.
- `SC-DATA-004`: the `counts` row covers the co-expression review; *Why* cites [Meier 2026] and states that, with SoupX
  or decontX, the corrected layer is created after doublet removal.
- §4 introduction: order of the QC steps with CellBender and with SoupX or decontX.
- `SC-QC-001`: *Why* adds that CellBender is the only compared method that models doublets; *Deviate when* adds decontX
  for data without the raw matrix, and both SoupX and decontX run after doublet removal.
- `SC-QC-010`: step 1 on `counts` with SoupX or decontX; step 3, co-expression measured on `counts`; *Note* corrected:
  [Neavin 2024] compared both combination rules.
- `SC-QC-011`: the recipe runs scDblFinder through a script, handles Scrublet thresholds and checks the input layer.
- `SC-REP-003`: report the ambient-correction method, the layer read by doublet detection and the Scrublet thresholds.
- §16: doublet tools updated.

### SAMPLE_IDENTITY_GUIDELINES.md [2.3.0]

#### Changed

- `ID-POOL-001`: with genotypes, Demuxalot with refined genotypes (combined with Dropulation when more than one tool can
  be run); without genotypes, vireo above 20 donors or as the only tool, souporcell with vireo otherwise; combined
  calls require more than half of the tools to call a singlet [Neavin 2024]; demuxlet as a deviation.
- `ID-POOL-009`: genotype-free demultiplexers are interchangeable for labelling, not for singlet and doublet calls.
- §9: demultiplexing tools updated.

### OMICS_REFERENCES.md

- Two references: [Meier 2026] (preprint) and [Yang 2020] (DecontX).

### recipes/ [0.4.0]

#### Changed

- `doublets`: `detect_doublets` runs the detectors named in `detectors` (default Scrublet and scDblFinder, the latter
  through the new `scdblfinder.R`; DoubletDetection optional, imported only when used); manual Scrublet thresholds per
  capture (`scrublet_thresholds`); a failed automatic threshold raises, or marks Scrublet as not evaluated for that
  capture (`on_scrublet_failure`), instead of being recorded as zero doublets; simulated Scrublet scores stored in
  `uns["doublets"]["scrublet"]`; new columns `doublet_n_evaluated` and `doublet_n_calls`; `combine_calls` takes two or
  more call vectors; `doublet_removal_mask` adds genotype doublets (`genotype_key`).
- `doublets`: new `ambient_method_of_layer`; `detect_doublets` and `doublet_cluster_review` reject
  `counts_ambient_corrected` when `uns["ambient_correction"]` names SoupX or decontX or names no method, and warn on
  layers other than `counts` and `counts_ambient_corrected`; `doublet_cluster_review` warns when a SoupX or decontX
  layer already exists (the first clustering must not use it). The per-capture summary reports `input_layer` and
  `ambient_method`.
- CI: R-backed job installs scDblFinder and SingleCellExperiment; the wheel check includes `scdblfinder.R`.
- `pyproject.toml`: R requirements list scDblFinder and SingleCellExperiment.

## 2026-10-26 (agentic guidelines)

### AGENTIC_GUIDELINES.md [1.0.0]

#### Added

- New guide (prefix `AGT`): how agents scope, verify, report, delegate and hand over work. Proportionality (§1);
  scope and main risks before work (§2); evidence before reasoning: read rather than recall, look up what cannot be
  placed, trace derived material to primary sources (§3); one honest refutation attempt before committing (§4);
  verification against the intended outcome, finished runs are not results, no loosening of criteria (§5); bounded
  universal and negative claims, regenerated numbers, kinds of evidence (§6); questions are not change requests,
  decisions put so they can be answered, corrections recorded at the narrowest scope (§7); unattended batches
  under a written brief (§8); delegation to subagents (§9); confidential and licensed material (§10); handoffs
  between sessions (§11); convergence of revision loops and idea triage (§12); what to record in
  `PROJECT_CONTEXT.md` (§13).
- Adapted in part from the Co-Scientist toolkit (Shocklab, University of Cape Town,
  <https://github.com/shocklab/co-scientist-toolkit>, CC BY 4.0): `guidance/CLAUDE-starter.md`,
  `guidance/paper-CLAUDE.md`, `guidance/handoff-CURRENT.md` and the skills `fable-mode`, `pre-register`,
  `convergence-check`, `decision-page` and `claim-fact-checker`.
- The working-process rules `CODE-AGENT-001` to `CODE-AGENT-007` stay in `CODING_GUIDELINES.md` §25; the new guide
  refers to them and does not repeat them.

### CODING_GUIDELINES.md [3.5.0]

#### Changed

- §2 (documentation set): the `AGENTS.md` template lists `guidelines/AGENTIC_GUIDELINES.md` after this file;
  handoff files in `handoffs/` added to the optional documents.
- `CODE-USE-008`: prefix `AGT` added to the list of guide prefixes.
- §25 (working process): points to `AGENTIC_GUIDELINES.md`. §26 (related guidelines): lists it, to be copied into
  every project developed with agents.

### STATISTICS_REFERENCES.md

- [Si 2025] (ideation-execution gap; arXiv preprint), section L; cited by `AGENTIC_GUIDELINES.md` §12.

### Repository

#### Changed

- `tools/rule_ids.py`: prefix `AGT` and its topics.
- `README.md`: contents table, copy example and licence note with the attribution of the adapted text.

## 2026-10-24 (repository layout)

### CODING_GUIDELINES.md [3.4.0]

#### Changed

- `CODE-USE-001`: in a project, the reference files the guides cite go to `guidelines/references/` (previously next
  to the guides in `guidelines/`); upstream, guides and site profiles are under `guidelines/` and reference files
  under `references/`. §2 (documentation set) and §26 (related guidelines) updated accordingly.

#### Fixed

- `CODE-DL-010` *Why*: line wrapping.

### Repository

#### Changed

- Guides moved to `guidelines/`, site profiles to `guidelines/hpc_sites/`, reference files (`.md`, `.bib`) to
  `references/`; contents unchanged, so their versions do not change. `README.md`, `CHANGELOG.md`,
  `RULE_IDS.tsv`, `tools/` and `recipes/` stay at the root; `RULE_IDS.tsv` still names guides by file name.
- `tools/rule_ids.py`, `tools/check_guides.py`, `tools/check_versions.py` read the new folders;
  `tools/check_versions.py` treats a guide moved without changes as unchanged.
- `README.md`: contents table and *Use in a project* (with a copy example) follow the new layout.

## 2026-10-24 (continuous integration; agent working process)

### CODING_GUIDELINES.md [3.3.0]

#### Changed

- `CODE-AGENT-001`: what counts as an ambiguous request (vague instruction, a term naming more than one thing,
  unstated scope); resolve first what the project answers; one focused question per decision, with options and a
  recommendation; a clear assumption is stated in the reply next to the result it affects, and in the decisions log
  if it outlives the task; no irreversible work under an unconfirmed assumption; *Why*, *Trade-off* and *Deviate
  when* (unattended runs) added.
- `CODE-AGENT-002`: refers to `CODE-AGENT-001`.

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
