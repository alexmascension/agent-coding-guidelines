# Recipes: reference implementations

Tested implementations of non-obvious procedures of the guidelines (`CODING_GUIDELINES.md` §26, recipes). Tests use
synthetic data: they check the contract of the code (inputs, outputs, edge cases, planted signals), not the
validity of a method on real data.

Status: **implemented** (module and tests here), **external** (run with the cited tool; the guide describes how),
**planned** (not yet available).

## Coverage by guideline

| Guideline | Procedure | Status | Module | Test scope and limits |
|---|---|---|---|---|
| SCRNASEQ §1.1 (`SC-DATA-004`) | Input layer per step (QC, doublets: `counts`; normalization, pseudobulk: `counts_ambient_corrected` if it exists, else `counts`); integer check for count-based DE; count layers protected from derived output | implemented | `scrnaseq/layers.py` | Objects with deliberately different layers: defaults follow the table, explicit layers honoured; a non-integer corrected layer refused with its name and the remedies (explicit rounding recorded, or `layer="counts"`); derived layers refused on count layers |
| SCRNASEQ §1.3 | Annotation from GTF, stable IDs, marker resolution | implemented | `common/gene_ids.py` | Ensembl/GENCODE versions, PAR suffix, duplicated and missing symbols |
| SCRNASEQ §3 | STARsolo mapping and counting | external | — | — |
| SCRNASEQ §3 | Library QC from STARsolo `Summary.csv` and `Barcodes.stats` (nf-core/scrnaseq layout), comparison across libraries, `uns["library_qc"]` | implemented | `scrnaseq/library_qc.py` | Feature names in rows (`GeneFull_Ex50pAS`), `NoMulti` kept missing, malformed lines, whitelist mismatch, cells vs loaded; MAD = 0 and fewer than 3 libraries not evaluated; nullable flags written to h5ad. Thresholds across libraries are a project default |
| SCRNASEQ §4.1 | CellBender, EmptyDrops_CR | external | — | — |
| SCRNASEQ §4.1–4.2 | Mito by chromosome, nuclear fraction, per-capture MAD flags; mitochondrial fraction on the logit scale with a project floor; top-20-gene fraction; flagged fraction per cluster; minimum MAD for fraction metrics; ribosomal and haemoglobin fractions from the default lists | implemented | `scrnaseq/qc.py` | minimum MAD applied and recorded (nuclei with MAD 0), not applied to log metrics, non-finite values, absolute mito limit, empty ribo/hb group missing (not 0), panels without mito genes; low mito not flagged, cells beyond the MAD threshold but below the floor kept; fractions outside [0, 1] (percentages) refused before the logit; top-20 fraction against a direct computation (sparse and dense), both tails; cluster that loses most cells marked. Inspection thresholds per cluster are a project default |
| SCRNASEQ §4.2 (`SC-QC-006`) | Default ribosomal and haemoglobin lists (human, mouse; Ensembl IDs) and their builder from HGNC gene groups and HGNC-curated MGI orthologues | implemented | `common/gene_groups.py`, `common/build_gene_groups.py`, `common/data/gene_groups/` | Shipped lists without RPS6K kinases, MRPL/MRPS, HBEGF, HBS1L; metadata present; malformed and duplicated files rejected. Builder, offline with hand-written responses: atomic writes, no overwrite without `--overwrite`, byte-for-byte source snapshot with checksums and identical rebuild with `--from-snapshot`. A networked build is not run in CI |
| SCRNASEQ §4.3 | Scrublet, scDblFinder (R) and optional DoubletDetection; consensus/union; Scrublet thresholds; cluster review; genotype doublets | implemented | `scrnaseq/doublets.py`, `scrnaseq/scdblfinder.R` | Planted doublets; mixed singlets vs true doublets; capture-adjusted rates; input layer refused when SoupX- or decontX-corrected (also by the direct scDblFinder wrapper), warning when such a layer exists before the review; failed automatic Scrublet threshold raised or marked not evaluated; manual thresholds; scDblFinder test skipped without R. The AND/OR combination is a project default |
| SCRNASEQ §5 | Pooled size factors (scranPY; R scran), median factors, log1p | implemented | `scrnaseq/normalization.py`, `scrnaseq/scran_size_factors.R` | Agreement scranPY–R scran on simulated composition differences; equivalence with `normalize_total(target_sum=None)`; invalid factors; default input layer per `SC-DATA-004`; output refused on count layers and on the input layer |
| SCRNASEQ §6 | Number of PCs and robustness to p ± 5 | implemented | `scrnaseq/pc_robustness.py` | Separated clusters vs pure noise; no truncation of the proposal; criteria evaluated recorded, warning when Gavish–Donoho is not evaluated; degenerate spectra and unknown clusters of interest refused |
| SCRNASEQ §7 | Harmony, scVI/scANVI | external | — | — |
| SCRNASEQ §8 | Resolution by subsampling stability | implemented | `scrnaseq/clustering_stability.py` | Separated clusters; does not test stability across donors |
| SCRNASEQ §9 | Marker-dictionary annotation (decoupler ULM) | implemented | `scrnaseq/marker_annotation.py` | Planted markers; ties, unscored cells, index alignment, contradictory signs |
| SCRNASEQ §9 | CellAssign, CellTypist | external | — | — |
| SCRNASEQ §10 | Pseudobulk aggregation and DE per cell type; pseudo-replicate guard; interval, rank and quartile of significance; detectability of genes of interest (equivalence test); global-null diagnostic; minimum-cells diagnostic; units, cells and UMIs per cell type | implemented | `scrnaseq/pseudobulk.py`, `scrnaseq/units.py` | Planted DE gene; numeric covariates kept numeric; absent combinations reported; samples split across units rejected, pooled captures accepted with a donor key; missing unit IDs and IDs that collide as strings (`1`, `"1"`) rejected; no gene passing the filter reported not testable; no finite p-value leaves the global adjustment missing; quartiles only for significant genes; changed, no-change-of-interest, inconclusive and not-tested statuses; few discoveries on null controls; drift of low-cell units. The null diagnostic does not establish FDR control |
| SCRNASEQ §10 | Cell types compared at matched depth (cells per unit, binomial thinning of UMIs) | implemented | `scrnaseq/depth_matching.py` | Cells per unit and mean UMIs per cell matched; counts integer and not above the originals; reproducible; non-integer counts rejected. Does not match composition, ambient RNA or capture efficiency |
| SCRNASEQ §2 | Empirical power of pseudobulk DE over donors × cells, against a held-out reference | implemented | `scrnaseq/empirical_power.py` | Large planted effects recovered; infeasible designs and splits without a testable reference marked. Recall and precision measure agreement with an imperfect reference, not true power or precision: recall tends to be optimistic and precision pessimistic, but neither is a guaranteed bound; every design is a full DE run |
| STATISTICS §10 | Held-out reference splits (disjoint, stratified, seeded) | implemented | `stats/heldout.py` | Disjointness, proportions per stratum, reproducibility, invalid input |
| SCRNASEQ §11 | Count table with explicit zeros (input of DA methods) | implemented | `scrnaseq/abundance.py` | Zeros, manifest with samples that lost all cells; a cell-type list that omits observed types refused unless declared a sub-composition (denominator recorded); manifest values that contradict `obs` refused; imports without pydeseq2. The DA tests (propeller, scCODA, sccomp, Milo) are external |
| SCRNASEQ §12–13 | Cell-cell communication, trajectories, velocity | external | — | — |
| RNASEQ_BULK §4 | Gene tables (nf-core) and tximport `lengthScaledTPM` | implemented | `bulk/counts_io.py`, `bulk/tximport_lengthscaled.R` | Version stripping; library-size preservation in R. No test of isoform-switch length correction yet |
| RNASEQ_BULK §5, §8 | Design-based filter, pydeseq2 with standard table, VST | implemented | `bulk/deseq2.py` | Planted DE genes, threshold test, unused categories, exact coefficients, count validation (exact integers, int64 range, no genes). Agreement with R DESeq2: see tests (skipped without R) |
| RNASEQ_BULK §5 | Gene-length bias diagnostics | implemented | `bulk/length_bias.py` | Planted slope vs none. Diagnostic only; correction (cqn, EDASeq) is external |
| RNASEQ_BULK §7 | Sensitivity to k hidden factors; residual PCs (exploratory) | implemented | `bulk/hidden_factors.py` | Recovery of an unconfounded factor; rank limits; factors actually used. Not qSVA/RUV/SVA (external) |
| RNASEQ_BULK §8 | Global-null permutation diagnostic | implemented | `bulk/fdr_calibration.py` | Null data, invalid p-values, Monte Carlo bound. Does not establish FDR control under alternatives |
| STATISTICS §10 | Permutation within strata, $p=(b+1)/(m+1)$, binomial bound | implemented | `stats/permutation.py` | Strata preserved, same permutations for a seed across `PYTHONHASHSEED` values, non-finite statistics rejected. Two-sided p doubles the smaller one-sided p (SciPy convention): comparing absolute values of $T$ is also valid under exchangeability but matches a two-sided alternative only when the null of $T$ is centred at 0 and roughly symmetric |
| SEX_INFERENCE | Sex calls, homolog ratio, in-silico mixtures | planned | — | — |
| SAMPLE_IDENTITY | Count hashing, correlation outliers, somalier parsing | planned | — | somalier, NGSCheckMate, Picard are external |
| DECONVOLUTION | Held-out pseudobulk validation, NNLS baseline | planned | — | Bisque, MuSiC, BayesPrism, hspe are external |
| META_ANALYSIS | Per-gene random effects (PM + HKSJ), prediction interval | planned | — | — |
| SPATIAL | Per-section statistics to donor-level test, negative-control FDR | planned | — | — |
| HPC §1–2, §4 | Role directories and links, staging, promotion with manifest; job template | implemented | `hpc/setup_storage.sh`, `hpc/stage.sh`, `hpc/promote.sh`, `hpc/job.slurm`, `hpc/sites/*.env` | Temporary directories as role roots: links, no replacement of existing directories, local site, checksummed read-only promotion, refusal to overwrite or to copy symbolic links; staging verified by checksum (same size and modification time, different content), path containment for staging and promotion, no `USER` needed, promotion manifest computed on the source and readable like the run, `--dry-run` without side effects, a differing staged replica replaced only with `--overwrite`, symbolic links and empty sources refused. No scheduler in tests; `job.slurm` is only linted |

## Installation and tests

```bash
pip install "agent-coding-guidelines-recipes @ git+https://github.com/alexmascension/agent-coding-guidelines.git@<commit>"
pip install -e ".[dev]"     # from a checkout, for development
pytest                        # R-dependent tests are skipped unless RECIPES_RSCRIPT or Rscript is available
ruff check recipes && ruff format --check recipes
```

R scripts (`*.R`) are shipped inside the package and located with `recipes.common.resources.package_script`. They
need R with scran, tximport, Matrix and optparse (and DESeq2 for the concordance test).

`scranPY` is installed from GitHub at a pinned commit; its distribution name collides with the Bioconductor
`scranpy`, so the two cannot be installed in the same environment.
