# Recipes: reference implementations

Tested implementations of non-obvious procedures of the guidelines (`CODING_GUIDELINES.md` §26, recipes). Tests use
synthetic data: they check the contract of the code (inputs, outputs, edge cases, planted signals), not the
validity of a method on real data.

Status: **implemented** (module and tests here), **external** (run with the cited tool; the guide describes how),
**planned** (not yet available).

## Coverage by guideline

| Guideline | Procedure | Status | Module | Test scope and limits |
|---|---|---|---|---|
| SCRNASEQ §1.3 | Annotation from GTF, stable IDs, marker resolution | implemented | `common/gene_ids.py` | Ensembl/GENCODE versions, PAR suffix, duplicated and missing symbols |
| SCRNASEQ §3 | STARsolo mapping and counting | external | — | — |
| SCRNASEQ §4.1 | CellBender, EmptyDrops_CR | external | — | — |
| SCRNASEQ §4.1–4.2 | Mito by chromosome, nuclear fraction, per-capture MAD flags | implemented | `scrnaseq/qc.py` | MAD = 0, non-finite values, absolute mito limit, panels without mito genes |
| SCRNASEQ §4.3 | Scrublet, scDblFinder (R) and optional DoubletDetection; consensus/union; Scrublet thresholds; cluster review; genotype doublets | implemented | `scrnaseq/doublets.py`, `scrnaseq/scdblfinder.R` | Planted doublets; mixed singlets vs true doublets; capture-adjusted rates; input layer refused when SoupX- or decontX-corrected, warning when such a layer exists before the review; failed automatic Scrublet threshold raised or marked not evaluated; manual thresholds; scDblFinder test skipped without R. The AND/OR combination is a project default |
| SCRNASEQ §5 | Pooled size factors (scranPY; R scran), median factors, log1p | implemented | `scrnaseq/normalization.py`, `scrnaseq/scran_size_factors.R` | Agreement scranPY–R scran on simulated composition differences; equivalence with `normalize_total(target_sum=None)`; invalid factors |
| SCRNASEQ §6 | Number of PCs and robustness to p ± 5 | implemented | `scrnaseq/pc_robustness.py` | Separated clusters vs pure noise; no truncation of the proposal |
| SCRNASEQ §7 | Harmony, scVI/scANVI | external | — | — |
| SCRNASEQ §8 | Resolution by subsampling stability | implemented | `scrnaseq/clustering_stability.py` | Separated clusters; does not test stability across donors |
| SCRNASEQ §9 | Marker-dictionary annotation (decoupler ULM) | implemented | `scrnaseq/marker_annotation.py` | Planted markers; ties, unscored cells, index alignment, contradictory signs |
| SCRNASEQ §9 | CellAssign, CellTypist | external | — | — |
| SCRNASEQ §10 | Pseudobulk aggregation and DE per cell type | implemented | `scrnaseq/pseudobulk.py` | Planted DE gene; numeric covariates kept numeric; absent combinations reported |
| SCRNASEQ §11 | Count table with explicit zeros (input of DA methods) | implemented | `scrnaseq/abundance.py` | Zeros, manifest with samples that lost all cells. The DA tests (propeller, scCODA, sccomp, Milo) are external |
| SCRNASEQ §12–13 | Cell-cell communication, trajectories, velocity | external | — | — |
| RNASEQ_BULK §4 | Gene tables (nf-core) and tximport `lengthScaledTPM` | implemented | `bulk/counts_io.py`, `bulk/tximport_lengthscaled.R` | Version stripping; library-size preservation in R. No test of isoform-switch length correction yet |
| RNASEQ_BULK §5, §8 | Design-based filter, pydeseq2 with standard table, VST | implemented | `bulk/deseq2.py` | Planted DE genes, threshold test, unused categories, exact coefficients, count validation. Agreement with R DESeq2: see tests (skipped without R) |
| RNASEQ_BULK §5 | Gene-length bias diagnostics | implemented | `bulk/length_bias.py` | Planted slope vs none. Diagnostic only; correction (cqn, EDASeq) is external |
| RNASEQ_BULK §7 | Sensitivity to k hidden factors; residual PCs (exploratory) | implemented | `bulk/hidden_factors.py` | Recovery of an unconfounded factor; rank limits; factors actually used. Not qSVA/RUV/SVA (external) |
| RNASEQ_BULK §8 | Global-null permutation diagnostic | implemented | `bulk/fdr_calibration.py` | Null data, invalid p-values, Monte Carlo bound. Does not establish FDR control under alternatives |
| STATISTICS §10 | Permutation within strata, $p=(b+1)/(m+1)$, binomial bound | implemented | `stats/permutation.py` | Strata preserved, non-finite statistics rejected, SciPy two-sided convention |
| SEX_INFERENCE | Sex calls, homolog ratio, in-silico mixtures | planned | — | — |
| SAMPLE_IDENTITY | Count hashing, correlation outliers, somalier parsing | planned | — | somalier, NGSCheckMate, Picard are external |
| DECONVOLUTION | Held-out pseudobulk validation, NNLS baseline | planned | — | Bisque, MuSiC, BayesPrism, hspe are external |
| META_ANALYSIS | Per-gene random effects (PM + HKSJ), prediction interval | planned | — | — |
| SPATIAL | Per-section statistics to donor-level test, negative-control FDR | planned | — | — |
| HPC §1–2, §4 | Role directories and links, staging, promotion with manifest; job template | implemented | `hpc/setup_storage.sh`, `hpc/stage.sh`, `hpc/promote.sh`, `hpc/job.slurm`, `hpc/sites/*.env` | Temporary directories as role roots: links, no replacement of existing directories, local site, checksummed read-only promotion, refusal to overwrite or to copy symbolic links. No scheduler in tests; `job.slurm` is only linted |

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
