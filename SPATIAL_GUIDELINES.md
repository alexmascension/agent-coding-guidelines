# Spatial transcriptomics guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 1.1.0 (2026-10-14)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys such as `[Palla 2022b]` point to
> `OMICS_REFERENCES.md`; keys marked with a dagger, e.g. `[Squair 2021†]`, point to `STATISTICS_REFERENCES.md`.

## 0. How to use this document

* Audience: coding agents and people with introductory statistics. Each rule is followed by its rationale (*Why*, *Trade-off*, *Deviate when*) where it applies.
* Copy this file with `SCRNASEQ_GUIDELINES.md`, `OMICS_REFERENCES.md`, `STATISTICS_GUIDELINES.md` and
  `STATISTICS_REFERENCES.md`. Once a cell- or spot-by-gene matrix exists, QC, normalization, clustering and annotation
  follow `SCRNASEQ_GUIDELINES.md` except where this file says otherwise.
* Platforms covered: 10x **Visium** (spot-based sequencing), 10x **Xenium** and **MERFISH** (Vizgen MERSCOPE)
  (imaging-based, single-cell resolution). Human and mouse.
* Ecosystem: `spatialdata` and `spatialdata-io` for data and readers [Marconato 2024], `squidpy` for spatial
  statistics [Palla 2022b], `scanpy` for expression analysis. R tools are called through scripts
  (`CODING_GUIDELINES.md` §3, languages). Field reviews: [Moses 2022; Rao 2021; Palla 2022a; Moffitt 2022].

## 1. Two kinds of data

| | Visium | Xenium, MERFISH |
|---|---|---|
| Unit of measurement | Spot of 55 µm containing several cells | Individual transcript with coordinates, assigned to segmented cells |
| Genes | Whole transcriptome (probe- or polyA-based) | Targeted panel (hundreds to a few thousand genes) |
| Key processing step | Deconvolution of cell types per spot (§7) | Cell segmentation (§5) |
| Main error sources | Mixed cell types per spot, tissue-dependent library size | Segmentation errors, transcript spillover, off-target signal |

* **For imaging platforms the gene panel is the universe**: genes not in the panel cannot be detected, and enrichment
  analyses MUST use the panel as background (`STATISTICS_GUIDELINES.md` §9, enrichment). Record the panel name and
  version, with the full gene list, in the project.
  * *Why:* an enrichment test against the whole genome compares the panel with genes that could never be measured,
    so any pathway over-represented in the panel design appears enriched.

## 2. Design

* **The experimental unit is the donor or animal.** Spots or cells of one section are not replicates, and several
  sections from one block are technical replicates [Squair 2021†; Zimmerman 2021†].
  * *Why:* spots and cells of one section share donor, tissue processing and slide; treating them as replicates
    turns every donor difference into a significant result. *Trade-off:* spatial studies often have few donors; the
    remedy is more donors or descriptive claims, not counting spots.
* Spread conditions across slides and runs; the slide or run is a batch variable and MUST be recorded.
  * *Why:* slide and run effects (permeabilization, imaging, chemistry lot) are large; if each condition sits on
    its own slide, condition and slide cannot be separated.
* Record tissue preparation (fresh frozen or FFPE), RNA quality (RIN or DV200), section thickness, and time from
  sectioning to run. MERSCOPE results are especially sensitive to sample quality and protocol deviations
  [Wang 2025].
* Where possible, profile matched scRNA-seq or snRNA-seq from the same tissue and species; it serves as reference for
  deconvolution, label transfer and contamination estimates.

## 3. Data structures and coordinates

* Load raw outputs with `spatialdata-io` readers into a `SpatialData` object that keeps images, transcript points,
  segmentation shapes and tables together, with explicit coordinate transformations [Marconato 2024].
* **State the units of every coordinate system** (micrometres or pixels) and the transformation between image and
  transcript coordinates. Use one internal convention and convert at load time (`CODING_GUIDELINES.md` §11,
  coordinate convention).
* For imaging platforms, keep the transcript-level table; it is needed to re-segment and to assess contamination.

## 4. Quality control

### 4.1 Imaging platforms (Xenium, MERFISH)

* Filter transcripts by the platform's quality score with the default threshold unless there is a reason to change it,
  and record the threshold.
* **Use negative controls** (negative-control probes, blank or unassigned codewords) to estimate the false-detection
  rate per section; platforms and runs differ widely in this rate [Ozirmak Lermi 2025].
  * *Why:* every imaging platform reports some signal that is not real; without an estimate of the false-detection
    rate per section, low-expressed genes cannot be interpreted.
* Per cell: transcripts, genes detected, cell area, and transcripts per µm². Thresholds are per section and adaptive
  (`SCRNASEQ_GUIDELINES.md` §4.2, cell-level filters); with panels of a few hundred genes, minimum counts are much lower
  than in scRNA-seq.
* Compare pseudobulk expression of each section with matched scRNA-seq or bulk RNA-seq for the panel genes; poor
  concordance indicates technical problems [Wang 2025; Ozirmak Lermi 2025].
* Practical QC and analysis recommendations for Xenium are benchmarked in [Marco Salas 2025].

### 4.2 Visium

* Keep only spots under tissue; inspect counts and genes per spot on the tissue image before filtering.
* **Do not apply single-cell thresholds blindly.** Library size correlates with tissue structure, so low-count spots
  can correspond to tissue regions with lower RNA content rather than to poor quality [Bhuva 2024]. Check where
  low-count spots lie on the image before removing them.

## 5. Cell segmentation (imaging platforms)

* Treat the platform segmentation as a baseline, not as ground truth. Alternatives: transcript-based Baysor
  [Petukhov 2021] or Proseg [Jones 2025]; image-based Cellpose [Stringer 2020] or Mesmer [Greenwald 2021] on nuclear and
  membrane stains.
  * *Why:* segmentation errors (merged neighbours, truncated cells) assign transcripts to the wrong cell, creating
    apparent co-expression and spurious cell states; the vendor default is optimized for general tissue, not for
    the tissue at hand.
* **Evaluate segmentation** before downstream analysis:
  * fraction of transcripts assigned to cells;
  * cells co-expressing markers of mutually exclusive types (e.g. a neuronal and an oligodendrocyte marker), which
    indicates merged cells or spillover [Ozirmak Lermi 2025];
  * cell size distribution against expected morphology.
* The expansion distance around nuclei trades sensitivity against contamination from neighbours; choose it
  deliberately and record it [Marco Salas 2025].
* For Xenium, transcript spillover between neighbouring cells can be estimated and corrected with snRNA-seq-informed
  methods (SPLIT) [Bilous 2026].
* Record the segmentation method, version, input stains and parameters; every downstream result depends on them.

## 6. Normalization

* **Library size in spatial data partly reflects biology**: it correlates with tissue structure, and standard
  scRNA-seq library-size normalization can impair the detection of spatial domains [Bhuva 2024].
  * *Trade-off:* not normalizing leaves technical depth differences between spots; normalizing can remove real
    differences between regions. Comparing both is the only way to know which dominates in a dataset.
* Keep raw counts in `layers["counts"]`. When domains or clusters are a goal, compare results with and without
  library-size normalization and report which was used and why.
* For imaging panels, total counts depend on panel composition; normalization by cell area MAY be compared with
  normalization by total counts.

## 7. Cell typing and deconvolution

* **Imaging platforms:** cluster and annotate as in `SCRNASEQ_GUIDELINES.md` §8–9 (clustering and annotation), or
  transfer labels from a matched single-cell reference. The panel limits the achievable resolution; MUST NOT report
  cell subtypes that the panel cannot distinguish.
* **Visium deconvolution:** use cell2location [Kleshchevnikov 2022] or RCTD [Cable 2021], the top performers in
  independent benchmarks [Li 2022b; Sang-aram 2024]. A simple regression baseline outperforms many dedicated methods
  [Sang-aram 2024], so a new method is justified only by a benchmark on comparable data.
  * *Why:* each 55 µm spot contains several cells, so expression mixes cell types; deconvolution separates them.
    General rules on references, scale and validation are in `DECONVOLUTION_GUIDELINES.md`.
* The reference MUST contain all cell types expected in the tissue and come from the same species and, ideally, the
  same tissue and condition.
* Deconvolution outputs are estimated abundances; report them with the reference used, and do not treat them as
  observed counts.

## 8. Spatially variable genes

* Methods: SpatialDE [Svensson 2018], SPARK [Sun 2020], nnSVG [Weber 2023]; Moran's I in `squidpy` as a fast screen.
* **Many spatially variable genes only reflect where cell types are.** Interpret them together with cell-type
  composition, or test them conditional on cell type.
  * *Why:* any marker of a cell type that is unevenly distributed is spatially variable; the interesting genes are
    those whose pattern is not explained by composition.
* With thousands of spots or cells, p-values are tiny for weak patterns; rank genes by effect size or statistic, and
  apply multiple-testing correction (`STATISTICS_GUIDELINES.md` §7–8, multiple testing and effect sizes).

## 9. Spatial domains

* **No method is best in all settings** [Yuan 2024]. Choose by data type and goal: BayesSpace for Visium
  [Zhao 2021a], graph-based methods such as STAGATE [Dong 2022] or SpaGCN [Hu 2021], or neighbourhood-aware clustering
  (e.g. Banksy).
  * *Trade-off:* without a single best method, results should be checked with a second method or across parameters;
    domains that depend on the method are not robust findings.
* Validate domains against histology or known anatomy, and check their stability across parameters and sections.

## 10. Neighbourhoods, co-localization and comparisons between conditions

* Neighbourhood enrichment and co-occurrence statistics (`squidpy`) use permutations of cell labels within a section.
  The permutation null assumes that labels are exchangeable inside the tissue; strong tissue structure makes many
  pairs "significant".
* **Inference across conditions uses sections or donors as units**: compute the statistic per section, then compare
  per-donor values between conditions (`STATISTICS_GUIDELINES.md` §4, hierarchical data).
  * *Why:* a neighbourhood statistic computed on pooled cells from several donors reflects the donors with most
    cells; per-section statistics weight donors correctly and give a valid between-condition test.
* **Differential expression between conditions:** pseudobulk per donor and cell type (imaging) or per donor and domain
  (Visium), then DE as in `SCRNASEQ_GUIDELINES.md` §10 (condition DE). Cells or spots MUST NOT be used as replicates.
* Neighbouring cells and spots are correlated (spatial autocorrelation), which violates the independence assumption
  of standard tests even within one section.

## 11. Several sections

* Each section keeps its own coordinate system; integration of expression across sections (§7 of
  `SCRNASEQ_GUIDELINES.md`, integration) is for clustering and annotation only.
* Batch variables for spatial data: slide, run, section, and processing date.

## 12. What to report

In addition to `STATISTICS_GUIDELINES.md` §12 (what to report):

1. Platform, chemistry, panel name and version (with gene list), instrument software version.
2. Tissue preparation, RNA quality, section thickness, number of sections per donor.
3. Transcript quality threshold, negative-control false-detection rate per section.
4. Segmentation method, stains, parameters, expansion distance, and evaluation results.
5. QC thresholds and cells or spots retained per section.
6. Normalization choice, deconvolution or label-transfer reference, and methods for spatially variable genes,
   domains and neighbourhood statistics, with parameters.

## 13. Python implementation

| Task | Python | R reference (via script) |
|---|---|---|
| Reading platform outputs | `spatialdata-io` | — |
| Data structure | `spatialdata` | `SpatialExperiment` [Righelli 2022] |
| Segmentation pipelines | `sopa` (Baysor, Proseg, Cellpose), `cellpose` | — |
| Expression QC, clustering | `scanpy` | — |
| Spatial statistics (Moran's I, neighbourhoods) | `squidpy` | — |
| Deconvolution | `cell2location` | `spacexr` (RCTD) |
| Spatially variable genes | `squidpy`, `spatialde` | `nnSVG`, `SPARK` |
| Spatial domains | `banksy-py` | `BayesSpace` |
| Label transfer, integration | `scvi-tools`, `tangram-sc` | — |
