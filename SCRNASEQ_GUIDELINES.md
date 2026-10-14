# Single-cell and single-nucleus RNA-seq guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.1.0 (2026-10-14)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys such as `[Heumos 2023]` point to
> `OMICS_REFERENCES.md`; keys with a dagger, e.g. `[Squair 2021†]`, point to `STATISTICS_REFERENCES.md`. Each rule is
> followed by ***Why*** (mechanism), ***Trade-off*** (cost) and ***Deviate when*** (justified exceptions) where they
> apply.

## 0. How to use this document

* Audience: coding agents and people with introductory statistics.
* Copy with `OMICS_REFERENCES.md`, `STATISTICS_GUIDELINES.md`, `STATISTICS_REFERENCES.md` and
  `RNASEQ_BULK_GUIDELINES.md` (pseudobulk differential expression, DE, follows its §7–8). Sex checks follow `SEX_INFERENCE_GUIDELINES.md`;
  identity and pooled designs follow `SAMPLE_IDENTITY_GUIDELINES.md`.
* Scope: droplet-based scRNA-seq (cells) and snRNA-seq (nuclei), human and mouse, from FASTQ to annotation,
  differential expression, abundance, communication and trajectories. Multimodal assays are out of scope.
* Ecosystem: scverse (`anndata`, `scanpy`, `scvi-tools`, `pertpy`, `decoupler`, `cellrank`, `liana`) [Virshup 2023;
  Wolf 2018; Gayoso 2022]. R tools are called through scripts (`CODING_GUIDELINES.md` §3, languages). Seurat
  [Hao 2023] is the main R alternative and is not used unless a project requires it. Workflow references:
  [Heumos 2023; Luecken 2019; Amezquita 2019].

## 1. Data conventions

### 1.1 Layers, not `.raw`

* Store each representation in a named layer: `layers["counts"]` (raw), `layers["counts_ambient_corrected"]`,
  `layers["spliced"]`, `layers["unspliced"]`, `layers["ambiguous"]`, `layers["log1p_norm"]`. `.X` holds one of them,
  and `adata.uns["X_content"]` states which.
* **Do not use `adata.raw`.**
  * *Why:* `.raw` is a frozen snapshot whose content differs between workflows (raw counts in some, log-normalized
    values in most), and several scanpy functions read it by default, so the same code can silently operate on
    different data. Named layers make the content explicit.
  * *Deviate when:* a tool requires `.raw`; fill it immediately before the call with a named layer and record which.
* Raw counts are never overwritten.

### 1.2 Observation metadata

* `adata.obs` MUST contain for every cell: `sample_id` (one capture), `donor_id` (the experimental unit), `condition`,
  and each technical variable (`batch`, `chemistry`, `capture_date`). Names are declared once
  (`CODING_GUIDELINES.md` §6, controlled vocabulary).

### 1.3 Gene identifiers and symbols

* **Index genes by the stable Ensembl ID without version** (`ENSG00000141510`), and store the versioned ID
  (`ENSG00000141510.18`) in `var["gene_id_version"]` and the symbol in `var["gene_symbol"]`, both from the same
  annotation release.
  * *Why:* symbols change between releases and are not unique. The stable ID identifies the gene across releases,
    so tables from different projects can be joined; the version records which gene model produced the counts.
* **Write marker dictionaries in symbols** (what people read and publish), and **resolve them to IDs when loading**
  with one project function that maps through the annotation table and reports:
  * symbols with no ID in this release (aliases or retired names): listed, never dropped silently;
  * symbols with several IDs (e.g. copies in the pseudoautosomal region (PAR), annotation duplicates): resolved by a fixed rule, stated once in
    `PROJECT_CONTEXT.md` (keep all IDs, or keep the most expressed).
* **Display symbols in plots and tables** through the same mapping; keep the ID in result tables next to the symbol.
* Cross-species marker lists use one-to-one orthologs from a stated Ensembl release; human symbols are not "converted"
  to mouse by changing case.
* Recipe: `recipes/common/gene_ids.py` (annotation table from the GTF, stable IDs, marker dictionary resolution
  with a report of unmatched, case-mismatched and multi-ID symbols).

## 2. Design

* **The experimental unit is the donor or animal, not the cell** [Squair 2021†; Zimmerman 2021†].
  * *Why:* cells of one donor share genotype, environment and processing; treating them as independent replicates
    makes every donor-level difference look significant. Thousands of cells from three mice are $n = 3$.
* Plan **at least 3, preferably 5 or more, donors per condition**.
* **Do not confound capture with condition**: multiplex samples (genetic or hashtag demultiplexing,
  `SAMPLE_IDENTITY_GUIDELINES.md` §6) or spread conditions across capture days.
* **Single-donor studies** support only descriptive conclusions about that donor:
  * comparisons between conditions within the donor (regions, treatments in vitro) use independent captures or
    dissociations as replicates; they measure technical and sampling variation, not biological variation;
  * pseudobulk per capture MAY be tested, and results are reported as valid for that donor only;
  * cell-level p-values are not evidence in any design.
* **Several samples per donor and several donors**: pseudobulk per sample with the donor in the model, i.e. paired
  design `~ donor + condition` when the condition varies within donors; mixed model or `limma::duplicateCorrelation`
  when it varies between donors; or aggregate to the donor.
* **Cells or nuclei**: nuclei suit frozen and hard-to-dissociate tissue (brain, spinal cord) and avoid dissociation
  artefacts; enzymatic dissociation of brain induces ex vivo activation signatures in microglia and other glia
  [Marsh 2022].
  * *Trade-off:* nuclei lose cytoplasmic transcripts and have more ambient RNA.
  * *Deviate when:* fresh tissue with a validated cold dissociation protocol with transcription and translation
    inhibitors [Marsh 2022].

## 3. Quantification

* **Default: STARsolo** [Kaminow 2021], gene-level counts, with these features in one run:
  `--soloFeatures Gene GeneFull Velocyto` (add `SJ` if splice junctions are needed).
  * `Gene` (exonic reads) is the main matrix for cells; `GeneFull` (exons and introns) is the main matrix for nuclei,
    where most reads are unspliced; `Velocyto` gives spliced, unspliced and ambiguous counts, stored as
    `layers["spliced"]`, `layers["unspliced"]`, `layers["ambiguous"]`.
  * To match Cell Ranger 4 or later, add `--clipAdapterType CellRanger4 --outFilterScoreMin 30
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR`
    (STARsolo documentation).
  * Keep the unique-gene matrix as the main one. Multi-gene reads (overlapping genes, paralog families) MAY be
    recovered with `--soloMultiMappers EM` as an additional matrix, for questions about those gene classes.
  * Cell calling: `--soloCellFilter EmptyDrops_CR` (EmptyDrops in the Cell Ranger flavour [Lun 2019]), combined with
    CellBender as described in §4.1.
  * *Why STARsolo:* it reproduces Cell Ranger results, is faster, was more accurate than pseudoalignment tools on
    realistic simulations, and handles multi-gene reads [Kaminow 2021]; in an independent comparison STARsolo, Cell
    Ranger and alevin-fry gave similar, reliable gene sets [Brüning 2022]. It also produces genome BAMs, needed for
    genotype checks (`SAMPLE_IDENTITY_GUIDELINES.md` §4).
  * *Why gene level:* in 3' or 5' tag protocols, reads concentrate at one transcript end, so most isoforms cannot be
    distinguished, and isoform analysis at single-cell resolution remains an open question [Arzalluz-Luque 2018].
  * *Why spliced/unspliced layers:* they give the nuclear fraction per cell for quality control (QC; §4.1) and RNA velocity without
    re-quantification (§13).
  * *Trade-off:* spliced/unspliced estimates differ between tools and feature definitions, and the differences change
    velocity results [Soneson 2021]; the STARsolo documentation also notes that `Velocyto` quantification is under
    active development. Fix one STAR version per project.
  * *Deviate when:* alevin-fry/simpleaf [He 2022] or kallisto|bustools in `nac` mode [Melsted 2021; Sullivan 2024] are
    acceptable alternatives when no BAM is needed.
* **Isoform-level questions** need full-length (e.g. Smart-seq) or long-read protocols [Lebrigand 2020], or analysis of
  transcript compatibility counts [Ntranos 2019]. For full-length per-cell isoform quantification, use Salmon rather
  than alevin. Isoforms are not quantified from 3' tag data.
* Record STAR version, all `--solo*` parameters, chemistry and whitelist, reference genome and annotation release.
* Default references can miss reads from genes with incomplete 3' annotation; optimized references recover them
  [Pool 2023]. Check the annotation of a gene of interest before concluding it is not expressed.

## 4. Quality control

All QC is done **per capture**, before merging or integration. Record, per capture, the number of barcodes removed at
each step.

### 4.1 Cell calling and ambient RNA

* **CellBender is the default** for ambient RNA removal and cell calling on the raw (unfiltered) STARsolo matrix
  [Fleming 2023]. Keep barcodes with cell probability above 0.5, store the corrected counts in
  `layers["counts_ambient_corrected"]`, and compare the calls with STARsolo `EmptyDrops_CR`; record barcodes called by
  only one of the two.
  * *Why:* it models ambient RNA and cell-containing droplets jointly. Ambient transcripts make every cell appear to
    express markers of abundant types (e.g. neuronal transcripts in glia in nuclei data).
  * *Trade-off:* needs a GPU to be practical.
  * *Deviate when:* without a GPU, SoupX (R) [Young 2020]; for statistical cell calling alone, EmptyDrops (R)
    [Lun 2019] or dropkick (Python) [Heiser 2021].
* **Nuclear fraction** per barcode is `unspliced / (spliced + unspliced)` from the STARsolo `Velocyto` layers.
  * *Why:* in scRNA-seq, barcodes with very low nuclear fraction are often cytoplasmic fragments without a nucleus,
    frequent even in reference atlases [Montserrat-Ayuso 2024]; damaged cells and empty droplets have characteristic
    nuclear fractions [Muskovic 2021]. In nuclei data, a low nuclear fraction indicates cytoplasmic contamination.

### 4.2 Cell-level filters

* **Thresholds are adaptive and per capture**: remove cells more than 3–5 median absolute deviations from the capture
  median in log total counts, log genes detected, and mitochondrial fraction, plus the nuclear-fraction filter above.
  * *Why:* depth and quality differ between captures; one fixed threshold removes good cells in one capture and keeps
    bad cells in another.
  * *Trade-off:* MAD (median absolute deviation) filters assume most cells are good; inspect distributions when a capture is mostly poor.
* **Mitochondrial genes are identified by chromosome** (`MT` in Ensembl, `chrM` in UCSC) from the GTF (gene transfer format) of the reference
  used, not by symbol prefix.
  * *Why:* with Ensembl IDs as index, prefixes depend on a symbol column that may be incomplete; a regex without the
    hyphen (`^MT`) also matches nuclear genes such as *MTOR* or *MTHFR*; capitalization differs between species.
* **Mitochondrial thresholds depend on species, tissue and assay**: the common 5% cut-off is often unsuitable for
  human tissues [Osorio 2020]; model-based alternatives exist [Hippen 2021]. In nuclei, mitochondrial reads should be
  near zero.
* Ribosomal and haemoglobin genes, which have no dedicated chromosome, are defined by explicit species-specific lists
  stored in the project.
* MUST NOT filter cells by the expression of markers of the hypothesis being tested.

* Recipe for §4.1–4.2: `recipes/scrnaseq/qc.py` (gene groups from the GTF chromosome, nuclear fraction, per-capture
  MAD flags with the thresholds used per capture).

### 4.3 Doublets

* **Detect doublets per capture, before integration, with two methods: Scrublet (`scanpy.pp.scrublet`)
  [Wolock 2019] and DoubletDetection** [Gayoso 2020]. Store both scores and calls, and derive `doublet_consensus` (both methods agree)
  and `doublet_any` (either method).
  * *Why two methods:* single methods are unstable across datasets, and combining methods improves accuracy and
    stability [Xiong 2022; Neavin 2024]. In a benchmark of nine methods, DoubletFinder (R) had the best accuracy
    [Xi 2021]; scDblFinder (R) ranked first in a later independent benchmark [Germain 2021]; they are the R
    alternatives. Solo (`scvi-tools`) is a further Python option [Bernstein 2020].
* **Flag first, remove later**:
  1. Keep all cells through a first clustering.
  2. Remove cells with `doublet_consensus`.
  3. Remove clusters enriched in `doublet_any` that co-express markers of two lineages; review borderline clusters
     manually.
  4. Run the final clustering and all downstream analyses without doublets.
  * *Why:* removing true singlets of rare types costs more than leaving isolated doublets, which the cluster-level
    evidence captures. Doublets MUST be removed before final clustering, DE and abundance analysis because they create
    spurious intermediate populations [Xi 2021].
  * *Note:* no benchmark compares union and intersection of these two methods; this rule is a precision-first default.

* Recipe: `recipes/scrnaseq/doublets.py` (both detectors per capture, consensus and union calls, per-cluster review
  table with an enrichment test and a lineage co-expression level).

## 5. Normalization and feature selection

* **Size factors from pooling (scran deconvolution) followed by `log1p`** [Lun 2016]; alternative:
  `scanpy.pp.normalize_total(target_sum=None)` (median total counts) followed by `log1p`.
  * In Python, pooled factors come from `scranPY` (`scranPY.compute_sum_factors`, a port of
    `scran::computeSumFactors`) [Fortmann 2023], installed from GitHub pinned to a commit. It is not the Bioconductor
    `scranpy` package (libscran bindings, which has no pooled size factors); the two share a distribution name and
    cannot be installed together. `scranPY` needs a dense matrix; when it does not fit in memory, use R scran through
    the recipe script. Pass clusters (e.g. a first Leiden clustering) so pools contain similar cells.
  * Recipe: `recipes/scrnaseq/normalization.py` (scranPY, R scran cross-check, median factors, `log1p`).
  * *Why:* "normalize every cell to 10,000 counts" is itself size-factor normalization with
    $s_i = \text{total}_i / 10^4$. The value 10,000 is arbitrary and changes the effective pseudocount of `log1p`; a
    shifted logarithm performs as well as more complex transformations when the reference total is reasonable
    (mean or median) [Ahlmann-Eltze 2023†]. Pooled size factors are more robust to many zeros and to differences in
    RNA content between cell types [Lun 2016]. Bulk normalization methods applied unchanged to single cells can
    mislead [Vallejos 2017].
  * *Deviate when:* analytic Pearson residuals (`scanpy.experimental.pp`) MAY be used for gene selection and PCA (principal component analysis); they
    outperform other transformations for these tasks [Lause 2021]. Regularized negative binomial residuals
    (sctransform) [Hafemeister 2019] and deviance-based gene selection [Townes 2019] are related alternatives.
* **Normalized values are for exploration, clustering and visualization**; condition DE uses raw counts (§10).
* **Select about 2,000 highly variable genes in a batch-aware way** (per batch, then combined). Feature selection affects
  integration and reference mapping, and standard highly variable gene selection is an effective default
  [Zappia 2025; Luecken 2021].

## 6. Principal components

1. Compute 50 PCs on the selected genes.
2. Propose a number $p$ with two criteria: the elbow of the explained variance (maximum distance to the chord between
   the first and last point of the curve), and the Gavish–Donoho hard threshold on the singular values
   [Gavish 2014]. If they disagree, take the larger.
3. **Check robustness over $p \pm 5$**:
   1. Build the neighbour graph and Leiden clusters at the working resolution for $p$ (the reference clustering).
   2. Repeat for every $n$ from $p-5$ to $p+5$ with the same number of neighbours, resolution and seed (use steps of 2
      when one clustering takes more than a few minutes).
   3. Compare each clustering with the reference: adjusted Rand index (ARI) overall, and, for each reference cluster,
      the Jaccard index with its best-matching cluster.
   4. $p$ is robust if ARI ≥ 0.9 for every $n$ and every reference cluster of interest keeps Jaccard ≥ 0.75. Otherwise
      increase $p$ (preferred) or lower the resolution, and repeat.
4. When in doubt, keep more PCs.

* Recipe: `recipes/scrnaseq/pc_robustness.py`.
* *Why:* the elbow is subjective; random-matrix thresholds assume homogeneous noise; permutation methods such as
  JackStraw are statistically sound but slow [Chung 2014]. A clustering that changes when one or two PCs are added
  depends on noise components. Extra PCs add little noise, whereas too few can merge rare populations.
* *Trade-off:* 11 extra clusterings. The thresholds (ARI 0.9, Jaccard 0.75) are project defaults; the Jaccard value
  follows the stability reading of §8 [Hennig 2007].
* Record $p$, both criteria and the robustness table.
* **UMAP (uniform manifold approximation and projection) and t-SNE (t-distributed stochastic neighbour embedding) are for display only** [Chari 2023†].
  * *Why:* 2-D embeddings distort distances and densities; gaps, cluster sizes and relative positions in the plot are
    not quantitative evidence, and their appearance changes with parameters. For t-SNE, use PCA initialization and
    settings adapted to dataset size [Kobak 2019].

## 7. Integration

* **Integrate only when samples separate by technical factors** that would otherwise dominate clustering; check first
  with PCA and UMAP coloured by `sample_id` and `batch`.
* **Default: Harmony** [Korsunsky 2019].
  * *Why:* among eight methods, Harmony was the only one that did not introduce detectable artefacts [Antonsson 2025];
    it is fast and among the best in an earlier benchmark [Tran 2020].
* **Second opinion: scVI or scANVI** [Lopez 2018], also the choice for complex atlas-level integration, where they
  performed best [Luecken 2021]; reference mapping: scArches [Lotfollahi 2021].
  * *Why a second opinion:* if a conclusion depends on the integration method, it is a property of the method.
* **Evaluate** with `scib-metrics` (batch mixing and conservation of biology) [Luecken 2021].
* **Integrated embeddings and corrected values are for clustering and visualization only**; DE uses raw counts with
  batch in the model [Nguyen 2023; Nygaard 2015†].
* If condition and batch are confounded, integration can remove the biological effect; report it as a limitation.

## 8. Clustering

* **Leiden** on the neighbour graph [Traag 2019].
* **Choose the resolution by stability**:
  1. Run a grid of resolutions (e.g. 0.1 to 2.0 in steps of 0.1).
  2. For each resolution, subsample 80% of cells $N$ times and re-cluster (neighbour graph recomputed on the subsample
     from the same embedding, same parameters, different seed). Default $N$ by dataset size: $N = 10$ up to 50,000
     cells, $N = 5$ up to 200,000, $N = 1$ above. With $N = 1$ the result is indicative only; say so.
  3. For each original cluster, compute the Jaccard index with its best-matching cluster in each subsample, and take
     the median [Tang 2020; Hennig 2007]. Reading [Hennig 2007]: below 0.5 the cluster dissolves; 0.6–0.75
     indicates a pattern but is doubtful; 0.75 or more is stable. The recipe labels the whole interval [0.5, 0.75) as
     `doubtful`.
  4. Choose the highest resolution at which all clusters of interest are stable.
  5. Inspect splits across resolutions with a clustering tree [Zappia 2018]; subsampling robustness scores
     [Patterson-Cross 2021] and stability-based estimates of the number of clusters [Yu 2022] are alternatives.
  6. Optionally, test whether neighbouring clusters differ significantly and merge those that do not [Grabski 2023].
  * *Why:* resolution is the main driver of how many "cell types" appear; without a criterion, over-clustering produces
    clusters that differ only by noise [Kiselev 2019; Grabski 2023].
  * *Trade-off:* $N$ clusterings per resolution; the size-dependent $N$ keeps the cost bounded on large datasets.
  * Recipe: `recipes/scrnaseq/clustering_stability.py`.
* Record resolution, PCs, number of neighbours, random seed and stability values.

## 9. Annotation and marker genes

* **Annotate in three steps**: automatic, manual review with known markers, verification [Clarke 2021].
* **With a reference atlas**: label transfer (scANVI) or CellTypist with a pretrained or custom-trained model
  [Domínguez Conde 2022]; a general classifier is also effective [Abdelaal 2019].
* **With a marker dictionary and no reference** (frequent case): CellTypist does not accept marker dictionaries (it
  needs a model trained on annotated cells). Use:
  * CellAssign (`scvi.external.CellAssign`), a probabilistic model on counts with a binary marker matrix
    [Zhang 2019]; or
  * marker-set scoring with `decoupler` (ULM: univariate linear model, or ORA: over-representation analysis, per cell), aggregated to clusters; positive and negative markers
    improve specificity [Ianevski 2022]. Hierarchical marker classifiers are an R alternative [Pliner 2019]; strategies
    are reviewed in [Pasquini 2021].
  * Optionally, train a CellTypist model on confidently annotated cells to transfer labels to **other** datasets; within
    the same dataset this is circular and does not validate the labels.
  * Recipe: `recipes/scrnaseq/marker_annotation.py` (decoupler ULM per cell, majority label and agreement per
    cluster).
* Dictionaries are species-specific and versioned in the project; record unmatched markers (§1.3).
* **Marker genes**: Wilcoxon, t-test and logistic regression perform best [Pullin 2024]. **Their p-values are not
  evidence of differences between clusters**, because clusters were defined on the same data
  (`STATISTICS_GUIDELINES.md` §6, circular analysis); report effect sizes (log fold change, fraction of expressing
  cells).

## 10. Differential expression between conditions

* **MUST use pseudobulk or mixed models; MUST NOT use cell-level tests** to compare conditions [Squair 2021†;
  Crowell 2020†; Junttila 2022†; Murphy 2023†].
  * *Why:* cell-level tests treat cells as replicates and report hundreds of false DE genes.
* **Pseudobulk procedure**:
  1. Sum raw counts per `sample_id` (or `donor_id`) and cell type.
  2. Drop pseudobulk samples with fewer than a minimum number of cells (e.g. 10–20) and record the threshold.
  3. Filter genes and run DE per cell type as in `RNASEQ_BULK_GUIDELINES.md` §5–8, with batch and other covariates in
     the design (§2 for several samples per donor).
* Recipe: `recipes/scrnaseq/pseudobulk.py` (aggregation with cell counts per pseudobulk sample; DE per cell type
  through `recipes/bulk/deseq2.py`; per-cell-type and global adjusted p-values).
* Mixed models with a random effect for donor are an alternative [Zimmerman 2021†].
* State whether multiple-testing correction is per cell type or across all tests (`STATISTICS_GUIDELINES.md` §7).

## 11. Differential abundance

* **Count cells per sample and cell type, with explicit zeros**, and use donors as replicates.
  * Recipe: `recipes/scrnaseq/abundance.py`.
* **Account for compositionality**: proportions sum to 1, so a real increase of one type lowers the proportions of all
  others; tests on separate proportions report spurious decreases. Compositional methods use a reference type or model
  the whole composition.
* **Choose the method by question and design** [Yi 2024]:
  * **Per cluster**: propeller (R; transformed proportions with limma; fast, needs a few samples per group)
    [Phipson 2022]; scCODA (`pertpy`; Bayesian compositional; suited to few samples; needs a reference cell type)
    [Büttner 2021]; sccomp (R; models variability and outliers) [Mangiola 2023].
  * **Without clusters** (changes inside or across clusters): Milo (`pertpy`) on k-NN (k-nearest neighbours) neighbourhoods [Dann 2021];
    DA-seq [Zhao 2021b].
  * In a benchmark, Meld was most accurate without batch effects and Milo with them; including batch in the model
    improved all methods, and all methods performed poorly with strong imbalance between cell types [Yi 2024].
* MUST NOT test proportions with a chi-square test on pooled cell counts: it treats cells as replicates.

## 12. Cell-cell communication

* Inference combines ligand–receptor databases with coordinated expression of ligands and receptors [Armingol 2020].
* Predictions depend strongly on the ligand–receptor resource and the method [Dimitrov 2022; Liu 2022]. Use a
  consensus across methods with LIANA+ [Dimitrov 2024], and state the resource and its version.
* **Comparisons between conditions** are done per sample, with samples as replicates (LIANA+ multi-sample workflows),
  not by comparing pooled scores.
* Results are hypotheses; support them with spatial co-localization (`SPATIAL_GUIDELINES.md` §10), expression of
  downstream targets, or experiments.

## 13. Trajectories, fate mapping and RNA velocity

* **Topology first**: PAGA summarizes connectivity between clusters and preserves global topology; use it to check that
  a trajectory is supported by the graph [Wolf 2019].
* **Pseudotime**: choose the method by the expected topology and dataset size [Saelens 2019]: diffusion pseudotime
  [Haghverdi 2016], Palantir [Setty 2019], Slingshot (R) [Street 2018]. A root cell or state MUST be justified by
  biology.
* **Fate mapping with CellRank 2** [Lange 2022; Weiler 2024]. Kernels provide directionality from different sources and
  can be combined with weights:
  * `ConnectivityKernel`: expression similarity, undirected; used in combination;
  * `PseudotimeKernel`: biases transitions along a pseudotime;
  * `CytoTRACEKernel`: developmental potential [Kang 2024];
  * `RealTimeKernel`: experimental time points;
  * `VelocityKernel`: RNA velocity (see caveats below).
  The GPCCA (generalized Perron cluster cluster analysis) estimator identifies macrostates and fate probabilities. Report the kernels, weights and parameters.
* **RNA velocity** [La Manno 2018; Bergen 2020] is a hypothesis, not evidence of direction:
  * results depend on spliced/unspliced quantification [Soneson 2021], on the k-NN graph [Zheng 2023] and on model
    assumptions [Gorin 2022; Bergen 2021];
  * use methods that report uncertainty (veloVI) to identify datasets where velocity is inappropriate [Gayoso 2023];
  * MUST NOT use velocity to validate the embedding it is projected on [Zheng 2023].
* **Comparisons across conditions** use methods that account for sample-to-sample variability [Hou 2023]; genes
  along trajectories are tested with dedicated models [Van den Berge 2020].

## 14. Human and mouse

* Species-specific annotation, marker lists, mitochondrial chromosome name and ribosomal lists.
* Cross-species comparisons use one-to-one orthologs from a stated release.

## 15. What to report

Following minSCe [Füllgrabe 2020] and `STATISTICS_GUIDELINES.md` §12:

1. Tissue handling: dissociation or nuclei protocol, chemistry, multiplexing.
2. STAR version and `--solo*` parameters (or the alternative tool), reference and annotation release.
3. Per capture: barcodes, cells after each QC step (cell calling, ambient, nuclear fraction, filters, doublets), final
   cells per cell type.
4. Normalization, gene selection, PCs and criteria, integration method and metrics, clustering resolution and stability.
5. Annotation method, references or dictionaries, unmatched markers.
6. DE, abundance, communication and trajectory methods with parameters, pseudobulk thresholds and designs.

## 16. Implementation

| Task | Python | R reference (via script) |
|---|---|---|
| Quantification | STARsolo (CLI); alternatives `simpleaf`/alevin-fry, `kb-python` | — |
| Data structure | `anndata` | — |
| QC, normalization, HVG (highly variable genes), PCA, Leiden, UMAP, PAGA | `scanpy` | — |
| Cell calling and ambient RNA | STARsolo `EmptyDrops_CR`, `cellbender` (CLI) | SoupX, DropletUtils (EmptyDrops) |
| Doublets | `scanpy.pp.scrublet`, `doubletdetection` | scDblFinder, DoubletFinder |
| Size factors | `scranPY` (GitHub, pinned commit) | scran |
| Integration | `harmonypy` (or `scanpy.external`), `scvi-tools` | — |
| Integration metrics | `scib-metrics` | — |
| Annotation | `celltypist`, `scvi-tools` (scANVI, CellAssign), `decoupler` | scType, Garnett |
| Pseudobulk and condition DE | `decoupler`, `pydeseq2` | muscat, DESeq2, edgeR, limma |
| Differential abundance | `pertpy` (scCODA, Milo) | speckle (propeller), sccomp, miloR |
| Cell-cell communication | `liana` | — |
| Pseudotime and fate | `scanpy` (DPT, diffusion pseudotime), `palantir`, `cellrank` | slingshot, tradeSeq, Lamian |
| RNA velocity | `scvelo`, `velovi` | — |
