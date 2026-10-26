# Sample identity and duplicate detection guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 2.3.0 (2026-10-26)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys point to `OMICS_REFERENCES.md`; keys with a
> dagger point to `STATISTICS_REFERENCES.md`. Each rule is followed by ***Why***, ***Trade-off*** and
> ***Deviate when*** where they apply.
> A rule without a key word is a default with the force of SHOULD. Rules carry stable identifiers (e.g.
> `ID-WHY-001`); overrides refer to them (`CODING_GUIDELINES.md` §0, rule identifiers).

## 0. Scope and purpose

* Verify that every sample comes from the individual and tissue it is labelled with, and find duplicated, swapped,
  mixed or contaminated samples, before any analysis that combines samples.
* Applies to bulk RNA-seq, single-cell and single-nucleus RNA-seq, and to comparisons of RNA with DNA (WGS: whole-genome sequencing, exome,
  arrays); human and mouse.
* Copy with `OMICS_REFERENCES.md` and `SEX_INFERENCE_GUIDELINES.md`, which provides the sex check used here.
* Vocabulary: a **duplicate** is the same library or sequencing data present twice (e.g. a re-deposit); a **replicate
  library** is a new library from the same RNA; a **swap** is a sample labelled with the wrong individual; a
  **mixture** or **contamination** is material from two individuals in one sample.

## 1. Why identity checks are mandatory

* `ID-WHY-001` **Every multi-sample dataset MUST be checked before analysis.**
  * *Why:* mix-ups are common and costly:
    * MixupMapper found on average 3% of samples mixed up across five genetics-of-expression datasets (up to 23%), and
      correcting them increased cis-eQTL discoveries by 15% on average [Westra 2011];
    * a mouse intercross had 18% genotype–expression mix-ups, mostly off-by-one plate pipetting errors [Broman 2015];
    * LD-based fingerprinting found about 1% mislabelled datasets in ENCODE [Javed 2020];
    * sex-annotation mismatches appear in 46% of public transcriptomics datasets [Toker 2016].
  * *Trade-off:* genotype extraction from BAMs adds minutes per sample; the fingerprints are small and kept forever.

## 2. Order of checks

1. `ID-ORDER-001` **Exact and near-exact duplicates** from the data themselves (§3).
2. `ID-ORDER-002` **Sex concordance** (`SEX_INFERENCE_GUIDELINES.md`).
3. `ID-ORDER-003` **Genotype fingerprints RNA↔RNA**: all samples of one donor must match each other and no other donor (§4).
4. `ID-ORDER-004` **Genotype RNA↔DNA**, when DNA exists: each sample must match its donor's DNA (§4).
5. `ID-ORDER-005` **Mixture and contamination** estimates (§5).
6. `ID-ORDER-006` **Resolution** of every flag (§7).

* `ID-ORDER-007` *Why this order:* the early checks are cheap and explain many later flags (a duplicate pair is also a "same
  individual" pair; a sex mismatch predicts a genotype mismatch).

## 3. Duplicates from the data

* `ID-DUP-001` **Exact duplicates:** hash the raw read files (checksum) and the count vectors; identical hashes are duplicates.
* `ID-DUP-002` **Near-duplicates** (same library re-sequenced or re-deposited, possibly with different read subsets): compare
  count vectors with a correlation of log counts; pairs whose correlation is an outlier **relative to the distribution
  of same-tissue pairs** are candidates.
  * *Why relative:* expression similarity is dominated by tissue and protocol, so different individuals of the same
    tissue are highly correlated; a fixed cut-off (e.g. $r > 0.99$) fails across tissues. Chance near-identical
    samples ("doppelgängers") also exist and inflate machine-learning evaluations when they fall in both training and
    test sets [Wang 2022].
  * *Trade-off / Deviate when:* expression **cannot establish identity**; it only nominates pairs. Identity is decided
    by genotype (§4).
* `ID-DUP-003` Check metadata in parallel: identical read counts, read names, flowcell or lane identifiers, or submission
  accessions pointing to the same run.
* `ID-DUP-004` **Re-deposited libraries** (the same library under several accessions, e.g. in more than one SubSeries of a GEO
  SuperSeries) are nominated by identical read counts with identical donor and tissue annotation. Process one copy
  and record the mapping between accessions. Keep a few re-deposited pairs as **positive controls**: both copies are
  processed in full and MUST come out as duplicates (identical counts, genotype relatedness close to 1); if they do
  not, the identity checks are faulty.
* `ID-DUP-005` **Compare metadata within duplicate pairs**: discordant annotations of the same library (e.g. tissue or region)
  are flagged and resolved in the decisions log, not silently harmonized.

## 4. Genotype fingerprints

* `ID-GENO-001` **Use a fingerprinting tool on the aligned reads of every sample**: somalier [Pedersen 2020] by default; NGSCheckMate
  [Lee 2017] or Picard CrosscheckFingerprints [Javed 2020] as alternatives; BAM-matcher [Wang 2016] for pairwise checks.
  * *Why somalier:* it extracts allele counts at a fixed panel of common coding SNPs into small files ("sketches")
    once per sample, compares all pairs in seconds, does not need joint variant calling, and was validated on GTEx
    WGS against RNA-seq [Pedersen 2020]. NGSCheckMate works across data types at depths above 0.5X [Lee 2017];
    Crosscheck uses linkage disequilibrium and remains effective on sparse data [Javed 2020].
  * *Trade-off:* fingerprints need the aligned reads; extract them before BAMs are deleted (or later from the archived CRAM,
    `HPC_GUIDELINES.md` §2, `HPC-FLOW-006`), and keep them.
* `ID-GENO-002` **Match site names to the reference**: rename contigs of the sites file (`chr1` versus `1`) to match the BAMs.
* `ID-GENO-003` **Metrics per pair**: relatedness (about 1 for the same individual, about 0.5 for first-degree relatives, about 0
  for unrelated), IBS0 (number of sites where the two samples are opposite homozygotes, i.e. share no allele; near 0 for
  the same individual and for parent–child), IBS2 (number of sites with identical genotypes), and the number of
  informative sites shared. somalier reports IBS0 and IBS2 as counts (`ibs0`, `ibs2`, out of `n` sites); divide by `n`
  to compare pairs with different numbers of shared sites.
* `ID-GENO-004` **Thresholds come from the cohort distribution, not from fixed values.**
  * *Why:* in RNA only expressed genes are covered, and allele-specific expression or RNA editing make heterozygous
    sites look homozygous in some tissues. Same-individual pairs across tissues therefore show relatedness below 1 and
    tissue-dependent values. The expected structure is bimodal (same individual versus different); place the cut in
    the gap and report pairs inside it.
  * *Deviate when:* with very few informative shared sites (e.g. low-expression tissues), treat the pair as
    undetermined rather than unrelated.
* `ID-GENO-005` **RNA↔DNA**: extract the same sites from the DNA (VCF or CRAM) and compare. Each RNA sample MUST match its donor's
  DNA.
  * *Why:* RNA↔RNA cannot detect a donor with a single sample, nor a donor whose samples were all swapped consistently
    with another donor's.
* `ID-GENO-006` **Expression-based genotype prediction** (eQTL-based, e.g. MixupMapper, R/lineup) MAY be used when no BAMs are
  available [Westra 2011; Broman 2015].

## 5. Mixtures and contamination

* `ID-MIX-001` **Male material in an XX sample** is detected by Y signal and homolog ratios (`SEX_INFERENCE_GUIDELINES.md` §5).
* `ID-MIX-002` **Genotype-based contamination** estimates (allele fractions inconsistent with one individual) exist for DNA:
  VerifyBamID detects contamination down to about 1% [Jun 2012]; Conpair down to 0.1% in tumour–normal pairs
  [Bergmann 2016].
  * *Trade-off:* in RNA, allelic imbalance from allele-specific expression mimics contamination; use DNA-derived
    methods on RNA only as screening, and confirm with DNA or with a second sample.
* `ID-MIX-003` A contaminated sample typically shows excess heterozygosity and elevated relatedness to the contaminating sample;
  these signals come from the fingerprint metrics of §4.

## 6. Single-cell pools

* `ID-POOL-001` **Pooled designs (cells of several donors in one capture) MUST be demultiplexed by genotype**, also when hashtags are used (hashtag calls are then checked against genotype):
  * with donor genotypes: Demuxalot with refined genotypes; Demuxalot and Dropulation combined when more than one tool
    can be run;
  * without genotypes: vireo [Huang 2019] with more than 20 donors per pool or when only one tool can be run; souporcell
    [Heaton 2020] and vireo combined with up to 20 donors, when more than one tool can be run;
  * combined calls: a droplet is a singlet only if more than half of the tools call it a singlet, and its donor is the
    one most tools assign; ties are left unassigned [Neavin 2024];
  * *Why:* on in silico pools of 2–128 donors with known doublets, Demuxalot and Dropulation reached a higher MCC than
    demuxlet, and Demuxalot with refined genotypes plus Dropulation was the recommended combination for large pools;
    without genotypes, souporcell lost accuracy above 32 donors, vireo alone performed best in large pools, and
    souporcell with vireo was part of the recommended combinations for smaller pools [Neavin 2024]. The authors place
    the limit at 16 donors in the text and at about 20 in their decision tree.
  * *Deviate when:* demuxlet [Kang 2017] MAY replace Demuxalot where it is already established; in its original
    evaluation about 50 SNPs per cell sufficed for 97% singlet assignment in pools of up to 64.
* `ID-POOL-002` **Match demultiplexed clusters to donors** with the fingerprints of §4 (pseudobulk BAM per cluster against donor DNA
  or other samples).
* `ID-POOL-003` Genotype demultiplexing also detects cross-donor doublets and estimates ambient RNA [Heaton 2020].
* `ID-POOL-005` **Label clusters with the strongest available evidence, in this order, and record per library which one was used**:
  (1) per-barcode donor labels from the data authors, checked against your own clusters (most labelled barcodes of a
  cluster agree on one donor); (2) genotype references of the same donors from other data (bulk RNA-seq, DNA genotypes,
  another tissue, or another pool containing the donor) matched to per-cluster fingerprints, RNA against RNA being
  adequate on common exonic sites; (3) sex, only in pools of one XX and one XY donor; (4) otherwise keep clusters as
  unnamed individuals of the pool.
  * *Why:* genotype-free demultiplexing returns anonymous clusters; donor attributes (condition, age, genotype) can
    only be attached through external evidence, and its strength differs between libraries.
* `ID-POOL-006` **Before deleting a pooled library's BAM, MUST store one genotype fingerprint per cluster** (reads split by cell
  barcode) together with the barcode-to-cluster assignment, so that clusters can be labelled when a reference becomes
  available without re-mapping.
* `ID-POOL-007` **A cluster-to-donor match MUST be one-to-one, MUST name a donor of that pool, and MUST pass stated thresholds** on
  relatedness, on the margin to the next-best donor and on the number of shared sites; report the thresholds and the
  full cluster-by-reference relatedness table.
* `ID-POOL-008` **Unlabelled clusters of a pool whose donors differ in condition MUST NOT be used for condition comparisons.**
  Unlabelled clusters of a pool whose donors share the condition MAY be used as individuals of that condition, with
  donor-level covariates set to missing.
* `ID-POOL-009` Genotype-free demultiplexers (vireo, souporcell, freemuxlet from popscle) are interchangeable for this purpose
  when run on common population variants (e.g. allele frequency 0.1–0.9) in transcribed regions; the labelling
  rules above do not depend on the demultiplexer. For singlet and doublet calls they are not interchangeable
  (`ID-POOL-001`).
* `ID-POOL-004` Reliability varies by species and depends on genetic diversity; inbred mouse strains cannot be demultiplexed by
  genotype [Cardiello 2023].

## 7. Resolving flags

| Finding | Likely explanation | Resolution |
|---|---|---|
| Two samples, same individual, labelled as different donors | Duplicate, or donor ID error | If exact duplicate, keep one; otherwise correct the donor ID, with evidence |
| Samples of one donor do not match each other (genotype) | Swap | Reassign if the true donor is identified by genotype; otherwise exclude |
| RNA does not match its DNA | Swap of RNA or of DNA | Check the other samples of the donor; exclude if unresolved |
| Sex mismatch, genotype matches the donor | Metadata error in sex; also mosaic loss of Y, tissue or pathology effects, low quality, sex-mismatched transplantation | Review; correct the sex annotation only with independent evidence (DNA karyotype, clinical record, concordant other tissues), recorded in the decisions log |
| Sex mismatch, genotype does not match | Sample swap | As for a swap |
| Excess heterozygosity, Y in XX | Mixture | Exclude, or keep with a flag only for analyses robust to it |

* `ID-RES-001` **Corrections MUST NOT be silent.** Every change and exclusion is recorded in the decisions log with the evidence and kept
  in a sample-level QC (quality control) table.

## 8. Outputs and report

* `ID-REP-001` Per sample: duplicate group, sex call and concordance, best-matching donor by genotype, relatedness to the labelled
  donor and to the best other match, contamination flags, final status (`ok`, `corrected`, `excluded`) and reason.
* `ID-REP-002` Per pair (for flagged pairs): relatedness, IBS0, IBS2, number of shared sites, expression correlation.
* `ID-REP-003` Figures: relatedness against IBS0 for all pairs coloured by expected relationship; heatmap of relatedness ordered by
  donor.
* `ID-REP-004` Tool versions, site panel, and thresholds with how they were chosen.

## 9. Implementation

| Task | Tool | Language |
|---|---|---|
| Fingerprints, relatedness | somalier | CLI |
| Alternatives | NGSCheckMate, Picard CrosscheckFingerprints, BAM-matcher | CLI |
| Expression-based genotype prediction | R/lineup, MixupMapper | R, Java |
| DNA contamination | VerifyBamID, Conpair | CLI |
| Single-cell demultiplexing | Demuxalot, `vireoSNP` (Python); Dropulation (Drop-seq tools), souporcell, demuxlet | Python, CLI |
| Hashing, correlations, tables | `hashlib`, `polars`, `numpy` | Python |
