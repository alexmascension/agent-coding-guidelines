# Sample identity and duplicate detection guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 1.0.0 (2026-10-14)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys point to `OMICS_REFERENCES.md`; keys with a
> dagger point to `STATISTICS_REFERENCES.md`. Each rule is followed by ***Why***, ***Trade-off*** and
> ***Deviate when*** where they apply.

## 0. Scope and purpose

* Verify that every sample comes from the individual and tissue it is labelled with, and find duplicated, swapped,
  mixed or contaminated samples, before any analysis that combines samples.
* Applies to bulk RNA-seq, single-cell and single-nucleus RNA-seq, and to comparisons of RNA with DNA (WGS, exome,
  arrays); human and mouse.
* Copy with `OMICS_REFERENCES.md` and `SEX_INFERENCE_GUIDELINES.md`, which provides the sex check used here.
* Vocabulary: a **duplicate** is the same library or sequencing data present twice (e.g. a re-deposit); a **replicate
  library** is a new library from the same RNA; a **swap** is a sample labelled with the wrong individual; a
  **mixture** or **contamination** is material from two individuals in one sample.

## 1. Why identity checks are mandatory

* **Every multi-sample dataset MUST be checked before analysis.**
  * *Why:* mix-ups are common and costly:
    * MixupMapper found on average 3% of samples mixed up across five genetics-of-expression datasets (up to 23%), and
      correcting them increased cis-eQTL discoveries by 15% on average [Westra 2011];
    * a mouse intercross had 18% genotype–expression mix-ups, mostly off-by-one plate pipetting errors [Broman 2015];
    * LD-based fingerprinting found about 1% mislabelled datasets in ENCODE [Javed 2020];
    * sex-annotation mismatches appear in 46% of public transcriptomics datasets [Toker 2016].
  * *Trade-off:* genotype extraction from BAMs adds minutes per sample; the fingerprints are small and kept forever.

## 2. Order of checks

1. **Exact and near-exact duplicates** from the data themselves (§3).
2. **Sex concordance** (`SEX_INFERENCE_GUIDELINES.md`).
3. **Genotype fingerprints RNA↔RNA**: all samples of one donor must match each other and no other donor (§4).
4. **Genotype RNA↔DNA**, when DNA exists: each sample must match its donor's DNA (§4).
5. **Mixture and contamination** estimates (§5).
6. **Resolution** of every flag (§7).

* *Why this order:* the early checks are cheap and explain many later flags (a duplicate pair is also a "same
  individual" pair; a sex mismatch predicts a genotype mismatch).

## 3. Duplicates from the data

* **Exact duplicates:** hash the raw read files (checksum) and the count vectors; identical hashes are duplicates.
* **Near-duplicates** (same library re-sequenced or re-deposited, possibly with different read subsets): compare
  count vectors with a correlation of log counts; pairs whose correlation is an outlier **relative to the distribution
  of same-tissue pairs** are candidates.
  * *Why relative:* expression similarity is dominated by tissue and protocol, so different individuals of the same
    tissue are highly correlated; a fixed cut-off (e.g. $r > 0.99$) fails across tissues. Chance near-identical
    samples ("doppelgängers") also exist and inflate machine-learning evaluations when they fall in both training and
    test sets [Wang 2022].
  * *Trade-off / Deviate when:* expression **cannot establish identity**; it only nominates pairs. Identity is decided
    by genotype (§4).
* Check metadata in parallel: identical read counts, read names, flowcell or lane identifiers, or submission
  accessions pointing to the same run.

## 4. Genotype fingerprints

* **Use a fingerprinting tool on the aligned reads of every sample**: somalier [Pedersen 2020] by default; NGSCheckMate
  [Lee 2017] or Picard CrosscheckFingerprints [Javed 2020] as alternatives; BAM-matcher [Wang 2016] for pairwise checks.
  * *Why somalier:* it extracts allele counts at a fixed panel of common coding SNPs into small files ("sketches")
    once per sample, compares all pairs in seconds, does not need joint variant calling, and was validated on GTEx
    WGS against RNA-seq [Pedersen 2020]. NGSCheckMate works across data types at depths above 0.5X [Lee 2017];
    Crosscheck uses linkage disequilibrium and remains effective on sparse data [Javed 2020].
  * *Trade-off:* fingerprints need the BAM; extract them before BAMs are deleted, and keep them.
* **Match site names to the reference**: rename contigs of the sites file (`chr1` versus `1`) to match the BAMs.
* **Metrics per pair**: relatedness (about 1 for the same individual, about 0.5 for first-degree relatives, about 0
  for unrelated), IBS0 (near 0 for the same individual and for parent–child), IBS2, and the number of informative sites
  shared.
* **Thresholds come from the cohort distribution, not from fixed values.**
  * *Why:* in RNA only expressed genes are covered, and allele-specific expression or RNA editing make heterozygous
    sites look homozygous in some tissues. Same-individual pairs across tissues therefore show relatedness below 1 and
    tissue-dependent values. The expected structure is bimodal (same individual versus different); place the cut in
    the gap and report pairs inside it.
  * *Deviate when:* with very few informative shared sites (e.g. low-expression tissues), treat the pair as
    undetermined rather than unrelated.
* **RNA↔DNA**: extract the same sites from the DNA (VCF or CRAM) and compare. Each RNA sample MUST match its donor's
  DNA.
  * *Why:* RNA↔RNA cannot detect a donor with a single sample, nor a donor whose samples were all swapped consistently
    with another donor's.
* **Expression-based genotype prediction** (eQTL-based, e.g. MixupMapper, R/lineup) MAY be used when no BAMs are
  available [Westra 2011; Broman 2015].

## 5. Mixtures and contamination

* **Male material in an XX sample** is detected by Y signal and homolog ratios (`SEX_INFERENCE_GUIDELINES.md` §5).
* **Genotype-based contamination** estimates (allele fractions inconsistent with one individual) exist for DNA:
  VerifyBamID detects contamination down to about 1% [Jun 2012]; Conpair down to 0.1% in tumour–normal pairs
  [Bergmann 2016].
  * *Trade-off:* in RNA, allelic imbalance from allele-specific expression mimics contamination; use DNA-derived
    methods on RNA only as screening, and confirm with DNA or with a second sample.
* A contaminated sample typically shows excess heterozygosity and elevated relatedness to the contaminating sample;
  these signals come from the fingerprint metrics of §4.

## 6. Single-cell pools

* **Pooled designs MUST be demultiplexed by genotype**:
  * with donor genotypes, demuxlet (about 50 SNPs per cell suffice for 97% singlet assignment in pools of up to 64)
    [Kang 2017];
  * without genotypes, souporcell [Heaton 2020] or vireo [Huang 2019];
  * combining several demultiplexing and doublet methods by consensus improves assignment [Neavin 2024].
* **Match demultiplexed clusters to donors** with the fingerprints of §4 (pseudobulk BAM per cluster against donor DNA
  or other samples).
* Genotype demultiplexing also detects cross-donor doublets and estimates ambient RNA [Heaton 2020].
* Reliability varies by species and depends on genetic diversity; inbred mouse strains cannot be demultiplexed by
  genotype [Cardiello 2023].

## 7. Resolving flags

| Finding | Likely explanation | Resolution |
|---|---|---|
| Two samples, same individual, labelled as different donors | Duplicate, or donor ID error | If exact duplicate, keep one; otherwise correct the donor ID, with evidence |
| Samples of one donor do not match each other | Swap | Reassign if the true donor is identified by genotype; otherwise exclude |
| RNA does not match its DNA | Swap of RNA or of DNA | Check the other samples of the donor; exclude if unresolved |
| Sex mismatch, genotype matches the donor | Metadata error in sex | Correct the sex annotation |
| Sex mismatch, genotype does not match | Sample swap | As for a swap |
| Excess heterozygosity, Y in XX | Mixture | Exclude, or keep with a flag only for analyses robust to it |

* **Never correct silently.** Every change and exclusion is recorded in the decisions log with the evidence and kept
  in a sample-level QC table.

## 8. Outputs and report

* Per sample: duplicate group, sex call and concordance, best-matching donor by genotype, relatedness to the labelled
  donor and to the best other match, contamination flags, final status (`ok`, `corrected`, `excluded`) and reason.
* Per pair (for flagged pairs): relatedness, IBS0, IBS2, number of shared sites, expression correlation.
* Figures: relatedness against IBS0 for all pairs coloured by expected relationship; heatmap of relatedness ordered by
  donor.
* Tool versions, site panel, and thresholds with how they were chosen.

## 9. Implementation

| Task | Tool | Language |
|---|---|---|
| Fingerprints, relatedness | somalier | CLI |
| Alternatives | NGSCheckMate, Picard CrosscheckFingerprints, BAM-matcher | CLI |
| Expression-based genotype prediction | R/lineup, MixupMapper | R, Java |
| DNA contamination | VerifyBamID, Conpair | CLI |
| Single-cell demultiplexing | `vireoSNP` (Python), souporcell, demuxlet | Python, CLI |
| Hashing, correlations, tables | `hashlib`, `polars`, `numpy` | Python |
