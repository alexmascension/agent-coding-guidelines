# Sex inference guidelines

> Source: <https://github.com/alexmascension/agent-coding-guidelines> · Version 1.0.1 (2026-10-14)
>
> Key words **MUST**, **SHOULD**, **MAY** follow RFC 2119. Bracketed keys such as `[Olney 2020]` point to
> `OMICS_REFERENCES.md`; keys with a dagger point to `STATISTICS_REFERENCES.md`. Each rule is followed by its
> rationale: ***Why*** (mechanism), ***Trade-off*** (cost) and ***Deviate when*** (justified exceptions), where they
> apply.

## 0. Scope and purpose

* Infer the sex-chromosome complement of each sample from transcriptomic data, compare it with the reported sex, and
  use it to (a) detect sample swaps and annotation errors, (b) choose sex-matched references, (c) flag anomalies
  (aneuploidy, mosaic loss of Y, mixtures).
* Applies to bulk RNA-seq and to single-cell and single-nucleus RNA-seq, human and mouse. The procedure is independent
  of aligner and counting tool; only the inputs (counts per marker gene and library size) are required.
* Copy with `OMICS_REFERENCES.md`. Related: `SAMPLE_IDENTITY_GUIDELINES.md` (genotype-based identity),
  `RNASEQ_BULK_GUIDELINES.md` §2 (reference genome).
* Terminology: this guide infers **chromosomal sex** (XX, XY and anomalies) from expression; it does not infer gender.
  Results are labelled `XX`, `XY`, etc., never "female" or "male", and are compared with the reported sex in a separate
  column.

## 1. Why run it on every dataset

* **Every dataset with sex metadata MUST be checked**, including public datasets and the user's own.
  * *Why:* when reported sex was compared with sex-specific gene expression, 46% of public human transcriptomics
    datasets contained at least one apparently mislabelled sample, and in a single cohort two of four datasets had
    mix-ups, pointing to laboratory errors rather than typos [Toker 2016]. A sex mismatch is the cheapest detectable
    symptom of a sample swap; genotype checks find more (`SAMPLE_IDENTITY_GUIDELINES.md`).
  * *Trade-off:* minutes of computation per dataset.
  * *Deviate when:* never for analyses that use sex as a variable or that combine samples across donors.

## 2. Marker genes

* **Female-specific signal: *XIST*** (mouse *Xist*), optionally with *TSIX*.
  * *Why:* *XIST* is expressed from the inactive X, so it is present in XX and essentially absent in XY cells in all
    adult tissues. It is nuclear-retained, which makes it especially strong in single-nucleus data.
  * *Deviate when:* do not add other X-linked genes as female markers. No other adult X gene is XX-specific; genes
    that escape X inactivation (*KDM5C*, *KDM6A*, *DDX3X*, *ZFX*, *EIF1AX*) are only modestly higher in XX and overlap
    between sexes [Tukiainen 2017; Oliva 2020].
* **Male-specific signal: broadly expressed male-specific region of Y (MSY) genes**, summed. Human: *RPS4Y1*,
  *DDX3Y*, *KDM5D*, *UTY*, *EIF1AY*, *USP9Y*, *ZFY* (and *NLGN4Y* in brain). Mouse: *Ddx3y*, *Eif2s3y*, *Kdm5d*, *Uty*.
  * *Why:* a five-gene panel (*RPS4Y1*, *EIF1AY*, *DDX3Y*, *KDM5D*, *XIST*) classified sex across 24 tissue collections
    and platforms [Staedtler 2013]; X-degenerate Y genes are expressed across most tissues [Godfrey 2020]. Summing
    several genes reduces dependence on any one gene's tissue-specific expression.
  * *Deviate when:* drop genes not expressed in the tissue at hand (check on samples of known sex).
* **Within-sample ratio from X/Y homolog pairs**: for each pair, $r = Y/(X+Y)$, then the median across pairs. Human
  pairs include *RPS4X/RPS4Y1*, *DDX3X/DDX3Y*, *KDM5C/KDM5D*, *KDM6A/UTY*, *EIF1AX/EIF1AY*, *USP9X/USP9Y*, *ZFX/ZFY*;
  mouse: *Ddx3x/Ddx3y*, *Eif2s3x/Eif2s3y*, *Kdm5c/Kdm5d*, *Kdm6a/Uty*.
  * *Why:* the X homolog is a built-in denominator. The ratio cancels library size, tissue and RNA quality, and
    separates XY ($r$ substantial) from XX ($r \approx 0$) without a cohort-level threshold. It is also the most
    sensitive signal for male material mixed into an XX sample.
  * *Trade-off:* homologs share sequence; with a reference that does not mask the Y PAR (pseudoautosomal region), or with multi-mapping reads,
    reads can be misassigned (§4).

## 3. Counting and normalization

* **Count marker genes over the whole gene body (exons and introns)** with uniquely mapped reads, when the library
  captures intronic reads (total RNA, nuclei).
  * *Why:* in total RNA and in nuclei about half of the reads are intronic, and several Y genes have long loci, so the
    signal multiplies several-fold. Introns of X-degenerate Y genes have diverged from their X homologs, so unique
    reads are unambiguous.
  * *Deviate when:* for poly(A) or 3'-tag libraries, exonic counts are sufficient.
* **Normalize to counts per million uniquely mapped reads (or fragments)** and analyse on the $\log_{10}(1+x)$ scale.
  * *Why:* the distributions are bimodal with many values near zero; the log scale separates the modes and stabilizes
    Poisson noise.
* **Screening depth:** for XX versus XY, about 1 million read pairs per sample is enough. Expected values are on the
  order of tens to hundreds of informative reads, and observing zero MSY reads in an XY sample is then negligible.
  Re-evaluate anomalies at full depth.
  * *Why:* detecting a 5–10% mixture or partial loss of Y needs the full library; a 10% male contribution to an XX
    sample adds only a few Y reads at 1 million pairs.
* **Minimum informative reads:** samples with fewer than a fixed number of informative reads (XIST + MSY + X homologs;
  about 30 as a starting value) are `LOW_SIGNAL` and are re-run at higher depth, not called.

## 4. Reference genome

* **Use sex-chromosome-complement-informed references** for the main alignment: Y fully masked for XX samples; Y PAR
  masked for XY samples (the PAR is then represented only on X) [Olney 2020].
  * *Why:* with the default reference, reads from X/Y homologous and PAR sequence are split or misassigned, which
    lowers PAR expression estimates and produces spurious Y expression in XX samples [Olney 2020; Webster 2019].
  * *Trade-off:* two references and a routing step based on the screening call (§6).
* The screening itself MUST use a reference that contains Y (XY reference); otherwise Y signal cannot be observed.
* The annotation MUST NOT contain features on masked sequence (e.g. genes annotated in the Y PAR for the XY reference);
  remove them, or the corresponding transcripts are all `N`.

## 5. Thresholds and categories

* **Thresholds are data-driven by default**: on the log scale, a two-group split (Otsu, or a two-component Gaussian
  mixture), with the cut placed between the medians of the two groups. Unsupervised clustering of Y-gene expression is
  the same idea applied to microarrays [Buckberry 2014].
  * *Why:* absolute levels depend on tissue, protocol and depth, so a fixed threshold that works in blood fails in
    brain.
  * *Trade-off / Deviate when:* a data-driven split needs both groups present in reasonable numbers. With a single-sex
    cohort, a strongly unbalanced cohort, or fewer than about 20 samples, use thresholds calibrated on samples of known
    sex processed the same way, and say so.
* **Y presence requires both** absolute Y signal above threshold **and** a median homolog ratio above a minimum (e.g.
  $r \ge 0.10$).
  * *Why:* requiring both avoids calling Y from mapping noise in XX samples.
* **Categories** (one per sample):

  | Call | XIST | Y | Meaning | Action |
  |---|---|---|---|---|
  | `XX` | high | absent | XX | XX reference |
  | `XY` | absent | present | XY | XY reference |
  | `XY_lowY` | absent | present but below a fraction (e.g. 0.5) of the XY median | Possible mosaic loss of Y, or mixture | XY reference; review (§7) |
  | `XIST+Y+` | high | present | XXY, or mixture of XX and XY material | XY reference; review with genotype |
  | `XIST-Y-` | absent | absent | X0, or technical failure | XY reference; review |
  | `LOW_SIGNAL` | — | — | Too few informative reads | Re-run deeper |

* **Flag `possible_mix`** in `XX` samples when the homolog ratio is above a trace level (e.g. 0.02) or Y signal is far
  above the XX median.
* **Report, for every call, the distance to each threshold** (in log units), not only the category.
  * *Why:* calls near a threshold are the ones that need review; a bare label hides them.
* **Calibrate the mixture and loss-of-Y flags with in-silico mixtures**: combine reads from an XX and an XY sample at
  known fractions (e.g. 1, 2, 5, 10, 20%) and record the fraction at which each flag fires.
  * *Why:* the trace thresholds above are starting values; their real sensitivity depends on depth and tissue.

## 6. Combining evidence and deciding

* **Order of evidence** for the karyotype used in analysis and for reference routing:
  1. WGS (whole-genome sequencing) or genotyping-array karyotype of the donor, when available.
  2. Expression call (this guide).
  3. Reported sex.
* **Genotype-based sex** from the same RNA BAM: heterozygosity on non-PAR X (XX and XXY have heterozygous sites; XY
  does not) and Y depth, e.g. with somalier [Pedersen 2020] or peddy (DNA) [Pedersen 2017].
  * *Why:* it is orthogonal to expression; a sample that is `XX` by expression but has no X heterozygosity is
    suspicious.
  * *Trade-off:* in RNA, only expressed genes are covered, and allele-specific expression (including X inactivation in
    clonal tissue) can hide heterozygosity; treat it as supporting evidence.
* **All discordances are flagged, never silently corrected.** Metadata changes are recorded in the decisions log with
  the evidence.
* **Within a donor, all samples MUST agree.** A donor with discordant tissues indicates a swap
  (`SAMPLE_IDENTITY_GUIDELINES.md`).

## 7. Loss of Y and other anomalies

* **Mosaic loss of chromosome Y (LOY) is common and biological.** About 20% of men in UK Biobank have detectable LOY in
  blood [Thompson 2019], and it is associated with mortality [Forsberg 2014] and disease [Sano 2022]. In brain, LOY is
  enriched in microglia and more frequent in Alzheimer's disease [Vermeulen 2022]. It can be scored from bulk and
  single-cell RNA-seq [Dumanski 2021].
  * *Consequence:* `XY_lowY` in aged tissue MUST NOT be interpreted as contamination by default. Interpret it with cell
    composition (e.g. microglial fraction from deconvolution) and, when single-cell data exist, with per-cell Y
    detection.
* `XIST+Y+` is resolved with genotype: XXY shows X heterozygosity and a single genotype; a mixture shows excess
  heterozygosity, allele fractions inconsistent with one individual, and high relatedness to another sample.
* `XIST-Y-` with good library quality and no X heterozygosity suggests X0; with poor quality, a failed library.

## 8. Single-cell and single-nucleus data

* **Call sex per sample (per donor), from pseudobulk counts**, using §2–5.
  * *Why:* per-cell detection is sparse (many cells have zero counts of every marker), so per-cell calls are noisy;
    the sample-level signal is robust.
* **Per-cell scores MAY be used** to (a) check that a sample is homogeneous, (b) detect cells from another donor in
  pooled designs, (c) estimate per-cell LOY in XY donors.
* **Ambient RNA carries marker transcripts**: in pools of XX and XY donors, ambient *XIST* and Y transcripts appear in
  every droplet. Use ambient-corrected counts and, for pools, genotype demultiplexing rather than sex markers
  (`SAMPLE_IDENTITY_GUIDELINES.md` §6).
* In nuclei, *XIST* is strong (nuclear RNA) and Y genes are often detected from intronic reads; count with introns.

## 9. Mouse

* Use *Xist* and the mouse MSY genes and homolog pairs of §2. Mouse Y carries few broadly expressed genes, so the Y sum
  is smaller than in human.
* Inbred strains have little heterozygosity, so genotype-based sex from X heterozygosity does not work; rely on
  expression (and on breeding records).
* Sex-chromosome complement affects autosomal gene expression independently of gonadal sex [Wijchers 2010]; models
  such as XY females are not identified by these markers alone.

## 10. Outputs and report

* One table per dataset with, per sample: reported sex, call, distance to each threshold, informative reads,
  normalized XIST and Y, median homolog ratio, `possible_mix`, chosen reference, concordance (`concordant`,
  `DISCORDANT`, `no_reported_sex`, `uncertain_call`), and review flag.
* A scatter plot of $\log_{10}(1+\text{XIST})$ against $\log_{10}(1+\text{Y})$ coloured by reported sex, with the
  thresholds, discordant samples circled, and a histogram of the homolog ratio.
* The thresholds used (and whether data-driven or calibrated), the marker lists and the reference versions.
* Counts of discordances and anomalies, and the actions taken for each.
