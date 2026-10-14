#!/usr/bin/env Rscript
# Recipe: Salmon transcript quantifications -> gene-level counts with tximport [Soneson 2016a].
# Guideline: RNASEQ_BULK_GUIDELINES.md §4 (quantification, transcript to gene).
#
# nf-core/rnaseq (star_salmon) already writes these tables (salmon.merged.gene_counts_length_scaled.tsv and
# salmon.merged.gene_lengths.tsv); read them with recipes/bulk/counts_io.py. Use this script for Salmon runs made
# outside nf-core.
#
# Inputs
#   --samples   TSV with columns sample, quant (path to quant.sf).
#   --tx2gene   TSV with columns transcript_id, gene_id (versions as in quant.sf; stripped on output).
#   --mode      lengthScaledTPM (for pydeseq2, default) or offsets (raw counts + length matrix for DESeq2/edgeR in R).
# Outputs (prefix --out)
#   <out>.counts.tsv   genes x samples; lengthScaledTPM counts or raw estimated counts (mode offsets).
#   <out>.length.tsv   genes x samples average transcript length (needed for offsets and for length-bias checks).
#   <out>.tpm.tsv      genes x samples abundance (TPM), for within-sample use only.

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
})

options <- list(
  make_option("--samples", type = "character"),
  make_option("--tx2gene", type = "character"),
  make_option("--mode", type = "character", default = "lengthScaledTPM"),
  make_option("--out", type = "character")
)
args <- parse_args(OptionParser(option_list = options))
stopifnot(args$mode %in% c("lengthScaledTPM", "offsets"))

samples <- read.delim(args$samples, stringsAsFactors = FALSE)
stopifnot(all(c("sample", "quant") %in% names(samples)), !anyDuplicated(samples$sample), all(file.exists(samples$quant)))
files <- setNames(samples$quant, samples$sample)
tx2gene <- read.delim(args$tx2gene, stringsAsFactors = FALSE)[, c("transcript_id", "gene_id")]

counts_from_abundance <- if (args$mode == "lengthScaledTPM") "lengthScaledTPM" else "no"
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = counts_from_abundance)

write_matrix <- function(m, path) {
  out <- data.frame(gene_id = sub("\\.[0-9]+", "", rownames(m)), m, check.names = FALSE)
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
write_matrix(txi$counts, paste0(args$out, ".counts.tsv"))
write_matrix(txi$length, paste0(args$out, ".length.tsv"))
write_matrix(txi$abundance, paste0(args$out, ".tpm.tsv"))
