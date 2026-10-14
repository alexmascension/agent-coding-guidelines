#!/usr/bin/env Rscript
# Recipe: pooled (deconvolution) size factors with scran [Lun 2016].
# Guideline: SCRNASEQ_GUIDELINES.md §5 (normalization). Called by recipes/scrnaseq/normalization.py through files;
# it can also be run on its own.
#
# Inputs
#   --counts    Matrix Market file, genes x cells, raw counts (integer).
#   --barcodes  Text file, one barcode per line, in column order.
#   --clusters  Optional text file, one cluster label per line (same order). Without it, scran::quickCluster is used.
#   --min-mean  Minimum mean count of genes used for pooling; 0.1 is the scran recommendation for UMI data.
#   --seed      Seed for quickCluster.
# Output
#   --out       TSV with columns barcode, size_factor (centred to mean 1), cluster (the clusters used).

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(scran)
})

options <- list(
  make_option("--counts", type = "character"),
  make_option("--barcodes", type = "character"),
  make_option("--clusters", type = "character", default = NA),
  make_option("--min-mean", type = "double", default = 0.1, dest = "min_mean"),
  make_option("--seed", type = "integer", default = 0L),
  make_option("--out", type = "character")
)
args <- parse_args(OptionParser(option_list = options))

counts <- as(readMM(args$counts), "CsparseMatrix")
barcodes <- readLines(args$barcodes)
stopifnot(ncol(counts) == length(barcodes))

if (is.na(args$clusters)) {
  set.seed(args$seed)
  clusters <- quickCluster(counts, min.mean = args$min_mean)
} else {
  clusters <- factor(readLines(args$clusters))
  stopifnot(length(clusters) == ncol(counts))
}

size_factors <- calculateSumFactors(counts, clusters = clusters, min.mean = args$min_mean)
if (any(size_factors <= 0)) {
  stop(sprintf("%d non-positive size factors; check low-quality cells", sum(size_factors <= 0)))
}
size_factors <- size_factors / mean(size_factors)

write.table(
  data.frame(barcode = barcodes, size_factor = size_factors, cluster = as.character(clusters)),
  file = args$out, sep = "\t", quote = FALSE, row.names = FALSE
)
