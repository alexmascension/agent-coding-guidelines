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
#   --overwrite Replace existing outputs (otherwise the script stops if they exist).
# Outputs
#   --out                    TSV with columns barcode, size_factor (centred to mean 1), cluster (the clusters used).
#   <out>.provenance.tsv     R and package versions, arguments and date.

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
  make_option("--out", type = "character"),
  make_option("--overwrite", action = "store_true", default = FALSE)
)
args <- parse_args(OptionParser(option_list = options))
script_name <- "scran_size_factors.R"

# Writes `df` to `path` through a temporary file in the same directory, renamed at the end, so that an interrupted run
# never leaves a partial output with the final name.
write_tsv_atomic <- function(df, path, overwrite) {
  if (file.exists(path) && !overwrite) stop(sprintf("%s exists; pass --overwrite to replace it", path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  write.table(df, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  if (!file.rename(tmp, path)) stop(sprintf("could not move %s to %s", tmp, path))
}

# Provenance next to the outputs: R and package versions, arguments, date.
write_provenance <- function(path, args, packages, overwrite) {
  versions <- vapply(packages, function(p) as.character(packageVersion(p)), character(1))
  arg_values <- vapply(args[setdiff(names(args), "help")], function(v) paste(v, collapse = ","), character(1))
  df <- data.frame(
    key = c("script", "date", "R", paste0("package:", packages), paste0("arg:", names(arg_values))),
    value = c(script_name, format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), R.version.string, versions, arg_values)
  )
  write_tsv_atomic(df, path, overwrite)
}

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

write_tsv_atomic(
  data.frame(barcode = barcodes, size_factor = size_factors, cluster = as.character(clusters)),
  args$out, args$overwrite
)
write_provenance(paste0(args$out, ".provenance.tsv"), args, c("scran", "Matrix"), args$overwrite)
