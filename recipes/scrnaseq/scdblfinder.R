#!/usr/bin/env Rscript
# Recipe: doublet scores and calls with scDblFinder [Germain 2021].
# Guideline: SCRNASEQ_GUIDELINES.md §4.3 (doublets, SC-QC-009). Called by recipes/scrnaseq/doublets.py through files;
# it can also be run on its own.
#
# Inputs
#   --counts    Matrix Market file, genes x cells, counts without ambient correction (SC-QC-012).
#   --barcodes  Text file, one barcode per line, in column order.
#   --genes     Text file, one gene identifier per line, in row order.
#   --samples   Text file, one capture label per line (same order as barcodes); doublets are simulated within each
#               capture.
#   --dbr       Optional expected doublet rate; without it the scDblFinder default is used.
#   --seed      Seed set before scDblFinder.
#   --overwrite Replace existing outputs (otherwise the script stops if they exist).
# Outputs
#   --out                    TSV with columns barcode, capture, score, class (singlet or doublet), in input order.
#   <out>.provenance.tsv     R and package versions, arguments and date.

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

options <- list(
  make_option("--counts", type = "character"),
  make_option("--barcodes", type = "character"),
  make_option("--genes", type = "character"),
  make_option("--samples", type = "character"),
  make_option("--dbr", type = "double", default = NA),
  make_option("--seed", type = "integer", default = 0L),
  make_option("--out", type = "character"),
  make_option("--overwrite", action = "store_true", default = FALSE)
)
args <- parse_args(OptionParser(option_list = options))
script_name <- "scdblfinder.R"

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
genes <- readLines(args$genes)
samples <- readLines(args$samples)
stopifnot(ncol(counts) == length(barcodes), nrow(counts) == length(genes), length(samples) == length(barcodes))
dimnames(counts) <- list(genes, barcodes)

sce <- SingleCellExperiment(assays = list(counts = counts))
dbr <- if (is.na(args$dbr)) NULL else args$dbr
set.seed(args$seed)
sce <- scDblFinder(sce, samples = samples, dbr = dbr)

write_tsv_atomic(
  data.frame(
    barcode = barcodes, capture = samples, score = sce$scDblFinder.score,
    class = as.character(sce$scDblFinder.class)
  ),
  args$out, args$overwrite
)
write_provenance(
  paste0(args$out, ".provenance.tsv"), args, c("scDblFinder", "SingleCellExperiment", "Matrix"), args$overwrite
)
