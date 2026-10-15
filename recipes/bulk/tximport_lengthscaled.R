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
#   <out>.provenance.tsv  R and package versions, arguments and date.
#   --overwrite        Replace existing outputs (otherwise the script stops if they exist).

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
})

options <- list(
  make_option("--samples", type = "character"),
  make_option("--tx2gene", type = "character"),
  make_option("--mode", type = "character", default = "lengthScaledTPM"),
  make_option("--out", type = "character"),
  make_option("--overwrite", action = "store_true", default = FALSE)
)
args <- parse_args(OptionParser(option_list = options))
script_name <- "tximport_lengthscaled.R"

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
stopifnot(args$mode %in% c("lengthScaledTPM", "offsets"))

samples <- read.delim(args$samples, stringsAsFactors = FALSE)
stopifnot(all(c("sample", "quant") %in% names(samples)), !anyDuplicated(samples$sample), all(file.exists(samples$quant)))
files <- setNames(samples$quant, samples$sample)
tx2gene <- read.delim(args$tx2gene, stringsAsFactors = FALSE)[, c("transcript_id", "gene_id")]

counts_from_abundance <- if (args$mode == "lengthScaledTPM") "lengthScaledTPM" else "no"
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = counts_from_abundance)

write_matrix <- function(m, path) {
  out <- data.frame(gene_id = sub("\\.[0-9]+", "", rownames(m)), m, check.names = FALSE)
  write_tsv_atomic(out, path, args$overwrite)
}
write_matrix(txi$counts, paste0(args$out, ".counts.tsv"))
write_matrix(txi$length, paste0(args$out, ".length.tsv"))
write_matrix(txi$abundance, paste0(args$out, ".tpm.tsv"))
write_provenance(paste0(args$out, ".provenance.tsv"), args, c("tximport"), args$overwrite)
