#!/usr/bin/env Rscript
# Test helper: DESeq2 reference results for the pydeseq2 concordance test (RNASEQ_BULK_GUIDELINES.md §8).
# Usage: deseq2_reference.R counts.tsv metadata.tsv out.tsv  (counts: samples x genes; metadata: sample, batch,
# condition). Design ~ batch + condition, contrast condition case vs control, no shrinkage.
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, stringsAsFactors = TRUE)
meta$condition <- relevel(meta$condition, ref = "control")
dds <- DESeqDataSetFromMatrix(t(counts), colData = meta[rownames(counts), ], design = ~ batch + condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "case", "control"))
out <- data.frame(gene_id = rownames(res), log2_fold_change = res$log2FoldChange, p_value = res$pvalue,
                  p_adjusted = res$padj)
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
