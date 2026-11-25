lib_path <- "//wsl.localhost/Ubuntu/home/iqra/R_libs"
if (dir.exists(lib_path)) {
  .libPaths(c(lib_path, .libPaths()))
}

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- "/home/iqra/Minor_Project_RNASeq"
counts_file <- if (length(args) >= 1) args[1] else file.path(project_dir, "data", "sample_counts.csv")
meta_file   <- if (length(args) >= 2) args[2] else file.path(project_dir, "data", "sample_metadata.csv")
out_dir     <- if (length(args) >= 3) args[3] else file.path(project_dir, "results")

counts_df <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta_df   <- read.csv(meta_file, row.names = 1)

meta_df$condition <- factor(meta_df$condition)
counts_mat <- as.matrix(counts_df[, rownames(meta_df)])
mode(counts_mat) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta_df, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- tryCatch(
  DESeq(dds),
  error = function(e) {
    dds_f <- estimateSizeFactors(dds)
    dds_f <- estimateDispersionsGeneEst(dds_f)
    dispersions(dds_f) <- mcols(dds_f)$dispGeneEst
    dds_f <- nbinomWaldTest(dds_f)
    message("Used gene-wise dispersion fallback due to narrow dispersion range.")
    dds_f
  }
)

res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res_df, file.path(out_dir, "deseq2_results.csv"), row.names = FALSE)
cat("DONE_DESEQ2 ->", file.path(out_dir, "deseq2_results.csv"), "\n")
