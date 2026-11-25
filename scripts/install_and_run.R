lib <- "//wsl.localhost/Ubuntu/home/iqra/R_libs"
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
}
BiocManager::install("DESeq2", lib = lib, update = FALSE, ask = FALSE)

suppressPackageStartupMessages(library(DESeq2))
counts <- read.csv("//wsl.localhost/Ubuntu/home/iqra/Minor_Project_RNASeq/data/sample_counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("//wsl.localhost/Ubuntu/home/iqra/Minor_Project_RNASeq/data/sample_metadata.csv", row.names = 1)
coldata$condition <- factor(coldata$condition, levels = c("Control", "Treated"))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- tryCatch(
  DESeq(dds),
  error = function(e) {
    dds_f <- estimateSizeFactors(dds)
    dds_f <- estimateDispersionsGeneEst(dds_f)
    dispersions(dds_f) <- mcols(dds_f)$dispGeneEst
    nbinomWaldTest(dds_f)
  }
)

res <- results(dds)
res_df <- as.data.frame(res[order(res$padj), ])
res_df$gene <- rownames(res_df)

write.csv(res_df, "//wsl.localhost/Ubuntu/home/iqra/Minor_Project_RNASeq/results/deseq2_results.csv", row.names = FALSE)
cat("DONE_DESEQ2\n")
