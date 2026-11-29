# OmicMind

**OmicMind** is an interactive web application for RNA-Seq differential gene expression analysis powered by **Bioconductor's `DESeq2`** via R, featuring strict input validation, a reproducible synthetic dataset generator, interactive visualizations, and plain-language result summaries.

## Key Features
- **Input Validation:** Checks for empty files, missing required columns, blank/duplicate gene IDs, non-numeric or negative count values, sample mismatches between counts and metadata, and insufficient condition groups.
- **Reproducible Synthetic Dataset:** Generates realistic count and metadata CSV files (`real_like_counts.csv` and `real_like_metadata.csv`) with a fixed random seed (`42`), including control/treated replicates and up/downregulated gene patterns.
- **DESeq2 Workflow & Visualizations:** Runs differential expression analysis and displays results alongside Volcano plots, MA plots, top DEGs, and a plain-language **Quick Summary**.

## Input File Formats
1. **Counts CSV:** First column as gene identifier (e.g., `gene_id` or `gene`) followed by numeric raw count columns matching sample IDs.
2. **Metadata CSV:** Must include `sample` and `condition` columns (with at least two condition groups).

## Notes & Limitations
- Built-in synthetic data (`seed=42`) is intended for testing, demonstration, and workflow verification rather than biological discovery.
- Always inspect count normalization, dispersion diagnostics, and batch/confounding variables when analyzing real experimental datasets.
