# Agentic RNA-Seq Differential Expression Assistant

An interactive Streamlit web application and DESeq2 analysis workflow for RNA-Seq differential expression analysis, volcano plot visualization, and automated interpretation reporting.

## Features
- Upload Counts CSV and Metadata CSV or run built-in preset datasets
- Controlled DESeq2 execution with a safe gene-wise dispersion fallback for narrow-dispersion datasets
- Configurable padj and |log2FoldChange| cutoffs
- Automated biological summary, top up/downregulated genes, volcano plots, and downloadable CSV/report outputs

## Run Locally
./run_app.sh
