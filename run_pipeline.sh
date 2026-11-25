#!/usr/bin/env bash
set -e
PROJECT_DIR="$HOME/Minor_Project_RNASeq"
MODE="${1:-default}"
if [ "$MODE" = "real_like" ]; then
  COUNTS="$PROJECT_DIR/data/real_like_counts.csv"
  META="$PROJECT_DIR/data/real_like_metadata.csv"
  OUT="$PROJECT_DIR/results/real_like"
else
  COUNTS="$PROJECT_DIR/data/sample_counts.csv"
  META="$PROJECT_DIR/data/sample_metadata.csv"
  OUT="$PROJECT_DIR/results"
fi
mkdir -p "$OUT"
/mnt/c/Program\ Files/R/R-4.6.2/bin/x64/Rscript.exe "$PROJECT_DIR/scripts/run_deseq2.R" "$COUNTS" "$META" "$OUT"
python3 "$PROJECT_DIR/scripts/ai_rnaseq_agent.py" "$OUT/deseq2_results.csv" "$OUT"
ls -lh "$OUT/"
