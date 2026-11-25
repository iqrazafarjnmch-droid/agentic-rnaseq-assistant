import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def build_agent_interpretation(df: pd.DataFrame, padj_cutoff: float = 0.05, lfc_cutoff: float = 1.0) -> str:
    total_genes = len(df)
    valid = df.dropna(subset=["log2FoldChange", "padj"]).copy()

    sig = valid[(valid["padj"] < padj_cutoff) & (valid["log2FoldChange"].abs() >= lfc_cutoff)].copy()
    up = sig[sig["log2FoldChange"] >= lfc_cutoff].sort_values(by=["padj", "log2FoldChange"], ascending=[True, False])
    down = sig[sig["log2FoldChange"] <= -lfc_cutoff].sort_values(by=["padj", "log2FoldChange"], ascending=[True, True])

    up_count = len(up)
    down_count = len(down)
    sig_count = len(sig)

    if up_count > down_count * 1.3:
        direction_pattern = "predominantly upregulated (activation-leaning response in Treated vs Control)"
    elif down_count > up_count * 1.3:
        direction_pattern = "predominantly downregulated (repression-leaning response in Treated vs Control)"
    else:
        direction_pattern = "balanced bidirectional response (comparable upregulation and downregulation)"

    lines = [
        "============================================================",
        "AGENTIC RNA-SEQ DIFFERENTIAL EXPRESSION INTERPRETATION REPORT",
        "============================================================",
        f"Total Genes Evaluated       : {total_genes}",
        f"Genes with Valid Stats      : {len(valid)}",
        f"Thresholds                  : padj < {padj_cutoff} and |log2FoldChange| >= {lfc_cutoff}",
        f"Significant DEGs Found      : {sig_count} (Upregulated: {up_count} | Downregulated: {down_count})",
        "",
        "1. EXPRESSION PROFILE SUMMARY",
        f"- The Treated vs. Control comparison shows a {direction_pattern}.",
        f"- Out of {total_genes} evaluated genes, {sig_count} ({(sig_count / max(total_genes, 1)) * 100:.1f}%) cross both significance and fold-change thresholds.",
        "",
        "2. TOP UPREGULATED CANDIDATES (Higher in Treated)",
    ]

    if up_count > 0:
        for _, row in up.head(5).iterrows():
            lines.append(f"  - {row['gene']}: log2FC = {row['log2FoldChange']:.3f}, padj = {row['padj']:.3e}")
    else:
        lines.append("  - None meeting threshold.")

    lines.append("")
    lines.append("3. TOP DOWNREGULATED CANDIDATES (Lower in Treated)")
    if down_count > 0:
        for _, row in down.head(5).iterrows():
            lines.append(f"  - {row['gene']}: log2FC = {row['log2FoldChange']:.3f}, padj = {row['padj']:.3e}")
    else:
        lines.append("  - None meeting threshold.")

    lines.extend([
        "",
        "4. AGENTIC BIOLOGICAL & NEXT-STEP HYPOTHESES",
        "- Prioritize top ranked genes with low padj and high |log2FC| for gene-set/pathway enrichment (GO / KEGG).",
        "- Check dispersion and count distributions before selecting biomarker/target candidates.",
        "- Validate key up/downregulated targets with independent assays (e.g., qRT-PCR or external cohort).",
        "============================================================",
    ])
    return "\n".join(lines)


def make_volcano_plot(df: pd.DataFrame, out_path: str, padj_cutoff: float = 0.05, lfc_cutoff: float = 1.0) -> None:
    plot_df = df.dropna(subset=["log2FoldChange", "padj"]).copy()
    plot_df["neg_log10_padj"] = -np.log10(np.clip(plot_df["padj"].astype(float), 1e-300, 1.0))

    is_up = (plot_df["padj"] < padj_cutoff) & (plot_df["log2FoldChange"] >= lfc_cutoff)
    is_down = (plot_df["padj"] < padj_cutoff) & (plot_df["log2FoldChange"] <= -lfc_cutoff)
    is_ns = ~(is_up | is_down)

    plt.figure(figsize=(8, 6))
    plt.scatter(plot_df.loc[is_ns, "log2FoldChange"], plot_df.loc[is_ns, "neg_log10_padj"],
                c="#9aa0a6", alpha=0.55, s=18, label=f"Non-significant ({is_ns.sum()})")
    plt.scatter(plot_df.loc[is_up, "log2FoldChange"], plot_df.loc[is_up, "neg_log10_padj"],
                c="#d93025", alpha=0.85, s=22, label=f"Upregulated ({is_up.sum()})")
    plt.scatter(plot_df.loc[is_down, "log2FoldChange"], plot_df.loc[is_down, "neg_log10_padj"],
                c="#1a73e8", alpha=0.85, s=22, label=f"Downregulated ({is_down.sum()})")

    plt.axvline(lfc_cutoff, color="#444444", linestyle="--", linewidth=1)
    plt.axvline(-lfc_cutoff, color="#444444", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(padj_cutoff), color="#444444", linestyle="--", linewidth=1)

    plt.title("RNA-Seq Volcano Plot (Treated vs. Control)")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 Adjusted P-value")
    plt.legend(loc="best", frameon=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def make_ma_plot(df: pd.DataFrame, out_path: str, padj_cutoff: float = 0.05, lfc_cutoff: float = 1.0) -> None:
    plot_df = df.dropna(subset=["baseMean", "log2FoldChange", "padj"]).copy()
    plot_df = plot_df[plot_df["baseMean"] > 0]
    plot_df["log10_baseMean"] = np.log10(plot_df["baseMean"].astype(float))

    is_sig = (plot_df["padj"] < padj_cutoff) & (plot_df["log2FoldChange"].abs() >= lfc_cutoff)

    plt.figure(figsize=(8, 6))
    plt.scatter(plot_df.loc[~is_sig, "log10_baseMean"], plot_df.loc[~is_sig, "log2FoldChange"],
                c="#9aa0a6", alpha=0.5, s=16, label=f"Non-significant ({(~is_sig).sum()})")
    plt.scatter(plot_df.loc[is_sig, "log10_baseMean"], plot_df.loc[is_sig, "log2FoldChange"],
                c="#e37400", alpha=0.85, s=22, label=f"Significant DEGs ({is_sig.sum()})")

    plt.axhline(0, color="#222222", linestyle="-", linewidth=0.9)
    plt.axhline(lfc_cutoff, color="#444444", linestyle="--", linewidth=1)
    plt.axhline(-lfc_cutoff, color="#444444", linestyle="--", linewidth=1)

    plt.title("RNA-Seq MA Plot (Mean Expression vs. Log2 Fold Change)")
    plt.xlabel("Log10 Mean Normalized Counts (baseMean)")
    plt.ylabel("Log2 Fold Change")
    plt.legend(loc="best", frameon=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def make_top_degs_plot(df: pd.DataFrame, out_path: str, padj_cutoff: float = 0.05, lfc_cutoff: float = 1.0) -> None:
    sig = df.dropna(subset=["log2FoldChange", "padj"]).copy()
    sig = sig[(sig["padj"] < padj_cutoff) & (sig["log2FoldChange"].abs() >= lfc_cutoff)]
    sig["abs_lfc"] = sig["log2FoldChange"].abs()
    top = sig.sort_values(by=["padj", "abs_lfc"], ascending=[True, False]).head(10)
    if len(top) == 0:
        return
    top = top.iloc[::-1]
    colors = ["#d93025" if x >= 0 else "#1a73e8" for x in top["log2FoldChange"]]

    plt.figure(figsize=(8, 5.5))
    plt.barh(top["gene"].astype(str), top["log2FoldChange"], color=colors)
    plt.axvline(0, color="#222222", linewidth=0.9)
    plt.title("Top 10 Significant Differentially Expressed Genes")
    plt.xlabel("Log2 Fold Change (Treated vs. Control)")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main() -> None:
    project_dir = os.path.expanduser("~/Minor_Project_RNASeq")
    results_csv = sys.argv[1] if len(sys.argv) >= 2 else os.path.join(project_dir, "results", "deseq2_results.csv")
    out_dir = sys.argv[2] if len(sys.argv) >= 3 else os.path.dirname(results_csv)
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(results_csv)
    report_text = build_agent_interpretation(df)

    report_path = os.path.join(out_dir, "automated_interpretation.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_text + "\n")

    make_volcano_plot(df, os.path.join(out_dir, "volcano_plot.png"))
    make_ma_plot(df, os.path.join(out_dir, "ma_plot.png"))
    make_top_degs_plot(df, os.path.join(out_dir, "top_degs_plot.png"))

    print(report_text)
    print(f"\nSaved report + plots -> {out_dir}")


if __name__ == "__main__":
    main()
