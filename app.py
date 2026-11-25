import re
import subprocess
from pathlib import Path

import pandas as pd
import streamlit as st

from scripts.ai_rnaseq_agent import (
    build_agent_interpretation,
    make_ma_plot,
    make_top_degs_plot,
    make_volcano_plot,
)

PROJECT_DIR = Path(__file__).resolve().parent

import subprocess, shutil
if shutil.which("Rscript") and (PROJECT_DIR / "setup.R").exists():
    check = subprocess.run(["Rscript", "-e", "quit(status = ifelse(requireNamespace('DESeq2', quietly=TRUE), 0, 1))"])
    if check.returncode != 0:
        subprocess.run(["Rscript", str(PROJECT_DIR / "setup.R")], check=False)

DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
RSCRIPT_EXE = "/mnt/c/Program Files/R/R-4.6.2/bin/x64/Rscript.exe"
DESEQ2_SCRIPT = PROJECT_DIR / "scripts" / "run_deseq2.R"

st.set_page_config(
    page_title="Agentic RNA-Seq Assistant",
    page_icon="🧬",
    layout="wide",
)


def validate_inputs(counts_df: pd.DataFrame, meta_df: pd.DataFrame) -> tuple[bool, str]:
    if counts_df.empty or meta_df.empty:
        return False, "Counts or metadata file is empty."
    if len(counts_df.columns) < 2:
        return False, "Counts matrix must contain a gene ID column and at least two sample columns."
    if not {"sample", "condition"}.issubset(set(meta_df.columns)):
        return False, "Metadata file must include 'sample' and 'condition' columns."

    sample_cols = list(counts_df.columns[1:])
    meta_samples = list(meta_df["sample"].astype(str))
    if set(sample_cols) != set(meta_samples):
        return False, f"Metadata samples ({meta_samples}) do not match count matrix columns ({sample_cols})."

    num_part = counts_df[sample_cols].apply(pd.to_numeric, errors="coerce")
    if num_part.isna().any().any():
        return False, "Counts matrix contains non-numeric or missing values."
    if (num_part.values < 0).any():
        return False, "Counts matrix cannot contain negative values."

    conds = meta_df["condition"].astype(str).unique().tolist()
    if len(conds) < 2:
        return False, "At least two experimental conditions (e.g., Control vs Treated) are required for differential expression."

    return True, "Validation passed: input files, sample IDs, and count values are valid."


def run_analysis(counts_path: Path, meta_path: Path, out_dir: Path, padj_cut: float, lfc_cut: float) -> tuple[bool, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [RSCRIPT_EXE, str(DESEQ2_SCRIPT), str(counts_path), str(meta_path), str(out_dir)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        err = res.stderr.strip() or res.stdout.strip() or "DESeq2 execution failed."
        return False, err

    csv_path = out_dir / "deseq2_results.csv"
    if not csv_path.exists():
        return False, "deseq2_results.csv was not generated."

    df = pd.read_csv(csv_path)
    report = build_agent_interpretation(df, padj_cutoff=padj_cut, lfc_cutoff=lfc_cut)
    (out_dir / "automated_interpretation.txt").write_text(report + "\n", encoding="utf-8")
    make_volcano_plot(df, str(out_dir / "volcano_plot.png"), padj_cutoff=padj_cut, lfc_cutoff=lfc_cut)
    make_ma_plot(df, str(out_dir / "ma_plot.png"), padj_cutoff=padj_cut, lfc_cutoff=lfc_cut)
    make_top_degs_plot(df, str(out_dir / "top_degs_plot.png"), padj_cutoff=padj_cut, lfc_cutoff=lfc_cut)
    return True, report


def get_dataset_paths(preset: str, uploaded_counts, uploaded_meta, out_dir: Path) -> tuple[Path, Path, str]:
    if preset == "Custom Upload" and uploaded_counts is not None and uploaded_meta is not None:
        c_path = out_dir / "uploaded_counts.csv"
        m_path = out_dir / "uploaded_metadata.csv"
        out_dir.mkdir(parents=True, exist_ok=True)
        c_path.write_bytes(uploaded_counts.getvalue())
        m_path.write_bytes(uploaded_meta.getvalue())
        return c_path, m_path, "custom_upload"
    if preset == "real_like (800 genes)":
        return DATA_DIR / "real_like_counts.csv", DATA_DIR / "real_like_metadata.csv", "real_like"
    if preset == "stress_test (1,500 genes)":
        return DATA_DIR / "stress_counts.csv", DATA_DIR / "stress_metadata.csv", "stress_test"
    return DATA_DIR / "sample_counts.csv", DATA_DIR / "sample_metadata.csv", "default"


if "messages" not in st.session_state:
    st.session_state.messages = [
        {
            "role": "assistant",
            "content": (
                "Hello! I am your **Agentic RNA-Seq Assistant**. "
                "You can ask me to run workflows or inspect results—for example: *'run real_like dataset'*, "
                "*'run stress test'*, *'show top genes summary'*, or *'run with padj < 0.01 and log2FC >= 1.5'*."
            ),
        }
    ]
if "last_out_dir" not in st.session_state:
    st.session_state.last_out_dir = RESULTS_DIR / "real_like"

st.title("🧬 AI-Assisted Agentic RNA-Seq Analysis Tool")
st.caption("Controlled DESeq2 Execution, Automated Result Interpretation & Interactive Assistant")

with st.sidebar:
    st.header("Workflow Controls")
    preset = st.selectbox(
        "Dataset Preset",
        ["real_like (800 genes)", "default (sample)", "stress_test (1,500 genes)", "Custom Upload"],
    )
    uploaded_counts = st.file_uploader("Upload Counts CSV", type=["csv"])
    uploaded_meta = st.file_uploader("Upload Metadata CSV", type=["csv"])

    padj_cut = st.number_input(
        "Adjusted p-value cutoff (padj)",
        min_value=0.0001,
        max_value=0.2,
        value=0.05,
        step=0.01,
        format="%.4f",
    )
    lfc_cut = st.number_input("Minimum |log2FoldChange| cutoff", min_value=0.2, max_value=5.0, value=1.0, step=0.2)

    run_btn = st.button("▶ Run Controlled Pipeline", type="primary", width="stretch")

if run_btn:
    out_dir = RESULTS_DIR / ("ui_custom" if preset == "Custom Upload" else preset.split()[0])
    c_path, m_path, _ = get_dataset_paths(preset, uploaded_counts, uploaded_meta, out_dir)
    if not c_path.exists() or not m_path.exists():
        st.error("Selected dataset files were not found. Please upload both CSV files when using Custom Upload.")
    else:
        ok_val, val_msg = validate_inputs(pd.read_csv(c_path), pd.read_csv(m_path))
        if not ok_val:
            st.error(val_msg)
        else:
            with st.spinner("Running DESeq2 and generating agentic interpretation..."):
                ok_run, info = run_analysis(c_path, m_path, out_dir, padj_cut, lfc_cut)
            if ok_run:
                st.session_state.last_out_dir = out_dir
                st.success(f"Analysis completed successfully. Outputs saved to `{out_dir}`.")
            else:
                st.error(info)

st.subheader("💬 Assistant Chat")
for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])

prompt = st.chat_input("Enter a command or question (e.g., 'run real_like' or 'top genes summary')...")
if prompt:
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    low = prompt.lower()
    use_padj = float(padj_cut)
    use_lfc = float(lfc_cut)

    m_padj = re.search(r"padj\s*[<>=]*\s*([0-9.]+)", low)
    if m_padj:
        use_padj = float(m_padj.group(1))
    m_lfc = re.search(r"(?:lfc|log2fc)\s*[<>=]*\s*([0-9.]+)", low)
    if m_lfc:
        use_lfc = float(m_lfc.group(1))

    if any(k in low for k in ["run", "execute", "start", "analysis"]):
        if "stress" in low:
            chosen = "stress_test"
            c_p, m_p = DATA_DIR / "stress_counts.csv", DATA_DIR / "stress_metadata.csv"
        elif "default" in low or "sample" in low:
            chosen = "default"
            c_p, m_p = DATA_DIR / "sample_counts.csv", DATA_DIR / "sample_metadata.csv"
        else:
            chosen = "real_like"
            c_p, m_p = DATA_DIR / "real_like_counts.csv", DATA_DIR / "real_like_metadata.csv"

        out_dir = RESULTS_DIR / chosen
        ok_val, val_msg = validate_inputs(pd.read_csv(c_p), pd.read_csv(m_p))
        if not ok_val:
            reply = f"Validation failed: {val_msg}"
        else:
            ok_run, info = run_analysis(c_p, m_p, out_dir, use_padj, use_lfc)
            if ok_run:
                st.session_state.last_out_dir = out_dir
                reply = (
                    f"Completed **`{chosen}`** workflow (`padj < {use_padj}`, `|log2FC| >= {use_lfc}`).\n\n"
                    f"```text\n{info}\n```"
                )
            else:
                reply = f"Execution error: {info}"
    else:
        res_csv = st.session_state.last_out_dir / "deseq2_results.csv"
        if res_csv.exists():
            df_cur = pd.read_csv(res_csv)
            reply = (
                f"Summary for current result (`{st.session_state.last_out_dir.name}`):\n\n"
                f"```text\n{build_agent_interpretation(df_cur, use_padj, use_lfc)}\n```"
            )
        else:
            reply = "No results loaded yet. Click **Run Controlled Pipeline** in the sidebar or type *'run real_like'*."

    st.session_state.messages.append({"role": "assistant", "content": reply})
    with st.chat_message("assistant"):
        st.markdown(reply)

st.divider()
st.subheader("📊 Results, Plots & Downloads")
cur_dir = st.session_state.last_out_dir
rep_file = cur_dir / "automated_interpretation.txt"
csv_file = cur_dir / "deseq2_results.csv"

if rep_file.exists() and csv_file.exists():
    st.text_area("Automated Interpretation Report", rep_file.read_text(encoding="utf-8"), height=260)

    c1, c2, c3 = st.columns(3)
    for col, fname, caption in [
        (c1, "volcano_plot.png", "Volcano Plot"),
        (c2, "ma_plot.png", "MA Plot"),
        (c3, "top_degs_plot.png", "Top DEGs"),
    ]:
        fpath = cur_dir / fname
        if fpath.exists():
            col.image(str(fpath), caption=caption, width="stretch")

    df_res = pd.read_csv(csv_file)
    st.dataframe(df_res.head(25), width="stretch")

    d1, d2 = st.columns(2)
    d1.download_button(
        "⬇ Download DESeq2 CSV",
        data=csv_file.read_bytes(),
        file_name="deseq2_results.csv",
        mime="text/csv",
    )
    d2.download_button(
        "⬇ Download Interpretation TXT",
        data=rep_file.read_bytes(),
        file_name="automated_interpretation.txt",
        mime="text/plain",
    )
else:
    st.info("Run the pipeline from the sidebar or type 'run real_like' in chat to view outputs.")
