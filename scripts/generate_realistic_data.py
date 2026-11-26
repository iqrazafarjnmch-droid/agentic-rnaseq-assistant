import os
from pathlib import Path
import numpy as np
import pandas as pd

def generate_realistic_dataset(
    output_dir: str = "data",
    n_genes: int = 500,
    n_per_group: int = 3,
    seed: int = 42,
) -> tuple[Path, Path]:
    rng = np.random.default_rng(seed)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    genes = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
    samples = [f"CTRL_{i+1}" for i in range(n_per_group)] + [f"TRT_{i+1}" for i in range(n_per_group)]
    conditions = ["control"] * n_per_group + ["treated"] * n_per_group

    base_means = rng.integers(25, 1500, size=n_genes).astype(float)
    matrix = np.zeros((n_genes, len(samples)), dtype=int)

    up_idx = np.arange(0, min(30, n_genes))
    down_idx = np.arange(min(30, n_genes), min(60, n_genes))

    for g in range(n_genes):
        mu = base_means[g]
        for s in range(len(samples)):
            sample_mu = mu
            if conditions[s] == "treated":
                if g in up_idx:
                    sample_mu = mu * rng.uniform(2.5, 6.0)
                elif g in down_idx:
                    sample_mu = mu * rng.uniform(0.15, 0.45)
            dispersion = 0.12
            shape = 1.0 / dispersion
            scale = sample_mu * dispersion
            count_val = int(rng.gamma(shape=shape, scale=scale))
            matrix[g, s] = max(0, count_val)

    counts_df = pd.DataFrame(matrix, columns=samples)
    counts_df.insert(0, "gene_id", genes)

    meta_df = pd.DataFrame({"sample": samples, "condition": conditions})

    counts_file = out_path / "real_like_counts.csv"
    meta_file = out_path / "real_like_metadata.csv"
    counts_df.to_csv(counts_file, index=False)
    meta_df.to_csv(meta_file, index=False)

    print(f"Generated realistic dataset: {counts_df.shape[0]} genes x {len(samples)} samples")
    return counts_file, meta_file

if __name__ == "__main__":
    generate_realistic_dataset()
