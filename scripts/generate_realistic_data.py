import os
import numpy as np
import pandas as pd

rng = np.random.default_rng(42)
project_dir = os.path.expanduser("~/Minor_Project_RNASeq")
data_dir = os.path.join(project_dir, "data")
os.makedirs(data_dir, exist_ok=True)

n_genes = 800
samples = ["CTRL_1", "CTRL_2", "CTRL_3", "TREAT_1", "TREAT_2", "TREAT_3"]
conditions = ["Control", "Control", "Control", "Treated", "Treated", "Treated"]
size_factors = np.array([0.75, 1.35, 0.95, 1.40, 0.80, 1.15])

genes = [f"ENSG{i:07d}" for i in range(1, n_genes + 1)]
base_means = rng.lognormal(mean=4.2, sigma=1.6, size=n_genes)
base_means[:40] = 0.0  # all-zero noise genes
base_means[40:120] = rng.uniform(0.1, 1.2, size=80)  # low-count genes

fc = np.ones(n_genes)
fc[120:155] = rng.uniform(2.2, 4.5, size=35)   # upregulated in Treated
fc[155:190] = rng.uniform(0.22, 0.45, size=35) # downregulated in Treated

counts = np.zeros((n_genes, len(samples)), dtype=int)
dispersion = 0.08

for j, sf in enumerate(size_factors):
    group_mult = fc if conditions[j] == "Treated" else np.ones(n_genes)
    mu = np.maximum(base_means * group_mult * sf, 1e-6)
    shape = 1.0 / dispersion
    scale = mu / shape
    counts[:, j] = rng.gamma(shape=shape, scale=scale)
    counts[:, j] = rng.poisson(counts[:, j])

counts_df = pd.DataFrame(counts, index=genes, columns=samples)
counts_df.index.name = "gene_id"
counts_df.to_csv(os.path.join(data_dir, "real_like_counts.csv"))

meta_df = pd.DataFrame({"sample": samples, "condition": conditions})
meta_df.to_csv(os.path.join(data_dir, "real_like_metadata.csv"), index=False)

print("Generated realistic dataset:", counts_df.shape)
