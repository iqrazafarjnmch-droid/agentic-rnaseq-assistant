import os
import numpy as np
import pandas as pd

rng = np.random.default_rng(108)
project_dir = os.path.expanduser("~/Minor_Project_RNASeq")
data_dir = os.path.join(project_dir, "data")
os.makedirs(data_dir, exist_ok=True)

# 1500-gene multi-replicate biological stress dataset (4 Control vs 4 Treated)
n_genes = 1500
samples = ["S1_CTRL", "S2_CTRL", "S3_CTRL", "S4_CTRL", "S5_TREAT", "S6_TREAT", "S7_TREAT", "S8_TREAT"]
conditions = ["Control", "Control", "Control", "Control", "Treated", "Treated", "Treated", "Treated"]
size_factors = np.array([0.70, 1.45, 0.90, 1.25, 0.78, 1.50, 1.05, 0.85])

genes = [f"GENE_BIO_{i:04d}" for i in range(1, n_genes + 1)]
base_means = rng.lognormal(mean=4.0, sigma=1.7, size=n_genes)
base_means[:100] = 0.0
base_means[100:260] = rng.uniform(0.05, 1.5, size=160)

fc = np.ones(n_genes)
fc[260:320] = rng.uniform(2.0, 5.0, size=60)
fc[320:380] = rng.uniform(0.20, 0.50, size=60)

counts = np.zeros((n_genes, len(samples)), dtype=int)
dispersion = 0.10

for j, sf in enumerate(size_factors):
    mult = fc if conditions[j] == "Treated" else np.ones(n_genes)
    mu = np.maximum(base_means * mult * sf, 1e-6)
    shape = 1.0 / dispersion
    gamma_rates = rng.gamma(shape=shape, scale=mu / shape)
    counts[:, j] = rng.poisson(gamma_rates)

counts_df = pd.DataFrame(counts, index=genes, columns=samples)
counts_df.index.name = "gene_id"
counts_df.to_csv(os.path.join(data_dir, "stress_counts.csv"))

meta_df = pd.DataFrame({"sample": samples, "condition": conditions})
meta_df.to_csv(os.path.join(data_dir, "stress_metadata.csv"), index=False)

print("Created 1500-gene stress test dataset:", counts_df.shape)
