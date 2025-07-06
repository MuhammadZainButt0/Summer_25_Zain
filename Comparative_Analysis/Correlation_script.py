import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys

# === File paths ===
snp_file = 'SNP_density.csv'  # Updated to match previous script output
expression_file = 'merged_gene_expression.csv'  # Replace with your file

# === Load data ===
try:
    snp_df = pd.read_csv(snp_file)
    expression_df = pd.read_csv(expression_file)
except FileNotFoundError as e:
    print(f"File not found: {e}")
    sys.exit(1)

# === Print column names (debugging) ===
print("SNP columns:", snp_df.columns.tolist())
print("Expression columns:", expression_df.columns.tolist())

# === Standardize column names if needed ===
snp_df.rename(columns={'gene': 'Gene'}, inplace=True)
expression_df.rename(columns={'GeneID': 'gene_id'}, inplace=True)

# === Merge datasets ===
merged_df = pd.merge(snp_df, expression_df, how='inner', left_on='Gene', right_on='gene_id')

if merged_df.empty:
    print("Merge failed: No matching genes found. Check gene ID consistency.")
    sys.exit(1)

# === Compute log-transformed average expression ===
merged_df["avg_expr"] = merged_df[["fpkm_Sample1", "fpkm_Sample2"]].mean(axis=1)
merged_df["log_avg_expr"] = np.log1p(merged_df["avg_expr"])

# === Pearson correlation for each sample ===
try:
    corr1, pval1 = pearsonr(merged_df['Density_SNPs_per_kb'], merged_df['fpkm_Sample1'])
    corr2, pval2 = pearsonr(merged_df['Density_SNPs_per_kb'], merged_df['fpkm_Sample2'])
except KeyError as e:
    print(f"Column missing: {e}")
    sys.exit(1)

# === Pearson correlation with log-transformed average expression ===
avg_corr, avg_pval = pearsonr(merged_df['Density_SNPs_per_kb'], merged_df['log_avg_expr'])

# === Statistic Summary of SNP Density ===
print("\nStatistic Summary of SNP Density")
print("SNP Density Summary:")
print(snp_df['Density_SNPs_per_kb'].describe())

# === Print all results ===
print("\nPearson Correlation with SNP Density")
print(f"  fpkm_Sample1:      r = {corr1:.4f}, p = {pval1:.4e}")
print(f"  fpkm_Sample2:      r = {corr2:.4f}, p = {pval2:.4e}")
print(f"  log(avg_expr):     r = {avg_corr:.4f}, p = {avg_pval:.4e}")

# === Save merged data ===
merged_df.to_csv("merged_snp_expression.csv", index=False)
print("Merged data saved as 'merged_snp_expression.csv'")

# === Scatter Plot ===
plt.figure(figsize=(8, 5))
plt.scatter(merged_df['Density_SNPs_per_kb'], merged_df['fpkm_Sample1'], label='Sample1', alpha=0.6)
plt.scatter(merged_df['Density_SNPs_per_kb'], merged_df['fpkm_Sample2'], label='Sample2', alpha=0.6)
plt.xlabel('SNP Density (per kb)')
plt.ylabel('FPKM')
plt.title('SNP Density vs Expression')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('snp_expression_scatter.png')
plt.show()
print("Scatter plot saved as 'snp_expression_scatter.png'")

# === Histogram of SNP Density ===
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.hist(merged_df["Density_SNPs_per_kb"], bins=30, color='skyblue', edgecolor='black')
plt.title("Histogram of SNP Density")
plt.xlabel("SNPs per kb")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig('snp_density_histogram.png')
print("Histogram saved as 'snp_density_histogram.png'")

# === Bar Plot: Top 20 Genes by SNP Density ===
top_genes = merged_df.sort_values(by="Density_SNPs_per_kb", ascending=False).head(20)
plt.figure(figsize=(12, 6))
plt.bar(top_genes["Gene"], top_genes["Density_SNPs_per_kb"], color='teal')
plt.xticks(rotation=90)
plt.title("Top 20 Genes with Highest SNP Density")
plt.xlabel("Gene")
plt.ylabel("SNPs per kb")
plt.tight_layout()
plt.savefig("top20_snp_density_barplot.png")
plt.show()
print("Bar plot saved as 'top20_snp_density_barplot.png'")

# === Save histogram ===
