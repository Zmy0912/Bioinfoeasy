# -*- coding: utf-8 -*-
"""Generate the example files used to demonstrate the KaKs analysis tool."""
import sys
import pandas as pd
from pathlib import Path

# Default to the directory containing this script so the example files are
# written next to it. Override OUT_DIR / sys.path here if your raw data
# lives somewhere else.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from kaks_analyzer import KaksProcessor

OUT_DIR = HERE  # Change this if your raw data is elsewhere.

# Load the raw data.
src = OUT_DIR / "merge.axt.kaks"
df = KaksProcessor.load_file(str(src))
print(f"Raw data: {len(df)} rows, {df['Gene'].nunique()} genes, "
      f"{df['Species2'].nunique()} species")

# Take the first 8 genes for the example (compact yet spans many
# different missing ratios).
sample_genes = sorted(df['Gene'].unique())[:8]
sample = df[df['Gene'].isin(sample_genes)].copy()
print(f"Example genes ({len(sample_genes)}): {sample_genes}")
print(f"Example rows: {len(sample)}")

# Simplify Species1/Species2: drop numeric ID suffixes and keep only
# genus + species (e.g. Asparagus_acutifolius).


def short_name(s):
    parts = s.split('_')
    # Keep only the first two tokens (genus + species).
    return '_'.join(parts[:2]) if len(parts) >= 2 else s


sample['Species1'] = sample['Species1'].apply(short_name)
sample['Species2'] = sample['Species2'].apply(short_name)

# Restructure to 6 columns: Species1, Species2, Gene, Ka, Ks, Ka/Ks
out = sample[['Species1', 'Species2', 'Gene', 'Ka', 'Ks', 'Ka/Ks']].copy()
# Sort species alphabetically within each gene for readability.
out = out.sort_values(['Gene', 'Species1', 'Species2']).reset_index(drop=True)

# Write xlsx.
out_xlsx = OUT_DIR / "example_kaks.xlsx"
out.to_excel(out_xlsx, index=False)
print(f"Wrote: {out_xlsx}")

# Also export the original .kaks format (the 22-column raw table).
out_kaks = OUT_DIR / "example_kaks.kaks"
raw_df = pd.read_csv(src, sep='\t')
# Keep only the example genes (matching by their original gene names
# ending in _mafft).
raw_sample_genes = sorted(raw_df['Sequence'].apply(
    lambda s: str(s).split('|')[-1].replace('_mafft', '').strip()
).unique())[:8]
sample_raw = raw_df[raw_df['Sequence'].apply(
    lambda s: str(s).split('|')[-1].replace('_mafft', '').strip()
).isin(raw_sample_genes)].copy()
sample_raw.to_csv(out_kaks, sep='\t', index=False)
print(f"Wrote: {out_kaks} (raw .kaks format)")

# Summary of the example dataset.
print()
print("Ka/Ks missing ratio per gene in the example dataset:")
sample["__miss__"] = sample['Ka/Ks'].isna() | (sample['Ka/Ks'] == 0)
miss_rate = sample.groupby('Gene')['__miss__'].mean()
print(miss_rate.describe())
print()
print("Ka/Ks == 0 count:", (sample['Ka/Ks'] == 0).sum())
print("Ka/Ks NaN count :", sample['Ka/Ks'].isna().sum())