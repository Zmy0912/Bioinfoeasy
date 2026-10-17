# Changelog

## v2.0 (2026-09-21)

### Major Improvements

#### 1. Smart Detection of Ka/Ks == 0 as Missing
**Problem:** In v1.0 the Ka/Ks median computation included 0 values, so
the imputed Ka/Ks remained 0.

**Root cause:** KaKs Calculator outputs 0 as a placeholder when Ka
cannot be computed. This value is not biologically meaningful.

**Fix:**
- Missing-value definition updated: Ka/Ks column treats **both NaN and
  == 0** as missing.
- Median computation: **only over Ka/Ks > 0 valid values**.
- The report now adds `MissKa` / `MissKs` columns to count missing
  values per column separately.

**Effect** (taking the `accD` gene as an example):
- Old median: ≈ 0.27 (zeros dragged it down)
- New median: **1.252** (valid values only)
- Imputed result: no longer 0

#### 2. Median Always Computed (Even When Imputation Is Disabled)
The report's `MedKa` / `MedKs` / `MedKa/Ks` are now shown regardless of
whether the user ticks the corresponding imputation box, so users can
inspect them.

#### 3. Stronger Gene-name Normalization
- Strips the `_mafft` suffix ✓
- Strips the `_copy\d+` suffix ✓ (e.g. `ndhB_copy2_mafft` → `ndhB`)
- Merges `_copy` variants (e.g. `ndhB` and `ndhB_copy2` are treated as
  the same gene)

**Effect:** the raw dataset's 86 "genes" are merged into 79 unique
genes, matching the organized xlsx.

#### 4. Bundled Example Files
- `example_kaks.xlsx`: 8 genes × 41 species organized sheet
- `example_kaks.kaks`: same content in raw KaKs output format
- The "Use Example" button now loads `example_kaks.xlsx` by default

### Other Improvements

- Report fields extended with `MissKa` / `MissKs` / `FilledCount` for
  more complete information
- Text report ends with field descriptions
- Added the "Open Output Folder" button
- GUI default title bumped to v2.0
- Added `requirements.txt`
- Polished `README.md` / `USAGE.md`

---

## v1.0 (2026-09-21)

### First Release

- Graphical interface (tkinter)
- Supports both `.xlsx` organized sheets and raw `.kaks` files
- 4 preset thresholds (20% / 30% / 50% / 70%)
- Median imputation (per gene)
- 3 output files (data + report + summary)
- Automatic deduplication of `(Sp1, Sp2, Gene)` duplicate rows