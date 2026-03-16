# Pi Sliding Window Gene Annotation Program

## Introduction

This program analyzes genomic nucleotide diversity (Pi) sliding window data and automatically annotates genes that are significantly elevated relative to their surrounding regions. The program identifies significant genes by calculating Z-score of each gene relative to its local region and generates high-quality visualization charts.

## Main Features

- Parse Pi sliding window data files
- Parse GFF3 gene annotation files
- Identify genes significantly elevated relative to surrounding regions
- Three significance criteria:
  - **Relative significance**: Gene Pi value is higher than surrounding region average by a certain multiple of standard deviation
  - **Absolute significance**: Gene Pi value is higher than set absolute threshold
  - **Lower threshold**: Gene Pi value is lower than set lower threshold (lowest Pi value genes)
- Smart label layout: Automatically avoid label overlaps, assign different heights to adjacent gene labels
- Multi-condition gene priority sorting: Prioritize display of genes meeting multiple conditions
- Generate high-quality PNG charts and detailed analysis reports

## System Requirements

- Python 3.6 or higher
- Operating system: Windows, Linux, or macOS

## Installation

### Method 1: Using pip

```bash
pip install pandas matplotlib numpy
```

### Method 2: Using requirements.txt

```bash
pip install -r requirements.txt
```

## Input File Format

### Pi File Format

The Pi file should be in text format and contain following columns:

| Column | Description |
|--------|-------------|
| Window | Sliding window range, format like "1-10000" |
| Midpoint | Window midpoint position |
| Pi | Pi value (nucleotide diversity) |
| Theta | Theta value (optional) |

Example:

```
Window	Midpoint	Pi	Theta
1-10000	5000	0.0234	0.0256
10001-20000	15000	0.0287	0.0301
20001-30000	25000	0.0198	0.0212
```

### GFF3 File Format

Standard GFF3 format gene annotation file. The program will extract rows with `feature type` as `gene`.

Example:

```
##gff-version 3
chr1	Ensembl	gene	1000	5000	.	+	.	ID=gene-001;Name=GENE1
chr1	Ensembl	gene	8000	12000	.	-	.	ID=gene-002;Name=GENE2
```

## Usage

### Basic Usage

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

### Complete Parameter Description

| Parameter | Short | Type | Required | Default | Description |
|-----------|-------|------|----------|---------|-------------|
| --pi-file | -p | File path | Yes | - | Pi value file path |
| --gff3-file | -g | File path | Yes | - | GFF3 gene annotation file path |
| --output | -o | File path | No | pi_sliding_window_plot.png | Output image file name |
| --gene-fontsize | - | Integer | No | 12 | Gene annotation font size |
| --max-genes | - | Integer | No | 15 | Maximum number of genes to annotate |
| --window-size | - | Integer | No | 10000 | Window size for calculating local average (bp) |
| --significance-std | - | Float | No | 2.0 | Relative significance standard (standard deviation multiplier) |
| --absolute-threshold | - | Float | No | None | Absolute Pi value threshold |
| --lower-threshold | - | Float | No | None | Lower Pi value threshold (lowest Pi value genes) |

### Usage Examples

#### 1. Relative significance only (default)

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

#### 2. Using absolute threshold

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --absolute-threshold 0.05
```

#### 3. Using relative significance and absolute threshold together

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 2.0 --absolute-threshold 0.05
```

#### 4. Using lower threshold (detect lowest Pi value genes)

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --lower-threshold 0.01
```

#### 5. Using all three thresholds together

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --significance-std 2.0 \
  --absolute-threshold 0.05 \
  --lower-threshold 0.01
```

#### 6. Custom output file and font size

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 -o my_plot.png --gene-fontsize 14
```

#### 7. Adjust significance parameters

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 15000 --significance-std 2.5
```

#### 8. Limit number of annotated genes

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --max-genes 10
```

## Output Files

After running the program, the following files will be generated:

1. **PNG Image File**: Shows Pi sliding window curve and significant gene annotations
   - File name: Default is `pi_sliding_window_plot.png` or user-specified name
   - Resolution: 300 DPI
   - Format: PNG

2. **Information Text File**: Contains detailed analysis results of significant genes
   - File name: Same name as image file with `_info.txt` suffix
   - Content: Detailed statistical information for each significant gene

## Output Chart Description

### Chart Elements

- **Blue line**: Pi value changes along genome position
- **Red dashed line**: Relative significance threshold (90th percentile of Pi values)
- **Purple dashed line**: Absolute significance threshold (if specified)
- **Green dashed line**: Lower threshold (if specified)
- **Colored dots**: Positions of significant genes
  - Red: Relative significance genes
  - Purple: Absolute significance genes
  - Green: Lower threshold genes (lowest Pi value)
  - Dark red: Meet both relative and absolute significance
  - Other colors: Combinations meeting multiple conditions

### Gene Annotation Colors

- **Yellow background, red border**: Relative only `[R]`
- **Light blue background, purple border**: Absolute only `[A]`
- **Light green background, green border**: Lower threshold only `[L]`
- **Orange background, red border**: Meet both relative and absolute significance `[R+A]`
- **Light yellow background, green border**: Meet lower and relative significance `[L+R]`
- **Light cyan background, green border**: Meet lower and absolute significance `[L+A]`
- **Light green background, green border**: Meet all three `[L+R+A]`

### Right Side Information Panel

The right side of the chart shows brief information of top 15 significant genes:
- `[L+R+A]`: Lower+Relative+Absolute (Highest priority)
- `[R+A]`: Relative+Absolute
- `[L+R]`: Lower+Relative
- `[L+A]`: Lower+Absolute
- `[R]`: Relative significant only
- `[A]`: Absolute significant only
- `[L]`: Lower threshold only (Lowest priority)

### Gene Sorting Priority

When multiple parameters are specified, the program displays genes in the following priority order (from high to low):
1. **L+R+A**: Meet lower, relative, and absolute thresholds (3 conditions)
2. **R+A**: Meet relative and absolute thresholds (2 conditions)
3. **L+R**: Meet lower and relative thresholds (2 conditions)
4. **L+A**: Meet lower and absolute thresholds (2 conditions)
5. **R**: Meet relative threshold only (1 condition)
6. **A**: Meet absolute threshold only (1 condition)
7. **L**: Meet lower threshold only (1 condition)

Within the same priority level, genes are sorted by Z-score from high to low.

## Parameter Tuning Suggestions

### Window Size (--window-size)

- **Small window (5000-10000 bp)**: Suitable for detecting local changes but may generate noise
- **Medium window (10000-20000 bp)**: Balance local and global features (recommended)
- **Large window (20000-50000 bp)**: Detect large-scale trends but may miss local peaks

### Significance Standard (--significance-std)

- **1.5**: Loose standard, will annotate more genes
- **2.0**: Standard balance (default)
- **2.5-3.0**: Strict standard, only annotate most significant genes

### Absolute Threshold (--absolute-threshold)

Choose based on your Pi value data distribution:
- Usually take 75-90 percentile of Pi value distribution
- Suggest to check data distribution first, then set threshold

### Number of Annotated Genes (--max-genes)

- **10-20**: Suitable for generating clear charts
- Increasing number may cause label crowding, program will automatically adjust layout

## FAQ

### Q1: Why are some high Pi value genes not annotated?

A: The program annotates based on significance relative to surrounding regions, not absolute Pi values. If a gene is in a high Pi value region of the entire genome, even if the absolute value is high, it may not be significant. You can:
- Lower the `--significance-std` parameter
- Add `--absolute-threshold` parameter
- Decrease the `--window-size` parameter

### Q2: What to do if gene labels overlap?

A: The program has built-in smart label layout algorithm, but if there are still overlaps, you can:
- Reduce the `--max-genes` parameter
- Increase the `--gene-fontsize` parameter
- Adjust the `label_spacing` parameter in the code

### Q3: How to get a clearer chart?

A: Suggest adjusting the following parameters:
- Reduce number of annotated genes: `--max-genes 10`
- Adjust font size: `--gene-fontsize 10`
- Adjust image size: Modify `figsize` parameter in code

### Q4: What to do if the program runs slowly?

A: You can optimize by:
- Reduce the `--window-size` parameter
- Use smaller input data files for testing
- Ensure the latest versions of pandas and numpy are installed

## Technical Details

### Significance Calculation Method

For each gene, the program performs the following steps:

1. Calculate the average Pi value of the gene region
2. Calculate the average and standard deviation of the surrounding region (±window_size)
3. Calculate Z-score: (Gene Pi - Surrounding Average) / Surrounding Standard Deviation
4. Determine significance:
   - Relative significant: Z-score ≥ significance_std
   - Absolute significant: Gene Pi ≥ absolute_threshold

### Gene Sorting

Significant genes are sorted by the following priority:
1. **Number of conditions met**: Prioritize display of genes meeting multiple conditions
2. **Z-score**: Within same priority, sort by Z-score from high to low
3. **Pi value**: When Z-score is same, sort by Pi value

For detailed priority, please refer to "Gene Sorting Priority" section.

### Label Layout Algorithm

- Group adjacent genes (within 15000bp) together
- Assign genes in group to different height layers
- Use collision detection to avoid label overlaps
- Maximum 50 iterations to adjust label positions

## Citation

If you use this program in your research, please cite:

```
PiPlot: A Python Tool for Visualizing Pi Values in Genome Regions
https://github.com/Zmy0912/Bioinfo/tree/master/PiPlot
```

## License

This program follows the MIT License. See LICENSE file for details.

## Changelog

### v1.1.0 (2026-01-04)
- Added `--lower-threshold` parameter to support detection of lowest Pi value genes
- Optimized gene sorting logic: prioritize display of genes meeting multiple conditions
- Support multiple significance type combinations (L+R+A, R+A, L+R, L+A, etc.)
- Updated color coding system to support more significance types
- Optimized parameter processing logic: enable relative threshold mode only when `--significance-std` is explicitly specified
- Updated documentation and examples

### v1.0.0 (2026-01-03)
- Initial version release
- Support relative and absolute significance detection
- Smart label layout algorithm
- Support custom parameters
- Generate detailed analysis reports

## Contact

For questions or suggestions, please contact via:
- Submit an Issue
- Send email to myzhang0726@foxmail.com

## Acknowledgments

Thanks to the following open source projects:
- pandas
- matplotlib
- numpy
