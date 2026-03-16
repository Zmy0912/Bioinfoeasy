# Pi Sliding Window Gene Labeling Program

## Introduction

This program analyzes sliding window nucleotide diversity (Pi) data and automatically labels genes that are significantly elevated relative to their surrounding regions. The program identifies significant genes by calculating the Z-score for each gene relative to its local region and generates high-quality visualization charts.

## Key Features

- Parse Pi sliding window data files
- Parse GFF3 gene annotation files
- Identify genes significantly elevated relative to surrounding regions
- Two significance criteria:
  - **Relative Significance**: Gene Pi value is a certain multiple of standard deviations above the surrounding region average
  - **Absolute Significance**: Gene Pi value is above a set absolute threshold
- Intelligent label layout: Automatically avoids label overlap, adjacent gene labels assigned different heights
- Generate high-quality PNG charts and detailed analysis reports

## System Requirements

- Python 3.6 or higher
- Operating System: Windows, Linux, or macOS

## Installing Dependencies

### Method 1: Install using pip

```bash
pip install pandas matplotlib numpy
```

### Method 2: Install using requirements.txt

```bash
pip install -r requirements.txt
```

## Input File Formats

### Pi File Format

The Pi file should be in text format containing the following columns:

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

Standard GFF3 format gene annotation file. The program extracts rows where `feature type` is `gene`.

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
|-----------|--------|------|----------|----------|-------------|
| --pi-file | -p | File path | Yes | - | Pi value file path |
| --gff3-file | -g | File path | Yes | - | GFF3 gene annotation file path |
| --output | -o | File path | No | pi_sliding_window_plot.png | Output image file name |
| --gene-fontsize | - | Integer | No | 12 | Gene label font size |
| --max-genes | - | Integer | No | 15 | Maximum number of genes to label |
| --window-size | - | Integer | No | 10000 | Window size for calculating local average (bp) |
| --significance-std | - | Float | No | 2.0 | Relative significance criterion (standard deviation multiplier) |
| --absolute-threshold | - | Float | No | None | Absolute Pi value threshold |

### Usage Examples

#### 1. Using relative significance only (default)

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

#### 2. Using absolute threshold

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --absolute-threshold 0.05
```

#### 3. Using both relative and absolute significance

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 2.0 --absolute-threshold 0.05
```

#### 4. Custom output file and font size

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 -o my_plot.png --gene-fontsize 14
```

#### 5. Adjusting significance parameters

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 15000 --significance-std 2.5
```

#### 6. Limiting the number of labeled genes

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --max-genes 10
```

## Output Files

After running, the program generates the following files:

1. **PNG Image File**: Displays Pi sliding window curve and significant gene labels
   - File name: Default `pi_sliding_window_plot.png` or user-specified name
   - Resolution: 300 DPI
   - Format: PNG

2. **Info Text File**: Contains detailed analysis results of significant genes
   - File name: Same name as image file with `_info.txt` suffix
   - Content: Detailed statistical information for each significant gene

## Output Chart Description

### Chart Elements

- **Blue Line**: Pi value changes along genome position
- **Red Dashed Line**: Relative significance threshold (90th percentile of Pi values)
- **Purple Dashed Line**: Absolute significance threshold (if specified)
- **Colored Dots**: Positions of significant genes
  - Red: Relative significance genes
  - Purple: Absolute significance genes
  - Dark red: Both criteria met

### Gene Label Colors

- **Yellow background, red border**: Relative significance only
- **Light blue background, purple border**: Absolute significance only
- **Orange background, red border**: Both relative and absolute significance

### Right Side Information Panel

The right side of the chart displays brief information for the top 15 significant genes:
- `[R]`: Relative significance
- `[A]`: Absolute significance
- `[R+A]`: Both criteria

## Parameter Tuning Recommendations

### Window Size (--window-size)

- **Small window (5000-10000 bp)**: Suitable for detecting local changes, but may produce noise
- **Medium window (10000-20000 bp)**: Balances local and global features (recommended)
- **Large window (20000-50000 bp)**: Detects large-scale trends, may miss local peaks

### Significance Criterion (--significance-std)

- **1.5**: Lenient criterion, labels more genes
- **2.0**: Standard balance (default)
- **2.5-3.0**: Strict criterion, labels only most significant genes

### Absolute Threshold (--absolute-threshold)

Choose based on your Pi value data distribution:
- Usually take the 75-90 percentile of Pi value distribution
- Recommend checking data distribution first, then setting threshold

### Number of Labeled Genes (--max-genes)

- **10-20**: Suitable for generating clear charts
- Increasing the number may lead to label crowding, program automatically adjusts layout

## FAQ

### Q1: Why are some high Pi value genes not labeled?

A: The program labels based on significance relative to surrounding regions, not absolute Pi values. If a gene is located in a high Pi value region of the entire genome, it may not be significant even if the absolute value is high. You can:
- Lower the `--significance-std` parameter
- Add the `--absolute-threshold` parameter
- Reduce the `--window-size` parameter

### Q2: What to do if gene labels overlap?

A: The program has built-in intelligent label layout algorithm, but if overlaps still occur, you can:
- Reduce the `--max-genes` parameter
- Adjust the `--gene-fontsize` parameter
- Modify the `label_spacing` parameter in the code

### Q3: How to get a clearer chart?

A: It is recommended to adjust the following parameters:
- Reduce the number of labeled genes: `--max-genes 10`
- Adjust font size: `--gene-fontsize 10`
- Adjust image size: Modify the `figsize` parameter in the code

### Q4: What if the program runs slowly?

A: You can optimize through the following methods:
- Reduce the `--window-size` parameter
- Use smaller input data files for testing
- Ensure the latest versions of pandas and numpy are installed

## Technical Details

### Significance Calculation Method

For each gene, the program performs the following steps:

1. Calculate the average Pi value for the gene region
2. Calculate the average and standard deviation of the surrounding region (±window_size)
3. Calculate Z-score: (Gene Pi - surrounding average) / surrounding standard deviation
4. Determine significance:
   - Relative significance: Z-score ≥ significance_std
   - Absolute significance: Gene Pi ≥ absolute_threshold

### Gene Sorting

Significant genes are sorted in the following order:
1. First by Z-score in descending order
2. When Z-scores are equal, by Pi value in descending order

### Label Layout Algorithm

- Group adjacent genes (within 15000bp) into the same group
- Assign genes within a group to different height layers
- Use collision detection to avoid label overlap
- Up to 50 iterations to adjust label positions

## Citation

If you use this program in your research, please cite:

```
PiPlot Program
https://github.com/Zmy0912/Bioinfo/tree/master/PiPlot
```

## License

This program follows the MIT License. See LICENSE file for details.

## Changelog

### v1.0.0 (2026-01-03)
- Initial release
- Support for relative and absolute significance detection
- Intelligent label layout algorithm
- Support for custom parameters
- Generate detailed analysis reports

## Contact

For questions or suggestions, please contact through:
- Submit an Issue
- Send email to myzhang0726@foxmail.com

## Acknowledgments

Thanks to the following open-source projects:
- pandas
- matplotlib
- numpy

## Author

Mingyuan Zhang
Email: myzhang0726@foxmail.com
