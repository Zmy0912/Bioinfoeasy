# Usage Guide

## Quick Start

### 1. Environment Setup

Ensure Python 3.6 or higher is installed:

```bash
python --version
```

### 2. Install Dependencies

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install pandas matplotlib numpy
```

### 3. Prepare Input Files

Ensure you have the following two files:
- `Pi.txt`: Pi sliding window data file
- `sequence.gff3`: Gene annotation file (GFF3 format)

### 4. Run the Program

The simplest way to run:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

## Detailed Usage Instructions

### Understanding Input Files

#### Pi File Structure

The Pi file should contain nucleotide diversity data for sliding windows:

```
Window	Midpoint	Pi	Theta
1-10000	5000	0.0234	0.0256
10001-20000	15000	0.0287	0.0301
20001-30000	25000	0.0198	0.0212
```

- **Window**: Window range, format "start_position-end_position"
- **Midpoint**: Window midpoint position
- **Pi**: Nucleotide diversity value of that window
- **Theta**: Theta value (optional, not used by program)

#### GFF3 File Structure

GFF3 file is the standard genomic annotation format:

```
##gff-version 3
chr1	Ensembl	gene	1000	5000	.	+	.	ID=gene-001;Name=GENE1
chr1	Ensembl	gene	8000	12000	.	-	.	ID=gene-002;Name=GENE2
chr1	Ensembl	gene	15000	18000	.	+	.	ID=gene-003;Name=GENE3
```

The program only processes rows with `feature type` as `gene`.

### Detailed Parameter Description

#### Required Parameters

##### `-p, --pi-file`
- **Purpose**: Specify Pi data file
- **Example**: `-p Pi.txt`
- **Note**: File must exist and have correct format

##### `-g, --gff3-file`
- **Purpose**: Specify GFF3 gene annotation file
- **Example**: `-g sequence.gff3`
- **Note**: File must contain gene features (feature type = "gene")

#### Optional Parameters

##### `-o, --output`
- **Purpose**: Specify output image file name
- **Default**: `pi_sliding_window_plot.png`
- **Example**: `-o my_result.png`
- **Note**: Program will automatically generate `_info.txt` file with same name

##### `--gene-fontsize`
- **Purpose**: Control font size of gene labels
- **Default**: `12`
- **Range**: Suggest 8-16
- **Example**: `--gene-fontsize 14`
- **Tip**: Reduce font size when there are too many labels to avoid overlap

##### `--max-genes`
- **Purpose**: Limit number of genes to annotate
- **Default**: `15`
- **Range**: Suggest 5-30
- **Example**: `--max-genes 10`
- **Note**: Program will select most significant genes to annotate

##### `--window-size`
- **Purpose**: Window size for calculating local average
- **Default**: `10000` (bp)
- **Range**: Suggest 5000-50000
- **Example**: `--window-size 15000`
- **Impact**:
  - Small window: More sensitive, may detect local peaks
  - Large window: More stable, may miss local changes

##### `--significance-std`
- **Purpose**: Relative significance standard (standard deviation multiplier)
- **Default**: `2.0`
- **Range**: Suggest 1.5-3.0
- **Example**: `--significance-std 2.5`
- **Note**:
  - Smaller value: More genes will be annotated
  - Larger value: Only most significant genes will be annotated

##### `--absolute-threshold`
- **Purpose**: Absolute Pi value threshold
- **Default**: None
- **Type**: Float
- **Example**: `--absolute-threshold 0.05`
- **Note**:
  - Genes with Pi value above this threshold are considered significant
  - Can be used together with relative significance standard
  - Can be used together with lower threshold

##### `--lower-threshold`
- **Purpose**: Lower Pi value threshold (lowest Pi value genes)
- **Default**: None
- **Type**: Float
- **Example**: `--lower-threshold 0.01`
- **Note**:
  - Genes with Pi value below this threshold are considered significant
  - Used to detect genes with lowest Pi values in genome
  - Can be used together with relative significance standard and absolute threshold

### Usage Scenarios

#### Scenario 1: Quick Data Overview

Use default parameters to quickly understand data:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

#### Scenario 2: Detect Local High Variation Regions

Use smaller window and lower significance standard:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --window-size 5000 \
  --significance-std 1.5
```

#### Scenario 3: Identify Global High Variation Regions

Use larger window:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --window-size 20000 \
  --significance-std 2.5
```

#### Scenario 4: Focus on Specific Pi Value Range Genes

Use absolute threshold:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --absolute-threshold 0.08
```

#### Scenario 5: Generate Clear Charts

Limit number of annotated genes and adjust font:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --max-genes 10 \
  --gene-fontsize 10
```

#### Scenario 6: Use Two Significance Standards Together

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --significance-std 2.0 \
  --absolute-threshold 0.06
```

#### Scenario 7: Use Lower Threshold to Detect Lowest Pi Value Genes

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --lower-threshold 0.01
```

#### Scenario 8: Use All Three Thresholds Together (Prioritize Genes Meeting Multiple Conditions)

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --significance-std 2.0 \
  --absolute-threshold 0.05 \
  --lower-threshold 0.01
```

### Output File Interpretation

#### PNG Image File

The image contains the following elements:

1. **Main curve**: Blue line shows Pi value changes along position
2. **Threshold lines**:
   - Red dashed line: Relative significance threshold (Pi 90th percentile)
   - Purple dashed line: Absolute significance threshold (if specified)
   - Green dashed line: Lower threshold (if specified)
3. **Gene markers**: Colored dots mark significant gene positions
4. **Gene labels**: Gene name labels with short pointers pointing to gene positions
5. **Right side information**: Brief list of top 15 significant genes, including significance type markers

#### Information Text File (_info.txt)

Contains detailed information for each annotated gene:

```
================================================================================
Significantly Elevated Genes Analysis
================================================================================

Analysis Parameters:
  Window size for local comparison: 10000 bp
  Relative significance threshold: 2.0 standard deviations
  Number of genes analyzed: 2500
  Number of significantly elevated genes found: 45

Pi Threshold (90th percentile): 0.0450

================================================================================
Top Significantly Elevated Genes:
================================================================================

1. Gene: GENE1
   Position: 125000-128500 bp
   Gene midpoint: 126750 bp
   Gene Pi value: 0.0850
   Region mean Pi: 0.0300
   Region std Pi: 0.0150
   Z-score: 3.6667
   Significance type: relative
   Significance: 3.67x above region average
```

### Gene Sorting Priority

When multiple threshold parameters are specified, the program will prioritize display of genes meeting multiple conditions. Priority from high to low is:

| Priority | Type Marker | Description | Number of Conditions Met |
|----------|-------------|-------------|--------------------------|
| 1 | `[L+R+A]` | Lower+Relative+Absolute | 3 |
| 2 | `[R+A]` | Relative+Absolute | 2 |
| 3 | `[L+R]` | Lower+Relative | 2 |
| 4 | `[L+A]` | Lower+Absolute | 2 |
| 5 | `[R]` | Relative only | 1 |
| 6 | `[A]` | Absolute only | 1 |
| 7 | `[L]` | Lower only | 1 |

Within same priority level, genes are sorted by Z-score from high to low.

### Parameter Tuning Strategy

#### Understanding Parameter Combination Effects

Different parameter combinations will produce different output modes:

- **Relative only**: `--significance-std` → Sort by relative significance
- **Absolute only**: `--absolute-threshold` → Sort by absolute threshold
- **Lower only**: `--lower-threshold` → Sort by lower threshold
- **Relative+Absolute**: `--significance-std` + `--absolute-threshold` → Prioritize `[R+A]`, then `[R]` or `[A]`
- **Relative+Lower**: `--significance-std` + `--lower-threshold` → Prioritize `[L+R]`, then `[R]` or `[L]`
- **Absolute+Lower** (no relative): `--absolute-threshold` + `--lower-threshold` → Evenly distribute `[A]` and `[L]`
- **All three**: All parameters → Prioritize `[L+R+A]`, then `[R+A]`, `[L+R]`, `[L+A]`, finally single conditions

#### Step 1: Initial Exploration

Run with default parameters and observe results:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

Check output image and information file to understand:
- Total number of significant genes
- Their distribution in the genome
- Approximate range of Pi values

#### Step 2: Adjust Window Size

Choose window size based on data characteristics:

```bash
# Detect local peaks
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 5000

# Detect large-scale trends
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 20000
```

#### Step 3: Adjust Significance Standard

Adjust sensitivity based on needs:

```bash
# More relaxed, annotate more genes
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 1.5

# More strict, only annotate most significant genes
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 2.5
```

#### Step 4: Set Absolute Threshold (Optional)

If filtering based on absolute Pi value is needed:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 --absolute-threshold 0.05
```

#### Step 5: Optimize Chart Display

Adjust display parameters for best visual effect:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --max-genes 12 \
  --gene-fontsize 11 \
  -o final_plot.png
```

### Common Errors and Solutions

#### Error 1: File Not Found

```
Error: Pi file does not exist: Pi.txt
```

**Solution**:
- Check if file path is correct
- Use absolute path or relative path
- Ensure file extension is correct

#### Error 2: Pi File Format Error

Error occurred while parsing Pi file.

**Solution**:
- Ensure file contains header line
- Check if column separator is tab or space
- Verify data columns are sufficient (at least 4 columns)

#### Error 3: GFF3 File Format Error

```
Error: GFF3 file does not exist: sequence.gff3
```

**Solution**:
- Ensure GFF3 file exists
- Check if file contains gene features
- Verify file format is standard GFF3

#### Error 4: Too Many Label Overlaps

Gene labels overlap severely in the chart.

**Solution**:
```bash
# Reduce number of annotated genes
python pi_plot.py -p Pi.txt -g sequence.gff3 --max-genes 8

# Reduce font size
python pi_plot.py -p Pi.txt -g sequence.gff3 --gene-fontsize 9
```

#### Error 5: No Significant Genes Found

```
Found 0 significant genes
```

**Solution**:
```bash
# Lower significance standard
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 1.5

# Add absolute threshold
python pi_plot.py -p Pi.txt -g sequence.gff3 --absolute-threshold 0.02

# Add lower threshold
python pi_plot.py -p Pi.txt -g sequence.gff3 --lower-threshold 0.01

# Reduce window size
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 5000
```

### Performance Optimization

#### Large Dataset Processing

When processing large genomes:

```bash
# Reduce window size to speed up calculation
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 5000

# Limit number of annotated genes
python pi_plot.py -p Pi.txt -g sequence.gff3 --max-genes 10
```

#### Batch Processing

Process multiple samples:

```bash
# Windows PowerShell
for ($i=1; $i -le 5; $i++) {
  python pi_plot.py -p sample${i}.txt -g genes.gff3 -o sample${i}_plot.png
}

# Linux/Mac Bash
for i in {1..5}; do
  python pi_plot.py -p sample${i}.txt -g genes.gff3 -o sample${i}_plot.png
done
```

### Advanced Tips

#### Custom Colors

Modify color settings in code:

```python
# In pi_plot.py, find color definitions
marker_color = 'red' if sig_type == 'relative' else 'purple'
label_color = 'yellow'
edge_color = 'red'
```

#### Adjust Chart Size

Modify figsize parameter in code:

```python
# In plot_pi_with_genes function
fig, ax = plt.subplots(figsize=(20, 12))
# Change to
fig, ax = plt.subplots(figsize=(24, 14))
```

#### Modify Label Layout Parameters

Adjust spacing parameters in code:

```python
# In pi_plot.py
label_spacing = max_pi * 0.15  # Vertical spacing
min_horizontal_distance = max_pi * 0.05  # Horizontal spacing
base_height = max_pi * 0.03  # Base height
```

### Data Preparation Suggestions

#### Pi File Generation

Use sliding window tools to generate Pi files, for example:

```bash
# Using vcftools
vcftools --vcf input.vcf --window-pi 10000 --out Pi

# Using popgenome (R)
library(PopGenome)
genome <- readVCF("input.vcf")
genome <- sliding.window(genome, win.size = 10000, step.size = 5000)
write.pPi(genome, "Pi.txt")
```

#### GFF3 File Preparation

Ensure GFF3 file contains necessary gene information:

```bash
# Use gffread to convert format
gffread input.gff3 -T -o genes.gff3

# Filter to keep only genes
awk '$3=="gene"' input.gff3 > genes.gff3
```

## Best Practices

1. **Start with small-scale testing**: Test with small files first to verify parameter settings
2. **Save intermediate results**: Keep output files with different parameters for comparison
3. **Check data quality**: Verify integrity of input files before running
4. **Adjust parameters step by step**: Only adjust one parameter at a time and observe effects
5. **Record parameter configurations**: Save successful parameter combinations for reproducibility

## Getting Help

View complete parameter descriptions:

```bash
python pi_plot.py --help
```

## Updates and Support

- View latest documentation: README.md
- Submit issues: GitHub Issues
- Contact author: myzhang0726@foxmail.com
