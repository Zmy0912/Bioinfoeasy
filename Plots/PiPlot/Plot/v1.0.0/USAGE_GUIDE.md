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

Pi file should contain nucleotide diversity data for sliding windows:

```
Window	Midpoint	Pi	Theta
1-10000	5000	0.0234	0.0256
10001-20000	15000	0.0287	0.0301
20001-30000	25000	0.0198	0.0212
```

- **Window**: Window range, format as "start-end"
- **Midpoint**: Window midpoint position
- **Pi**: Nucleotide diversity value for that window
- **Theta**: Theta value (optional, not used by program)

#### GFF3 File Structure

GFF3 file is the standard genome annotation format:

```
##gff-version 3
chr1	Ensembl	gene	1000	5000	.	+	.	ID=gene-001;Name=GENE1
chr1	Ensembl	gene	8000	12000	.	-	.	ID=gene-002;Name=GENE2
chr1	Ensembl	gene	15000	18000	.	+	.	ID=gene-003;Name=GENE3
```

The program only processes rows where `feature type` is `gene`.

### Parameter Details

#### Required Parameters

##### `-p, --pi-file`
- **Purpose**: Specify Pi data file
- **Example**: `-p Pi.txt`
- **Note**: File must exist and be in correct format

##### `-g, --gff3-file`
- **Purpose**: Specify GFF3 gene annotation file
- **Example**: `-g sequence.gff3`
- **Note**: File must contain gene features (feature type = "gene")

#### Optional Parameters

##### `-o, --output`
- **Purpose**: Specify output image file name
- **Default**: `pi_sliding_window_plot.png`
- **Example**: `-o my_result.png`
- **Description**: Program will automatically generate a `_info.txt` file with the same name

##### `--gene-fontsize`
- **Purpose**: Control font size of gene labels
- **Default**: `12`
- **Range**: Recommend 8-16
- **Example**: `--gene-fontsize 14`
- **Tip**: Reduce font size when there are many labels to avoid overlap

##### `--max-genes`
- **Purpose**: Limit number of genes to label
- **Default**: `15`
- **Range**: Recommend 5-30
- **Example**: `--max-genes 10`
- **Description**: Program will select the most significant genes to label

##### `--window-size`
- **Purpose**: Window size for calculating local average
- **Default**: `10000` (bp)
- **Range**: Recommend 5000-50000
- **Example**: `--window-size 15000`
- **Impact**:
  - Small window: More sensitive, may detect local peaks
  - Large window: More stable, may miss local changes

##### `--significance-std`
- **Purpose**: Relative significance criterion (standard deviation multiplier)
- **Default**: `2.0`
- **Range**: Recommend 1.5-3.0
- **Example**: `--significance-std 2.5`
- **Description**:
  - Smaller value: More genes labeled
  - Larger value: Only most significant genes labeled

##### `--absolute-threshold`
- **Purpose**: Absolute Pi value threshold
- **Default**: None
- **Type**: Float
- **Example**: `--absolute-threshold 0.05`
- **Description**:
  - Genes with Pi value above this threshold are considered significant
  - Can be used together with relative significance criterion

### Usage Scenarios

#### Scenario 1: Quick Data Overview

Use default parameters to quickly understand data:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

#### Scenario 2: Detect Local High Variation Regions

Use smaller window and lower significance criterion:

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

#### Scenario 4: Focus on Specific Pi Value Range

Use absolute threshold:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --absolute-threshold 0.08
```

#### Scenario 5: Generate Clear Charts

Limit labeled gene count and adjust font:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --max-genes 10 \
  --gene-fontsize 10
```

#### Scenario 6: Use Both Significance Criteria

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3 \
  --significance-std 2.0 \
  --absolute-threshold 0.06
```

### Interpreting Output Files

#### PNG Image File

The image contains the following elements:

1. **Main Curve**: Blue line showing Pi value changes along position
2. **Threshold Lines**:
   - Red dashed line: Relative significance threshold (90th percentile of Pi)
   - Purple dashed line: Absolute significance threshold (if specified)
3. **Gene Markers**: Colored dots marking significant gene positions
4. **Gene Labels**: Gene name labels with short pointers pointing to gene positions
5. **Right Side Info**: Brief list of top 15 significant genes

#### Info Text File (_info.txt)

Contains detailed information for each labeled gene:

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

### Parameter Tuning Strategy

#### Step 1: Initial Exploration

Run with default parameters to observe results:

```bash
python pi_plot.py -p Pi.txt -g sequence.gff3
```

Check output image and info file to understand:
- How many significant genes in total
- Their distribution across genome
- Approximate range of Pi values

#### Step 2: Adjust Window Size

Choose window size based on data characteristics:

```bash
# Detect local peaks
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 5000

# Detect large-scale trends
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 20000
```

#### Step 3: Adjust Significance Criterion

Adjust sensitivity as needed:

```bash
# More lenient, label more genes
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 1.5

# More strict, label only most significant genes
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 2.5
```

#### Step 4: Set Absolute Threshold (Optional)

If screening based on absolute Pi values is needed:

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
- Use absolute or relative paths
- Ensure file extension is correct

#### Error 2: Pi File Format Error

Program errors while parsing Pi file.

**Solution**:
- Ensure file contains a header row
- Check if column separator is tab or space
- Verify enough data columns (at least 4 columns)

#### Error 3: GFF3 File Format Error

```
Error: GFF3 file does not exist: sequence.gff3
```

**Solution**:
- Ensure GFF3 file exists
- Check if file contains gene features
- Verify file is in standard GFF3 format

#### Error 4: Excessive Label Overlap

Gene labels in chart overlap severely.

**Solution**:
```bash
# Reduce number of labeled genes
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
# Lower significance criterion
python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 1.5

# Add absolute threshold
python pi_plot.py -p Pi.txt -g sequence.gff3 --absolute-threshold 0.02

# Reduce window size
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 5000
```

### Performance Optimization

#### Processing Large Datasets

If working with large genomes:

```bash
# Reduce window size to speed up calculation
python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 5000

# Limit number of labeled genes
python pi_plot.py -p Pi.txt -g sequence.gff3 --max-genes 10
```

#### Batch Processing

Processing multiple samples:

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

#### Customizing Colors

Modify color settings in code:

```python
# Find color definitions in pi_plot.py
marker_color = 'red' if sig_type == 'relative' else 'purple'
label_color = 'yellow'
edge_color = 'red'
```

#### Adjusting Chart Size

Modify figsize parameter in code:

```python
# In plot_pi_with_genes function
fig, ax = plt.subplots(figsize=(20, 12))
# Change to
fig, ax = plt.subplots(figsize=(24, 14))
```

#### Modifying Label Layout Parameters

Adjust spacing parameters in code:

```python
# In pi_plot.py
label_spacing = max_pi * 0.15  # Vertical spacing
min_horizontal_distance = max_pi * 0.05  # Horizontal spacing
base_height = max_pi * 0.03  # Base height
```

### Data Preparation Recommendations

#### Pi File Generation

Generate Pi file using sliding window tools, for example:

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
# Convert format using gffread
gffread input.gff3 -T -o genes.gff3

# Filter to keep only genes
awk '$3=="gene"' input.gff3 > genes.gff3
```

## Best Practices

1. **Start with small-scale testing**: Test with small files first to verify parameter settings
2. **Save intermediate results**: Keep output files from different parameters for comparison
3. **Check data quality**: Verify integrity of input files before running
4. **Adjust parameters gradually**: Change one parameter at a time and observe effects
5. **Record parameter configurations**: Save successful parameter combinations for reproducibility

## Getting Help

View complete parameter description:

```bash
python pi_plot.py --help
```

## Updates and Support

- Check latest documentation: README.md
- Submit issues: GitHub Issues
- Contact author: myzhang0726@foxmail.com
