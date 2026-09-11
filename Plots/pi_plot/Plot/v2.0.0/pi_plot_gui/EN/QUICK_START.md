# Quick Start Guide

## I. Installation Steps

### 1. Check Python Environment

Open Command Prompt (CMD) or PowerShell and type:

```bash
python --version
```

Ensure Python version is 3.6 or higher.

### 2. Install Dependencies

```bash
pip install pandas matplotlib numpy
```

Or use requirements.txt:

```bash
pip install -r requirements.txt
```

### 3. Start Program

**Method 1: Double-click to start**
- Double-click `start_program.bat` file

**Method 2: Command line start**
```bash
python pi_plot_gui.py
```

## II. Quick Usage Examples

### Example 1: First Time Use (Recommended)

1. Start the program
2. Select Pi file
3. Select GFF3 file
4. Select output directory
5. Keep default parameters
6. Select "PNG Format" (quick preview)
7. Click "Start Analysis"

### Example 2: Academic Publishing Purpose

1. Start the program
2. Select files and output directory
3. Select "All Formats" (PNG+SVG+PDF)
4. Set font size to 10
5. Set max genes to 12
6. Click "Start Analysis"
7. Use SVG or PDF format from output directory for publication

### Example 3: Detect Lowest Pi Value Genes

1. Start the program
2. Select files and output directory
3. Set lower threshold to 0.01 (adjust according to actual situation)
4. Select SVG format
5. Click "Start Analysis"

## III. Output Files Explanation

After running, the following files will be generated in the output directory:

```
output_directory/
├── pi_sliding_window_plot.png    # PNG format chart (if PNG selected)
├── pi_sliding_window_plot.svg    # SVG format chart (if SVG selected)
├── pi_sliding_window_plot.pdf    # PDF format chart (if PDF selected)
└── pi_sliding_window_plot_info.txt  # Detailed analysis results
```

### Chart File Explanation

| Format | Extension | Purpose | Editing Software |
|--------|-----------|---------|------------------|
| PNG | .png | Screen viewing | Any image viewer |
| SVG | .svg | Academic publishing, printing | Inkscape, Adobe Illustrator |
| PDF | .pdf | Document embedding | Adobe Acrobat, Inkscape |

### Info File Explanation

The `_info.txt` file contains:
- Analysis parameter settings
- Total number of significant genes found
- Detailed information for each significant gene:
  - Gene name
  - Gene position
  - Pi value
  - Z-score
  - Significance type

## IV. Common Parameter Settings

### Quick Preview
```
Output Format: PNG
Font Size: 12
Max Genes: 15
```

### Academic Publishing
```
Output Format: SVG or PDF
Font Size: 10
Max Genes: 10-12
```

### Large Genomes
```
Output Format: SVG (vector, smaller file)
Font Size: 8-10
Max Genes: 20-30
Window Size: 15000-20000
```

### Local Region Analysis
```
Window Size: 5000-8000
Relative Significance: 1.5-2.0
Max Genes: 10
```

## V. Troubleshooting

### Problem 1: Program cannot start

**Solution**:
1. Confirm Python is installed: `python --version`
2. Confirm dependencies are installed: `pip list | findstr pandas matplotlib numpy`
3. Check if tkinter is available: `python -c "import tkinter"`

### Problem 2: Text displays as boxes

**Solution**:
- Windows: Built-in support, should display normally
- Linux: Install Chinese font `sudo apt-get install fonts-wqy-zenhei`
- macOS: Built-in support, should display normally

### Problem 3: Output files are very large

**Solution**:
- Reduce number of annotated genes (set max_genes to 10-15)
- Reduce font size (set gene_fontsize to 8-10)
- Larger SVG/PDF files are normal

### Problem 4: Analysis is slow

**Solution**:
- Increase window size (set window_size to 15000-20000)
- Reduce number of annotated genes
- Close other CPU-intensive programs

## VI. Getting Help

The program includes detailed usage instructions:
1. Start the program
2. Click "Help" tab
3. Read the complete help documentation

For more help:
- View `README.md` file
- Send email to myzhang0726@foxmail.com

## VII. File Format Reference

### Pi File Example

```
Window    Midpoint    Pi    Theta
1-10000   5000        0.0234 0.0256
10001-20000  15000    0.0287 0.0301
20001-30000  25000    0.0198 0.0212
```

### GFF3 File Example

```
##gff-version 3
chr1    Ensembl    gene    1000    5000    .    +    .    ID=gene-001;Name=GENE1
chr1    Ensembl    gene    8000    12000    .    -    .    ID=gene-002;Name=GENE2
chr1    Ensembl    gene    15000    18000    .    +    .    ID=gene-003;Name=GENE3
```

## VIII. Quick Operation Tips

1. **Batch Processing**: If you need to process multiple files, you can use the command-line version (v1.1.0) to write batch scripts

2. **Save Parameters**: Record successful parameter configurations for reproducing results

3. **Preview and Adjust**: Use PNG format for quick preview first, then generate SVG/PDF after confirming results

4. **Multi-format Output**: Select "All Formats" to generate all three formats at once, meeting various needs

5. **Check Results**: Check the `_info.txt` file first to understand analysis results, then view the charts

---

Happy using!
