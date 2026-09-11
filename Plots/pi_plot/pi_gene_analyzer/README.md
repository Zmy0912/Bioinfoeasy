# Gene Nucleotide Polymorphism (Pi) Analyzer

A graphical tool for analyzing genome nucleotide polymorphism, combining Pi value data from DNAsp software with GFF3 annotation files to calculate nucleotide polymorphism indices for individual genes.

## Features

- 📊 **Graphical Interface** - User-friendly GUI
- 📁 **Auto File Detection** - Automatically recognizes and preloads Pi.txt and GFF3 files from the test folder
- 🔍 **Multi-dimensional Analysis** - Calculates average Pi values and Pi at gene start, end, and midpoint positions
- 💾 **Multiple Export Formats** - Supports export to Excel (.xlsx) and CSV formats
- 📈 **Real-time Preview** - Results preview window for easy data inspection

## Input Files

### Pi.txt File
Nucleotide polymorphism data file output by DNAsp software, format example:

```
Window    Midpoint    Pi    Theta    S
   1-677       377  0.00257  0.00470  12
  278-877       577  0.00125  0.00235   6
```

### GFF3 File
Genome annotation file containing gene position information, standard GFF3 format:

```
##gff-version 3
OQ851342.1    Genbank    gene    76    1137    .    -    .    ID=gene-psbA;Name=psbA;...
OQ851342.1    Genbank    CDS     76    1137    .    -    0    ID=cds-...;...
```

## Output Results

The program generates a table with the following fields:

| Field | Description |
|-------|-------------|
| Gene Name | Gene identifier |
| Start | Start position of the gene in the genome |
| End | End position of the gene in the genome |
| Midpoint | Midpoint of the gene interval |
| Pi at Start | Pi value at the gene start position |
| Pi at End | Pi value at the gene end position |
| Pi at Midpoint | Pi value at the gene midpoint |
| Average Pi | Weighted average Pi value for the entire gene region |

### Calculation Methods

- **Average Pi**: Weighted average based on overlap length between gene region and Pi windows
- **Position Pi**: Find the Pi window containing the position and return the corresponding Pi value
- **N/A Marker**: If no Pi data is available for a position, it displays as "N/A"

## Installation

Before running the program, ensure the following Python packages are installed:

```bash
pip install pandas openpyxl
```

Or install using requirements.txt:

```bash
pip install -r requirements.txt
```

## Usage

### Method 1: Direct Execution

```bash
python pi_gene_analyzer.py
```

### Method 2: Using Python Interpreter

```bash
python3 pi_gene_analyzer.py
```

## Operation Steps

1. **Launch Program**
   - Run the script and the program will automatically preload Pi.txt and sequence.gff3 from the test folder

2. **Select Files (Optional)**
   - Click "Browse" buttons to select Pi file and GFF3 file if needed

3. **Run Analysis**
   - Click "Analyze" button to start data processing
   - The program will parse files and calculate Pi values for each gene

4. **View Results**
   - View analysis results in the results preview window
   - Scroll to see data for all genes

5. **Export Results**
   - Click "Export Results" button
   - Choose save location and file format (Excel or CSV)
   - Results will be saved to the specified file

## Output Example

Exported Excel or CSV file content example:

```
Gene Name,Start,End,Midpoint,Pi at Start,Pi at End,Pi at Midpoint,Average Pi
psbA,76,1137,606,0.00257,0.00343,0.00170,0.00235
trnK-UUU,1369,4025,2697,0.00387,0.00573,0.00439,0.00468
matK,1662,3218,2440,0.00387,0.00457,0.00370,0.00423
```

## File Structure

```
test/
├── pi_gene_analyzer.py    # Main program file
├── Pi.txt                # Pi value data from DNAsp
├── sequence.gff3         # Genome annotation file
├── README.md             # Documentation (this file)
└── requirements.txt      # Python dependencies list
```

## Notes

1. **File Encoding**: Input files should use UTF-8 encoding
2. **Data Coverage**: Ensure Pi data covers all gene positions, otherwise some genes cannot calculate Pi values
3. **Gene Names**: Genes in GFF3 file are extracted using `gene=` or `Name=` attributes
4. **Window Overlap**: Overlap region between gene and Pi windows is used to calculate weighted average Pi

## FAQ

### Q: The program won't start?
A: Ensure Python version is 3.6 or higher and all required dependencies are installed.

### Q: Some genes have no Pi data?
A: This may be because Pi windows do not cover that gene region. Please check the coverage range of Pi.txt file.

### Q: Exported Excel file shows garbled characters?
A: Try using the latest version of Excel or WPS, or choose CSV format for export.

### Q: How to handle multiple genes with the same name in GFF3?
A: The program automatically merges genes with the same name, taking the minimum and maximum positions.

## Technical Support

For questions or suggestions, please contact the developer or submit an issue.

## License

This tool is for academic research use only.

## Version History

### v1.1.0 (2024)
- Added Pi at Start, Pi at End, and Pi at Midpoint fields
- Optimized interface display and export functionality

### v1.0.0 (2024)
- Initial release
- Basic Pi value analysis and export functionality
