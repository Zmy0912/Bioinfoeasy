# Genome Annotation File Analysis Tool

## Overview

This tool can extract the following information from genome annotation files in GB (GenBank), GFF, and other formats:

### Extracted Information
- **Basic Information**: Sequence ID, species name, description, genome sequence length
- **Gene Statistics**: Gene types and counts, total gene count, CDS count
- **RNA Statistics**: rRNA count, tRNA count, other ncRNA count
- **Other Information**: GC content, topology, taxonomy, annotation date, etc.

### Supported File Formats
- GenBank (.gb, .gbk)
- GFF/GFF3 (.gff, .gff3)
- FASTA (.fasta, .fa)

## Installation

```bash
pip install -r requirements.txt
```

Or install individually:

```bash
pip install biopython pandas openpyxl BCBio-GFF
```

## Usage

### Basic Usage

```bash
python genome_annotation_analyzer.py <input_file>
```

### Examples

1. **Analyze GenBank file**
   ```bash
   python genome_annotation_analyzer.py chloroplast.gb
   ```

2. **Analyze GFF file**
   ```bash
   python genome_annotation_analyzer.py annotation.gff
   ```

3. **Analyze FASTA file**
   ```bash
   python genome_annotation_analyzer.py genome.fasta
   ```

4. **Specify output filename**
   ```bash
   python genome_annotation_analyzer.py annotation.gb result.xlsx
   ```

5. **Debug mode (view all features)**
   ```bash
   python genome_annotation_analyzer.py annotation.gb --debug
   ```

## Output

### Excel Table Fields

| Field Name | Description |
|------------|-------------|
| Sequence ID | Unique identifier for the genome |
| Sequence Name | Genome name |
| Species Name | Scientific name of the species |
| Description | Sequence description |
| Genome Length | Total length of genome (bp) |
| Topology | Linear or circular |
| Molecule Type | DNA/RNA type |
| Taxonomy | Taxonomic classification |
| GC Content (%) | GC base pair percentage |
| Total Genes | Total number of genes |
| CDS Count | Number of protein coding sequences |
| rRNA Count | Number of ribosomal RNA |
| tRNA Count | Number of transfer RNA |
| Other ncRNA Count | Number of other non-coding RNA |
| Annotation Date | Data annotation date |
| Data Version | Database version |
| Annotation Source | Data source |

### Output Filename

If no output filename is specified, the program will use the sequence ID as the filename, for example:
- Input: `NC_001665.gb`
- Output: `NC_001665.xlsx`

## FAQ

### Q1: Why do some fields show "N/A"?

Possible reasons:
- The file format does not contain that information (e.g., FASTA files)
- The file lacks corresponding annotations
- Field naming is non-standard

### Q2: Does it support batch processing?

The current version does not support batch processing, but you can use shell scripts:
```bash
for file in *.gb; do
    python genome_annotation_analyzer.py "$file"
done
```

## Features

### 1. Multi-format Support
- GenBank: Most complete annotation format
- GFF: General-purpose genome feature format
- FASTA: Basic sequence format

### 2. Smart Recognition
- Automatic file format detection
- Multi-dimensional information extraction

### 3. Error Tolerance
- Missing information marked as N/A
- Automatic exception handling and prompts
- Debug mode for troubleshooting

### 4. Efficient Processing
- Efficient parsing with BioPython
- Support for large files
- Fast Excel export

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| biopython | >=1.81 | Genome file parsing |
| pandas | >=1.5.0 | Data processing |
| openpyxl | >=3.0.0 | Excel file writing |
| BCBio-GFF | >=0.6.0 | GFF file parsing (optional) |

## Example Output

```
============================================================
Genome Annotation File Analysis Tool
============================================================

Analyzing GenBank file: chloroplast.gb

============================================================
Analysis Results Summary:
============================================================

Sequence: NC_001665
------------------------------------------------------------
  CDS Count: 78
  GC Content (%): 38.50
  Other ncRNA Count: 4
  rRNA Count: 4
  tRNA Count: 30
  Topology: circular
  Data Version: PLN
  Sequence ID: NC_001665
  Sequence Name: NC_001665
  Description: Nicotiana tabacum chloroplast, complete genome
  Species Name: Nicotiana tabacum
  Annotation Source: chloroplast
  Molecule Type: DNA
  Total Genes: 112
  Genome Length: 155500
  Annotation Date: 15-MAR-2024
  Taxonomy: Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Core eudicots; Gunneridae; Pentapetalae; asterids; lamiids; Solanales; Solanaceae; Nicotiana

Results successfully exported to: NC_001665.xlsx

Number of records extracted: 1
Number of fields: 18
```

## Version History

### v2.0
- Removed chloroplast region length statistics
- Retained core genome annotation information extraction
- Improved error handling

### v1.0
- Initial version
- Basic format support
- Gene statistics functionality

## License

This tool is for learning and research purposes only.

## Contact

For questions or suggestions, please feel free to provide feedback.
email: <EMAIL>myzhang0726@foxmail.com
