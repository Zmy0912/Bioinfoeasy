# Genome Annotation File Analysis Tool

A Python tool for extracting genome annotation information and generating an Excel summary table.

## Overview

This tool parses genome annotation files in **GenBank (.gb, .gbk)**, **GFF/GFF3 (.gff, .gff3)** and **FASTA (.fasta, .fa)** formats, extracts sequence and annotation statistics, and exports the results to a single Excel workbook.

### Extracted Information

- **Basic Information**: Source file, sequence ID, species name, description, genome length
- **Gene Statistics**: Total gene count, CDS count
- **RNA Statistics**: rRNA count, tRNA count, other ncRNA count
- **Other Information**: GC content, topology, molecule type, taxonomy, annotation date, etc.

## Features

- **Single-file analysis** — analyze one annotation file from the command line
- **Folder batch analysis** — automatically scan all annotation files in a folder (including subfolders) and merge the results into one summary Excel file
- **Multi-file analysis** — select several annotation files at once and merge them into one summary Excel file
- **Double-click to run** — no argument launches an input chooser that lets you pick a folder, a single file, or multiple files
- **Drag & drop support** — drag one or more files, or a folder, onto the script (or a shortcut) to process it directly
- **Error tolerance** — a single failed file does not stop the whole batch; a failure report is printed at the end
- **Auto file naming** — if the target Excel file is locked (e.g. opened in Excel), a timestamp-suffixed copy is saved automatically

### Supported File Formats

| Format | Extensions | Annotation Level |
|--------|-----------|------------------|
| GenBank | `.gb`, `.gbk` | Most complete annotation information |
| GFF/GFF3 | `.gff`, `.gff3` | Feature table (gene/CDS/rRNA/tRNA/ncRNA counts) |
| FASTA | `.fasta`, `.fa` | Sequence length + GC content only |

## Installation

Requires **Python 3.7+**.

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install biopython pandas openpyxl bcbio-gff
```

> `bcbio-gff` is only needed for full GFF/GFF3 parsing. If it is not installed, the tool automatically falls back to a simplified built-in parser.

## Usage

### 1. Double-click (no command line)

Double-click `genome_annotation_analyzer.py`. An input chooser dialog appears with three options:

| Button | Behavior |
|--------|----------|
| Select a Folder | Analyze all supported files inside that folder (including subfolders) |
| Select a Single File | Analyze just one file and export its result |
| Select Multiple Files | Merge the selected files into one summary Excel file |

After processing, the console stays open so you can review the results.

### 2. Drag & drop

Drag one or more files, or a folder, onto the script icon (or a desktop shortcut to the script):

- a dropped **folder** is analyzed as a batch (all supported files inside, subfolders included);
- dropped **files** are merged into one summary table;
- dropping a single file is treated as single-file analysis.

### 3. Command line

Analyze a single file:

```bash
python genome_annotation_analyzer.py annotation.gb
```

Analyze a single file with a custom output path:

```bash
python genome_annotation_analyzer.py annotation.gb result.xlsx
```

Analyze every annotation file in a folder:

```bash
python genome_annotation_analyzer.py /path/to/folder
```

Analyze several files at once (merged into one summary):

```bash
python genome_annotation_analyzer.py file1.gb file2.gb file3.gff
```

Batch mode with a custom output file:

```bash
python genome_annotation_analyzer.py /path/to/folder summary.xlsx
```

Debug mode (print every feature of a GenBank file):

```bash
python genome_annotation_analyzer.py annotation.gb --debug
```

## Output

The Excel output contains one row per sequence record. All rows from all files in a folder batch are merged into the same workbook.

### Excel Columns

| Field Name | Description |
|------------|-------------|
| Source File | File the record came from (folder-relative path in batch mode) |
| Sequence ID | Unique identifier for the sequence |
| Sequence Name | Sequence name |
| Species Name | Scientific name of the species |
| Description | Sequence description |
| Genome Length | Total length of the sequence (bp) |
| Topology | Linear or circular |
| Molecule Type | DNA/RNA type |
| Taxonomy | Taxonomic classification |
| GC Content (%) | GC base-pair percentage |
| Total Genes | Total number of `gene` features |
| CDS Count | Number of protein-coding sequences |
| rRNA Count | Number of ribosomal RNA features |
| tRNA Count | Number of transfer RNA features |
| Other ncRNA Count | Number of other non-coding RNA features |
| Annotation Date | Date of annotation |
| Data Version | Database / data-file division |
| Annotation Source | Data source annotation |

> Additional gene-type columns (e.g. protein-coding gene counts) are appended automatically when present in a GenBank file.

### Default Output Filename

| Mode | Default Output |
|------|----------------|
| Single file | `<first sequence ID>.xlsx` in the current working directory |
| Folder batch | `<folder name>_annotation_summary.xlsx` inside the analyzed folder |
| Multiple files (merged) | `annotation_summary.xlsx` next to the first selected file |

If the default filename is already open in Excel, a timestamp-suffixed file is saved instead (e.g. `<folder name>_annotation_summary_20260907_153000.xlsx`).

## Notes & FAQ

### Q1: Why do some fields show "N/A"?

Possible reasons:
- The format does not contain that information (e.g. FASTA files have no gene annotations).
- The file lacks the corresponding annotation.
- The qualifier/field naming is non-standard.

### Q2: Does the tool support batch processing?

Yes. Since v3.0 the tool supports batch processing: double-click the script and choose a folder or multiple files, or drag files/folders onto the script. Folders are scanned recursively for annotation files.

### Q3: How are duplicate sequence IDs handled in batch mode?

Each record key is prefixed internally with its file index, and the `Source File` column records where each row came from, so identical sequence IDs from different files will not overwrite each other.

### Q4: What happens if one file fails to parse in a folder batch?

Only the failed file is skipped and recorded in the failure report at the end; the remaining files are still processed and exported.

### Q5: Can I select only specific files instead of a whole folder?

Yes. In the double-click chooser, use "Select Multiple Files" (or "Select a Single File"), or pass several file paths on the command line. The selected files are merged into one summary Excel file.

## Example Output

### Single-file mode

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
  Source File: chloroplast.gb
  CDS Count: 78
  GC Content (%): 38.50
  ...

Results successfully exported to: NC_001665.xlsx

Number of records extracted: 1
Number of fields: 18
```

### Folder batch mode (final summary)

```
Found 20 annotation file(s). Starting batch analysis...

[1/20] Analyzing: sp1/NC_001665.gb
...
[20/20] Analyzing: sp20/NC_020000.gb

Results successfully exported to: /path/to/folder/chloroplasts_annotation_summary.xlsx

Number of records extracted: 20
Number of fields: 18

============================================================
Batch Analysis Summary:
============================================================
  Files processed   : 20
  Files succeeded   : 20
  Files failed      : 0
  Records extracted : 20
```

## Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| biopython | >=1.81 | Genome/annotation file parsing |
| pandas | >=1.5.0 | Data processing and table construction |
| openpyxl | >=3.0.0 | Excel file writing |
| bcbio-gff | >=0.6.0 | GFF/GFF3 parsing (optional) |

## Version History

### v3.0
- Added folder batch analysis (recursive scan of `.gb/.gbk/.gff/.gff3/.fasta/.fa`)
- Added double-click launch with an input chooser (folder / single file / multiple files)
- Added multi-file analysis that merges selected files into one summary Excel file
- Added drag & drop support for folders and files
- Added the `Source File` column to track which file each record came from
- All user-facing messages converted to English

### v2.0
- Removed chloroplast region length statistics
- Retained core genome annotation information extraction
- Improved error handling

### v1.0
- Initial version
- Basic format support
- Gene statistics functionality

## License

This tool is provided for learning and research purposes only.

## Contact

For questions or suggestions, feel free to provide feedback.
Email: myzhang0726@foxmail.com
