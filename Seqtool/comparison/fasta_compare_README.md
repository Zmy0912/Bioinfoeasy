# FASTA Sequence Comparison Tool

A desktop application for comparing sequences across multiple FASTA files, capable of identifying common and unique sequences and providing detailed visual statistical information.

## Features

- **Multi-file Comparison**: Compare sequences from multiple FASTA files simultaneously
- **Common Sequence Identification**: Automatically identify sequences present in all files
- **Unique Sequence Identification**: Identify sequences unique to each file
- **Sequence Viewing**: View detailed information and complete content of any sequence
- **Statistics**: Provide detailed statistics and frequency distribution
- **Export**: Support exporting common or unique sequences to FASTA files
- **Graphical Interface**: Friendly GUI based on Tkinter

## System Requirements

- Python 3.6 or higher
- Windows / macOS / Linux
- No additional dependencies (uses Python standard library)

## Installation

1. Ensure Python 3.x is installed
2. Place `fasta_compare.py` and `fasta_compare.bat` in the same folder
3. Double-click `fasta_compare.bat` to launch the program

## Quick Start

### Basic Usage Flow

1. **Launch Program**: Double-click `fasta_compare.bat`
2. **Select Folder**: Click "Select Folder" button, select folder containing FASTA files
3. **Start Comparison**: Click "Start Comparison" button to compare sequences
4. **View Results**: View common sequences, unique sequences, and statistics in different tabs
5. **Export Results**: Click export button to save results as FASTA files

### Tab Description

#### 1. Common Sequences Tab
Display sequences present in all loaded files.
- **Sequence Name**: Gene/sequence name extracted from FASTA header
- **Count**: Number of times this sequence appears in files
- **Files**: List of files containing this sequence
- **Double-click Action**: Double-click any sequence to view details in "Sequence Viewer" tab

#### 2. Unique Sequences Tab
Display sequences unique to each file (sequences present in only one file).
- **Sequence Name**: Name of the sequence
- **Source File**: File containing this sequence
- **Double-click Action**: Double-click any sequence to view details

#### 3. Statistics Tab
Provide detailed statistical information:
- Number of files compared
- Number of sequences in each file
- Total number of distinct sequences
- Number of common sequences
- Number of unique sequences for each file
- Sequence frequency distribution (number of sequences appearing in N files)

#### 4. Sequence Viewer Tab
View detailed information of specific sequence:
- Enter sequence name and source file to query
- Display complete header and metadata of sequence
- Display sequence length and complete sequence content
- Support copying sequence to clipboard
- Support exporting current sequence as FASTA file

## Sequence Recognition Rules

Program uses following rules to extract sequence names from FASTA headers:

1. **Contains `_join{`**: Extract `rps12_join{position}` → `rps12`
2. **Contains `[_`**: Extract `psbA_[start:end](strand)` → `psbA`
3. **Contains `[`**: Extract `gene_name[info]` → `gene_name`
4. **Other cases**: Use complete header as identifier

Examples:
```
>rps12_join{...}|A10.fasta     →  rps12
>psbA_[76:1137](-)             →  psbA
>cemA_[...] |A11.fasta         →  cemA
```

## Export Features

### Export Common Sequences
Export sequences common to all files as a FASTA file:
1. Click "Export Common Sequences" button
2. Select save location and file name
3. Program adds source file marker for each sequence

### Export Unique Sequences
Create separate FASTA file for each file with unique sequences:
1. Click "Export Unique Sequences" button
2. Select save folder
3. Program creates `filename_unique.fasta` file for each file

### Export Single Sequence
In "Sequence Viewer" tab:
1. View the sequence to export
2. Click "Export" button
3. Select save location

## Workflow Example

### Scenario: Compare CDS sequences from multiple species

1. Place all CDS FASTA files in same folder (e.g., `CDS/`)
2. Launch program, select `CDS/` folder
3. Click "Start Comparison"
4. In "Common Sequences" tab, view genes common to all species
5. In "Unique Sequences" tab, view genes unique to each species
6. In "Statistics" tab, view detailed statistical distribution
7. Click "Export Common Sequences" to save common genes
8. Click "Export Unique Sequences" to save species-specific genes

## Notes

1. **File Format**: FASTA files must consist of header starting with `>` and subsequent sequence content
2. **File Encoding**: Recommend using UTF-8 encoding
3. **Sequence Comparison**: Based on header name matching, not sequence content
4. **Large Files**: Program supports processing files with large number of sequences, but may take longer
5. **File Overwriting**: Export will overwrite same-name files, please confirm before proceeding

## Troubleshooting

### Problem: Program won't start
- Ensure Python 3.6 or higher is installed
- Check if you can run `python --version` in command line

### Problem: Can't find FASTA files
- Ensure selected folder contains `.fasta` or `.fa` files
- Check if file extensions are correct

### Problem: Display garbled text
- Ensure FASTA file uses UTF-8 encoding
- Can use text editor (e.g., Notepad++) to view and convert encoding

### Problem: Sequence recognition incorrect
- Check if FASTA header format matches recognition rules
- Can use "Sequence Viewer" feature to view complete header information

## Technical Details

### Core Classes

- **FastaParser**: FASTA file parser, responsible for reading and parsing FASTA files
- **SequenceComparator**: Sequence comparator, responsible for comparing sequences across multiple files
- **FastaCompareApp**: Graphical interface application, providing user interaction features

### Data Structure

```python
{
    'all_keys': set(),           # All distinct sequence keys
    'common': {},                # Common sequences {key: [filenames]}
    'unique': {},                # Unique sequences {filename: set(keys)}
    'file_keys': {}              # Sequence keys per file {filename: set(keys)}
}
```

## Version History

- **v1.0.0**: Initial version
  - Support multi-file comparison
  - Common and unique sequence identification
  - Graphical interface and export features

## License

This tool is open source software, free to use and modify.

## Contact

For questions or suggestions, please contact developer.

---

**Tip**: Recommend using with "FASTA Sample Classification Tool", first compare then classify to get sample-organized sequence files.
