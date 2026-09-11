# FASTA Sample Classification Tool

A desktop application that classifies sequences in mixed FASTA files by sample names, supporting selective export of sequences for one or multiple samples.

## Features

- **Automatic Classification**: Automatically recognize sample names from FASTA sequence headers and classify
- **Flexible Export**: Support exporting all samples or selectively exporting some samples
- **Visual Interface**: Provide intuitive sample list and sequence viewing features
- **Batch Processing**: Process large files with thousands of sequences at once
- **Sequence Preview**: Preview sequence content of any sample before export
- **Statistics**: Display sequence count and total length for each sample

## System Requirements

- Python 3.6 or higher
- Windows / macOS / Linux
- No additional dependencies (uses Python standard library)

## Installation

1. Ensure Python 3.x is installed
2. Place `fasta_classifier.py` and `fasta_classifier.bat` in the same folder
3. Double-click `fasta_classifier.bat` to launch the program

## Quick Start

### Basic Usage Flow

1. **Launch Program**: Double-click `fasta_classifier.bat`
2. **Select File**: Click "Select File" button, select FASTA file to process
3. **View Samples**: View all recognized samples in the left sample list
4. **Select Export**: Check samples to export (or export all)
5. **Export Sequences**: Click export button, select save location
6. **Preview Content**: Click sample to preview all sequences of that sample

### Interface Description

#### Top Toolbar
- **Select File**: Select FASTA file to process
- **Export All Samples**: Export all samples as separate FASTA files
- **Export Selected Samples**: Export only checked samples
- **Select All**: Select all samples
- **Deselect All**: Deselect all selections
- **Clear**: Clear all data, restart

#### File Information Area
Displays information of currently loaded file.

#### Left: Sample List
- Displays all recognized sample names
- Displays sequence count for each sample
- Displays total length of sequences for each sample
- Click checkbox in first column (☐/☑) to select/deselect sample
- Click sample name to view sequence details of that sample on right

#### Right: Sequence Details
- **Sequence List**: Display all sequences of selected sample
- **Sequence Content**: Display complete content of selected sequence

#### Bottom Status Bar
Displays current operation status and prompt information.

## Sample Name Recognition Rules

Program supports multiple header formats, tries in priority order:

### Format 1: Pipe separator + .fasta suffix
```
>rps12_join{...}|A10.fasta     →  Sample: A10
>cemA_[...] |A11.fasta         →  Sample: A11
```

### Format 2: Pipe separator, no suffix
```
>gene_name|A10                →  Sample: A10
>sequence_id|sample1          →  Sample: sample1
```

### Format 3: Space separator + .fasta suffix
```
>gene_name A10.fasta          →  Sample: A10
```

### Format 4: End directly with .fasta
```
>rps12_A10.fasta              →  Sample: rps12_A10
```

### Format 5: Other cases
If above formats don't match, try to extract last part containing letters and numbers combination. If still can't recognize, mark as `unknown`.

**Note**: Program automatically removes `.fasta` and `.fa` suffixes.

## Usage Examples

### Example 1: Process CDS sequences from multiple species

Assume there's a `redtr.fasta` file containing CDS sequences from multiple samples:

```
>rps12_join{...}|A10.fasta
ATGCGTACG...
>cemA_[...] |A11.fasta
ATGCGTACG...
>psbA_[76:1137](-)|A12.fasta
ATGCGTACG...
```

**Steps**:
1. Launch program, select `redtr.fasta`
2. Program automatically recognizes samples: A10, A11, A12, etc.
3. Click "Select All" or check needed samples
4. Click "Export Selected Samples"
5. Select save folder
6. Program generates following files:
   - `A10.fasta` (contains all A10 sample sequences)
   - `A11.fasta` (contains all A11 sample sequences)
   - `A12.fasta` (contains all A12 sample sequences)

### Example 2: Selective export of specific samples

1. Load FASTA file
2. Check only needed samples in sample list (e.g., A10 and A11)
3. Click "Export Selected Samples"
4. Only checked samples will be exported

### Example 3: Preview sample sequences

1. Load FASTA file
2. Click any sample in left list (e.g., A10)
3. Right side displays all sequence list of that sample
4. Click sequence to view complete sequence content

## Export Features Explained

### Export All Samples
- Create separate FASTA file for each recognized sample
- File name format: `sample_name.fasta`
- Preserve original FASTA header format

### Export Selected Samples
- Only export checked samples
- Unchecked samples won't be exported
- Other features same as "Export All Samples"

## Workflow Recommendations

### Recommended Workflow: Compare + Classify

1. **Step 1: Use FASTA Sequence Comparison Tool**
   - Compare multiple FASTA files
   - Export common or unique sequences

2. **Step 2: Use Sample Classification Tool**
   - Load exported sequence files
   - Classify by sample name
   - Export sample-organized sequence files

This gives you:
- Sequences common to all samples (e.g., common genes)
- Sequences unique to each sample
- Sample-organized FASTA files

## Notes

1. **Sample Name Uniqueness**: Ensure different samples have unique identifiers
2. **File Overwriting**: Export will overwrite same-name files, please confirm before proceeding
3. **Encoding Format**: Recommend using UTF-8 encoded FASTA files
4. **Large File Processing**: Processing large files may take longer, please be patient
5. **unknown samples**: If header can't be recognized, sample will be marked as "unknown"
6. **Spaces and Separators**: Maintain consistency in header format, helps correct sample recognition

## Troubleshooting

### Problem: Can't recognize sample name
- **Cause**: Header format doesn't match recognition rules
- **Solution**: Check header format, ensure sample identifier is present

### Problem: Wrong number of recognized samples
- **Cause**: Some sequence header formats not uniform
- **Solution**: Unify header formats for all sequences, or manually adjust

### Problem: Exported files are empty
- **Cause**: Selected samples have no sequences, or recognized as "unknown"
- **Solution**: Check sample recognition results, confirm there are valid samples

### Problem: Program won't start
- **Solution**: Ensure Python 3.6 or higher is installed

### Problem: Display garbled text
- **Solution**: Ensure FASTA file uses UTF-8 encoding

## Technical Details

### Core Classes

- **FastaSampleClassifier**: Sample sequence classifier
  - Parse FASTA files
  - Extract sample names from headers
  - Organize sequences by sample
  - Provide sample information queries

- **FastaClassifierApp**: Graphical interface application
  - Provide user interaction interface
  - Handle file selection and export
  - Display sample list and sequence details

### Sample Name Extraction Algorithm

Program uses multi-stage matching algorithm:

```python
def extract_sample_name(header):
    # Stage 1: Check pipe separator + .fasta suffix
    if '|' in header and header.endswith(('.fasta', '.fa')):
        return header.split('|')[-1].replace('.fasta', '')

    # Stage 2: Check pipe separator
    if '|' in header:
        return header.split('|')[-1].strip()

    # Stage 3: Check space separator + .fasta suffix
    # Stage 4: Check ending .fasta
    # Stage 5: Try to extract alphanumeric combination part

    return "unknown"
```

## Advanced Features

### Checkbox Operations

- **Click Checkbox**: Toggle sample selection state
- **Select All Button**: Select all samples with one click
- **Deselect All Button**: Deselect all selections with one click

### Sequence Viewing

- **Sequence List**: Display all sequence headers of current sample
- **Sequence Details**: Display complete sequence and metadata (length, etc.)
- **Real-time Preview**: Display content immediately after selecting sequence

## Best Practices

1. **Unified Header Format**: Use unified header format when generating FASTA files
2. **Clear Sample Naming**: Use clear, meaningful sample names
3. **Regular Backup**: Backup original files before export
4. **Batch Processing**: Use "Export All Samples" feature for batch processing
5. **Verify Results**: Spot-check some files after export to ensure correctness

## Common Application Scenarios

1. **Multi-sample Genomic Analysis**: Classify mixed gene sequences by sample
2. **Species Comparison Studies**: Organize sequence data from different species
3. **Population Genetics**: Organize sequences by population or population group
4. **Data Preprocessing**: Prepare sample-organized data for downstream analysis
5. **Quality Control**: Check sequence completeness for each sample

## Version History

- **v1.0.0**: Initial version
  - Automatic sample recognition and classification
  - Visual interface
  - Selective export feature
  - Sequence preview feature

## License

This tool is open source software, free to use and modify.

## Contact

For questions or suggestions, please contact developer.

---

**Tip**: Recommend using with "FASTA Sequence Comparison Tool", first compare then classify by sample to improve work efficiency.
