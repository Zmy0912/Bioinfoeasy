# GenBank/GFF3 to mVISTA Converter

A graphical tool for batch converting GenBank or GFF3 format annotation files to mVISTA annotation files.

## Features

- ✅ Graphical interface, simple and intuitive operation
- ✅ Supports GenBank format (.gb, .gbk, .gbff)
- ✅ Supports GFF3 format (.gff, .gff3)
- ✅ Supports batch conversion of multiple files
- ✅ Supports importing all files from a folder
- ✅ Automatically identifies gene strand direction (forward strand uses `>`, reverse strand uses `<`)
- ✅ Automatically groups exon, CDS and UTR features by gene
- ✅ Intelligent GFF3 file handling: prioritizes exon, uses CDS as exon if exon is not available
- ✅ Supports extracting gene names from multiple attribute fields (gene, gene_id, Name, ID, Parent, etc.)
- ✅ Real-time display of conversion progress and detailed logs
- ✅ Shows statistics after conversion is complete
- ✅ Cross-platform support (Windows, Linux, macOS)

## Installation

```bash
pip install -r requirements.txt
```

Or install individually:

```bash
pip install biopython PyQt5 bcbio-gff
```

## Usage

### Starting the Program

```bash
python genbank_to_mvista.py
```

### Operation Steps

1. **Select File Format**
   - Select the input file format at the top of the interface (GenBank or GFF3)
   - Switching formats will automatically clear the file list

2. **Add Input Files**
   - Click the "Add Files" button to select single or multiple files
   - Or click the "Add Folder" button to select a folder containing files
   - GenBank format: automatically recognizes files with extensions `.gb`, `.gbk`, `.gbff`
   - GFF3 format: automatically recognizes files with extensions `.gff`, `.gff3`

3. **Select Output Directory**
   - Click the "Browse..." button to select the output directory
   - Converted mVISTA format files will be saved in this directory

4. **Start Conversion**
   - Click the "Start Conversion" button
   - View the conversion progress in the log panel on the right
   - Wait for conversion to complete

## mVISTA Format Description

The converted files use the standard mVISTA annotation file format, displayed in groups by gene.

### Format Rules

1. **Gene Line**: Starts with `>` (forward strand) or `<` (reverse strand), contains the gene's start position, end position, and gene name
   ```
   > start_position end_position gene_name
   ```
   or
   ```
   < start_position end_position gene_name
   ```

2. **UTR Line**: Contains UTR's start position, end position, type is "utr"
   ```
   start_position end_position utr
   ```

3. **exon Line**: Contains exon's start position, end position, type is "exon"
   ```
   start_position end_position exon
   ```

4. **Gene Separator**: Each gene block is separated by a blank line

### File Example

```
< 106481 116661 gene1
106481 106497 utr
107983 108069 exon
109884 110033 exon
111865 112023 exon

> 39424 42368 gene2
39424 39820 exon
41401 42368 exon

> 77817 81088 gene3
77817 78820 utr
79538 80107 exon
```

### Feature Type Description

mVISTA format only contains the following feature types:

- **gene**: Gene (use `>` at the beginning of the line to indicate forward strand, `<` to indicate reverse strand)
- **exon**: Exon
- **utr**: Untranslated region (including 5'UTR and 3'UTR)

### Supported Input Feature Types

**GenBank Format**:
- gene, exon, 5'UTR, 3'UTR, utr

**GFF3 Format**:
- gene, exon, CDS, 5'UTR, 3'UTR, UTR, five_prime_utr, three_prime_utr

**Note**: For GFF3 format, the program prioritizes exon features; if no exon is available, it will automatically use CDS regions as exons.

## Output File Naming

Converted file naming rule:
```
{input_filename}_mvista.txt
```

For example:
- GenBank input: `chloroplast.gb` → output: `chloroplast_mvista.txt`
- GFF3 input: `annotation.gff3` → output: `annotation_mvista.txt`

## Coordinate System

- Uses **1-based** coordinate system (consistent with GenBank format)
- Start position counts from 1
- Start position ≤ End position

## Notes

1. **File Format Requirements**:
   - GenBank files must be in standard GenBank format (containing LOCUS, FEATURES, ORIGIN tags, etc.)
   - GFF3 files must comply with GFF3 specifications (9 columns tab-separated)

2. **Coordinate System**:
   - GenBank: Uses 1-based coordinate system
   - GFF3: Already 1-based coordinate system (consistent with mVISTA)

3. **Output Overwrite**: If the output file already exists, it will be overwritten

4. **Gene Grouping**:
   - The program automatically groups by gene name and classifies related exon, UTR and other features under the corresponding gene
   - GenBank: Uses gene, locus_tag, protein_id and other fields as grouping identifiers
   - GFF3: Prioritizes gene, gene_id, Name fields, or uses ID or Parent fields if not available

5. **Strand Direction**:
   - GenBank: Extracts strand direction from the feature's strand attribute
   - GFF3: Extracts from the 7th column strand field (+, -, .)
   - Uses `>` to indicate forward strand, `<` to indicate reverse strand

6. **Gene Line**: If there is no explicit gene type feature, the program will infer the gene's start and end positions based on the position ranges of all related features

7. **GFF3 Special Handling**:
   - Prioritizes exon features as exons
   - If no exon but has CDS, automatically uses CDS regions as exons
   - UTR features support multiple naming conventions: 5UTR, 3UTR, UTR, five_prime_utr, three_prime_utr

## Example Data

### GenBank Input Example

```
LOCUS       NC_001665             155500 bp    DNA     circular PLN 15-MAR-2024
DEFINITION  Nicotiana tabacum chloroplast, complete genome.
SOURCE      Nicotiana tabacum
  ORGANISM  Nicotiana tabacum
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta;
            Tracheophyta; Spermatophyta; Magnoliopsida;
            eudicotyledons; Core eudicots; Gunneridae; Pentapetalae;
            asterids; lamiids; Solanales; Solanaceae; Nicotiana.
FEATURES             Location/Qualifiers
     gene            complement(106481..116661)
                     /gene="gene1"
                     /locus_tag="gene1"
     mRNA            complement(106481..116661)
                     /gene="gene1"
     exon            complement(107983..108069)
                     /gene="gene1"
     exon            complement(109884..110033)
                     /gene="gene1"
     5'UTR           complement(106481..106497)
                     /gene="gene1"
     gene            39424..42368
                     /gene="gene2"
                     /locus_tag="gene2"
     exon            39424..39820
                     /gene="gene2"
     exon            41401..42368
                     /gene="gene2"
```

### GFF3 Input Example

```
chr1	NCBI	gene	39424	42368	.	+	.	ID=gene2;Name=gene2
chr1	NCBI	mRNA	39424	42368	.	+	.	ID=mRNA2;Parent=gene2
chr1	NCBI	exon	39424	39820	.	+	.	ID=exon1;Parent=mRNA2
chr1	NCBI	exon	41401	42368	.	+	.	ID=exon2;Parent=mRNA2
chr1	NCBI	gene	106481	116661	.	-	.	ID=gene1;Name=gene1
chr1	NCBI	mRNA	106481	116661	.	-	.	ID=mRNA1;Parent=gene1
chr1	NCBI	exon	107983	108069	.	-	.	ID=exon3;Parent=mRNA1
chr1	NCBI	exon	109884	110033	.	-	.	ID=exon4;Parent=mRNA1
chr1	NCBI	5'UTR	106481	106497	.	-	.	ID=utr1;Parent=mRNA1
```

### mVISTA Output Example

```
< 106481 116661 gene1
106481 106497 utr
107983 108069 exon
109884 110033 exon

> 39424 42368 gene2
39424 39820 exon
41401 42368 exon
```

## Troubleshooting

### Q1: Conversion failed, message "No valid GenBank record found"

**Cause**: The file is not in valid GenBank format

**Solution**:
- Confirm that the file extension is `.gb`, `.gbk` or `.gbff`
- Check if the file content contains a complete GenBank format (LOCUS, ORIGIN tags, etc.)
- You can use a text editor to open the file and check if the format is correct

### Q2: GFF3 file conversion failed

**Cause**: The file is not in valid GFF3 format

**Solution**:
- Confirm that the file extension is `.gff` or `.gff3`
- Check if the file contains 9 columns of tab-separated data
- Skip comment lines starting with `#`
- Ensure the 3rd column (feature type) contains gene, exon, CDS or UTR related features

### Q3: Some genes do not appear in the output file

**Cause**: The file may not have explicit gene, exon or UTR features

**Solution**:
- Check if the file has these types of features
- GenBank: Ensure feature type naming follows standards (e.g., "gene", "exon", "5'UTR", "3'UTR", etc.)
- GFF3: Ensure the 3rd column feature type is gene, exon, CDS, 5UTR, 3UTR, etc.

### Q4: Gene strand direction is incorrect

**Cause**: The gene feature in the file does not have explicit strand direction information

**Solution**:
- GenBank: Check the strand attribute of the gene feature
- GFF3: Check the 7th column strand field (should be `+`, `-` or `.`)
- If there is no strand direction information, the program defaults to forward strand (`>`)

### Q5: GFF3 file has no exon, only CDS

**Cause**: Some GFF3 files only contain CDS features and do not contain exons

**Solution**: The program will automatically treat CDS regions as exons, no manual intervention required

### Q6: Gene name displays as feature_xxx

**Cause**: Unable to extract a valid gene name from the feature

**Solution**:
- Check if the 9th column attributes of the GFF3 file contain gene, gene_id, Name, ID fields
- Check if the GenBank file feature qualifiers contain gene, locus_tag, protein_id fields
- You can manually edit the file to add appropriate attribute fields

### Q7: Program startup failed

**Cause**: Missing dependency libraries

**Solution**:
```bash
pip install -r requirements.txt
```

## Technical Implementation

### Dependencies

| Library Name | Purpose |
|--------------|---------|
| biopython | GenBank file parsing |
| PyQt5 | Graphical interface |
| bcbio-gff | GFF3 file processing (optional) |

### Core Algorithms

#### GenBank Format Processing
1. **File Reading**: Use BioPython's SeqIO module to read GenBank files
2. **Feature Extraction**: Iterate through all features, only extract gene, exon and UTR types
3. **Gene Grouping**: Classify exons and UTRs under corresponding genes based on gene name
4. **Strand Direction Judgment**: Get strand direction from the gene feature's strand attribute
5. **Coordinate Conversion**: Convert 0-based coordinates to 1-based coordinates
6. **Data Sorting**: Sort output by gene start position and feature position
7. **File Writing**: Write to text file in mVISTA format

#### GFF3 Format Processing
1. **File Reading**: Read GFF3 file line by line, skip comment lines
2. **Line Parsing**: Parse 9 columns of tab-separated data
3. **Feature Filtering**: Only process gene, exon, CDS, UTR related features
4. **Attribute Parsing**: Parse the 9th column attribute fields (key=value pairs, semicolon separated)
5. **Gene Grouping**: Classify features based on attribute fields (gene, ID, Parent, etc.)
6. **Smart Exon**: Prioritize exon, use CDS as exon if exon is not available
7. **Strand Direction Judgment**: Get from the 7th column strand field (+, -, .)
8. **Data Sorting**: Sort output by gene start position and feature position
9. **File Writing**: Write to text file in mVISTA format

### GFF3 Format Specification

GFF3 files are text files with 9 columns separated by tabs:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | seqid | Sequence ID |
| 2 | source | Source |
| 3 | type | Feature type |
| 4 | start | Start position (1-based) |
| 5 | end | End position (1-based) |
| 6 | score | Score (can be `.`) |
| 7 | strand | Strand direction (`+`, `-`, `.`) |
| 8 | phase | Phase (0, 1, 2, `.`) |
| 9 | attributes | Attributes (semicolon-separated key-value pairs) |

Attribute field format: `key1=value1;key2=value2;...`

## mVISTA Format Specification

According to mVISTA official documentation:

> If a gene annotation of the sequence is available you can submit it in a simple plain text format to be displayed on the plot
>
> Each gene is defined by its start and end coordinates on the sequence, and the name listed on one line. A greater than (>) or less than (<) sign should be placed before this line to indicate plus or minus strand, although the numbering should be according to the plus strand. The exons are listed individually with the word "exon", after the start and end coordinates of each exon. UTRs are annotated the same way exons, with the word "utr" replacing "exon".

## Version History

### v3.0 (2026-02-01)
- Added GFF3 file format support
- Added file format selection feature (GenBank/GFF3)
- Intelligent handling of GFF3 file exon and CDS features
- Support for extracting gene names from multiple GFF3 attribute fields
- Updated file selection dialog to match corresponding format
- Optimized file selection logic (clear list when switching formats)
- Updated dependencies (added bcbio-gff)

### v2.0 (2026-02-01)
- Updated to correct mVISTA format specification
- Gene line uses `>` (forward strand) or `<` (reverse strand) to indicate direction
- Only outputs three feature types: gene, exon and UTR
- Removed feature type filter option (format is now fixed)
- Fixed missing strand attribute error

### v1.0 (2026-02-01)
- Initial version
- Supports basic GenBank to mVISTA format conversion
- Graphical interface
- Batch processing functionality

## License

This tool is for learning and research purposes only.

## Contact

If you have any questions or suggestions, please feel free to provide feedback.
