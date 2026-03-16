# GenBank/GFF3 to mVISTA Converter

A graphical tool for batch converting GenBank or GFF3 format annotation files to mVISTA annotation file format.

## Features

- ✅ Graphical interface, simple and intuitive operation
- ✅ Supports GenBank format (.gb, .gbk, .gbff)
- ✅ Supports GFF3 format (.gff, .gff3)
- ✅ Batch conversion of multiple files
- ✅ Import all files from a folder
- ✅ Automatic gene strand direction recognition (use `>` for plus strand, `<` for minus strand)
- ✅ Automatic grouping of exon, CDS, and UTR features by gene
- ✅ Intelligent GFF3 processing: prioritize exon, use CDS as exon if no exon available
- ✅ Support multiple attribute fields for gene name extraction (gene, gene_id, Name, ID, Parent, etc.)
- ✅ Comprehensive RNA feature support: rRNA, tRNA, mRNA, snRNA, scRNA, ncRNA, miRNA, lnc_RNA
- ✅ RNA features treated as independent genes with automatic RNA type annotation
- ✅ Real-time conversion progress and detailed logs
- ✅ Statistics display after conversion completion
- ✅ Cross-platform support (Windows, Linux, macOS)

## Installation

```bash
pip install -r requirements.txt
```

Or install individually:

```bash
pip install biopython PyQt5
```

## Usage

### Launch the Program

```bash
python genbank_to_mvista.py
```

### Operation Steps

1. **Select File Format**
   - Choose the input file format (GenBank or GFF3) at the top of the interface
   - Switching formats will automatically clear the file list

2. **Add Input Files**
   - Click "Add Files" button to select single or multiple files
   - Or click "Add Folder" button to select a folder containing files
   - GenBank format: automatically recognizes `.gb`, `.gbk`, `.gbff` format files
   - GFF3 format: automatically recognizes `.gff`, `.gff3` format files

3. **Select Output Directory**
   - Click "Browse..." button to select output directory
   - Converted mVISTA format files will be saved in this directory

4. **Start Conversion**
   - Click "Start Conversion" button
   - Check the log panel on the right for conversion progress
   - Wait for conversion to complete

## mVISTA Format Specification

The converted file uses the standard mVISTA annotation file format, organized by gene groups.

### Format Rules

1. **Gene Line**: Starts with `>` (plus strand) or `<` (minus strand), containing gene start position, end position, and gene name
   ```
   > start_position end_position gene_name
   ```
   or
   ```
   < start_position end_position gene_name
   ```

2. **UTR Line**: Contains UTR start position, end position, type is "utr"
   ```
   start_position end_position utr
   ```

3. **Exon Line**: Contains exon start position, end position, type is "exon"
   ```
   start_position end_position exon
   ```

4. **Gene Separation**: Empty line between each gene block

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

> 1234 1567 5S_RRNA (RRNA)

< 2345 3210 tRNA_Tyr (TRNA)
```

### Feature Type Description

mVISTA format only includes the following feature types:

- **gene**: Gene (use `>` at line start for plus strand, `<` for minus strand)
- **exon**: Exon
- **utr**: Untranslated region (including 5'UTR and 3'UTR)
- **RNA features**: RNA types are displayed in parentheses after gene name, e.g., `gene_name (RRNA)`

### Supported Input Feature Types

**GenBank format**:
- Basic features: gene, exon, 5'UTR, 3'UTR, utr
- RNA features: rna, rRNA, tRNA, snRNA, scRNA, mRNA, ncRNA, miRNA, lnc_RNA

**GFF3 format**:
- Basic features: gene, exon, CDS, 5UTR, 3UTR, UTR, five_prime_utr, three_prime_utr
- RNA features: rna, rrna, trna, snrna, scrna, mrna, ncrna, mirna, lnc_rna

**Note**:
- For GFF3 format, the program prioritizes exon features; if no exon, CDS regions are automatically used as exons
- RNA features are treated as independent genes and RNA types are annotated in the output file

## Output File Naming

Converted file naming convention:
```
{input_filename}_mvista.txt
```

For example:
- GenBank input: `chloroplast.gb` → Output: `chloroplast_mvista.txt`
- GFF3 input: `annotation.gff3` → Output: `annotation_mvista.txt`

## Coordinate System

- Uses **1-based** coordinate system (consistent with GenBank format)
- Start position counts from 1
- Start position ≤ End position

## Notes

1. **File Format Requirements**:
   - GenBank files must be in standard GenBank format (containing LOCUS, FEATURES, ORIGIN tags, etc.)
   - GFF3 files must comply with GFF3 specification (9 tab-separated columns)

2. **Coordinate System**:
   - GenBank: Uses 1-based coordinate system
   - GFF3: Already uses 1-based coordinate system (consistent with mVISTA)

3. **Output Overwrite**: If output file already exists, it will be overwritten

4. **Gene Grouping**:
   - The program automatically groups by gene name and classifies related exon, UTR and other features under the corresponding gene
   - GenBank: Uses gene, locus_tag, protein_id fields as grouping identifiers
   - GFF3: Prioritizes gene, gene_id, Name fields; if not available, uses ID or Parent fields

5. **Strand Direction**:
   - GenBank: Extracted from feature's strand attribute
   - GFF3: Extracted from column 7 strand field (+, -, .)
   - Use `>` for plus strand, `<` for minus strand

6. **Gene Line**: If there is no explicit gene type feature, the program infers gene start and end positions based on the position range of all related features

7. **GFF3 Special Handling**:
   - Prioritize using exon features as exons
   - If no exon but CDS exists, automatically use CDS regions as exons
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

> 1234 1567 5S_RRNA (RRNA)

< 2345 3210 tRNA_Tyr (TRNA)
```

## Troubleshooting

### Q1: Conversion failed with error "No valid GenBank records found"

**Cause**: File is not in valid GenBank format

**Solution**:
- Confirm file extension is `.gb`, `.gbk` or `.gbff`
- Check if file contains complete GenBank format (LOCUS, ORIGIN tags, etc.)
- Use text editor to open file and verify correct format

### Q2: GFF3 file conversion failed

**Cause**: File is not in valid GFF3 format

**Solution**:
- Confirm file extension is `.gff` or `.gff3`
- Check if file contains 9-column tab-separated data
- Skip lines starting with `#`
- Ensure column 3 (feature type) contains gene, exon, CDS or UTR related features

### Q3: Some genes do not appear in output file

**Cause**: File may not contain explicit gene, exon, UTR or RNA features

**Solution**:
- Check if file contains these feature types
- GenBank: Ensure feature type naming follows conventions (e.g., "gene", "exon", "5'UTR", "3'UTR", "rRNA", "tRNA")
- GFF3: Ensure column 3 feature type is gene, exon, CDS, 5UTR, 3UTR or related RNA types
- Program requires at least gene, exon/CDS or RNA features to output

### Q4: Gene strand direction is incorrect

**Cause**: File gene features do not have explicit strand direction information

**Solution**:
- GenBank: Check gene feature's strand attribute
- GFF3: Check column 7 strand field (should be `+`, `-` or `.`)
- If no strand direction information, program defaults to plus strand (`>`)

### Q5: GFF3 file has no exon, only CDS

**Cause**: Some GFF3 files only contain CDS features, not exon

**Solution**: Program automatically uses CDS regions as exons, no manual intervention needed

### Q6: Gene name displays as feature_xxx

**Cause**: Unable to extract valid gene name from features

**Solution**:
- Check if GFF3 file column 9 attributes contain gene, gene_id, Name, ID fields
- Check if GenBank file feature qualifiers contain gene, locus_tag, protein_id fields
- Can manually edit file to add appropriate attribute fields

### Q7: Program failed to start

**Cause**: May be missing dependency libraries

**Solution**:
```bash
pip install -r requirements.txt
```

## Technical Implementation

### Dependencies

| Library | Usage |
|---------|-------|
| biopython | GenBank file parsing |
| PyQt5 | Graphical interface |

### Core Algorithms

#### GenBank Format Processing
1. **File Reading**: Use BioPython's SeqIO module to read GenBank files
2. **Feature Extraction**: Traverse all features, extract gene, exon, UTR and RNA type features
3. **Gene Grouping**:
   - Normal genes: Group exon and UTR under corresponding gene based on gene name
   - RNA features: Treated as independent genes with automatic RNA type annotation
4. **Strand Direction**: Extracted from feature's strand attribute
5. **Coordinate Conversion**: Convert 0-based coordinates to 1-based coordinates
6. **Data Sorting**: Sort output by gene start position and feature position
7. **File Writing**: Write to text file in mVISTA format

#### GFF3 Format Processing
1. **File Reading**: Read GFF3 file line by line, skip comment lines
2. **Line Parsing**: Parse 9-column tab-separated data
3. **Feature Filtering**: Only process gene, exon, CDS, UTR and RNA related features
4. **Attribute Parsing**: Parse column 9 attribute fields (key=value pairs, semicolon separated)
5. **Gene Grouping**:
   - Normal genes: Group features based on attribute fields (gene, ID, Parent, etc.)
   - RNA features: Treated as independent genes with automatic RNA type annotation
6. **Intelligent Exon**: Prioritize using exon, use CDS as exon if no exon available
7. **Strand Direction**: Extracted from column 7 strand field (+, -, .)
8. **Data Sorting**: Sort output by gene start position and feature position
9. **File Writing**: Write to text file in mVISTA format

### GFF3 Format Specification

GFF3 file is a 9-column tab-separated text file:

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

### v4.0 (2026-02-06)
- Comprehensive RNA feature type support: rRNA, tRNA, mRNA, snRNA, scRNA, ncRNA, miRNA, lnc_RNA
- RNA features treated as independent genes with automatic RNA type annotation (e.g., gene_name (RRNA))
- Optimized feature extraction logic, supports more RNA-related naming variants
- Updated documentation with RNA feature processing description
- Removed unnecessary dependencies (bcbio-gff), simplified installation process

### v3.0 (2026-02-01)
- Added GFF3 file format support
- Added file format selection feature (GenBank/GFF3)
- Intelligent GFF3 exon and CDS feature processing
- Support multiple GFF3 attribute fields for gene name extraction
- Updated file selection dialogs to match corresponding formats
- Optimized file selection logic (clear list when format switches)

### v2.0 (2026-02-01)
- Updated to correct mVISTA format specification
- Gene lines use `>` (plus strand) or `<` (minus strand) to indicate direction
- Only output three feature types: gene, exon, UTR
- Removed feature type filter option (format fixed)
- Fixed strand attribute missing error

### v1.0 (2026-02-01)
- Initial version
- Basic GenBank to mVISTA format conversion support
- Graphical interface
- Batch processing feature

## License

This tool is for learning and research purposes only.

## Contact

For questions or suggestions, feedback is welcome.
