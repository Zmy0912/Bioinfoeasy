# MISA SSR Parser Tool

## Overview

This program parses MISA website output files for genomic simple sequence repeats (SSRs), with the following features:

- Read and parse MISA output files (tab-separated format)
- Count SSRs of different types
- Separately count mono-, di-, tri-, tetra-, penta-, and hexanucleotide repeats for each sequence
- Classify and count by repeat unit size
- Generate detailed tables, summary tables, and per-sequence statistics tables
- Export to Excel and CSV formats

## MISA File Format

The program supports standard MISA output format (tab-separated):

```
ID        SSR nr.  SSR type  SSR      size  start  end
NC1       1        p2        (TA)8    16    26     41
NC1       2        p1        (T)13    13    286    298
```

**Field descriptions:**
- `ID`: Sequence name
- `SSR nr.`: SSR number
- `SSR type`: SSR type (p1-p6 for simple SSRs, c/c* for compound SSRs)
- `SSR`: SSR sequence, e.g., `(TA)8` means TA repeats 8 times
- `size`: Total SSR length
- `start`: Start position
- `end`: End position

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```bash
python misa_parser.py -i input.misa
```

### Specify Output File Prefix

```bash
python misa_parser.py -i input.misa -o result
```

### Display Statistics Only, Do Not Export

```bash
python misa_parser.py -i input.misa --no-export
```

### Command Line Arguments

| Argument | Type | Description | Default |
|----------|------|-------------|---------|
| `-i, --input` | Required | Path to MISA output file | - |
| `-o, --output` | Optional | Output file prefix | `misa_result` |
| `--no-export` | Optional | Display statistics only, do not export files | False |

## Output Files

The program will generate the following files:

1. `{prefix}_detailed.xlsx` - Detailed table (Excel format)
   - Contains detailed information for each SSR
   - Columns: Sequence_ID, SSR_Number, SSR_Type, SSR, SSR_Length, Unit_Size, Start, End

2. `{prefix}_by_sequence.xlsx` - Per-sequence statistics table (Excel format)
   - SSR statistics for each sequence
   - Columns: Sequence_Name, Mononucleotide, Dinucleotide, Trinucleotide, Tetranucleotide, Pentanucleotide, Hexanucleotide, Total_SSR

3. `{prefix}_summary.xlsx` - Summary table (Excel format)
   - Counts and percentages by repeat unit size
   - Columns: Category, Count, Percentage(%)

4. `{prefix}_detailed.csv` - Detailed table (CSV format)

5. `{prefix}_by_sequence.csv` - Per-sequence statistics table (CSV format)

## Output Statistics

The program will display the following statistics in the console:

- Total SSR count
- Number of sequences analyzed
- Global statistics by repeat unit size
- Per-sequence statistics for mono-, di-, tri-, tetra-, penta-, and hexanucleotide repeats

## Statistics Description

### Per-Sequence Statistics Table Description

This table separately counts the following information for each sequence:

| Column | Description |
|--------|-------------|
| **Mononucleotide** | Number of SSRs with 1 bp repeat unit (e.g., A/T/C/G repeats, corresponding to p1 type) |
| **Dinucleotide** | Number of SSRs with 2 bp repeat unit (e.g., AT/AG/AC repeats, corresponding to p2 type) |
| **Trinucleotide** | Number of SSRs with 3 bp repeat unit (corresponding to p3 type) |
| **Tetranucleotide** | Number of SSRs with 4 bp repeat unit (corresponding to p4 type) |
| **Pentanucleotide** | Number of SSRs with 5 bp repeat unit (corresponding to p5 type) |
| **Hexanucleotide** | Number of SSRs with 6 bp repeat unit (corresponding to p6 type) |
| **Total_SSR** | Total number of SSRs for the sequence (including simple and compound SSRs) |

### SSR Type Description

| Type | Description |
|------|-------------|
| **p1-p6** | Simple SSRs, where the number indicates the repeat unit size in bp (e.g., p2 = 2 bp repeat unit) |
| **c** | Compound SSR, containing multiple different repeat unit types |
| **c*** | Compound SSR with interruptions |

> **Note:** For compound SSRs (c/c*), the program will identify all repeat units within them and classify based on the smallest repeat unit size.

## Notes

1. Input file must be in standard MISA format (tab-separated)
2. First line is the header, data starts from the second line
3. Recommended file encoding is UTF-8
4. Input file can contain merged results from multiple sequences; the program will automatically group and count by sequence name (ID column)
5. For compound SSRs, all repeat units will be counted, but classification uses the smallest repeat unit size
6. Ensure Python 3.6+ and required dependencies are installed before running the program

## Example

```
==========================================================================================
MISA SSR Statistics
==========================================================================================
Total SSR Count: 514
Number of Sequences: 8

Global Statistics by Unit Size:
--------------------------------------------------
  1bp: 385 (74.90%)
  2bp: 78 (15.18%)
  3bp: 9 (1.75%)
  4bp: 28 (5.45%)
  5bp: 14 (2.72%)

Statistics by Sequence:
------------------------------------------------------------------------------------------
Sequence_Name     Mono     Di       Tri      Tetra    Penta    Hexa     Total    
------------------------------------------------------------------------------------------
NC1              49       11       1        3        1        0        65    
NC2              35       10       0        6        1        0        52    
NC3              33       8        1        2        0        0        44    
NC4              32       8        1        2        0        0        43    
NC5              33       8        2        5        2        0        50    
NC6              32       8        2        5        2        0        49    
NC7              30       9        3        4        3        0        49    
NC8              33       8        2        5        2        0        50    
==========================================================================================
```
