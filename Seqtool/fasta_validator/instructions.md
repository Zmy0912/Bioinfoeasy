# FASTA Format Validator - Detailed User Guide

## Table of Contents

- [Features](#features)
- [FASTA Format Specification](#fasta-format-specification)
- [Common Errors and Fixes](#common-errors-and-fixes)
- [Usage Steps](#usage-steps)
- [Interface Description](#interface-description)
- [Error Position Highlighting](#error-position-highlighting)
- [Advanced Features](#advanced-features)
- [FAQ](#faq)

## Features

1. Identify FASTA format files
2. Validate FASTA format correctness
3. Detailed error reporting (including line numbers and specific issues)
4. **Mark error positions in sequences (with `^` marking invalid characters)**
5. Auto-fix common FASTA format errors
6. Intuitive graphical user interface
7. Support for strict ATCG validation mode

## FASTA Format Specification

### Basic Format

- First line must start with `>` symbol, followed by sequence ID and description
- Sequence ID and description are typically separated by space
- Subsequent lines contain sequence data
- Case-insensitive (A, C, G, T are valid)
- Empty lines allowed (will be ignored)

### Standard Mode Valid Characters

- A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N, -
- Case-insensitive
- Empty lines ignored

**Meanings**:
- A = Adenine
- C = Cytosine
- G = Guanine
- T = Thymine
- U = Uracil
- R = A or G (puRine)
- Y = C or T (pYrimidine)
- And so on...

### Strict ATCG Mode Valid Characters

- A, C, G, T, -
- Only four base letters and gap symbol
- Suitable for standard DNA sequence validation only

## Common Errors and Fixes

The program can detect and fix these error types:

### 1. Line 1 must start with '>' (header line)

**Error Example**:
```
sequence without header
ATCGATCG
```

**Fix**: Automatically add header

### 2. Sequence data without header

**Error Example**:
```
>seq1
ATCGATCG
wrong_line_without_gt
```

**Fix**: Automatically add header line

### 3. Header line missing sequence ID

**Error Example**:
```
>
ATCGATCG
```

**Fix**: Auto-generate default ID (e.g., sequence_1, sequence_2...)

### 4. Contains invalid characters

**Standard mode error example**:
```
>seq1
ATCG12345
```
- **Error**: Sequence contains numbers, special symbols
- **Fix**: Remove non-sequence characters (keep ACGTURYSWKMBDHVN-)
- **Display**: `Line X: Contains invalid characters 'XXX' (expected ACGTURYSWKMBDHVN)`

**Strict mode error example**:
```
>seq1
ATCGNNNNURYS
```
- **Error**: Sequence contains N, U, R, etc. (non-ATCG)
- **Fix**: Remove non-sequence characters based on current mode
- **Display**: `Line X: Contains invalid characters 'XXX' (expected ATCG)`

### 5. No sequence headers found

**Error Example**:
```
ATCGATCG
GCTAGCTA
```

**Fix**: Add headers for each sequence segment

### 6. File is empty

- **Error message**: File is empty
- **Cannot be fixed automatically**

## Usage Steps

### 1. Select File

Click "Select File" button, choose FASTA file

**Supported formats**: .fasta, .fa, .fas, .fna, .ffn, .faa, .frn, .txt

### 2. Choose Validation Mode (Optional)

Select validation mode in "Validation Options" area:

- **Standard Mode**: All IUPAC nucleotide codes (ACGTURYSWKMBDHVN)
- **Strict ATCG Mode**: Only A, C, G, T
- Toggle "Strict ATCG Mode" checkbox to switch

### 3. Enable Error Highlighting (Optional)

- Check "Highlight Error Positions" checkbox
- Error markers displayed in "Error Position Markers" text box

### 4. Validate Format

Click "Validate Format" button

**Display**:
- Error information in bottom "Validation Results/Error Information" text box
- Error position markers in "Error Position Markers" text box

### 5. Fix File

If errors found, click "Fix File"

**Result**:
- Fixed content shown in "Fixed Content" text box
- Fix based on current selected mode

### 6. Save Result

Review and click "Save Fixed"

**Save options**:
- Choose save path and filename
- Default: add `_fixed` suffix to original filename

### 7. Clear and Reset

Click "Clear" to process another file

## Interface Description

### Layout (top to bottom)

```
┌─────────────────────────────────────────┐
│ [Select File][Validate][Fix][Save][Clear] │  ← Toolbar
├─────────────────────────────────────────┤
│ Validation Options                     │
│ ☑ Strict ATCG Mode                    │
│ ☑ Highlight Error Positions (mark with ^)│
│ Current Mode: Standard (ACGTURYSWKMBDHVN) │
├─────────────────────────────────────────┤
│ Current File: xxx.fasta                │  ← File Info
├─────────────────────────────────────────┤
│ ┌─ Original File Content ──────────────┐ │  ← Box 1
│ │ >seq1                              │ │
│ │ ATCGATCG                            │ │
│ └─────────────────────────────────────┘ │
│ ┌─ Error Position Markers (^) ─────┐ │  ← Box 2
│ │ >seq1                              │ │
│ │ Line 2: ATCG12345                 │ │
│ │         ^^^^^                      │ │
│ └─────────────────────────────────────┘ │
│ ┌─ Fixed Content ───────────────────────┐ │  ← Box 3
│ │ >seq1                              │ │
│ │ ATCG                               │ │
│ └─────────────────────────────────────┘ │
│ ┌─ Validation Results / Errors ───────┐ │  ← Box 4
│ │ === FASTA Format Validation ===    │ │
│ │ [X] Found 1 error                 │ │
│ │ 1. Line 2: Contains invalid...  │ │
│ └─────────────────────────────────────┘ │
├─────────────────────────────────────────┤
│ Validation complete: 1 error found     │  ← Status Bar
└─────────────────────────────────────────┘
```

### Window Features

- ✅ Resizable: Drag separators to adjust heights
- ✅ Auto-scroll: Auto-scroll to error area after validation
- ✅ Color coding:
  - 🔵 Blue bold: Headers
  - 🔴 Red: Errors
  - 🟢 Green: Success
  - 🟠 Orange: Warnings
  - ⚫ Black: Normal

## Error Position Highlighting

### Feature Overview

Mark invalid characters in sequences with `^` symbols for clear visualization.

### Enable

1. Check "Highlight Error Positions" checkbox
2. Click "Validate Format"
3. Error markers displayed in second text box

### Format

```
>seq1_with_errors
Line 2: ATCGATCG12345
            ^^^^^
            ^ Invalid characters: 12345

Line 3: ATCG@#$%^&*()
            ^^^^^^^^^
            ^ Invalid characters: #$%&()*@^
```

### Color Coding

- 🔵 **Blue**: Sequence headers
- ⚫ **Black**: Original sequence lines
- 🔴 **Red**: Error marker lines (`^`)
- ⚪ **Gray**: Informational text

### Symbol Meaning

- **`^`**: Marks invalid character position
- **Space**: Character is valid

### Examples

#### Example 1: Standard Mode

**Input**:
```
>seq1
ATCGATCG12345
ATCG@#$%^&*()
ATCGNNNNURYS
```

**Markers**:
```
>seq1
Line 2: ATCGATCG12345
            ^^^^^
            ^ Invalid characters: 12345

Line 3: ATCG@#$%^&*()
            ^^^^^^^^^
            ^ Invalid characters: #$%&()*@^
```

#### Example 2: Strict ATCG Mode

**Markers**:
```
>seq1
Line 2: ATCGATCG12345
            ^^^^^
            ^ Invalid characters: 12345

Line 3: ATCG@#$%^&*()
            ^^^^^^^^^
            ^ Invalid characters: #$%&()*@^

Line 4: ATCGNNNNURYS
            ^^^^^^^^
            ^ Invalid characters: NRSUY
```

## Advanced Features

### Batch Processing Suggestions

- Use scripts for multiple files
- Manually process important files
- Prioritize files with most errors

### Custom Validation Rules

Current versions provide two modes:
1. **Standard Mode**: Most FASTA files
2. **Strict ATCG Mode**: Pure DNA sequences

Modify `valid_chars` set in source code for custom rules.

### Performance Tips

- Large files: Disable error highlighting for speed
- Frequent validation: Review all errors before fixing

## FAQ

### Q: Why is error marker box empty?

**A**: Possible reasons:
1. "Highlight Error Positions" not checked
2. File has no errors
3. File is empty

### Q: Are `^` markers aligned?

**A**: Yes, auto-aligned to mark exact positions.

### Q: Do different modes show different errors?

**A**: Yes, strict ATCG mode marks more errors.

### Q: Can I resize error marker box?

**A**: Yes, drag separators to adjust heights.

### Q: Is format preserved after fixing?

**A**: Valid sequence info preserved, invalid chars removed.

### Q: Max file size?

**A**: No theoretical limit, but <100MB recommended.

### Q: Can I undo fixes?

**A**: Not currently. Backup before fixing.

### Q: Why are some chars valid in standard but invalid in strict?

**A**: Standard supports IUPAC codes (including degenerate bases), strict only A/C/G/T.

### Q: Which mode should I use?

**A**:
- **Standard**: Most FASTA files, especially with degenerate bases
- **Strict ATCG**: Sequences requiring only standard DNA bases

## Test Files

| File | Description | Purpose |
|------|-------------|---------|
| `correct.fasta` | Correct FASTA format | Test correct file handling |
| `error1.fasta` | Missing `>` | Test missing header detection |
| `error2.fasta` | Invalid characters | Test invalid char detection/fix |
| `error3.fasta` | Empty header | Test empty header detection |
| `error_markers_test.fasta` | Error markers | Test error highlighting |
| `non_atcg.fasta` | Non-ATCG chars | Test strict ATCG mode |

### Using Test Files

1. Launch program
2. Click "Select File"
3. Choose any test file
4. Click "Validate Format"
5. Try fixing and compare results

## Important Notes

1. Backup before fixing
2. Manually review results before saving
3. Manually edit complex format errors
4. Program preserves valid sequence info
5. Large files may take longer

## System Requirements

- Python 3.x
- Windows OS
- UTF-8, GBK, GB2312, Latin-1 encodings

---

**Last Updated**: 2026-02-03
**Version**: 3.0
