# FASTA Sequence Joiner

A Tkinter + Biopython based desktop tool that concatenates multiple sequences in FASTA files into one long sequence.

Two working modes are supported:

- **Single-file mode**: concatenate all sequences of one FASTA file into a single sequence and save to a specified output file.
- **Batch (folder) mode**: read a folder and concatenate the sequences inside **every** FASTA file in it. Output files keep the **same file names** as the inputs (e.g. input `1.fasta`, `2.fasta`, `3.fasta` produce `1.fasta`, `2.fasta`, `3.fasta`, but each now contains ONE long joined sequence).

---

## Requirements

- Python 3.6 or higher
- Biopython (see installation below)
- Works on Windows / macOS / Linux (Tkinter ships with the standard Python installation)

## Install Dependencies

```bash
pip install -r requirements.txt
```

## Run the Program

```bash
python fasta_joiner.py
```

After launch:

1. Choose a mode in the **Working Mode** section:
   - `Single FASTA File` – single-file mode
   - `Batch Folder (all FASTA files)` – batch folder mode
2. Follow Step 1 → Step 2 → Step 3 in the interface.

---

## Usage

### Single-file Mode

| Step | Action |
|------|--------|
| Step 1 | Browse to select the FASTA input file to process |
| Step 2 | Set the new sequence name (default `joined_sequence`, click `Custom` to change it); browse to choose the output file path |
| Step 3 | Click `Start Joining Sequences` |

### Batch (folder) Mode

| Step | Action |
|------|--------|
| Step 1 | Browse to select the **input folder** (FASTA files inside are detected automatically) |
| Step 2 | Browse to choose the **output folder** (created automatically if missing); optionally set a *Name suffix* — e.g. `_joined` makes the sequence ID inside `1.fasta` become `1_joined`, leave it empty to use the plain file name |
| Step 3 | Click `Start Joining Sequences` |

> ⚠ If the output folder is the same as the input folder, a confirmation dialog is shown; after confirming, the original files will be **overwritten** by the joined results.

---

## Example Result

In batch mode, suppose the folder `seqs/` contains:

```
seqs/1.fasta
>frag1
ATGC
>frag2
GCTA
>frag3
CGAT
```

After processing, `1.fasta` in the output folder becomes:

```
>1
ATGCGCTACGAT
```

i.e. each output file holds a single long sequence; fragments are concatenated directly (no separators) in the order they appear in the file.

---

## Supported File Formats

- Input: `.fasta`, `.fa`, `.fna`, `.ffn`, `.faa`, `.frn`
- Output: FASTA format (`.fasta` / `.fa`)

---

## File Structure

```
0908/
├── README.md           # This documentation
├── requirements.txt    # Dependency list
└── fasta_joiner.py     # Main program (single-file + batch modes)
```

---

## Features

- Concatenates sequences in file order without any separators
- Automatically detects and processes every FASTA file in the selected folder
- Output keeps original file names for easy downstream per-file comparison
- Optional name suffix for the joined sequence ID (batch mode)
- Real-time processing log; a single failing file does not stop the others
- Input/output path validation throughout the interface to prevent mistakes

---

## FAQ

**Q1: I get "No module named Bio" when running.**
Biopython is not installed — run `pip install -r requirements.txt` first.

**Q2: Are the original files deleted in batch mode?**
No. However, if the output folder is the same as the input folder, the original contents are overwritten by the joined results (a confirmation dialog appears first).

**Q3: Are separators inserted between sequences?**
No. Sequences are joined head-to-tail directly in their file order.

**Q4: Does folder mode process FASTA files in subfolders?**
No. Only FASTA files located **directly** in the selected folder are processed.

---

## Version History

- v1.1.0: Added batch (folder) mode (2026)
- v1.0.0: Initial single-file joining feature

---

## License

Open-source software, free to use and modify.
