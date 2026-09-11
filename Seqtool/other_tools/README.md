# Sequence Tools (Other Tools)

A set of small command-line / GUI utilities for everyday FASTA handling.

## Contents

| Directory | Description |
|-----------|-------------|
| [`EN/`](EN/) | English versions of the tools |
| [`CN/`](CN/) | Chinese versions of the tools |

Both directories contain the same eight scripts:

| Script | Description |
|--------|-------------|
| `extract_sequences.py` | Extract sequences by ID / pattern |
| `extract_sequence_info.py` | Extract summary information (length, GC, etc.) |
| `rename_sequences.py` | Rename sequence headers |
| `remove_sequences.py` | Remove sequences by ID / pattern |
| `merge_sequences.py` | Merge sequences |
| `complement_sequences.py` | Generate reverse-complement sequences |
| `fasta_cleaner_gui.py` | GUI FASTA validator and cleaner |
| `text_to_uppercase.py` | Convert sequence text to uppercase |

## Requirements

- Python 3.6+
- `biopython`
- `tkinter` (bundled with CPython) for the GUI tool

```bash
pip install biopython
```

## Usage

```bash
cd EN
python extract_sequences.py
```

Back to [Seqtool](../README.md) · [project root](../../README.md).
