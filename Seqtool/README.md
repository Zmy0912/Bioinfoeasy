# Seqtool

Sequence processing, validation, comparison and format-conversion tools.

## Contents

| Module | Description |
|--------|-------------|
| [`bed_converter/`](bed_converter/) | BED/GFF3 coordinate converters (`gene_to_0`, `gff3_to_bed_converter`), CN/EN |
| [`comparison/`](comparison/) | FASTA sample classification and sequence comparison |
| [`fasta_header_renamer/`](fasta_header_renamer/) | Batch FASTA header ⇄ filename renaming with a multi-sequence editor (CN/EN); standard library only |
| [`fasta_joiner/`](fasta_joiner/) | Merge multiple FASTA files (`v1.0.0`, `v2.0.0`) |
| [`fasta_validator/`](fasta_validator/) | FASTA format validation and fixing |
| [`genbank_to_mvista_converter/`](genbank_to_mvista_converter/) | GenBank/GFF3 to mVISTA converter (`v1.0.0`, `v1.1.0`) |
| [`other_tools/`](other_tools/) | Sequence extraction, renaming, removal, complement, cleaner GUI, etc. (CN/EN) |

## Common dependencies

Most tools are built on Biopython; some also use `pandas` / `openpyxl`.
`fasta_header_renamer` needs no third-party packages at all — it runs on the
standard library alone.

```bash
pip install biopython pandas openpyxl
```

See each sub-module's `README.md` and `requirements.txt` for details.

Back to the [project root](../README.md).
