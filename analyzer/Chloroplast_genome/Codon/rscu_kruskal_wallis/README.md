# RSCU Kruskal-Wallis Analysis Tool

Statistical comparison of Relative Synonymous Codon Usage (RSCU) between
groups using the Kruskal-Wallis test.

## Contents

| Directory | Description | Files |
|-----------|-------------|-------|
| [`EN/`](EN/) | English version | `rscu_kruskal_wallis.py`, `README.md`, `requirements.txt`, `example_groups.csv` |
| [`CN/`](CN/) | Chinese version | `rscu_kruskal_wallis.py`, `README.md`, `requirements.txt`, `example_groups.csv` |

## Requirements

```bash
pip install -r EN/requirements.txt
```

- `pandas>=1.3.0`
- `numpy>=1.20.0`
- `scipy>=1.7.0`
- `openpyxl>=3.0.0`
- `matplotlib>=3.5.0`

## Usage

```bash
cd EN
python rscu_kruskal_wallis.py
```

See `EN/README.md` / `CN/README.md` for the full input format description and
an example grouping file (`example_groups.csv`).

Back to [Codon analysis](../README.md) · [project root](../../../../README.md).
