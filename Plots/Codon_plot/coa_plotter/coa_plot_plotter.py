# -*- coding: utf-8 -*-
"""COA Plot Plotter - Correspondence Analysis (CodonW) supplementary figure.

Reads the correspondence-analysis output files produced by CodonW
(`genes.coa`, `eigen.coa`, `coa_raw` or a CodonW bulk `.out` file) and draws
the supplementary figure requested by the reviewer:

  * projection of the genes on the first two COA axes
    (axes are orthogonal; axis 1 explains the largest proportion of the
     total inertia, subsequent axes explain decreasing amounts);
  * a scree plot showing the fraction of variation explained by each axis;
  * scatter plots correlating the gene positions on axes 1/2 with
    nucleotide-composition parameters (GC3s, GC, GC12, ...);
  * an exported correlation table (Pearson + Spearman, with P values)
    that can be used for the discussion in section 4.2.

The program follows the structure of `pr2_plot_plotter.py` (tkinter GUI,
single-file / folder input modes, colouring by gene functional category,
labelling of outlying genes, adjustable style settings, scrollable canvas).
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import os
import glob
import math
import re


# ---------------------------------------------------------------------------
# Gene functional classification (standard chloroplast genome categories)
# ---------------------------------------------------------------------------
CATEGORY_ORDER = [
    'Photosystem I',
    'Photosystem II',
    'Cytochrome b/f complex',
    'ATP synthase',
    'NADH dehydrogenase',
    'Rubisco (rbcL)',
    'Ribosomal proteins (SSU)',
    'Ribosomal proteins (LSU)',
    'RNA polymerase',
    'Maturase (matK)',
    'Protease (clpP)',
    'Envelope membrane (cemA)',
    'C-type cytochrome synthesis (ccsA)',
    'Translation initiation factor (infA)',
    'Acetyl-CoA carboxylase (accD)',
    'Conserved ORFs (ycf)',
    'Transfer RNA',
    'Ribosomal RNA',
    'Other / Unclassified',
]

DEFAULT_CATEGORY_COLORS = {
    'Photosystem I': '#1f77b4',
    'Photosystem II': '#d62728',
    'Cytochrome b/f complex': '#2ca02c',
    'ATP synthase': '#ff7f0e',
    'NADH dehydrogenase': '#9467bd',
    'Rubisco (rbcL)': '#8c564b',
    'Ribosomal proteins (SSU)': '#17becf',
    'Ribosomal proteins (LSU)': '#e377c2',
    'RNA polymerase': '#bcbd22',
    'Maturase (matK)': '#7f7f7f',
    'Protease (clpP)': '#aec7e8',
    'Envelope membrane (cemA)': '#ffbb78',
    'C-type cytochrome synthesis (ccsA)': '#98df8a',
    'Translation initiation factor (infA)': '#ff9896',
    'Acetyl-CoA carboxylase (accD)': '#c49c94',
    'Conserved ORFs (ycf)': '#c5b0d5',
    'Transfer RNA': '#f7b6d2',
    'Ribosomal RNA': '#c7c7c7',
    'Other / Unclassified': '#404040',
}

SPECIFIC_GENE_FUNCTIONS = {
    'ycf3': 'Photosystem I',
    'ycf4': 'Photosystem I',
    'ycf1': 'Conserved ORFs (ycf)',
    'ycf2': 'Conserved ORFs (ycf)',
    'rbcl': 'Rubisco (rbcL)',
    'matk': 'Maturase (matK)',
    'clpp': 'Protease (clpP)',
    'cema': 'Envelope membrane (cemA)',
    'ccsa': 'C-type cytochrome synthesis (ccsA)',
    'infa': 'Translation initiation factor (infA)',
    'accd': 'Acetyl-CoA carboxylase (accD)',
}

FAMILY_GENE_FUNCTIONS = [
    ('psa', 'Photosystem I'),
    ('psb', 'Photosystem II'),
    ('pet', 'Cytochrome b/f complex'),
    ('atp', 'ATP synthase'),
    ('ndh', 'NADH dehydrogenase'),
    ('rps', 'Ribosomal proteins (SSU)'),
    ('rpl', 'Ribosomal proteins (LSU)'),
    ('rpo', 'RNA polymerase'),
    ('rbc', 'Rubisco (rbcL)'),
    ('ycf', 'Conserved ORFs (ycf)'),
    ('rrn', 'Ribosomal RNA'),
    ('trn', 'Transfer RNA'),
    ('orf', 'Conserved ORFs (ycf)'),
]

# First character that terminates the gene symbol inside a raw identifier
GENE_SYMBOL_TERMINATORS = '_( [{:|,;.'

# Nucleotide-composition parameters that can be correlated with the axes.
# The first three are the ones requested by the reviewer.
CORR_PARAM_ORDER = ['GC3s', 'GC', 'GC12', 'CAI', 'ENC', 'Fop', 'L_aa']
CORR_PARAM_DEFAULT = {'GC3s': True, 'GC': True, 'GC12': True,
                      'CAI': False, 'ENC': False, 'Fop': False, 'L_aa': False}

# CodonW COA works on the 59 sense codons (AUG, UGG and stops excluded).
# This is the exact column order used by CodonW in `coa_raw` / `codon.coa`.
CODONW_CODONS = [
    'UUU', 'UCU', 'UAU', 'UGU', 'UUC', 'UCC', 'UAC', 'UGC',
    'UUA', 'UCA', 'UUG', 'UCG',
    'CUU', 'CCU', 'CAU', 'CGU', 'CUC', 'CCC', 'CAC', 'CGC',
    'CUA', 'CCA', 'CAA', 'CGA', 'CUG', 'CCG', 'CAG', 'CGG',
    'AUU', 'ACU', 'AAU', 'AGU', 'AUC', 'ACC', 'AAC', 'AGC',
    'AUA', 'ACA', 'AAA', 'AGA', 'ACG', 'AAG', 'AGG',
    'GUU', 'GCU', 'GAU', 'GGU', 'GUC', 'GCC', 'GAC', 'GGC',
    'GUA', 'GCA', 'GAA', 'GGA', 'GUG', 'GCG', 'GAG', 'GGG',
]

# Columns of a CodonW bulk `.out` file that we know how to interpret
OUT_COLUMN_ALIASES = {
    'gc3s': 'GC3s', 'gc': 'GC', 'gc12': 'GC12', 'gc3': 'GC3s',
    'cai': 'CAI', 'enc': 'ENC', 'nc': 'ENC', 'fop': 'Fop',
    'l_aa': 'L_aa', 'laa': 'L_aa', 'laa': 'L_aa',
    'kind': None, 'title': None, 'species': None, 'name': None,
    't3s': None, 'c3s': None, 'a3s': None, 'g3s': None,
    'cbi': None, 'gravy': None, 'aromo': None, 'seq_hash': None,
}

LABEL_METHOD_ORIGIN = 'Farthest from origin (0, 0)'
LABEL_METHOD_CENTER = 'Farthest from data center'
LABEL_METHODS = (LABEL_METHOD_ORIGIN, LABEL_METHOD_CENTER)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _is_float(token):
    try:
        float(token)
        return True
    except (TypeError, ValueError):
        return False


def extract_gene_symbol(gene_name):
    """Return the bare gene symbol of a raw identifier.

    ``rpl20_[69725:70078](-)|01``  ->  ``rpl20``
    ``trnL-UAA_[107000:107063]``   ->  ``trnL-UAA``
    """
    raw = str(gene_name or '').strip()
    name = raw
    if '>' in name:
        name = name.rsplit('>', 1)[-1].strip()

    symbol = []
    for ch in name:
        if ch in GENE_SYMBOL_TERMINATORS or ch.isspace():
            break
        symbol.append(ch)

    symbol = ''.join(symbol).strip().strip('.')
    return symbol or raw


def coa_gene_symbol(label):
    """Gene symbol of a `genes.coa` label such as ``11_psbM_[29928:30032](...``.

    CodonW prefixes every label with the running gene number, which is
    stripped before the symbol is extracted.
    """
    s = str(label or '').strip()
    m = re.match(r'^\d+[_\-](.*)$', s)
    if m:
        s = m.group(1)
    return extract_gene_symbol(s)


def classify_gene_symbol(symbol):
    """Functional category of an already-extracted gene symbol."""
    gene_id = str(symbol).split('_')[0].split('|')[0].strip().lower()
    if not gene_id:
        return 'Other / Unclassified'
    if gene_id in SPECIFIC_GENE_FUNCTIONS:
        return SPECIFIC_GENE_FUNCTIONS[gene_id]
    for prefix, category in FAMILY_GENE_FUNCTIONS:
        if gene_id.startswith(prefix):
            return category
    return 'Other / Unclassified'


def gene_key(label):
    """Merge key of a gene label: the leading running number if present."""
    m = re.match(r'^\s*[>#]*(\d+)', str(label or ''))
    if m:
        return ('idx', int(m.group(1)))
    return ('name', str(label).strip().lower())


# ---------------------------------------------------------------------------
# CodonW file parsers
# ---------------------------------------------------------------------------
def _read_lines(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.readlines()


def parse_genes_coa(filepath):
    """Parse `genes.coa` (gene coordinates on the first COA axes).

    Returns ``{'labels': [...], 'coords': [np.array, ...], 'n_axes': int}``
    or ``None`` when no data could be extracted.
    """
    lines = _read_lines(filepath)

    n_axes = 0
    labels, rows = [], []

    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.lower().startswith('label'):
            toks = line.split()
            if len(toks) > 1:
                n_axes = len(toks) - 1
            continue

        tokens = line.split()

        # Number of trailing numeric tokens (the axis coordinates)
        k = 0
        while k < len(tokens) and _is_float(tokens[-1 - k]):
            k += 1
        if k < 2:
            continue
        if n_axes and k >= n_axes:
            k = n_axes
        if len(tokens) - k < 1:
            continue

        labels.append(' '.join(tokens[:-k]))
        rows.append([float(v) for v in tokens[-k:]])

    if len(labels) < 3:
        return None

    coords = [np.array([row[i] for row in rows], dtype=float)
              for i in range(len(rows[0]))]
    return {'labels': labels, 'coords': coords,
            'n_axes': len(rows[0]), 'n_genes': len(labels)}


def parse_eigen_coa(filepath):
    """Parse `eigen.coa` (axis number, eigenvalue, fraction, cumulative)."""
    rows = []
    for raw in _read_lines(filepath):
        parts = raw.split()
        if len(parts) < 4:
            continue
        try:
            n = int(parts[0])
            eig = float(parts[1])
            frac = float(parts[2])
            cum = float(parts[3])
        except ValueError:
            continue
        rows.append((n, eig, frac, cum))

    if not rows:
        return None
    return {'axis': [r[0] for r in rows],
            'eigen': [r[1] for r in rows],
            'frac': [r[2] for r in rows],
            'cum': [r[3] for r in rows]}


def compute_params_from_counts(counts):
    """Nucleotide-composition parameters from 59 codon counts (`coa_raw`).

    Follows the CodonW conventions (fractions between 0 and 1).
    """
    counts = np.asarray(counts, dtype=float)

    t3 = c3 = a3 = g3 = 0.0
    gc12_count = 0.0
    total_bases = 0.0
    for cnt, cod in zip(counts, CODONW_CODONS):
        b3 = cod[2]
        if b3 == 'T':
            t3 += cnt
        elif b3 == 'C':
            c3 += cnt
        elif b3 == 'A':
            a3 += cnt
        else:
            g3 += cnt
        if cod[0] in 'GC':
            gc12_count += cnt
        if cod[1] in 'GC':
            gc12_count += cnt
        total_bases += cnt * 3.0

    n_codons = counts.sum()
    n3 = t3 + c3 + a3 + g3

    params = {}
    params['GC3s'] = (g3 + c3) / n3 if n3 > 0 else float('nan')
    params['GC'] = (gc12_count + g3 + c3) / total_bases if total_bases > 0 else float('nan')
    params['GC12'] = gc12_count / (2.0 * n_codons) if n_codons > 0 else float('nan')
    params['L_aa'] = n_codons
    return params


def parse_coa_raw(filepath):
    """Parse `coa_raw` (raw codon usage, 59 counts per gene)."""
    labels, all_counts = [], []
    for raw in _read_lines(filepath):
        line = raw.strip()
        if not line:
            continue
        tokens = line.split()
        k = 0
        while k < len(tokens) and _is_float(tokens[-1 - k]):
            k += 1
        if k < 20:            # need (almost) all 59 codon columns
            continue
        if len(tokens) - k < 1:
            continue
        labels.append(' '.join(tokens[:-k]))
        all_counts.append([float(v) for v in tokens[-k:]])

    if not labels:
        return None

    # Keep exactly the 59 CodonW codon columns when more numbers are present
    param_rows = []
    for counts in all_counts:
        if len(counts) >= len(CODONW_CODONS):
            counts = counts[-len(CODONW_CODONS):]
        else:
            counts = (counts + [0.0] * len(CODONW_CODONS))[:len(CODONW_CODONS)]
        param_rows.append(compute_params_from_counts(counts))

    params = {name: np.array([row[name] for row in param_rows], dtype=float)
              for name in param_rows[0]}
    return {'labels': labels, 'params': params}


def parse_codonw_out(filepath):
    """Parse a CodonW bulk output file (`.out`).

    The first (commented) line holds the column names.  Only the columns
    listed in ``OUT_COLUMN_ALIASES`` are kept.  Returns ``{'labels': ...,
    'params': {...}}`` or ``None`` if the file cannot be interpreted.
    """
    lines = _read_lines(filepath)

    header = None
    records = []
    for raw in lines:
        s = raw.strip()
        if not s:
            continue
        if raw.lstrip().startswith('#'):
            if header is None:
                header = s.lstrip('#').replace(',', ' ').split()
            continue
        tokens = s.replace(',', ' ').split()
        if header is None:
            # A data line before any header -> this is not a bulk .out file
            if _is_float(tokens[-1]):
                return None
            header = tokens
            continue
        ncols = len(header) - 1
        if ncols <= 0 or len(tokens) <= ncols:
            continue
        name = ' '.join(tokens[:-ncols])
        try:
            values = [float(v) if _is_float(v) else float('nan')
                      for v in tokens[-ncols:]]
        except ValueError:
            continue
        records.append((name, values))

    if header is None or not records:
        return None

    wanted = []
    for i, col in enumerate(header[1:], start=1):
        key = OUT_COLUMN_ALIASES.get(col.strip().lower())
        if key:
            wanted.append((i, key))
    if not wanted:
        return None

    params = {}
    for idx, key in wanted:
        if key in params:
            continue
        params[key] = np.array([rec[1][idx] for rec in records], dtype=float)

    return {'labels': [rec[0] for rec in records], 'params': params}


def load_composition(filepath):
    """Load nucleotide-composition parameters from `coa_raw` or a `.out` file."""
    base = os.path.basename(filepath).lower()
    if base == 'coa_raw':
        return parse_coa_raw(filepath)
    if base.endswith('.out'):
        result = parse_codonw_out(filepath)
        if result and result['params']:
            return result
        return parse_coa_raw(filepath)
    # unknown extension: try both
    result = parse_codonw_out(filepath)
    if result and result['params']:
        return result
    return parse_coa_raw(filepath)


def merge_composition(gene_labels, comp):
    """Row index in ``comp`` matching each label of ``genes.coa`` (or None)."""
    by_key = {}
    prefixes = []
    for i, lab in enumerate(comp['labels']):
        by_key[gene_key(lab)] = i
        prefixes.append(str(lab)[:12])

    mapping = []
    for lab in gene_labels:
        k = gene_key(lab)
        if k in by_key:
            mapping.append(by_key[k])
            continue
        pref = str(lab)[:12]
        hit = None
        for i, p in enumerate(prefixes):
            if p == pref:
                hit = i
                break
        mapping.append(hit)
    return mapping


def build_dataset(genes_path, eigen_path=None, comp_path=None, name=None):
    """Assemble one dataset from the CodonW output files of one species."""
    genes = parse_genes_coa(genes_path)
    if genes is None:
        return None

    eigen = None
    if eigen_path and os.path.isfile(eigen_path):
        eigen = parse_eigen_coa(eigen_path)

    params = None
    if comp_path and os.path.isfile(comp_path):
        comp = load_composition(comp_path)
        if comp and comp.get('params'):
            mapping = merge_composition(genes['labels'], comp)
            n = len(genes['labels'])
            params = {}
            for pname, arr in comp['params'].items():
                params[pname] = np.array(
                    [arr[m] if m is not None else float('nan') for m in mapping],
                    dtype=float)
            matched = sum(1 for m in mapping if m is not None)
            print(f"Composition match: {matched}/{n} genes "
                  f"(parameters: {', '.join(sorted(params))})")

    if len(genes['coords']) < 2:
        return None

    symbols = [coa_gene_symbol(lab) for lab in genes['labels']]
    return {
        'name': name or os.path.basename(genes_path),
        'filename': os.path.basename(genes_path),
        'x': genes['coords'][0],           # Axis 1
        'y': genes['coords'][1],           # Axis 2
        'coords': genes['coords'],
        'n_axes': genes['n_axes'],
        'genes': genes['labels'],
        'symbols': symbols,
        'categories': [classify_gene_symbol(s) for s in symbols],
        'eigen': eigen,
        'params': params,
    }


# ---------------------------------------------------------------------------
# Statistics (Pearson & Spearman with P values, no SciPy required)
# ---------------------------------------------------------------------------
def _betacf(a, b, x):
    """Continued fraction for the incomplete beta function (Numerical Recipes)."""
    MAXIT, EPS, FPMIN = 200, 3e-12, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS:
            break
    return h


def _betai(a, b, x):
    """Regularised incomplete beta function I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_bt = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log(1.0 - x))
    bt = math.exp(ln_bt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_two_tailed_p(t, df):
    """Two-tailed P value of Student's t with df degrees of freedom."""
    if df <= 0:
        return float('nan')
    return _betai(0.5 * df, 0.5, df / (df + t * t))


def _average_ranks(v):
    order = np.argsort(v, kind='mergesort')
    ranks = np.empty(len(v), dtype=float)
    sv = v[order]
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def correlations(x, y):
    """Pearson and Spearman correlation between two arrays (NaN-safe)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 3:
        return None

    xm, ym = x.mean(), y.mean()
    ssx = float(((x - xm) ** 2).sum())
    ssy = float(((y - ym) ** 2).sum())
    if ssx <= 0 or ssy <= 0:
        return None

    r = float(((x - xm) * (y - ym)).sum() / math.sqrt(ssx * ssy))
    r = max(-1.0, min(1.0, r))
    if abs(r) >= 1.0:
        p = 0.0
    else:
        t = r * math.sqrt((n - 2) / (1.0 - r * r))
        p = t_two_tailed_p(t, n - 2)

    rx, ry = _average_ranks(x), _average_ranks(y)
    rxm, rym = rx.mean(), ry.mean()
    srx = float(((rx - rxm) ** 2).sum())
    sry = float(((ry - rym) ** 2).sum())
    if srx <= 0 or sry <= 0:
        return {'n': n, 'r': r, 'p': p, 'rho': float('nan'),
                'p_rho': float('nan')}
    rho = float(((rx - rxm) * (ry - rym)).sum() / math.sqrt(srx * sry))
    rho = max(-1.0, min(1.0, rho))
    if abs(rho) >= 1.0:
        p_rho = 0.0
    else:
        t = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
        p_rho = t_two_tailed_p(t, n - 2)

    return {'n': n, 'r': r, 'p': p, 'rho': rho, 'p_rho': p_rho}


def fmt_p(p):
    if p != p:  # NaN
        return 'P = NA'
    if p < 0.001:
        return 'P < 0.001'
    return f'P = {p:.3f}'


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------
class COAPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("COA Plot Plotter (CodonW Correspondence Analysis)")
        self.root.geometry("1150x950")

        self.style_settings = {
            'title': {'fontsize': 13, 'fontweight': 'bold', 'color': 'black', 'visible': True},
            'xlabel': {'fontsize': 11, 'fontweight': 'bold', 'color': 'black', 'visible': True},
            'ylabel': {'fontsize': 11, 'fontweight': 'bold', 'color': 'black', 'visible': True},
            'xticklabels': {'fontsize': 9, 'color': 'black'},
            'yticklabels': {'fontsize': 9, 'color': 'black'},
            'refline': {'color': 'gray', 'linestyle': '--', 'linewidth': 1.0, 'visible': True},
            'scatter': {'color': '#1f77b4', 'edgecolor': 'darkblue', 'size': 40,
                        'alpha': 0.7, 'linewidth': 0.5, 'visible': True},
            'legend': {'fontsize': 9, 'color': 'black', 'visible': True},
        }

        # Gene labels for the outlying points
        self.label_outliers = tk.BooleanVar(value=True)
        self.label_method = tk.StringVar(value=LABEL_METHOD_ORIGIN)
        self.label_count = tk.IntVar(value=5)
        self.label_fontsize = tk.IntVar(value=8)
        self.label_color = 'darkgreen'

        # Options
        self.color_by_function = tk.BooleanVar(value=False)
        self.show_regression = tk.BooleanVar(value=True)
        self.show_corr_value = tk.BooleanVar(value=True)
        self.input_mode = tk.StringVar(value='file')
        self.folder_path = tk.StringVar(value='')
        self.grid_columns = tk.StringVar(value='0')
        self.scree_count = tk.StringVar(value='10')

        # Correlation parameter selection (populated after loading)
        self.corr_vars = {}

        # Datasets: each {'name', 'filename', 'x', 'y', 'coords', 'genes',
        #                  'symbols', 'categories', 'eigen', 'params'}
        self.datasets = []

        self.category_colors = dict(DEFAULT_CATEGORY_COLORS)
        self.category_visibility = {cat: True for cat in CATEGORY_ORDER}

        self.create_widgets()

    # ------------------------------------------------------------------ UI
    def create_widgets(self):
        control_frame = ttk.LabelFrame(self.root, text="Input Selection (CodonW COA output)", padding=10)
        control_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(control_frame, text="Input mode:").grid(row=0, column=0, sticky=tk.W, pady=5)
        mode_frame = ttk.Frame(control_frame)
        mode_frame.grid(row=0, column=1, columnspan=3, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(mode_frame, text="Single species (full supplementary figure)",
                        variable=self.input_mode, value='file',
                        command=self.on_mode_change).pack(side=tk.LEFT)
        ttk.Radiobutton(mode_frame, text="Folder (multiple species, COA panel grid)",
                        variable=self.input_mode, value='folder',
                        command=self.on_mode_change).pack(side=tk.LEFT, padx=(15, 0))

        # genes.coa
        ttk.Label(control_frame, text="genes.coa (required):").grid(row=1, column=0, sticky=tk.W, pady=4)
        self.genes_path = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.genes_path, width=52).grid(
            row=1, column=1, padx=5, pady=4, sticky=tk.EW)
        ttk.Button(control_frame, text="Browse...", command=self.browse_genes).grid(
            row=1, column=2, padx=5, pady=4)

        # eigen.coa
        ttk.Label(control_frame, text="eigen.coa (optional):").grid(row=2, column=0, sticky=tk.W, pady=4)
        self.eigen_path = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.eigen_path, width=52).grid(
            row=2, column=1, padx=5, pady=4, sticky=tk.EW)
        ttk.Button(control_frame, text="Browse...", command=self.browse_eigen).grid(
            row=2, column=2, padx=5, pady=4)

        # composition source
        ttk.Label(control_frame, text="Composition file\n(coa_raw or .out, optional):").grid(
            row=3, column=0, sticky=tk.W, pady=4)
        self.comp_path = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.comp_path, width=52).grid(
            row=3, column=1, padx=5, pady=4, sticky=tk.EW)
        ttk.Button(control_frame, text="Browse...", command=self.browse_comp).grid(
            row=3, column=2, padx=5, pady=4)

        # folder mode
        ttk.Label(control_frame, text="Input folder\n(each subfolder = one species):").grid(
            row=4, column=0, sticky=tk.W, pady=4)
        self.folder_entry = ttk.Entry(control_frame, textvariable=self.folder_path, width=52)
        self.folder_entry.grid(row=4, column=1, padx=5, pady=4, sticky=tk.EW)
        ttk.Button(control_frame, text="Browse Folder...", command=self.browse_folder).grid(
            row=4, column=2, padx=5, pady=4)

        # species name / actions
        ttk.Label(control_frame, text="Species name:").grid(row=5, column=0, sticky=tk.W, pady=6)
        self.species_name = tk.StringVar(value="Species")
        ttk.Entry(control_frame, textvariable=self.species_name, width=28).grid(
            row=5, column=1, padx=5, pady=6, sticky=tk.W)

        action_frame = ttk.Frame(control_frame)
        action_frame.grid(row=5, column=2, columnspan=2, sticky=tk.E, padx=5)
        ttk.Button(action_frame, text="Plot Graph", command=self.plot_graph).pack(side=tk.LEFT, padx=2)
        ttk.Button(action_frame, text="Style Settings", command=self.open_style_settings).pack(side=tk.LEFT, padx=2)
        ttk.Button(action_frame, text="Save Image", command=self.save_image).pack(side=tk.LEFT, padx=2)
        ttk.Button(action_frame, text="Export Stats (CSV)", command=self.export_stats).pack(side=tk.LEFT, padx=2)

        control_frame.columnconfigure(1, weight=1)

        # ---- Options -----------------------------------------------------
        options_frame = ttk.LabelFrame(self.root, text="Plot Options", padding=8)
        options_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Checkbutton(options_frame, text="Color points by gene function",
                        variable=self.color_by_function,
                        command=self.redraw).pack(side=tk.LEFT)
        ttk.Checkbutton(options_frame, text="Show regression line",
                        variable=self.show_regression,
                        command=self.redraw).pack(side=tk.LEFT, padx=(15, 0))
        ttk.Checkbutton(options_frame, text="Show r / P on correlation panels",
                        variable=self.show_corr_value,
                        command=self.redraw).pack(side=tk.LEFT, padx=(15, 0))

        ttk.Label(options_frame, text="  Scree axes:").pack(side=tk.LEFT)
        ttk.Spinbox(options_frame, from_=2, to=58, width=4,
                    textvariable=self.scree_count,
                    command=self.redraw).pack(side=tk.LEFT, padx=(3, 0))

        ttk.Label(options_frame, text="  Subplot columns (folder mode, 0 = auto):").pack(side=tk.LEFT)
        ttk.Spinbox(options_frame, from_=0, to=12, width=4,
                    textvariable=self.grid_columns).pack(side=tk.LEFT, padx=(3, 0))

        # ---- Correlation parameters --------------------------------------
        self.corr_frame = ttk.LabelFrame(self.root, text="Composition Parameters Correlated with Axis 1 / Axis 2", padding=8)
        self.corr_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(self.corr_frame, text="(available parameters appear here after a file has been plotted)",
                  foreground='gray').pack(anchor=tk.W)

        # ---- Gene labels ---------------------------------------------------
        label_frame = ttk.LabelFrame(self.root, text="Gene Labels (outliers)", padding=5)
        label_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Checkbutton(label_frame, text="Label outlying genes",
                        variable=self.label_outliers,
                        command=self.redraw).pack(side=tk.LEFT)

        ttk.Label(label_frame, text="  Rank by:").pack(side=tk.LEFT)
        method_box = ttk.Combobox(label_frame, textvariable=self.label_method,
                                  values=list(LABEL_METHODS), width=30, state='readonly')
        method_box.pack(side=tk.LEFT, padx=(4, 0))
        method_box.bind('<<ComboboxSelected>>', lambda event: self.redraw())

        ttk.Label(label_frame, text="  Number:").pack(side=tk.LEFT)
        ttk.Spinbox(label_frame, from_=1, to=30, width=4,
                    textvariable=self.label_count,
                    command=self.redraw).pack(side=tk.LEFT, padx=(4, 0))

        ttk.Label(label_frame, text="  Font size:").pack(side=tk.LEFT)
        ttk.Spinbox(label_frame, from_=4, to=16, width=4,
                    textvariable=self.label_fontsize,
                    command=self.redraw).pack(side=tk.LEFT, padx=(4, 0))

        ttk.Label(label_frame, text="  Color:").pack(side=tk.LEFT)
        self.label_color_btn = ttk.Button(label_frame, text="", width=8,
                                          command=self.choose_label_color)
        self.label_color_btn.pack(side=tk.LEFT, padx=(4, 0))
        self.update_color_button(self.label_color_btn, self.label_color)

        # ---- Plot area (scrollable) ---------------------------------------
        plot_frame = ttk.LabelFrame(self.root, text="Graph Display", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.canvas_viewport = tk.Canvas(plot_frame, bg='white', highlightthickness=0)
        vbar = ttk.Scrollbar(plot_frame, orient=tk.VERTICAL,
                             command=self.canvas_viewport.yview)
        hbar = ttk.Scrollbar(plot_frame, orient=tk.HORIZONTAL,
                             command=self.canvas_viewport.xview)
        self.canvas_viewport.configure(yscrollcommand=vbar.set,
                                       xscrollcommand=hbar.set)
        self.canvas_viewport.grid(row=0, column=0, sticky='nsew')
        vbar.grid(row=0, column=1, sticky='ns')
        hbar.grid(row=1, column=0, sticky='ew')
        plot_frame.rowconfigure(0, weight=1)
        plot_frame.columnconfigure(0, weight=1)

        self.canvas_container = tk.Frame(self.canvas_viewport, width=600, height=600,
                                         bg='white')
        self.canvas_container.pack_propagate(False)
        self.canvas_viewport.create_window((0, 0), window=self.canvas_container,
                                           anchor='nw')

        self.fig, self.ax = plt.subplots(figsize=(7, 7))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_container)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(fill=tk.X, padx=5, pady=5)

        self.on_mode_change()

    # ------------------------------------------------------------ callbacks
    def on_mode_change(self):
        if self.input_mode.get() == 'folder':
            self.status_bar.config(
                text="Folder mode: every subfolder containing genes.coa becomes one COA panel "
                     "(eigen.coa / coa_raw / *.out in the same subfolder are used automatically)")
        else:
            self.status_bar.config(
                text="Single-species mode: the supplementary figure contains the COA plot, "
                     "scree plot and correlation panels")

    def browse_genes(self):
        filename = filedialog.askopenfilename(
            title="Select genes.coa",
            filetypes=[("COA gene coordinates", "genes.coa"),
                       ("All files", "*.*")])
        if filename:
            self.input_mode.set('file')
            self.genes_path.set(filename)
            self.autodetect_siblings(filename)
            self.status_bar.config(text=f"Selected: {os.path.basename(filename)}")

    def browse_eigen(self):
        filename = filedialog.askopenfilename(
            title="Select eigen.coa",
            filetypes=[("COA eigenvalues", "eigen.coa"), ("All files", "*.*")])
        if filename:
            self.eigen_path.set(filename)

    def browse_comp(self):
        filename = filedialog.askopenfilename(
            title="Select coa_raw or CodonW .out file",
            filetypes=[("CodonW output", "*.out"), ("coa_raw", "coa_raw"),
                       ("All files", "*.*")])
        if filename:
            self.comp_path.set(filename)

    def browse_folder(self):
        folder = filedialog.askdirectory(title="Select folder containing species subfolders")
        if folder:
            self.input_mode.set('folder')
            self.folder_path.set(folder)
            n = len(self.find_species_subfolders(folder))
            self.status_bar.config(text=f"Selected folder: {folder}  ({n} species folder(s) found)")

    def autodetect_siblings(self, genes_file):
        """Fill the eigen / composition entries with files next to genes.coa."""
        d = os.path.dirname(genes_file) or '.'
        if not self.eigen_path.get():
            cand = os.path.join(d, 'eigen.coa')
            if os.path.isfile(cand):
                self.eigen_path.set(cand)
        if not self.comp_path.get():
            cand = os.path.join(d, 'coa_raw')
            if os.path.isfile(cand):
                self.comp_path.set(cand)
            else:
                for f in sorted(glob.glob(os.path.join(d, '*.out'))):
                    if os.path.getsize(f) > 0:
                        self.comp_path.set(f)
                        break

    @staticmethod
    def find_species_subfolders(root):
        """Subfolders of `root` that contain a `genes.coa` (or `root` itself)."""
        found = []
        for d in sorted(glob.glob(os.path.join(root, '*'))):
            if os.path.isdir(d) and glob.glob(os.path.join(d, 'genes.coa')):
                found.append(d)
        if not found and glob.glob(os.path.join(root, 'genes.coa')):
            found = [root]
        return found

    # ------------------------------------------------------------- loading
    def load_data(self):
        genes_file = self.genes_path.get()
        if not genes_file or not os.path.exists(genes_file):
            messagebox.showerror("Error", "Please select a valid genes.coa file!")
            return False
        try:
            dataset = build_dataset(genes_file,
                                    self.eigen_path.get(),
                                    self.comp_path.get(),
                                    self.species_name.get().strip() or os.path.basename(genes_file))
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file: {e}")
            return False
        if dataset is None:
            messagebox.showerror("Error", "Could not extract valid data from genes.coa!")
            return False

        self.datasets = [dataset]
        self.rebuild_param_checkboxes()
        self.report_dataset(dataset)
        return True

    def load_folder_data(self):
        folder = self.folder_path.get()
        if not folder or not os.path.isdir(folder):
            messagebox.showerror("Error", "Please select a valid input folder!")
            return False

        subfolders = self.find_species_subfolders(folder)
        if not subfolders:
            messagebox.showerror("Error",
                                 "No subfolder containing genes.coa was found!")
            return False

        datasets, failed = [], []
        for d in subfolders:
            g = os.path.join(d, 'genes.coa')
            e = os.path.join(d, 'eigen.coa')
            c = os.path.join(d, 'coa_raw')
            if not os.path.isfile(c):
                c = None
                for f in sorted(glob.glob(os.path.join(d, '*.out'))):
                    if os.path.getsize(f) > 0:
                        c = f
                        break
            try:
                dataset = build_dataset(g, e if os.path.isfile(e) else None, c,
                                        name=os.path.basename(d))
            except Exception as exc:
                print(f"Failed to read {g}: {exc}")
                dataset = None
            if dataset is None:
                failed.append(os.path.basename(d))
            else:
                datasets.append(dataset)

        if not datasets:
            messagebox.showerror("Error", "No valid data could be extracted!")
            return False

        self.datasets = datasets
        if failed:
            print(f"Skipped folders without valid data: {', '.join(failed)}")
        self.rebuild_param_checkboxes()
        total = sum(len(d['x']) for d in datasets)
        print(f"Loaded {len(datasets)} species, {total} genes in total")
        self.status_bar.config(
            text=f"Loaded {len(datasets)} species, {total} genes in total")
        return True

    def report_dataset(self, dataset):
        n = len(dataset['x'])
        eig = dataset['eigen']
        if eig:
            print(f"Axis 1 explains {eig['frac'][0]*100:.2f}% of the total inertia, "
                  f"axis 2 {eig['frac'][1]*100:.2f}% (cumulative "
                  f"{eig['cum'][1]*100:.2f}%)")
        self.status_bar.config(
            text=f"Loaded {n} genes"
                 + (f" - Axis1 {eig['frac'][0]*100:.2f}%, Axis2 {eig['frac'][1]*100:.2f}%"
                    if eig else ""))

    def rebuild_param_checkboxes(self):
        """(Re)create the correlation-parameter checkboxes for the loaded data."""
        for w in self.corr_frame.winfo_children():
            w.destroy()

        available = []
        for d in self.datasets:
            if d.get('params'):
                for p in d['params']:
                    if p not in available:
                        available.append(p)
        ordered = [p for p in CORR_PARAM_ORDER if p in available]
        ordered += sorted(p for p in available if p not in CORR_PARAM_ORDER)

        if not ordered:
            ttk.Label(self.corr_frame,
                      text="No composition data loaded - correlation panels are skipped "
                           "(select a coa_raw or .out file)",
                      foreground='gray').pack(anchor=tk.W)
            return

        ttk.Label(self.corr_frame, text="Correlate axes with:").pack(side=tk.LEFT)
        for p in ordered:
            var = self.corr_vars.get(p)
            if var is None:
                var = tk.BooleanVar(value=CORR_PARAM_DEFAULT.get(p, False))
                self.corr_vars[p] = var
            ttk.Checkbutton(self.corr_frame, text=p, variable=var,
                            command=self.redraw).pack(side=tk.LEFT, padx=(8, 0))

    def selected_params(self):
        return [p for p in CORR_PARAM_ORDER
                if self.corr_vars.get(p) is not None and self.corr_vars[p].get()]

    # -------------------------------------------------------------- drawing
    def fit_canvas_to_figure(self):
        try:
            w_in, h_in = self.fig.get_size_inches()
            if w_in <= 0 or h_in <= 0:
                return
            display_dpi = max(int(self.fig.dpi), 100)
            width = max(600, int(round(w_in * display_dpi)))
            height = max(600, int(round(h_in * display_dpi)))
            self.canvas_container.config(width=width, height=height)
            self.canvas_viewport.update_idletasks()
            self.canvas_viewport.config(scrollregion=(0, 0, width, height))
        except Exception:
            pass

    def get_active_categories(self, categories):
        present = set(categories)
        ordered = [cat for cat in CATEGORY_ORDER if cat in present]
        ordered += sorted(cat for cat in present if cat not in CATEGORY_ORDER)
        return ordered

    def get_category_color(self, category):
        if category not in self.category_colors:
            self.category_colors[category] = DEFAULT_CATEGORY_COLORS.get(category, '#404040')
        return self.category_colors[category]

    def axis_label(self, dataset, i):
        eig = dataset.get('eigen')
        if eig and i - 1 < len(eig['frac']):
            return f'Axis {i} ({eig["frac"][i - 1] * 100:.2f}%)'
        return f'Axis {i}'

    # -- outliers -----------------------------------------------------------
    def _label_count(self):
        try:
            count = int(float(self.label_count.get()))
        except (TypeError, ValueError):
            count = 5
        return max(0, count)

    def find_outlier_genes(self, dataset, xlim=None, ylim=None):
        """Genes farthest from the origin (0,0) or from the data centroid."""
        count = self._label_count()
        if count <= 0 or len(dataset['x']) == 0:
            return []

        xs = np.asarray(dataset['x'], dtype=float)
        ys = np.asarray(dataset['y'], dtype=float)

        if xlim is None:
            xlim = (float(np.nanmin(xs)), float(np.nanmax(xs)))
        if ylim is None:
            ylim = (float(np.nanmin(ys)), float(np.nanmax(ys)))
        span = max((xlim[1] - xlim[0]), (ylim[1] - ylim[0]), 1e-9)
        span_x = max(xlim[1] - xlim[0], 1e-9)
        span_y = max(ylim[1] - ylim[0], 1e-9)

        nx = xs / span_x
        ny = ys / span_y
        if self.label_method.get() == LABEL_METHOD_CENTER:
            nx = nx - np.nanmean(nx)
            ny = ny - np.nanmean(ny)
        scores = np.hypot(nx, ny)
        scores = np.where(np.isfinite(scores), scores, -np.inf)

        order = np.argsort(-scores)[:count]
        return [(int(i), float(scores[i]), dataset['symbols'][i],
                 float(xs[i]), float(ys[i])) for i in order]

    def annotate_outliers(self, ax, dataset, xlim=None, ylim=None):
        if not self.label_outliers.get():
            return
        genes = self.find_outlier_genes(dataset, xlim, ylim)
        if not genes:
            return
        try:
            fontsize = max(4, int(float(self.label_fontsize.get())))
        except (TypeError, ValueError):
            fontsize = 8
        color = self.label_color

        px_per_point = ax.figure.dpi / 72.0
        offsets = ((7, 7), (-7, 7), (7, -7), (-7, -7),
                   (7, 20), (7, -20), (-7, 20), (-7, -20))
        placed = []
        for _, _, symbol, x, y in genes:
            px, py = ax.transData.transform((x, y))
            width = (len(symbol) * 0.62 + 0.7) * fontsize * px_per_point + 6
            height = fontsize * 1.35 * px_per_point + 4
            dx, dy = offsets[0]
            for dx, dy in offsets:
                x0 = px + dx if dx > 0 else px + dx - width
                y0 = py + dy if dy > 0 else py + dy - height
                box = (x0, y0, x0 + width, y0 + height)
                if not any(self._boxes_overlap(box, other) for other in placed):
                    break
            placed.append(box)
            ax.annotate(symbol, xy=(x, y),
                        xytext=(dx / px_per_point, dy / px_per_point),
                        textcoords='offset points',
                        ha='left' if dx > 0 else 'right',
                        va='bottom' if dy > 0 else 'top',
                        fontsize=fontsize, color=color, fontweight='bold',
                        zorder=6,
                        arrowprops=dict(arrowstyle='->', color=color,
                                        lw=0.5, shrinkA=0, shrinkB=3))

    @staticmethod
    def _boxes_overlap(a, b):
        return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])

    def choose_label_color(self):
        color = colorchooser.askcolor(title="Choose gene label color",
                                      color=self.label_color)[1]
        if color:
            self.label_color = color
            self.update_color_button(self.label_color_btn, color)
            self.redraw()

    # -- COA panel ----------------------------------------------------------
    def draw_scatter_by_function(self, ax, dataset):
        style = self.style_settings['scatter']
        groups = {}
        for x, y, cat in zip(dataset['x'], dataset['y'], dataset['categories']):
            if not self.category_visibility.get(cat, True):
                continue
            groups.setdefault(cat, {'x': [], 'y': []})
            groups[cat]['x'].append(x)
            groups[cat]['y'].append(y)
        for cat in self.get_active_categories(dataset['categories']):
            if cat not in groups:
                continue
            pts = groups[cat]
            ax.scatter(pts['x'], pts['y'],
                       c=self.get_category_color(cat),
                       alpha=style['alpha'], s=style['size'],
                       edgecolors=style['edgecolor'],
                       linewidths=style['linewidth'], zorder=2,
                       label=f'{cat} (n={len(pts["x"])})')

    def draw_coa_panel(self, ax, dataset, title_text, add_legend=True,
                       title_fontsize=None, equal_aspect=False):
        style = self.style_settings

        if style['refline']['visible']:
            ax.axvline(0.0, color=style['refline']['color'],
                       linestyle=style['refline']['linestyle'],
                       linewidth=style['refline']['linewidth'], zorder=1)
            ax.axhline(0.0, color=style['refline']['color'],
                       linestyle=style['refline']['linestyle'],
                       linewidth=style['refline']['linewidth'], zorder=1)

        if style['scatter']['visible']:
            if self.color_by_function.get():
                self.draw_scatter_by_function(ax, dataset)
            else:
                ax.scatter(dataset['x'], dataset['y'],
                           c=style['scatter']['color'],
                           alpha=style['scatter']['alpha'],
                           s=style['scatter']['size'],
                           edgecolors=style['scatter']['edgecolor'],
                           linewidths=style['scatter']['linewidth'],
                           zorder=2,
                           label=f'Genes (n={len(dataset["x"])})')

        if style['xlabel']['visible']:
            ax.set_xlabel(self.axis_label(dataset, 1),
                          fontsize=style['xlabel']['fontsize'],
                          fontweight=style['xlabel']['fontweight'],
                          color=style['xlabel']['color'])
        if style['ylabel']['visible']:
            ax.set_ylabel(self.axis_label(dataset, 2),
                          fontsize=style['ylabel']['fontsize'],
                          fontweight=style['ylabel']['fontweight'],
                          color=style['ylabel']['color'])

        ax.tick_params(axis='x', labelsize=style['xticklabels']['fontsize'],
                       labelcolor=style['xticklabels']['color'])
        ax.tick_params(axis='y', labelsize=style['yticklabels']['fontsize'],
                       labelcolor=style['yticklabels']['color'])

        if style['title']['visible'] and title_text:
            size = title_fontsize if title_fontsize is not None else style['title']['fontsize']
            ax.set_title(title_text, fontsize=size,
                         fontweight=style['title']['fontweight'],
                         color=style['title']['color'], pad=12)

        if equal_aspect:
            ax.set_aspect('equal', adjustable='datalim')
        ax.grid(True, alpha=0.3, linestyle=':')

        if add_legend and style['legend']['visible']:
            ax.legend(loc='best', fontsize=style['legend']['fontsize'],
                      labelcolor=style['legend']['color'])

        self.annotate_outliers(ax, dataset, ax.get_xlim(), ax.get_ylim())

    # -- scree plot ---------------------------------------------------------
    def draw_scree(self, ax, dataset):
        eig = dataset.get('eigen')
        if not eig:
            return False
        try:
            n_show = max(2, int(float(self.scree_count.get())))
        except (TypeError, ValueError):
            n_show = 10
        n_show = min(n_show, len(eig['frac']))

        xs = np.arange(1, n_show + 1)
        frac = np.array(eig['frac'][:n_show]) * 100.0
        cum = np.array(eig['cum'][:n_show]) * 100.0

        ax.bar(xs, frac, color='#4c72b0', alpha=0.85,
               label='Fraction of inertia')
        ax2 = ax.twinx()
        ax2.plot(xs, cum, color='#c44e52', marker='o', markersize=3,
                 linewidth=1.5, label='Cumulative')
        ax2.set_ylim(0, 105)
        ax2.set_ylabel('Cumulative (%)', fontsize=9, color='#c44e52')
        ax2.tick_params(axis='y', labelsize=8, labelcolor='#c44e52')

        ax.set_xlabel('Axis number', fontsize=9)
        ax.set_ylabel('Variation explained (%)', fontsize=9)
        ax.set_title('Variance explained by COA axes',
                     fontsize=10, fontweight='bold')
        ax.tick_params(labelsize=8)
        ax.grid(True, axis='y', alpha=0.3, linestyle=':')
        ax.set_xticks(xs[::max(1, len(xs) // 10)])

        handles = [Patch(facecolor='#4c72b0', label='Fraction of inertia'),
                   Line2D([0], [0], color='#c44e52', marker='o', label='Cumulative')]
        ax.legend(handles=handles, loc='center right', fontsize=8, frameon=False)
        return True

    # -- correlation panel --------------------------------------------------
    def draw_corr_panel(self, ax, dataset, param, axis_index):
        if not dataset.get('params') or param not in dataset['params']:
            ax.set_visible(False)
            return

        x = np.asarray(dataset['params'][param], dtype=float)
        y = np.asarray(dataset['coords'][axis_index], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if len(x) < 3:
            ax.set_visible(False)
            return

        style = self.style_settings['scatter']
        ax.scatter(x, y, c=style['color'], alpha=style['alpha'],
                   s=max(15, style['size'] * 0.6),
                   edgecolors=style['edgecolor'],
                   linewidths=style['linewidth'], zorder=2)

        if self.show_regression.get():
            coef = np.polyfit(x, y, 1)
            xs = np.linspace(float(x.min()), float(x.max()), 50)
            ax.plot(xs, np.polyval(coef, xs), color='#d62728',
                    linewidth=1.5, zorder=3)

        ax.set_xlabel(param, fontsize=9)
        ax.set_ylabel(self.axis_label(dataset, axis_index + 1), fontsize=9)
        ax.set_title(f'{self.axis_label(dataset, axis_index + 1)} vs {param}',
                     fontsize=10)
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3, linestyle=':')

        if self.show_corr_value.get():
            stats = correlations(dataset['params'][param],
                                 dataset['coords'][axis_index])
            if stats:
                text = (f"Pearson r = {stats['r']:.3f}\n{fmt_p(stats['p'])}\n"
                        f"Spearman rho = {stats['rho']:.3f}")
                ax.text(0.03, 0.97, text, transform=ax.transAxes,
                        fontsize=8, va='top', ha='left',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  edgecolor='gray', alpha=0.85))

    # -- figure assembly ------------------------------------------------------
    def plot_graph(self):
        if self.input_mode.get() == 'folder':
            self.plot_combined()
        else:
            self.plot_single()

    def plot_single(self):
        if not self.load_data():
            return
        self.render_single()

    def render_single(self):
        if not self.datasets:
            return
        dataset = self.datasets[0]

        params = [p for p in self.selected_params()
                  if dataset.get('params') and p in dataset['params']]
        has_scree = dataset.get('eigen') is not None

        self.fig.clear()

        if not params and not has_scree:
            # COA plot only
            self.fig.set_size_inches(7, 7)
            ax = self.fig.add_subplot(111)
            self.draw_coa_panel(ax, dataset, dataset['name'], add_legend=True,
                                 equal_aspect=True)
        else:
            ncols = max(3, len(params))
            self.fig.set_size_inches(2.9 * ncols + 1.5, 13.0)
            gs = self.fig.add_gridspec(4, ncols,
                                       height_ratios=[1.15, 1.15, 1.0, 1.0],
                                       left=0.06, right=0.97,
                                       top=0.95, bottom=0.06,
                                       wspace=0.32, hspace=0.42)

            ax_coa = self.fig.add_subplot(gs[0:2, 0:2])
            self.draw_coa_panel(ax_coa, dataset, dataset['name'],
                                add_legend=True, equal_aspect=True)

            if has_scree:
                ax_scree = self.fig.add_subplot(gs[0:2, 2:ncols])
                self.draw_scree(ax_scree, dataset)

            for r, axis_index in enumerate((0, 1)):        # Axis 1, Axis 2
                for c, param in enumerate(params[:ncols]):
                    ax = self.fig.add_subplot(gs[2 + r, c])
                    self.draw_corr_panel(ax, dataset, param, axis_index)

        self.fit_canvas_to_figure()
        self.canvas.draw()

        eig = dataset.get('eigen')
        extra = (f" - Axis1 {eig['frac'][0]*100:.2f}%, "
                 f"Axis2 {eig['frac'][1]*100:.2f}%") if eig else ""
        self.status_bar.config(
            text=f"Supplementary figure: {len(dataset['x'])} genes{extra}")

    def plot_combined(self):
        if not self.load_folder_data():
            return
        self.render_combined()

    def render_combined(self):
        import math as _math
        datasets = self.datasets
        if not datasets:
            return
        n = len(datasets)

        try:
            user_cols = int(float(self.grid_columns.get()))
        except (TypeError, ValueError):
            user_cols = 0
        ncols = min(user_cols, n) if user_cols and user_cols > 0 else min(4, n)
        ncols = max(1, ncols)
        nrows = _math.ceil(n / ncols)

        panel, side_width = 4.0, 3.0
        fig_w = ncols * panel + side_width
        fig_h = max(nrows * panel, 5.0)

        self.fig.clear()
        self.fig.set_size_inches(fig_w, fig_h)

        right = max(0.45, 1.0 - side_width / fig_w)
        gs = self.fig.add_gridspec(nrows, ncols, left=0.07, right=right,
                                   top=0.95, bottom=0.06,
                                   wspace=0.30, hspace=0.35)

        title_fontsize = max(7, self.style_settings['title']['fontsize'] - 3)
        for i, dataset in enumerate(datasets):
            r, c = divmod(i, ncols)
            ax = self.fig.add_subplot(gs[r, c])
            self.draw_coa_panel(ax, dataset, dataset['name'],
                                add_legend=False,
                                title_fontsize=title_fontsize)

        self.add_side_legend(right)
        self.fit_canvas_to_figure()
        self.canvas.draw()

        total = sum(len(d['x']) for d in datasets)
        self.status_bar.config(
            text=f"Combined COA figure: {n} species, {total} genes, "
                 f"{nrows}x{ncols} grid")

    def add_side_legend(self, right):
        if not self.style_settings['legend']['visible']:
            return
        handles, labels = [], []
        style = self.style_settings['scatter']

        if self.color_by_function.get():
            present = []
            for dataset in self.datasets:
                for cat in self.get_active_categories(dataset['categories']):
                    if self.category_visibility.get(cat, True) and cat not in present:
                        present.append(cat)
            present = self.get_active_categories(present)
            for cat in present:
                handles.append(Patch(facecolor=self.get_category_color(cat),
                                     edgecolor=style['edgecolor'],
                                     alpha=style['alpha']))
                labels.append(cat)
        else:
            handles.append(Patch(facecolor=style['color'],
                                 edgecolor=style['edgecolor'],
                                 alpha=style['alpha']))
            labels.append('Genes')

        if not handles:
            return
        x_side = min(0.985, right + 0.03)
        self.fig.legend(handles, labels, loc='center left',
                        bbox_to_anchor=(x_side, 0.5),
                        fontsize=self.style_settings['legend']['fontsize'],
                        labelcolor=self.style_settings['legend']['color'],
                        frameon=False, borderaxespad=0.0)

    def redraw(self):
        if not self.datasets:
            return
        if self.input_mode.get() == 'folder':
            self.render_combined()
        else:
            self.render_single()

    # -- statistics export ---------------------------------------------------
    def collect_stat_rows(self):
        rows = []

        def add_rows(name, x_params, coords, eigen=None):
            for axis_index in (0, 1):
                for param, x in x_params.items():
                    stats = correlations(x, coords[axis_index])
                    if not stats:
                        continue
                    rows.append({
                        'dataset': name,
                        'axis': axis_index + 1,
                        'parameter': param,
                        'n': stats['n'],
                        'pearson_r': f"{stats['r']:.4f}",
                        'pearson_P': f"{stats['p']:.4g}",
                        'spearman_rho': f"{stats['rho']:.4f}",
                        'spearman_P': f"{stats['p_rho']:.4g}",
                    })

        params = self.selected_params()
        if not params:
            return rows

        for dataset in self.datasets:
            if dataset.get('params'):
                x_params = {p: dataset['params'][p] for p in params
                            if p in dataset['params']}
                if x_params:
                    add_rows(dataset['name'], x_params, dataset['coords'])

        # pooled rows when several species are loaded
        if len(self.datasets) > 1:
            pooled_params, pooled_axes = {}, {0: [], 1: []}
            for p in params:
                xs = [d['params'][p] for d in self.datasets
                      if d.get('params') and p in d['params']]
                if xs:
                    pooled_params[p] = np.concatenate(xs)
            for ai in (0, 1):
                pooled_axes[ai] = np.concatenate(
                    [d['coords'][ai] for d in self.datasets
                     if len(d['coords']) > ai])
            if pooled_params:
                add_rows('All species (pooled)', pooled_params, pooled_axes)
        return rows

    def export_stats(self):
        if not self.datasets:
            messagebox.showwarning("Warning", "Please plot the graph first!")
            return
        rows = self.collect_stat_rows()
        if not rows:
            messagebox.showwarning(
                "Warning", "No correlation could be computed.\n"
                           "Load a composition file (coa_raw or .out) and "
                           "select at least one parameter.")
            return

        filename = filedialog.asksaveasfilename(
            title="Export Correlation Statistics",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not filename:
            return

        header = ['dataset', 'axis', 'parameter', 'n',
                  'pearson_r', 'pearson_P', 'spearman_rho', 'spearman_P']
        try:
            with open(filename, 'w', encoding='utf-8-sig') as f:
                f.write('# Correlations between gene positions on the first two '
                        'COA axes (CodonW) and nucleotide composition parameters\n')
                f.write(','.join(header) + '\n')
                for row in rows:
                    f.write(','.join(str(row[k]) for k in header) + '\n')
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")
            return

        self.status_bar.config(text=f"Statistics exported: {os.path.basename(filename)} "
                                    f"({len(rows)} rows)")
        messagebox.showinfo("Success",
                            f"{len(rows)} correlation rows exported to:\n{filename}")

    # -- saving ---------------------------------------------------------------
    def save_image(self):
        if not self.datasets:
            messagebox.showwarning("Warning", "Please plot the graph first!")
            return
        filename = filedialog.asksaveasfilename(
            title="Save Image",
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("SVG files", "*.svg"),
                       ("PNG files", "*.png"), ("EPS files", "*.eps"),
                       ("All files", "*.*")])
        if not filename:
            return
        try:
            self.fig.savefig(filename, dpi=300, bbox_inches='tight')
        except Exception as e:
            messagebox.showerror("Error", f"Save failed: {e}")
            return
        self.status_bar.config(text=f"Image saved: {os.path.basename(filename)}")
        messagebox.showinfo("Success", f"Image saved to:\n{filename}")

    # -- style settings ---------------------------------------------------------
    def _add_font_scale(self, parent, label, section, lo, hi):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text=label, width=22).pack(side=tk.LEFT)
        var = tk.DoubleVar(value=self.style_settings[section]['fontsize'])
        scale = ttk.Scale(row, from_=lo, to=hi, orient=tk.HORIZONTAL,
                          variable=var,
                          command=lambda _v, s=section, v=var:
                          self._set_fontsize(s, v))
        scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

    def _set_fontsize(self, section, var):
        try:
            self.style_settings[section]['fontsize'] = float(var.get())
        except (TypeError, ValueError):
            return

    def _add_color_row(self, parent, label, section, key='color'):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text=label, width=22).pack(side=tk.LEFT)
        btn = ttk.Button(row, text=self.style_settings[section][key], width=12,
                         command=lambda: self.choose_color(section, btn))
        btn.pack(side=tk.LEFT)

    def _add_visibility(self, parent, label, section):
        var = tk.BooleanVar(value=self.style_settings[section]['visible'])
        ttk.Checkbutton(parent, text=label, variable=var,
                        command=lambda: self._set_visible(section, var)).pack(
            anchor=tk.W, padx=6, pady=2)

    def _set_visible(self, section, var):
        self.style_settings[section]['visible'] = var.get()
        self.redraw()

    def choose_color(self, element, button):
        color = colorchooser.askcolor(title="Choose Color")[1]
        if color:
            self.style_settings[element]['color'] = color
            self.update_color_button(button, color)
            self.redraw()

    def update_color_button(self, button, color):
        button.configure(text=color)

    def open_style_settings(self):
        win = tk.Toplevel(self.root)
        win.title("Style Settings")
        win.geometry("420x640")

        notebook = ttk.Notebook(win)
        notebook.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        tab_labels = ttk.Frame(notebook)
        notebook.add(tab_labels, text="Title & Labels")
        self._add_visibility(tab_labels, "Show title", 'title')
        self._add_font_scale(tab_labels, "Title font size", 'title', 8, 24)
        self._add_color_row(tab_labels, "Title color", 'title')
        self._add_font_scale(tab_labels, "X label font size", 'xlabel', 8, 20)
        self._add_font_scale(tab_labels, "Y label font size", 'ylabel', 8, 20)
        self._add_font_scale(tab_labels, "X tick font size", 'xticklabels', 6, 18)
        self._add_font_scale(tab_labels, "Y tick font size", 'yticklabels', 6, 18)

        tab_ref = ttk.Frame(notebook)
        notebook.add(tab_ref, text="Reference Lines")
        self._add_visibility(tab_ref, "Show zero reference lines", 'refline')
        row = ttk.Frame(tab_ref)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text="Line style", width=22).pack(side=tk.LEFT)
        combo = ttk.Combobox(row, values=['-', '--', ':', '-.'], width=10,
                             state='readonly')
        combo.set(self.style_settings['refline']['linestyle'])
        combo.pack(side=tk.LEFT)
        combo.bind('<<ComboboxSelected>>',
                   lambda e: self._set_linestyle(combo.get()))
        self._add_color_row(tab_ref, "Line color", 'refline')

        tab_scatter = ttk.Frame(notebook)
        notebook.add(tab_scatter, text="Scatter Points")
        self._add_visibility(tab_scatter, "Show scatter points", 'scatter')
        self._add_color_row(tab_scatter, "Fill color", 'scatter')
        self._add_color_row(tab_scatter, "Edge color", 'scatter', 'edgecolor')
        row = ttk.Frame(tab_scatter)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text="Point size", width=22).pack(side=tk.LEFT)
        var = tk.DoubleVar(value=self.style_settings['scatter']['size'])
        ttk.Scale(row, from_=10, to=200, orient=tk.HORIZONTAL, variable=var,
                  command=lambda _v: self._set_scatter('size', var)).pack(
            side=tk.LEFT, fill=tk.X, expand=True)
        row = ttk.Frame(tab_scatter)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text="Transparency", width=22).pack(side=tk.LEFT)
        var = tk.DoubleVar(value=self.style_settings['scatter']['alpha'])
        ttk.Scale(row, from_=0.1, to=1.0, orient=tk.HORIZONTAL, variable=var,
                  command=lambda _v: self._set_scatter('alpha', var)).pack(
            side=tk.LEFT, fill=tk.X, expand=True)

        tab_legend = ttk.Frame(notebook)
        notebook.add(tab_legend, text="Legend")
        self._add_visibility(tab_legend, "Show legend", 'legend')
        self._add_font_scale(tab_legend, "Legend font size", 'legend', 6, 16)
        self._add_color_row(tab_legend, "Text color", 'legend')

        tab_function = ttk.Frame(notebook)
        notebook.add(tab_function, text="Gene Function")
        ttk.Checkbutton(tab_function, text="Color scatter points by gene function",
                        variable=self.color_by_function,
                        command=self.redraw).pack(anchor=tk.W, padx=6, pady=6)
        canvas = tk.Canvas(tab_function, highlightthickness=0)
        scrollbar = ttk.Scrollbar(tab_function, orient=tk.VERTICAL,
                                  command=canvas.yview)
        inner = ttk.Frame(canvas)
        inner.bind('<Configure>',
                   lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=inner, anchor='nw')
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        for category in CATEGORY_ORDER:
            rowf = ttk.Frame(inner)
            rowf.pack(fill=tk.X, pady=2)
            vis_var = tk.BooleanVar(value=self.category_visibility.get(category, True))
            ttk.Checkbutton(rowf, text="Show", variable=vis_var,
                            command=lambda c=category, v=vis_var:
                            self.toggle_category(c, v)).pack(side=tk.LEFT)
            text = category if len(category) <= 32 else category[:30] + '...'
            ttk.Label(rowf, text=text, width=34, anchor=tk.W).pack(
                side=tk.LEFT, padx=(5, 5))
            color_btn = ttk.Button(rowf, text="", width=10)
            color_btn.configure(command=lambda c=category, b=color_btn:
                                self.choose_category_color(c, b))
            color_btn.pack(side=tk.LEFT)
            self.update_color_button(color_btn, self.get_category_color(category))

    def _set_linestyle(self, style):
        self.style_settings['refline']['linestyle'] = style
        self.redraw()

    def _set_scatter(self, key, var):
        try:
            self.style_settings['scatter'][key] = float(var.get())
        except (TypeError, ValueError):
            return

    def toggle_category(self, category, var):
        self.category_visibility[category] = var.get()

    def choose_category_color(self, category, button):
        color = colorchooser.askcolor(title=f"Choose color: {category}")[1]
        if color:
            self.category_colors[category] = color
            self.update_color_button(button, color)
            self.redraw()


def main():
    root = tk.Tk()
    app = COAPlotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
