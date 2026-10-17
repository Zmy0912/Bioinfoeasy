# -*- coding: utf-8 -*-
"""COA Simple Plotter - only the orthogonal-axes main plot (Axis 1 x Axis 2).

A stripped-down companion of `coa_plot_plotter.py`: the scree plot and the
correlation panels are removed, and only the projection of the genes on the
first two (orthogonal) COA axes is drawn.

Supported features (same input modes as the previous programs):

  * single-species mode  : one genes.coa  -> one square main plot
  * folder mode          : every subfolder containing genes.coa becomes one
                           panel of a combined grid figure (shared legend)
  * optional per-file export of individual images (folder mode)
  * axis labels annotated with the fraction of inertia explained
    (when eigen.coa is available): Axis 1 (9.57%) / Axis 2 (8.31%)
  * colouring by gene functional category, outlying-gene labels,
    adjustable style settings, scrollable canvas, image export

Input files (CodonW COA output, same as coa_plot_plotter.py):

  genes.coa   required   gene coordinates on the COA axes
  eigen.coa   optional   eigenvalues / fraction of variation per axis
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
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

GENE_SYMBOL_TERMINATORS = '_( [{:|,;.'

LABEL_METHOD_ORIGIN = 'Farthest from origin (0, 0)'
LABEL_METHOD_CENTER = 'Farthest from data center'
LABEL_METHODS = (LABEL_METHOD_ORIGIN, LABEL_METHOD_CENTER)

# CodonW COA outputs that are NOT gene-coordinate tables; they are ignored
# when several .coa files share one folder
NONGENE_COA_FILES = {
    'eigen.coa', 'codon.coa', 'cusort.coa', 'hilo.coa', 'fop.coa',
    'cai.coa', 'cbi.coa', 'summary.coa', 'aop.coa', 'inertia.coa',
    'base.coa', 'coa.coa',
}


def default_style_settings():
    """Fresh copy of the default style settings (used on reset)."""
    return {
        'title': {'fontsize': 13, 'fontweight': 'bold', 'color': 'black',
                  'visible': True},
        'xlabel': {'fontsize': 11, 'fontweight': 'bold', 'color': 'black',
                   'visible': True},
        'ylabel': {'fontsize': 11, 'fontweight': 'bold', 'color': 'black',
                   'visible': True},
        'xticklabels': {'fontsize': 9, 'color': 'black'},
        'yticklabels': {'fontsize': 9, 'color': 'black'},
        'refline': {'color': 'gray', 'linestyle': '--', 'linewidth': 1.0,
                    'visible': True},
        'scatter': {'color': '#1f77b4', 'edgecolor': 'darkblue', 'size': 40,
                    'alpha': 0.7, 'linewidth': 0.5, 'visible': True},
        'legend': {'fontsize': 9, 'color': 'black', 'visible': True},
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _is_float(token):
    try:
        float(token)
        return True
    except (TypeError, ValueError):
        return False


def extract_gene_symbol(gene_name):
    """Bare gene symbol of a raw identifier (``rpl20_[...](-)|01`` -> ``rpl20``)."""
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
    """Gene symbol of a `genes.coa` label (leading running number stripped)."""
    s = str(label or '').strip()
    m = re.match(r'^\d+[_\-](.*)$', s)
    if m:
        s = m.group(1)
    return extract_gene_symbol(s)


def classify_gene_symbol(symbol):
    gene_id = str(symbol).split('_')[0].split('|')[0].strip().lower()
    if not gene_id:
        return 'Other / Unclassified'
    if gene_id in SPECIFIC_GENE_FUNCTIONS:
        return SPECIFIC_GENE_FUNCTIONS[gene_id]
    for prefix, category in FAMILY_GENE_FUNCTIONS:
        if gene_id.startswith(prefix):
            return category
    return 'Other / Unclassified'


def _read_lines(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.readlines()


# ---------------------------------------------------------------------------
# CodonW file parsers
# ---------------------------------------------------------------------------
def parse_genes_coa(filepath):
    """Parse `genes.coa`: labels + gene coordinates on the COA axes."""
    labels, rows, n_axes = [], [], 0
    for raw in _read_lines(filepath):
        line = raw.strip()
        if not line:
            continue
        if line.lower().startswith('label'):
            toks = line.split()
            if len(toks) > 1:
                n_axes = len(toks) - 1
            continue

        tokens = line.split()
        k = 0
        while k < len(tokens) and _is_float(tokens[-1 - k]):
            k += 1
        if k < 2 or len(tokens) - k < 1:
            continue
        if n_axes and k >= n_axes:
            k = n_axes
        labels.append(' '.join(tokens[:-k]))
        rows.append([float(v) for v in tokens[-k:]])

    if len(labels) < 3:
        return None
    coords = [np.array([row[i] for row in rows], dtype=float)
              for i in range(len(rows[0]))]
    if len(coords) < 2:
        return None
    return {'labels': labels, 'coords': coords, 'n_genes': len(labels)}


def parse_eigen_coa(filepath):
    """Parse `eigen.coa`: axis number, eigenvalue, fraction, cumulative."""
    rows = []
    for raw in _read_lines(filepath):
        parts = raw.split()
        if len(parts) < 4:
            continue
        try:
            rows.append((int(parts[0]), float(parts[1]),
                         float(parts[2]), float(parts[3])))
        except ValueError:
            continue
    if not rows:
        return None
    return {'axis': [r[0] for r in rows], 'eigen': [r[1] for r in rows],
            'frac': [r[2] for r in rows], 'cum': [r[3] for r in rows]}


def build_dataset(genes_path, eigen_path=None, name=None):
    """Assemble one dataset from the CodonW output files of one species."""
    genes = parse_genes_coa(genes_path)
    if genes is None:
        return None

    eigen = None
    if eigen_path and os.path.isfile(eigen_path):
        eigen = parse_eigen_coa(eigen_path)

    symbols = [coa_gene_symbol(lab) for lab in genes['labels']]
    return {
        'name': name or os.path.basename(genes_path),
        'filename': os.path.basename(genes_path),
        'x': genes['coords'][0],      # Axis 1
        'y': genes['coords'][1],      # Axis 2
        'genes': genes['labels'],
        'symbols': symbols,
        'categories': [classify_gene_symbol(s) for s in symbols],
        'eigen': eigen,
    }


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------
class COASimplePlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("COA Simple Plotter - Orthogonal Axes Main Plot (CodonW)")
        self.root.geometry("1100x900")

        self.style_settings = default_style_settings()

        # Gene labels for the outlying points
        self.label_outliers = tk.BooleanVar(value=True)
        self.label_method = tk.StringVar(value=LABEL_METHOD_ORIGIN)
        self.label_count = tk.IntVar(value=5)
        self.label_fontsize = tk.IntVar(value=8)
        self.label_color = 'darkgreen'

        # Input mode / options
        self.input_mode = tk.StringVar(value='file')
        self.folder_path = tk.StringVar(value='')
        self.grid_columns = tk.StringVar(value='0')
        self.export_individual = tk.BooleanVar(value=False)
        self.color_by_function = tk.BooleanVar(value=False)

        # Datasets: each {'name', 'filename', 'x', 'y', 'genes', 'symbols',
        #                  'categories', 'eigen'}
        self.datasets = []

        self.category_colors = dict(DEFAULT_CATEGORY_COLORS)
        self.category_visibility = {cat: True for cat in CATEGORY_ORDER}

        # Scale variables of the currently open style window (for Apply)
        self._style_scale_vars = {}

        self.create_widgets()

    # ------------------------------------------------------------------ UI
    def create_widgets(self):
        control_frame = ttk.LabelFrame(
            self.root, text="Input Selection (CodonW COA output)", padding=10)
        control_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(control_frame, text="Input mode:").grid(row=0, column=0, sticky=tk.W, pady=5)
        mode_frame = ttk.Frame(control_frame)
        mode_frame.grid(row=0, column=1, columnspan=3, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(mode_frame, text="Single species (one main plot)",
                        variable=self.input_mode, value='file',
                        command=self.on_mode_change).pack(side=tk.LEFT)
        ttk.Radiobutton(mode_frame, text="Folder (all species, combined figure)",
                        variable=self.input_mode, value='folder',
                        command=self.on_mode_change).pack(side=tk.LEFT, padx=(15, 0))

        # genes.coa
        ttk.Label(control_frame, text="genes.coa (required):").grid(
            row=1, column=0, sticky=tk.W, pady=4)
        self.genes_path = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.genes_path, width=52).grid(
            row=1, column=1, padx=5, pady=4, sticky=tk.EW)
        ttk.Button(control_frame, text="Browse...",
                   command=self.browse_genes).grid(row=1, column=2, padx=5, pady=4)

        # eigen.coa
        ttk.Label(control_frame, text="eigen.coa (optional):").grid(
            row=2, column=0, sticky=tk.W, pady=4)
        self.eigen_path = tk.StringVar()
        ttk.Entry(control_frame, textvariable=self.eigen_path, width=52).grid(
            row=2, column=1, padx=5, pady=4, sticky=tk.EW)
        ttk.Button(control_frame, text="Browse...",
                   command=self.browse_eigen).grid(row=2, column=2, padx=5, pady=4)

        # folder mode
        ttk.Label(control_frame, text="Input folder\n(.coa files and/or subfolders):").grid(
            row=3, column=0, sticky=tk.W, pady=4)
        ttk.Entry(control_frame, textvariable=self.folder_path, width=52).grid(
            row=3, column=1, padx=5, pady=4, sticky=tk.EW)
        ttk.Button(control_frame, text="Browse Folder...",
                   command=self.browse_folder).grid(row=3, column=2, padx=5, pady=4)

        # species name / grid columns
        ttk.Label(control_frame, text="Species name:").grid(row=4, column=0, sticky=tk.W, pady=6)
        self.species_name = tk.StringVar(value="Species")
        ttk.Entry(control_frame, textvariable=self.species_name, width=28).grid(
            row=4, column=1, padx=5, pady=6, sticky=tk.W)

        ttk.Label(control_frame, text="Subplot columns (0 = auto):").grid(
            row=5, column=0, sticky=tk.W, pady=5)
        ttk.Spinbox(control_frame, from_=0, to=12, width=6,
                    textvariable=self.grid_columns).grid(
            row=5, column=1, padx=5, pady=5, sticky=tk.W)

        # actions
        action_frame = ttk.Frame(control_frame)
        action_frame.grid(row=4, column=2, rowspan=2, columnspan=2,
                          sticky=tk.E, padx=5)
        ttk.Button(action_frame, text="Plot Graph",
                   command=self.plot_graph).pack(side=tk.LEFT, padx=2)
        ttk.Button(action_frame, text="Style Settings",
                   command=self.open_style_settings).pack(side=tk.LEFT, padx=2)
        ttk.Button(action_frame, text="Save Image",
                   command=self.save_image).pack(side=tk.LEFT, padx=2)
        ttk.Button(action_frame, text="Export Axis Stats (CSV)",
                   command=self.export_axis_stats).pack(side=tk.LEFT, padx=2)

        control_frame.columnconfigure(1, weight=1)

        # ---- Options -------------------------------------------------------
        options_frame = ttk.LabelFrame(self.root, text="Plot Options", padding=8)
        options_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Checkbutton(options_frame, text="Color points by gene function",
                        variable=self.color_by_function,
                        command=self.redraw).pack(side=tk.LEFT)
        ttk.Checkbutton(options_frame,
                        text="Also export one image per file (folder mode)",
                        variable=self.export_individual).pack(side=tk.LEFT, padx=(20, 0))

        # ---- Gene labels ----------------------------------------------------
        label_frame = ttk.LabelFrame(self.root, text="Gene Labels (outliers)", padding=5)
        label_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Checkbutton(label_frame, text="Label outlying genes",
                        variable=self.label_outliers,
                        command=self.redraw).pack(side=tk.LEFT)
        ttk.Label(label_frame, text="  Rank by:").pack(side=tk.LEFT)
        method_box = ttk.Combobox(label_frame, textvariable=self.label_method,
                                  values=list(LABEL_METHODS), width=30,
                                  state='readonly')
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

        # ---- Plot area (scrollable) -----------------------------------------
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

        self.canvas_container = tk.Frame(self.canvas_viewport, width=600,
                                         height=600, bg='white')
        self.canvas_container.pack_propagate(False)
        self.canvas_viewport.create_window((0, 0), window=self.canvas_container,
                                           anchor='nw')

        self.fig, self.ax = plt.subplots(figsize=(7, 7))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_container)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN,
                                    anchor=tk.W)
        self.status_bar.pack(fill=tk.X, padx=5, pady=5)

        self.on_mode_change()

    # ------------------------------------------------------------ callbacks
    def on_mode_change(self):
        if self.input_mode.get() == 'folder':
            self.status_bar.config(
                text="Folder mode: every .coa gene-coordinate file and every subfolder "
                     "containing genes.coa becomes one panel of the combined figure "
                     "(eigen.coa is paired automatically when available)")
        else:
            self.status_bar.config(
                text="Single-species mode: one main plot of the genes on COA axes 1 and 2")

    def browse_genes(self):
        filename = filedialog.askopenfilename(
            title="Select gene-coordinate .coa file (e.g. genes.coa)",
            filetypes=[("COA gene coordinates", "*.coa"),
                       ("All files", "*.*")])
        if filename:
            self.input_mode.set('file')
            self.genes_path.set(filename)
            if not self.eigen_path.get():
                eigen = self.autodetect_eigen(filename)
                if eigen:
                    self.eigen_path.set(eigen)
            self.status_bar.config(text=f"Selected: {os.path.basename(filename)}")

    def browse_eigen(self):
        filename = filedialog.askopenfilename(
            title="Select eigen.coa",
            filetypes=[("COA eigenvalues", "eigen.coa"), ("All files", "*.*")])
        if filename:
            self.eigen_path.set(filename)

    def browse_folder(self):
        folder = filedialog.askdirectory(
            title="Select folder containing species subfolders")
        if folder:
            self.input_mode.set('folder')
            self.folder_path.set(folder)
            n = len(self.find_species_datasets(folder))
            self.status_bar.config(
                text=f"Selected folder: {folder}  ({n} dataset(s) found)")

    @staticmethod
    def eigen_candidates(folder, stem):
        return [os.path.join(folder, f'{stem}_eigen.coa'),
                os.path.join(folder, f'eigen_{stem}.coa')]

    def autodetect_eigen(self, genes_file):
        """Find an eigen.coa next to the selected gene-coordinate file."""
        folder = os.path.dirname(genes_file) or '.'
        stem = os.path.splitext(os.path.basename(genes_file))[0]
        for cand in self.eigen_candidates(folder, stem):
            if os.path.isfile(cand):
                return cand
        cand = os.path.join(folder, 'eigen.coa')
        return cand if os.path.isfile(cand) else None

    def find_species_datasets(self, folder):
        """Datasets in `folder`: loose .coa gene tables + genes.coa subfolders."""
        datasets = []
        for d in sorted(glob.glob(os.path.join(folder, '*'))):
            if not os.path.isdir(d):
                continue
            g = os.path.join(d, 'genes.coa')
            if not os.path.isfile(g):
                continue
            e = os.path.join(d, 'eigen.coa')
            datasets.append({'kind': 'dir', 'genes': g,
                             'eigen': e if os.path.isfile(e) else None,
                             'name': os.path.basename(d)})
        return self._add_loose_coa_files(folder, datasets)

    @staticmethod
    def _add_loose_coa_files(folder, datasets):
        """Append loose *.coa gene tables found directly in `folder`."""
        loose = [f for f in sorted(glob.glob(os.path.join(folder, '*.coa')))
                 if os.path.basename(f).lower() not in NONGENE_COA_FILES]
        if not loose:
            return datasets
        root_eigen = os.path.join(folder, 'eigen.coa')
        root_eigen_ok = os.path.isfile(root_eigen) and (
            any(os.path.basename(f) == 'genes.coa' for f in loose)
            or len(loose) == 1)
        for f in loose:
            stem = os.path.splitext(os.path.basename(f))[0]
            eigen = None
            for cand in COASimplePlotter.eigen_candidates(folder, stem):
                if os.path.isfile(cand):
                    eigen = cand
                    break
            if eigen is None and root_eigen_ok:
                eigen = root_eigen
            datasets.append({'kind': 'file', 'genes': f,
                             'eigen': eigen, 'name': stem})
        return datasets

    # ------------------------------------------------------------- loading
    def load_data(self):
        genes_file = self.genes_path.get()
        if not genes_file or not os.path.exists(genes_file):
            messagebox.showerror("Error", "Please select a valid genes.coa file!")
            return False
        try:
            dataset = build_dataset(genes_file, self.eigen_path.get(),
                                    self.species_name.get().strip()
                                    or os.path.basename(genes_file))
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file: {e}")
            return False
        if dataset is None:
            messagebox.showerror("Error",
                                 "Could not extract valid data from genes.coa!")
            return False

        self.datasets = [dataset]
        self.report_dataset(dataset)
        return True

    def load_folder_data(self):
        folder = self.folder_path.get()
        if not folder or not os.path.isdir(folder):
            messagebox.showerror("Error", "Please select a valid input folder!")
            return False

        items = self.find_species_datasets(folder)
        if not items:
            messagebox.showerror(
                "Error", "No .coa gene-coordinate file or genes.coa subfolder "
                         "was found in the selected folder!")
            return False

        datasets, failed = [], []
        for item in items:
            g, e, name = item['genes'], item['eigen'], item['name']
            try:
                dataset = build_dataset(g, e, name=name)
            except Exception as exc:
                print(f"Failed to read {g}: {exc}")
                dataset = None
            if dataset is None:
                failed.append(name)
            else:
                datasets.append(dataset)

        if not datasets:
            messagebox.showerror("Error", "No valid data could be extracted!")
            return False

        self.datasets = datasets
        if failed:
            print(f"Skipped datasets without valid data: {', '.join(failed)}")
        total = sum(len(d['x']) for d in datasets)
        print(f"Loaded {len(datasets)} species, {total} genes in total")
        self.status_bar.config(
            text=f"Loaded {len(datasets)} species, {total} genes in total")
        return True

    def report_dataset(self, dataset):
        n = len(dataset['x'])
        eig = dataset['eigen']
        if eig:
            print(f"Axis 1 explains {eig['frac'][0]*100:.2f}% of the total "
                  f"inertia, axis 2 {eig['frac'][1]*100:.2f}% (cumulative "
                  f"{eig['cum'][1]*100:.2f}%)")
        self.status_bar.config(
            text=f"Loaded {n} genes"
                 + (f" - Axis1 {eig['frac'][0]*100:.2f}%, "
                    f"Axis2 {eig['frac'][1]*100:.2f}%" if eig else ""))

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
            self.category_colors[category] = DEFAULT_CATEGORY_COLORS.get(
                category, '#404040')
        return self.category_colors[category]

    def axis_label(self, dataset, i):
        eig = dataset.get('eigen')
        if eig and i - 1 < len(eig['frac']):
            return f'Axis {i} ({eig["frac"][i - 1] * 100:.2f}%)'
        return f'Axis {i}'

    def get_symmetric_range(self, x_values, y_values):
        """ONE square range centred on the COA origin (0, 0).

        Both axes get the same limits, so the X and Y numeric spans are
        equal and the orthogonal axes stay centred.  This avoids the empty
        band ("gap") that appears when ``adjustable='datalim'`` silently
        stretches one axis to satisfy the equal-aspect constraint.
        """
        xs = np.abs(np.asarray(x_values, dtype=float))
        ys = np.abs(np.asarray(y_values, dtype=float))
        finite = np.concatenate([xs[np.isfinite(xs)], ys[np.isfinite(ys)]])
        half = (float(np.nanmax(finite)) * 1.1) if finite.size else 1.0
        half = max(half, 0.1)
        return -half, half

    def compute_shared_range(self, datasets):
        """Axis range shared by every panel of a combined figure."""
        xs = [v for d in datasets for v in d['x']]
        ys = [v for d in datasets for v in d['y']]
        if not xs or not ys:
            return -1.0, 1.0
        return self.get_symmetric_range(xs, ys)

    # -- outliers -----------------------------------------------------------
    def _label_count(self):
        try:
            count = int(float(self.label_count.get()))
        except (TypeError, ValueError):
            count = 5
        return max(0, count)

    def find_outlier_genes(self, dataset, xlim=None, ylim=None):
        count = self._label_count()
        if count <= 0 or len(dataset['x']) == 0:
            return []

        xs = np.asarray(dataset['x'], dtype=float)
        ys = np.asarray(dataset['y'], dtype=float)

        if xlim is None:
            xlim = (float(np.nanmin(xs)), float(np.nanmax(xs)))
        if ylim is None:
            ylim = (float(np.nanmin(ys)), float(np.nanmax(ys)))
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

    # -- main panel -----------------------------------------------------------
    def draw_scatter_by_function(self, ax, dataset):
        style = self.style_settings['scatter']
        groups = {}
        for x, y, cat in zip(dataset['x'], dataset['y'],
                             dataset['categories']):
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

    def draw_panel(self, ax, dataset, title_text, add_legend=True,
                   title_fontsize=None, axis_range=None):
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

        ax.tick_params(axis='x',
                       labelsize=style['xticklabels']['fontsize'],
                       labelcolor=style['xticklabels']['color'])
        ax.tick_params(axis='y',
                       labelsize=style['yticklabels']['fontsize'],
                       labelcolor=style['yticklabels']['color'])

        if style['title']['visible'] and title_text:
            size = title_fontsize if title_fontsize is not None \
                else style['title']['fontsize']
            ax.set_title(title_text, fontsize=size,
                         fontweight=style['title']['fontweight'],
                         color=style['title']['color'], pad=12)

        # ONE square, origin-centred range keeps the orthogonal axes
        # centred and prevents empty bands at the plot borders
        if axis_range is None:
            axis_range = self.get_symmetric_range(dataset['x'], dataset['y'])
        ax.set_xlim(*axis_range)
        ax.set_ylim(*axis_range)

        # Keep the plot square (adjustable='box' shrinks the frame instead
        # of stretching the data limits)
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, alpha=0.3, linestyle=':')

        if add_legend and style['legend']['visible']:
            ax.legend(loc='best', fontsize=style['legend']['fontsize'],
                      labelcolor=style['legend']['color'])

        self.annotate_outliers(ax, dataset, ax.get_xlim(), ax.get_ylim())

    # -- figure assembly -------------------------------------------------------
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

        self.fig.clear()
        self.fig.set_size_inches(7, 7)
        self.ax = self.fig.add_subplot(111)
        self.draw_panel(self.ax, dataset, dataset['name'], add_legend=True)

        self.fig.tight_layout()
        self.fit_canvas_to_figure()
        self.canvas.draw()

        eig = dataset.get('eigen')
        extra = (f" - Axis1 {eig['frac'][0]*100:.2f}%, "
                 f"Axis2 {eig['frac'][1]*100:.2f}%") if eig else ""
        self.status_bar.config(
            text=f"Orthogonal-axes plot: {len(dataset['x'])} genes{extra}")

    def plot_combined(self):
        if not self.load_folder_data():
            return
        self.render_combined()

    def render_combined(self):
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
        nrows = math.ceil(n / ncols)

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
        axis_range = self.compute_shared_range(datasets)
        for i, dataset in enumerate(datasets):
            r, c = divmod(i, ncols)
            ax = self.fig.add_subplot(gs[r, c])
            self.draw_panel(ax, dataset, dataset['name'],
                            add_legend=False,
                            title_fontsize=title_fontsize,
                            axis_range=axis_range)

        self.add_side_legend(right)
        self.fit_canvas_to_figure()
        self.canvas.draw()

        total = sum(len(d['x']) for d in datasets)
        self.status_bar.config(
            text=f"Combined orthogonal-axes figure: {n} species, "
                 f"{total} genes, {nrows}x{ncols} grid")

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

    # -- saving -----------------------------------------------------------------
    def save_image(self):
        if not self.datasets:
            messagebox.showwarning("Warning", "Please plot the graph first!")
            return
        filename = filedialog.asksaveasfilename(
            title="Save Image", defaultextension=".pdf",
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

        saved_msg = f"Image saved to:\n{filename}"
        self.status_bar.config(text=f"Image saved: {os.path.basename(filename)}")

        # Optionally also export one image per input file (folder mode only)
        if self.input_mode.get() == 'folder' and self.export_individual.get():
            individual_dir = os.path.join(
                os.path.dirname(filename), 'individual_plots')
            ext = os.path.splitext(filename)[1] or '.pdf'
            try:
                os.makedirs(individual_dir, exist_ok=True)
                count = self.export_individual_images(individual_dir, ext)
                saved_msg += (f"\n\nAlso exported {count} individual image(s) "
                              f"to:\n{individual_dir}")
            except Exception as e:
                messagebox.showwarning(
                    "Partial success",
                    f"{saved_msg}\n\nBut exporting individual images failed:\n{e}")
                return

        messagebox.showinfo("Success", saved_msg)

    # -- statistics export -------------------------------------------------
    def export_axis_stats(self):
        """Export a paper-ready table of the variation explained by each axis."""
        if not self.datasets:
            messagebox.showwarning("Warning", "Please plot the graph first!")
            return
        with_eigen = [d for d in self.datasets if d.get('eigen')]
        if not with_eigen:
            messagebox.showwarning(
                "Warning", "No eigen.coa data loaded.\n"
                           "The variation explained by each axis comes from "
                           "CodonW's eigen.coa file.\n"
                           "Place eigen.coa next to the .coa files "
                           "(or select it in single-file mode) and plot again.")
            return
        filename = filedialog.asksaveasfilename(
            title="Export Axis Statistics", defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not filename:
            return
        try:
            rows = self._write_axis_stats(filename, with_eigen)
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")
            return
        self.status_bar.config(
            text=f"Axis statistics exported: {os.path.basename(filename)}")
        messagebox.showinfo(
            "Success",
            f"Axis statistics exported to:\n{filename}\n\n"
            f"{len(rows)} data rows (Table 1: per-axis eigenvalues; "
            f"Table 2: summary per species)")

    @staticmethod
    def _write_axis_stats(path, datasets):
        """Write the CSV tables and return the flat list of data rows."""
        rows = []
        datasets = [d for d in datasets if d.get('eigen')]
        with open(path, 'w', encoding='utf-8-sig', newline='') as f:
            f.write('# Variation explained by the orthogonal axes of the '
                    'Correspondence Analysis (CodonW eigen.coa)\n')
            f.write('Table 1. Variation explained by each COA axis per '
                    'species\n')
            f.write('species,n_genes,axis,eigenvalue,'
                    'variation_explained_%,cumulative_%,total_inertia\n')
            for d in datasets:
                eig = d['eigen']
                inertia = (eig['eigen'][0] / eig['frac'][0]
                           if eig['frac'][0] > 0 else float('nan'))
                for i, (ev, frac, cum) in enumerate(
                        zip(eig['eigen'], eig['frac'], eig['cum']), 1):
                    row = [d['name'], len(d['x']), i, f'{ev:.6f}',
                           f'{frac*100:.2f}', f'{cum*100:.2f}',
                           f'{inertia:.6f}']
                    f.write(','.join(str(v) for v in row) + '\n')
                    rows.append(row)
            f.write('\n')
            f.write('Table 2. Summary per species (first four orthogonal '
                    'axes)\n')
            f.write('species,n_genes,n_axes,total_inertia,'
                    'axis1_%,axis2_%,axis3_%,axis4_%,axis1+2_cumulative_%\n')
            for d in datasets:
                eig = d['eigen']
                inertia = (eig['eigen'][0] / eig['frac'][0]
                           if eig['frac'][0] > 0 else float('nan'))
                fr = list(eig['frac']) + [float('nan')] * 4
                row = [d['name'], len(d['x']), len(eig['frac']),
                       f'{inertia:.6f}']
                for i in range(4):
                    row.append(f'{fr[i]*100:.2f}')
                row.append(f'{eig["cum"][1]*100:.2f}'
                           if len(eig['cum']) > 1 else '')
                f.write(','.join(str(v) for v in row) + '\n')
                rows.append(row)
        return rows

    def export_individual_images(self, out_dir, extension):
        """Draw and save one orthogonal-axes plot per dataset."""
        count = 0
        for dataset in self.datasets:
            fig = plt.figure(figsize=(7, 7))
            ax = fig.add_subplot(111)
            self.draw_panel(ax, dataset, dataset['name'], add_legend=True)
            fig.tight_layout()
            safe_name = os.path.splitext(
                dataset.get('filename', dataset['name']))[0]
            out_path = os.path.join(out_dir, f'{safe_name}{extension}')
            fig.savefig(out_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            count += 1
        return count

    # -- style settings -----------------------------------------------------------
    def _add_font_scale(self, parent, label, section, lo, hi):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text=label, width=22).pack(side=tk.LEFT)
        var = tk.DoubleVar(value=self.style_settings[section]['fontsize'])
        ttk.Scale(row, from_=lo, to=hi, orient=tk.HORIZONTAL, variable=var,
                  command=lambda _v, s=section, v=var:
                  self._set_fontsize(s, v)).pack(
            side=tk.LEFT, fill=tk.X, expand=True)
        self._style_scale_vars[f'{section}.fontsize'] = var

    def _set_fontsize(self, section, var):
        try:
            self.style_settings[section]['fontsize'] = float(var.get())
        except (TypeError, ValueError):
            return

    def _add_value_scale(self, parent, label, section, key, lo, hi):
        """Generic numeric style control registered for the Apply button."""
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text=label, width=22).pack(side=tk.LEFT)
        var = tk.DoubleVar(value=float(self.style_settings[section][key]))
        ttk.Scale(row, from_=lo, to=hi, orient=tk.HORIZONTAL,
                  variable=var).pack(side=tk.LEFT, fill=tk.X, expand=True)
        self._style_scale_vars[f'{section}.{key}'] = var

    def _add_color_row(self, parent, label, section):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, padx=6, pady=2)
        ttk.Label(row, text=label, width=22).pack(side=tk.LEFT)
        btn = ttk.Button(row, text=self.style_settings[section]['color'],
                         width=12,
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
        self._add_visibility(tab_labels, "Show X-axis label", 'xlabel')
        self._add_font_scale(tab_labels, "X label font size", 'xlabel', 8, 20)
        self._add_color_row(tab_labels, "X label color", 'xlabel')
        self._add_visibility(tab_labels, "Show Y-axis label", 'ylabel')
        self._add_font_scale(tab_labels, "Y label font size", 'ylabel', 8, 20)
        self._add_color_row(tab_labels, "Y label color", 'ylabel')
        self._add_font_scale(tab_labels, "X tick font size", 'xticklabels', 6, 18)
        self._add_color_row(tab_labels, "X tick color", 'xticklabels')
        self._add_font_scale(tab_labels, "Y tick font size", 'yticklabels', 6, 18)
        self._add_color_row(tab_labels, "Y tick color", 'yticklabels')

        tab_ref = ttk.Frame(notebook)
        notebook.add(tab_ref, text="Reference Lines")
        self._add_visibility(tab_ref, "Show zero reference lines", 'refline')
        self._add_value_scale(tab_ref, "Line width", 'refline',
                              'linewidth', 0.5, 5.0)
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
        self._add_color_row(tab_scatter, "Edge color", 'scatter')
        self._add_value_scale(tab_scatter, "Edge width", 'scatter',
                              'linewidth', 0.0, 3.0)
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
        ttk.Checkbutton(tab_function,
                        text="Color scatter points by gene function",
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
            vis_var = tk.BooleanVar(
                value=self.category_visibility.get(category, True))
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
            self.update_color_button(color_btn,
                                     self.get_category_color(category))

        # ---- Apply / Reset buttons (like pr2_plot_plotter.py) ----
        button_frame = ttk.Frame(win)
        button_frame.pack(fill=tk.X, padx=8, pady=8)
        ttk.Button(button_frame, text="Apply",
                   command=lambda: self.apply_settings(win)).pack(
            side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Reset to Default",
                   command=lambda: self.reset_settings(win)).pack(
            side=tk.RIGHT, padx=5)

    def apply_settings(self, window):
        """Apply every scale control of the style window and repaint."""
        for key, var in self._style_scale_vars.items():
            section, attr = key.split('.', 1)
            try:
                self.style_settings[section][attr] = float(var.get())
            except (TypeError, ValueError):
                continue
        self.redraw()
        messagebox.showinfo("Success", "Style settings applied.")
        window.destroy()

    def reset_settings(self, window):
        """Restore all style settings to their defaults."""
        self.style_settings = default_style_settings()
        self.color_by_function.set(False)
        self.category_colors = dict(DEFAULT_CATEGORY_COLORS)
        self.category_visibility = {cat: True for cat in CATEGORY_ORDER}
        self.label_outliers.set(True)
        self.label_method.set(LABEL_METHOD_ORIGIN)
        self.label_count.set(5)
        self.label_fontsize.set(8)
        self.label_color = 'darkgreen'
        self.update_color_button(self.label_color_btn, self.label_color)
        self.grid_columns.set('0')
        self.export_individual.set(False)
        self.redraw()
        messagebox.showinfo("Success", "Settings reset to default.")
        window.destroy()

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
    app = COASimplePlotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
