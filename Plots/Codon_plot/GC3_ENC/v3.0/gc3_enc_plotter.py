import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import os
import glob
import math


# ---------------------------------------------------------------------------
# Gene functional classification (standard chloroplast genome categories)
# ---------------------------------------------------------------------------
# Canonical display order of the functional categories
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

# Default color for every functional category
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

# Exact gene names (lowercase) with a dedicated category
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

# Gene families resolved by name prefix (checked after SPECIFIC_GENE_FUNCTIONS)
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


# ---------------------------------------------------------------------------
# Outlier gene labels
# ---------------------------------------------------------------------------
# Fixed axis ranges (GC3 proportion / ENC score). They are also used to
# normalise the distances when looking for the most "outlying" genes.
AXIS_XLIM = (0.0, 1.0)
AXIS_YLIM = (20.0, 64.05)

# How the genes that get a label are ranked
LABEL_METHOD_CENTER = 'Farthest from data center'
LABEL_METHOD_CURVE = 'Farthest from theoretical curve'
LABEL_METHODS = (LABEL_METHOD_CENTER, LABEL_METHOD_CURVE)

# First character that terminates the gene symbol inside a raw identifier
GENE_SYMBOL_TERMINATORS = '_( [{:|,;.'


def extract_gene_symbol(gene_name):
    """Return only the bare gene symbol of a raw ``.out`` identifier.

    The identifiers produced by the annotation pipeline are FASTA headers
    (with or without the leading ``>``) followed by the gene location and
    orientation, for example::

        >rpl20_[69725:70078](-)|01   ->  rpl20
        psbB_[73610:75136](+)|01.    ->  psbB
        atpF_join{[13161:13305](-    ->  atpF
        psbK_[7404:7589](+)|01.fa    ->  psbK
        trnL-UAA_[107000:107063]|01. ->  trnL-UAA

    Everything in front of the last ``>`` is dropped, then the symbol is taken
    up to the first special character (``_``, ``(``, ``[``, ``{``, ``:``,
    ``|``, ``,``, ``;``, ``.`` or a space). Characters that really belong to a
    gene name (letters, digits and ``-``) are preserved.
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


def classify_gene_function(gene_name):
    """Return the functional category of a gene based on its name.

    Handles names such as ``rpl20_[69725:70078](-)|01`` by keeping only the
    leading gene identifier (``rpl20``).
    """
    gene_id = str(gene_name).split('_')[0].split('|')[0].strip().lower()
    if not gene_id:
        return 'Other / Unclassified'
    if gene_id in SPECIFIC_GENE_FUNCTIONS:
        return SPECIFIC_GENE_FUNCTIONS[gene_id]
    for prefix, category in FAMILY_GENE_FUNCTIONS:
        if gene_id.startswith(prefix):
            return category
    return 'Other / Unclassified'


class GC3ENCPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("GC3-ENC Plotter")
        self.root.geometry("1100x900")
        
        # Style settings - default values
        self.style_settings = {
            'title': {
                'fontsize': 14,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xlabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'ylabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'yticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'curve': {
                'color': 'red',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'scatter': {
                'color': 'blue',
                'edgecolor': 'darkblue',
                'size': 50,
                'alpha': 0.6,
                'linewidth': 0.5,
                'visible': True
            },
            'legend': {
                'fontsize': 10,
                'color': 'black',
                'visible': True
            }
        }
        
        # Data storage (mirrors the single dataset for legacy helpers)
        self.gc3_data = []         # GC3 (0-1)  - X axis
        self.enc_data = []         # ENC        - Y axis
        self.gene_names = []
        self.gene_categories = []  # functional category of each point

        # Input mode: 'file' (single .out file) or 'folder' (all .out files)
        self.input_mode = tk.StringVar(value='file')
        self.folder_path = tk.StringVar(value=self._default_out_folder())
        # Number of subplot columns for the combined figure (0 = auto)
        self.grid_columns = tk.StringVar(value='0')
        # Also export one image per input file when saving the combined figure
        self.export_individual = tk.BooleanVar(value=False)

        # Datasets currently used for drawing.
        # Each item: {'name', 'filename', 'gc3', 'enc', 'genes', 'categories'}
        self.datasets = []

        # Gene function coloring
        self.color_by_function = tk.BooleanVar(value=False)
        self.category_colors = dict(DEFAULT_CATEGORY_COLORS)
        self.category_visibility = {cat: True for cat in CATEGORY_ORDER}

        # Labels of the genes that lie farthest from the bulk of the data
        self.label_outliers = tk.BooleanVar(value=False)
        self.label_method = tk.StringVar(value=LABEL_METHOD_CENTER)
        self.label_count = tk.IntVar(value=5)
        self.label_fontsize = tk.IntVar(value=7)
        self.label_color = 'darkgreen'
        
        # Create interface
        self.create_widgets()

    def _default_out_folder(self):
        """Guess the ``out`` folder located next to the project directory."""
        try:
            project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            candidate = os.path.join(project_dir, 'out')
            if os.path.isdir(candidate):
                return candidate
        except Exception:
            pass
        return ''
        
    def create_widgets(self):
        # Top control panel
        control_frame = ttk.LabelFrame(self.root, text="Input Selection", padding=10)
        control_frame.pack(fill=tk.X, padx=10, pady=5)

        # ---- Input mode ----
        ttk.Label(control_frame, text="Input mode:").grid(row=0, column=0, sticky=tk.W, pady=5)
        mode_frame = ttk.Frame(control_frame)
        mode_frame.grid(row=0, column=1, columnspan=2, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(mode_frame, text="Single file",
                        variable=self.input_mode, value='file',
                        command=self.on_mode_change).pack(side=tk.LEFT)
        ttk.Radiobutton(mode_frame, text="Folder (all .out files, combined figure)",
                        variable=self.input_mode, value='folder',
                        command=self.on_mode_change).pack(side=tk.LEFT, padx=(15, 0))
        
        # File path display
        ttk.Label(control_frame, text="Input file (.out):").grid(row=1, column=0, sticky=tk.W, pady=5)

        self.file_path = tk.StringVar()
        file_entry = ttk.Entry(control_frame, textvariable=self.file_path, width=50)
        file_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        browse_btn = ttk.Button(control_frame, text="Browse File...", command=self.browse_file)
        browse_btn.grid(row=1, column=2, padx=5, pady=5)

        # ---- Folder path ----
        ttk.Label(control_frame, text="Input folder:").grid(row=2, column=0, sticky=tk.W, pady=5)

        folder_entry = ttk.Entry(control_frame, textvariable=self.folder_path, width=50)
        folder_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)

        browse_folder_btn = ttk.Button(control_frame, text="Browse Folder...",
                                       command=self.browse_folder)
        browse_folder_btn.grid(row=2, column=2, padx=5, pady=5)

        # ---- Grid columns for the combined figure ----
        ttk.Label(control_frame, text="Subplot columns (0 = auto):").grid(
            row=3, column=0, sticky=tk.W, pady=5)
        col_spin = ttk.Spinbox(control_frame, from_=0, to=12, width=6,
                               textvariable=self.grid_columns)
        col_spin.grid(row=3, column=1, sticky=tk.W, padx=5, pady=5)

        # ---- Species name / action buttons ----
        ttk.Label(control_frame, text="Species name:").grid(row=4, column=0, sticky=tk.W, pady=5)
        self.species_name = tk.StringVar(value="Species")
        name_entry = ttk.Entry(control_frame, textvariable=self.species_name, width=30)
        name_entry.grid(row=4, column=1, sticky=tk.W, padx=5, pady=5)

        plot_btn = ttk.Button(control_frame, text="Plot Graph", command=self.plot_graph)
        plot_btn.grid(row=4, column=2, padx=5, pady=5, sticky=tk.W)

        style_btn = ttk.Button(control_frame, text="Style Settings",
                               command=self.open_style_settings)
        style_btn.grid(row=4, column=3, padx=5, pady=5)

        save_btn = ttk.Button(control_frame, text="Save Image", command=self.save_image)
        save_btn.grid(row=4, column=4, padx=5, pady=5)

        # ---- Options ----
        options_frame = ttk.Frame(control_frame)
        options_frame.grid(row=5, column=0, columnspan=5, sticky=tk.W, padx=5, pady=5)

        ttk.Checkbutton(options_frame,
                        text="Color points by gene function",
                        variable=self.color_by_function,
                        command=self.on_function_toggle).pack(side=tk.LEFT)

        ttk.Checkbutton(options_frame,
                        text="Also export one image per file (folder mode)",
                        variable=self.export_individual).pack(side=tk.LEFT, padx=(20, 0))

        control_frame.columnconfigure(1, weight=1)

        # ---- Gene labels for the outlying points ----
        label_frame = ttk.LabelFrame(control_frame, text="Gene Labels (outliers)",
                                     padding=5)
        label_frame.grid(row=6, column=0, columnspan=5, sticky=tk.EW, padx=5, pady=(5, 0))

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
        label_color_btn = ttk.Button(
            label_frame, text="", width=8,
            command=lambda: self.choose_label_color(label_color_btn))
        label_color_btn.pack(side=tk.LEFT, padx=(4, 0))
        self.label_color_btn = label_color_btn
        self.update_color_button(label_color_btn, self.label_color)
        
        # Plot area: a scrollable viewport so that large (combined) figures
        # can still be displayed at their true size.
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

        # Fixed-size container inside the scrollable viewport
        self.canvas_container = tk.Frame(self.canvas_viewport, width=600, height=600,
                                         bg='white')
        self.canvas_container.pack_propagate(False)  # Prevent the container from shrinking
        self.canvas_window = self.canvas_viewport.create_window(
            (0, 0), window=self.canvas_container, anchor='nw')

        # Create matplotlib figure (square size, fixed)
        self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_container)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Status bar
        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(fill=tk.X, padx=5, pady=5)

        self.on_mode_change()

    def on_mode_change(self):
        """Update the status bar / hint when the input mode changes."""
        if self.input_mode.get() == 'folder':
            self.status_bar.config(
                text="Folder mode: all .out files in the folder will be plotted into one combined figure")
        else:
            self.status_bar.config(text="Single-file mode")

    def fit_canvas_to_figure(self):
        """Size the canvas container to the figure's natural pixel size.

        The container lives inside a scrollable viewport, so even the large
        combined figures are rendered at their true (full) size instead of
        being squeezed into the window.
        """
        try:
            w_in, h_in = self.fig.get_size_inches()
            if w_in <= 0 or h_in <= 0:
                return

            # A minimum preview size keeps very small figures readable
            display_dpi = max(int(self.fig.dpi), 100)
            width = max(600, int(round(w_in * display_dpi)))
            height = max(600, int(round(h_in * display_dpi)))

            self.canvas_container.config(width=width, height=height)
            self.canvas_viewport.update_idletasks()
            self.canvas_viewport.config(scrollregion=(0, 0, width, height))
        except Exception:
            pass
        
    def browse_file(self):
        """Browse and select a single .out file"""
        filename = filedialog.askopenfilename(
            title="Select .out file",
            filetypes=[("OUT files", "*.out"), ("All files", "*.*")]
        )
        if filename:
            self.input_mode.set('file')
            self.file_path.set(filename)
            self.status_bar.config(text=f"Selected: {os.path.basename(filename)}")

    def browse_folder(self):
        """Browse and select a folder containing .out files"""
        folder = filedialog.askdirectory(title="Select folder containing .out files")
        if folder:
            self.input_mode.set('folder')
            self.folder_path.set(folder)
            n = len(glob.glob(os.path.join(folder, '*.out')))
            self.status_bar.config(
                text=f"Selected folder: {folder}  ({n} .out file(s) found)")
            
    def parse_out_file(self, filepath):
        """Parse a single .out file and return a dataset dict or ``None``."""
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        # Skip header line
        data_lines = lines[1:]

        gc3_data, enc_data = [], []
        gene_names, gene_categories = [], []

        for line in data_lines:
            line = line.strip()
            if not line:
                continue

            parts = line.split()

            # Need at least 10 columns (indices 0-9): name ... Nc (8), GC3s (9)
            if len(parts) < 10:
                continue

            try:
                gene_name = parts[0]
                enc_str = parts[8]   # Nc   -> ENC
                gc3_str = parts[9]   # GC3s -> GC3

                # Skip invalid data
                if enc_str == '*****' or gc3_str == '*****':
                    continue

                gc3 = float(gc3_str)
                enc = float(enc_str)

                gc3_data.append(gc3)
                enc_data.append(enc)
                gene_names.append(gene_name)
                gene_categories.append(classify_gene_function(gene_name))

            except (ValueError, IndexError) as e:
                print(f"Error parsing line: {e}")
                continue

        if len(gc3_data) == 0:
            return None

        return {
            'name': os.path.basename(filepath),
            'gc3': gc3_data,
            'enc': enc_data,
            'genes': gene_names,
            'categories': gene_categories,
        }

    def load_data(self):
        """Load a single .out file into ``self.datasets``."""
        filepath = self.file_path.get()

        if not filepath or not os.path.exists(filepath):
            messagebox.showerror("Error", "Please select a valid input file!")
            return False

        try:
            dataset = self.parse_out_file(filepath)
            if dataset is None:
                messagebox.showerror("Error", "Could not extract valid data from file!")
                return False

            # Use the species name as the display name/title when provided
            dataset['name'] = self.species_name.get().strip() or os.path.basename(filepath)
            dataset['filename'] = os.path.basename(filepath)

            self.datasets = [dataset]
            self._sync_single_arrays(dataset)

            print(f"Loaded {len(dataset['gc3'])} data points")
            print(f"GC3 range: {min(dataset['gc3']):.3f} to {max(dataset['gc3']):.3f}")
            print(f"ENC range: {min(dataset['enc']):.2f} to {max(dataset['enc']):.2f}")

            self.status_bar.config(
                text=f"Successfully loaded {len(dataset['gc3'])} data points")
            return True

        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file: {str(e)}")
            print(f"Error: {e}")
            return False

    def load_folder_data(self):
        """Load every .out file inside the selected folder into ``self.datasets``."""
        folder = self.folder_path.get()

        if not folder or not os.path.isdir(folder):
            messagebox.showerror("Error", "Please select a valid input folder!")
            return False

        files = sorted(glob.glob(os.path.join(folder, '*.out')))
        if not files:
            messagebox.showerror("Error", "No .out files found in the selected folder!")
            return False

        datasets = []
        failed = []
        for filepath in files:
            try:
                dataset = self.parse_out_file(filepath)
            except Exception as e:
                print(f"Failed to read {filepath}: {e}")
                dataset = None
            if dataset is None:
                failed.append(os.path.basename(filepath))
                continue
            datasets.append(dataset)

        if not datasets:
            messagebox.showerror("Error", "No valid data could be extracted from any file!")
            return False

        self.datasets = datasets
        if failed:
            print(f"Skipped files without valid data: {', '.join(failed)}")

        total = sum(len(d['gc3']) for d in datasets)
        print(f"Loaded {len(datasets)} files, {total} data points in total")
        self.status_bar.config(
            text=f"Loaded {len(datasets)} file(s), {total} data points in total")
        return True

    def _sync_single_arrays(self, dataset):
        """Mirror a single dataset into the legacy per-point arrays."""
        self.gc3_data = dataset['gc3']
        self.enc_data = dataset['enc']
        self.gene_names = dataset['genes']
        self.gene_categories = dataset['categories']

    def calculate_enc_curve(self):
        """Calculate theoretical ENC curve (Wright 1990)"""
        gc3_values = []
        enc_values = []

        # Generate GC3 values from 0.05 to 0.95
        for i in range(5, 96):
            gc3 = i / 100.0

            # Theoretical ENC calculation formula (Wright 1990)
            # ENC = 2 + GC3 + 29 / (GC3^2 + (1-GC3)^2)
            enc = 2 + gc3 + 29 / (gc3**2 + (1-gc3)**2)

            gc3_values.append(gc3)
            enc_values.append(enc)

        return np.array(gc3_values), np.array(enc_values)

    def get_active_categories(self, categories=None):
        """Categories present in the data, sorted by the canonical order."""
        if categories is None:
            categories = self.gene_categories
        present = set(categories)
        ordered = [cat for cat in CATEGORY_ORDER if cat in present]
        ordered += sorted(cat for cat in present if cat not in CATEGORY_ORDER)
        return ordered

    def get_category_color(self, category):
        """Return the color assigned to a functional category."""
        if category not in self.category_colors:
            self.category_colors[category] = DEFAULT_CATEGORY_COLORS.get(
                category, '#404040')
        return self.category_colors[category]

    def draw_scatter_by_function(self, ax, dataset):
        """Draw scatter points, one color per gene functional category."""
        style = self.style_settings['scatter']

        groups = {}
        for gc3, enc, category in zip(dataset['gc3'], dataset['enc'],
                                      dataset['categories']):
            if not self.category_visibility.get(category, True):
                continue
            groups.setdefault(category, {'x': [], 'y': []})
            groups[category]['x'].append(gc3)
            groups[category]['y'].append(enc)

        for category in self.get_active_categories(dataset['categories']):
            if category not in groups:
                continue
            points = groups[category]
            ax.scatter(points['x'], points['y'],
                       c=self.get_category_color(category),
                       alpha=style['alpha'],
                       s=style['size'],
                       edgecolors=style['edgecolor'],
                       linewidths=style['linewidth'],
                       zorder=3,
                       label=f'{category} (n={len(points["x"])})')

    def draw_panel(self, ax, dataset, title_text, add_legend=True,
                   title_fontsize=None):
        """Draw one GC3-ENC panel for ``dataset`` on ``ax``.

        When ``add_legend`` is False the legend is omitted, so that a combined
        figure can place a single shared legend on the side instead.
        """
        style_curve = self.style_settings['curve']
        style_scatter = self.style_settings['scatter']

        # Theoretical ENC curve FIRST (drawn below the data points)
        if style_curve['visible']:
            gc3_theory, enc_theory = self.calculate_enc_curve()
            ax.plot(gc3_theory, enc_theory,
                    color=style_curve['color'],
                    linestyle=style_curve['linestyle'],
                    linewidth=style_curve['linewidth'],
                    zorder=1,
                    label='Theoretical ENC curve (no selection)')

        # Observed data points SECOND (on top of the curve)
        if style_scatter['visible']:
            if self.color_by_function.get():
                self.draw_scatter_by_function(ax, dataset)
            else:
                ax.scatter(dataset['gc3'], dataset['enc'],
                           c=style_scatter['color'],
                           alpha=style_scatter['alpha'],
                           s=style_scatter['size'],
                           edgecolors=style_scatter['edgecolor'],
                           linewidths=style_scatter['linewidth'],
                           zorder=2,
                           label=f'Observed data (n={len(dataset["gc3"])})')

        # Axis labels
        if self.style_settings['xlabel']['visible']:
            ax.set_xlabel('GC3 (GC content at third codon position)',
                          fontsize=self.style_settings['xlabel']['fontsize'],
                          fontweight=self.style_settings['xlabel']['fontweight'],
                          color=self.style_settings['xlabel']['color'])

        if self.style_settings['ylabel']['visible']:
            ax.set_ylabel('ENC (Effective Number of Codons)',
                          fontsize=self.style_settings['ylabel']['fontsize'],
                          fontweight=self.style_settings['ylabel']['fontweight'],
                          color=self.style_settings['ylabel']['color'])

        # Tick labels
        ax.tick_params(axis='x',
                       labelsize=self.style_settings['xticklabels']['fontsize'],
                       labelcolor=self.style_settings['xticklabels']['color'])
        ax.tick_params(axis='y',
                       labelsize=self.style_settings['yticklabels']['fontsize'],
                       labelcolor=self.style_settings['yticklabels']['color'])

        # Title
        if self.style_settings['title']['visible'] and title_text:
            size = title_fontsize if title_fontsize is not None \
                else self.style_settings['title']['fontsize']
            ax.set_title(title_text,
                         fontsize=size,
                         fontweight=self.style_settings['title']['fontweight'],
                         color=self.style_settings['title']['color'],
                         pad=10)

        # Axis ranges: GC3 0-1, ENC 20 to 61 plus 5% head-room
        ax.set_xlim(*AXIS_XLIM)
        ax.set_ylim(*AXIS_YLIM)
        ax.grid(True, alpha=0.3, linestyle=':')

        # In-panel legend
        if add_legend and self.style_settings['legend']['visible']:
            ax.legend(loc='best',
                      fontsize=self.style_settings['legend']['fontsize'],
                      labelcolor=self.style_settings['legend']['color'])

        # Gene symbols of the most outlying points (drawn last, on top)
        self.annotate_outliers(ax, dataset)

    def find_outlier_genes(self, dataset):
        """Return the genes lying farthest away from the rest of the data.

        The result is a list of ``(index, score, symbol, gc3, enc)`` tuples
        sorted by decreasing "outlyingness". ``score`` is expressed as a
        fraction of the visible axis span, so both ranking methods are
        directly comparable.
        """
        count = self._label_count()
        if count <= 0 or not dataset['gc3']:
            return []

        xs = np.asarray(dataset['gc3'], dtype=float)
        ys = np.asarray(dataset['enc'], dtype=float)

        x_span = max(AXIS_XLIM[1] - AXIS_XLIM[0], 1e-9)
        y_span = max(AXIS_YLIM[1] - AXIS_YLIM[0], 1e-9)

        if self.label_method.get() == LABEL_METHOD_CURVE:
            # Vertical distance to the theoretical ENC curve (Wright 1990)
            theory = 2 + xs + 29 / (xs ** 2 + (1 - xs) ** 2)
            scores = np.abs(ys - theory) / y_span
        else:
            # Distance to the center of the cloud, in normalised axes units
            nx = xs / x_span
            ny = (ys - AXIS_YLIM[0]) / y_span
            scores = np.hypot(nx - nx.mean(), ny - ny.mean())

        order = np.argsort(-scores)[:count]
        return [(int(i), float(scores[i]), extract_gene_symbol(dataset['genes'][i]),
                 float(xs[i]), float(ys[i])) for i in order]

    def annotate_outliers(self, ax, dataset):
        """Write the gene symbol next to the most outlying points."""
        if not self.label_outliers.get():
            return

        genes = self.find_outlier_genes(dataset)
        if not genes:
            return

        try:
            fontsize = max(4, int(float(self.label_fontsize.get())))
        except (TypeError, ValueError):
            fontsize = 7
        color = self.label_color

        # Offsets are given in points; the collision test needs pixels
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

    def _label_count(self):
        """Number of gene labels requested by the user."""
        try:
            count = int(float(self.label_count.get()))
        except (TypeError, ValueError):
            count = 5
        return max(0, count)

    @staticmethod
    def _boxes_overlap(a, b):
        """True when two (x0, y0, x1, y1) rectangles overlap."""
        return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])

    def choose_label_color(self, button):
        """Choose the color used for the gene labels."""
        color = colorchooser.askcolor(title="Choose gene label color",
                                      color=self.label_color)[1]
        if color:
            self.label_color = color
            self.update_color_button(button, color)
            self.redraw()

    def build_shared_legend_handles(self):
        """Build legend handles shared by every panel of the combined figure."""
        handles, labels = [], []

        if self.style_settings['curve']['visible']:
            handles.append(Line2D([0], [0],
                                  color=self.style_settings['curve']['color'],
                                  linestyle=self.style_settings['curve']['linestyle'],
                                  linewidth=self.style_settings['curve']['linewidth']))
            labels.append('Theoretical ENC curve (no selection)')

        if self.style_settings['scatter']['visible']:
            style = self.style_settings['scatter']
            if self.color_by_function.get():
                # Union of the visible categories across all datasets
                present = []
                for dataset in self.datasets:
                    for category in self.get_active_categories(dataset['categories']):
                        if self.category_visibility.get(category, True) and category not in present:
                            present.append(category)
                present = self.get_active_categories(present)
                for category in present:
                    handles.append(Patch(facecolor=self.get_category_color(category),
                                         edgecolor=style['edgecolor'],
                                         alpha=style['alpha']))
                    labels.append(category)
            else:
                handles.append(Patch(facecolor=style['color'],
                                     edgecolor=style['edgecolor'],
                                     alpha=style['alpha']))
                labels.append('Observed data')

        return handles, labels

    def redraw(self):
        """Re-render the current figure from the already loaded datasets."""
        if not self.datasets:
            return
        if self.input_mode.get() == 'folder':
            self.render_combined()
        else:
            self.render_single()

    def on_function_toggle(self):
        """Repaint the graph when the function coloring option changes."""
        self.redraw()

    def toggle_category(self, category, var):
        """Show/hide a functional category."""
        self.category_visibility[category] = var.get()

    def choose_category_color(self, category, button):
        """Choose the color of a functional category."""
        color = colorchooser.askcolor(title=f"Choose color: {category}")[1]
        if color:
            self.category_colors[category] = color
            self.update_color_button(button, color)

    def plot_graph(self):
        """Dispatch to single-file or folder (combined) plotting."""
        if self.input_mode.get() == 'folder':
            self.plot_combined()
        else:
            self.plot_single()

    def plot_single(self):
        """Load a single .out file and render it."""
        if not self.load_data():
            return
        self.render_single()

    def render_single(self):
        """Render the current single dataset (without re-reading the file)."""
        if not self.datasets:
            return

        # Rebuild the figure from scratch
        self.fig.clear()
        self.fig.set_size_inches(6, 6)
        self.ax = self.fig.add_subplot(111)

        dataset = self.datasets[0]
        self.draw_panel(self.ax, dataset,
                        f'{dataset["name"]} GC3-ENC Plot',
                        add_legend=True)

        self.fig.tight_layout()
        self.fit_canvas_to_figure()
        self.canvas.draw()

        self.status_bar.config(
            text=f"Graph plotted - {len(dataset['gc3'])} data points")

    def plot_combined(self):
        """Load a folder of .out files and render the combined grid figure."""
        if not self.load_folder_data():
            return
        self.render_combined()

    def render_combined(self):
        """Render the combined grid figure from the loaded datasets."""
        datasets = self.datasets
        if not datasets:
            return

        n = len(datasets)

        # Determine the grid layout
        try:
            user_cols = int(float(self.grid_columns.get()))
        except (TypeError, ValueError):
            user_cols = 0

        if user_cols and user_cols > 0:
            ncols = min(user_cols, n)
        else:
            ncols = min(4, n)
        ncols = max(1, ncols)
        nrows = math.ceil(n / ncols)

        # Figure size: ~4.0 in per panel, plus room on the right for the legend
        panel = 4.0
        side_width = 3.0
        fig_w = ncols * panel + side_width
        fig_h = max(nrows * panel, 5.0)

        self.fig.clear()
        self.fig.set_size_inches(fig_w, fig_h)

        # Reserve the right margin for the shared legend
        right = max(0.45, 1.0 - side_width / fig_w)
        gs = self.fig.add_gridspec(nrows, ncols,
                                   left=0.06, right=right,
                                   top=0.95, bottom=0.07,
                                   wspace=0.30, hspace=0.40)

        # Smaller title font when many panels are shown
        title_fontsize = max(7, self.style_settings['title']['fontsize'] - 3)

        for i, dataset in enumerate(datasets):
            r, c = divmod(i, ncols)
            ax = self.fig.add_subplot(gs[r, c])
            # Panel title is the file name; only the legend moves to the side
            self.draw_panel(ax, dataset, dataset['name'],
                            add_legend=False,
                            title_fontsize=title_fontsize)

        self.add_side_legend(right)

        self.fit_canvas_to_figure()
        self.canvas.draw()

        total = sum(len(d['gc3']) for d in datasets)
        self.status_bar.config(
            text=f"Combined figure: {n} file(s), {total} data points, {nrows}×{ncols} grid")

    def add_side_legend(self, right):
        """Place the unified legend on the side of the combined figure."""
        if not self.style_settings['legend']['visible']:
            return

        handles, labels = self.build_shared_legend_handles()
        if not handles:
            return

        x_side = min(0.985, right + 0.03)
        self.fig.legend(handles, labels,
                        loc='center left',
                        bbox_to_anchor=(x_side, 0.5),
                        fontsize=self.style_settings['legend']['fontsize'],
                        labelcolor=self.style_settings['legend']['color'],
                        frameon=False,
                        borderaxespad=0.0)

    def save_image(self):
        """Save the currently displayed image (combined or single)."""
        if not self.datasets:
            messagebox.showwarning("Warning", "Please plot the graph first!")
            return

        filename = filedialog.asksaveasfilename(
            title="Save Image",
            defaultextension=".pdf",
            filetypes=[
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("PNG files", "*.png"),
                ("EPS files", "*.eps"),
                ("All files", "*.*")
            ]
        )

        if not filename:
            return

        try:
            self.fig.savefig(filename, dpi=300, bbox_inches='tight')
        except Exception as e:
            messagebox.showerror("Error", f"Save failed: {str(e)}")
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
                saved_msg += f"\n\nAlso exported {count} individual image(s) to:\n{individual_dir}"
            except Exception as e:
                messagebox.showwarning(
                    "Partial success",
                    f"{saved_msg}\n\nBut exporting individual images failed:\n{e}")
                return

        messagebox.showinfo("Success", saved_msg)

    def export_individual_images(self, out_dir, extension):
        """Draw and save one GC3-ENC figure per dataset into ``out_dir``."""
        os.makedirs(out_dir, exist_ok=True)
        count = 0
        for dataset in self.datasets:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
            self.draw_panel(ax, dataset,
                            f'{dataset["name"]} GC3-ENC Plot',
                            add_legend=True)
            fig.tight_layout()
            safe_name = os.path.splitext(dataset['name'])[0]
            out_path = os.path.join(out_dir, f'{safe_name}{extension}')
            fig.savefig(out_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            count += 1
        return count

    def open_style_settings(self):
        """Open style settings window"""
        style_window = tk.Toplevel(self.root)
        style_window.title("Style Settings")
        style_window.geometry("500x650")
        
        # Create notebook for tabs
        notebook = ttk.Notebook(style_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Tab 1: Title and Labels
        tab_labels = ttk.Frame(notebook)
        notebook.add(tab_labels, text="Title & Labels")

        # Title settings
        ttk.Label(tab_labels, text="Title", font=('Arial', 10, 'bold')).pack(pady=(10, 5))

        title_frame = ttk.Frame(tab_labels)
        title_frame.pack(fill=tk.X, padx=10)

        show_title_var = tk.BooleanVar(value=self.style_settings['title']['visible'])
        ttk.Checkbutton(title_frame, text="Show title",
                       variable=show_title_var,
                       command=lambda: self.toggle_visibility('title', show_title_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(title_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        title_size = ttk.Scale(title_frame, from_=8, to=24, orient=tk.HORIZONTAL)
        title_size.set(self.style_settings['title']['fontsize'])
        title_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(title_frame, text="Color",
                  command=lambda: self.choose_color('title', title_color_btn)).grid(row=2, column=0, pady=5)
        title_color_btn = ttk.Button(title_frame, text="", width=15,
                                    command=lambda: self.choose_color('title', title_color_btn))
        title_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(title_color_btn, self.style_settings['title']['color'])

        # X-axis label settings
        ttk.Label(tab_labels, text="X-axis Label", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        xlabel_frame = ttk.Frame(tab_labels)
        xlabel_frame.pack(fill=tk.X, padx=10)

        show_xlabel_var = tk.BooleanVar(value=self.style_settings['xlabel']['visible'])
        ttk.Checkbutton(xlabel_frame, text="Show X-axis label",
                       variable=show_xlabel_var,
                       command=lambda: self.toggle_visibility('xlabel', show_xlabel_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(xlabel_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        xlabel_size = ttk.Scale(xlabel_frame, from_=8, to=20, orient=tk.HORIZONTAL)
        xlabel_size.set(self.style_settings['xlabel']['fontsize'])
        xlabel_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(xlabel_frame, text="Color",
                  command=lambda: self.choose_color('xlabel', xlabel_color_btn)).grid(row=2, column=0, pady=5)
        xlabel_color_btn = ttk.Button(xlabel_frame, text="", width=15,
                                    command=lambda: self.choose_color('xlabel', xlabel_color_btn))
        xlabel_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(xlabel_color_btn, self.style_settings['xlabel']['color'])

        # X-axis tick labels settings
        ttk.Label(tab_labels, text="X-axis Tick Labels", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        xtick_frame = ttk.Frame(tab_labels)
        xtick_frame.pack(fill=tk.X, padx=10)

        ttk.Label(xtick_frame, text="Font size:").grid(row=0, column=0, sticky=tk.W, pady=5)
        xtick_size = ttk.Scale(xtick_frame, from_=6, to=18, orient=tk.HORIZONTAL)
        xtick_size.set(self.style_settings['xticklabels']['fontsize'])
        xtick_size.grid(row=0, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(xtick_frame, text="Color",
                  command=lambda: self.choose_color('xticklabels', xtick_color_btn)).grid(row=1, column=0, pady=5)
        xtick_color_btn = ttk.Button(xtick_frame, text="", width=15,
                                    command=lambda: self.choose_color('xticklabels', xtick_color_btn))
        xtick_color_btn.grid(row=1, column=1, padx=5, pady=5)
        self.update_color_button(xtick_color_btn, self.style_settings['xticklabels']['color'])

        # Y-axis label settings
        ttk.Label(tab_labels, text="Y-axis Label", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        ylabel_frame = ttk.Frame(tab_labels)
        ylabel_frame.pack(fill=tk.X, padx=10)

        show_ylabel_var = tk.BooleanVar(value=self.style_settings['ylabel']['visible'])
        ttk.Checkbutton(ylabel_frame, text="Show Y-axis label",
                       variable=show_ylabel_var,
                       command=lambda: self.toggle_visibility('ylabel', show_ylabel_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(ylabel_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ylabel_size = ttk.Scale(ylabel_frame, from_=8, to=20, orient=tk.HORIZONTAL)
        ylabel_size.set(self.style_settings['ylabel']['fontsize'])
        ylabel_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(ylabel_frame, text="Color",
                  command=lambda: self.choose_color('ylabel', ylabel_color_btn)).grid(row=2, column=0, pady=5)
        ylabel_color_btn = ttk.Button(ylabel_frame, text="", width=15,
                                    command=lambda: self.choose_color('ylabel', ylabel_color_btn))
        ylabel_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(ylabel_color_btn, self.style_settings['ylabel']['color'])

        # Y-axis tick labels settings
        ttk.Label(tab_labels, text="Y-axis Tick Labels", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        ytick_frame = ttk.Frame(tab_labels)
        ytick_frame.pack(fill=tk.X, padx=10)

        ttk.Label(ytick_frame, text="Font size:").grid(row=0, column=0, sticky=tk.W, pady=5)
        ytick_size = ttk.Scale(ytick_frame, from_=6, to=18, orient=tk.HORIZONTAL)
        ytick_size.set(self.style_settings['yticklabels']['fontsize'])
        ytick_size.grid(row=0, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(ytick_frame, text="Color",
                  command=lambda: self.choose_color('yticklabels', ytick_color_btn)).grid(row=1, column=0, pady=5)
        ytick_color_btn = ttk.Button(ytick_frame, text="", width=15,
                                    command=lambda: self.choose_color('yticklabels', ytick_color_btn))
        ytick_color_btn.grid(row=1, column=1, padx=5, pady=5)
        self.update_color_button(ytick_color_btn, self.style_settings['yticklabels']['color'])
        
        # Tab 2: Curve
        tab_curve = ttk.Frame(notebook)
        notebook.add(tab_curve, text="Theoretical Curve")
        
        ttk.Label(tab_curve, text="Theoretical ENC Curve", font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        curve_frame = ttk.Frame(tab_curve)
        curve_frame.pack(fill=tk.X, padx=10)
        
        show_curve_var = tk.BooleanVar(value=self.style_settings['curve']['visible'])
        ttk.Checkbutton(curve_frame, text="Show theoretical curve",
                       variable=show_curve_var,
                       command=lambda: self.toggle_visibility('curve', show_curve_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        ttk.Label(curve_frame, text="Line width:").grid(row=1, column=0, sticky=tk.W, pady=5)
        curve_width = ttk.Scale(curve_frame, from_=0.5, to=5, orient=tk.HORIZONTAL)
        curve_width.set(self.style_settings['curve']['linewidth'])
        curve_width.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Label(curve_frame, text="Line style:").grid(row=2, column=0, sticky=tk.W, pady=5)
        curve_style = ttk.Combobox(curve_frame, values=['-', '--', ':', '-.'], state='readonly', width=18)
        curve_style.set(self.style_settings['curve']['linestyle'])
        curve_style.grid(row=2, column=1, padx=5, pady=5)
        
        ttk.Button(curve_frame, text="Color",
                  command=lambda: self.choose_color('curve', curve_color_btn)).grid(row=3, column=0, pady=5)
        curve_color_btn = ttk.Button(curve_frame, text="", width=15,
                                   command=lambda: self.choose_color('curve', curve_color_btn))
        curve_color_btn.grid(row=3, column=1, padx=5, pady=5)
        self.update_color_button(curve_color_btn, self.style_settings['curve']['color'])
        
        # Tab 3: Scatter
        tab_scatter = ttk.Frame(notebook)
        notebook.add(tab_scatter, text="Scatter Points")
        
        ttk.Label(tab_scatter, text="Scatter Points", font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        scatter_frame = ttk.Frame(tab_scatter)
        scatter_frame.pack(fill=tk.X, padx=10)
        
        show_scatter_var = tk.BooleanVar(value=self.style_settings['scatter']['visible'])
        ttk.Checkbutton(scatter_frame, text="Show scatter points",
                       variable=show_scatter_var,
                       command=lambda: self.toggle_visibility('scatter', show_scatter_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        ttk.Label(scatter_frame, text="Point size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        scatter_size = ttk.Scale(scatter_frame, from_=10, to=200, orient=tk.HORIZONTAL)
        scatter_size.set(self.style_settings['scatter']['size'])
        scatter_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Label(scatter_frame, text="Transparency:").grid(row=2, column=0, sticky=tk.W, pady=5)
        scatter_alpha = ttk.Scale(scatter_frame, from_=0.1, to=1.0, orient=tk.HORIZONTAL)
        scatter_alpha.set(self.style_settings['scatter']['alpha'])
        scatter_alpha.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Label(scatter_frame, text="Edge width:").grid(row=3, column=0, sticky=tk.W, pady=5)
        scatter_edge = ttk.Scale(scatter_frame, from_=0, to=3, orient=tk.HORIZONTAL)
        scatter_edge.set(self.style_settings['scatter']['linewidth'])
        scatter_edge.grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Button(scatter_frame, text="Fill color",
                  command=lambda: self.choose_color('scatter', scatter_fill_btn)).grid(row=4, column=0, pady=5)
        scatter_fill_btn = ttk.Button(scatter_frame, text="", width=15,
                                    command=lambda: self.choose_color('scatter', scatter_fill_btn))
        scatter_fill_btn.grid(row=4, column=1, padx=5, pady=5)
        self.update_color_button(scatter_fill_btn, self.style_settings['scatter']['color'])
        
        ttk.Button(scatter_frame, text="Edge color",
                  command=lambda: self.choose_color('scatter_edge', scatter_edge_btn)).grid(row=5, column=0, pady=5)
        scatter_edge_btn = ttk.Button(scatter_frame, text="", width=15,
                                    command=lambda: self.choose_color('scatter_edge', scatter_edge_btn))
        scatter_edge_btn.grid(row=5, column=1, padx=5, pady=5)
        self.update_color_button(scatter_edge_btn, self.style_settings['scatter']['edgecolor'])
        
        # Tab 4: Legend
        tab_legend = ttk.Frame(notebook)
        notebook.add(tab_legend, text="Legend")
        
        ttk.Label(tab_legend, text="Legend", font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        legend_frame = ttk.Frame(tab_legend)
        legend_frame.pack(fill=tk.X, padx=10)
        
        show_legend_var = tk.BooleanVar(value=self.style_settings['legend']['visible'])
        ttk.Checkbutton(legend_frame, text="Show legend",
                       variable=show_legend_var,
                       command=lambda: self.toggle_visibility('legend', show_legend_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        ttk.Label(legend_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        legend_size = ttk.Scale(legend_frame, from_=6, to=16, orient=tk.HORIZONTAL)
        legend_size.set(self.style_settings['legend']['fontsize'])
        legend_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Button(legend_frame, text="Text color",
                  command=lambda: self.choose_color('legend', legend_color_btn)).grid(row=2, column=0, pady=5)
        legend_color_btn = ttk.Button(legend_frame, text="", width=15,
                                    command=lambda: self.choose_color('legend', legend_color_btn))
        legend_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(legend_color_btn, self.style_settings['legend']['color'])
        
        # Tab 5: Gene Function
        tab_function = ttk.Frame(notebook)
        notebook.add(tab_function, text="Gene Function")
        
        ttk.Label(tab_function, text="Color Points by Gene Function",
                 font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        func_top = ttk.Frame(tab_function)
        func_top.pack(fill=tk.X, padx=10)
        
        ttk.Checkbutton(func_top, text="Color scatter points by gene function",
                       variable=self.color_by_function,
                       command=self.on_function_toggle).pack(anchor=tk.W, pady=5)
        ttk.Label(func_top, text="Uncheck 'Show' to hide a category from the graph.",
                 foreground='gray').pack(anchor=tk.W, pady=(0, 5))
        
        # Scrollable category list
        list_container = ttk.Frame(tab_function)
        list_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        list_canvas = tk.Canvas(list_container, highlightthickness=0)
        list_scrollbar = ttk.Scrollbar(list_container, orient=tk.VERTICAL,
                                       command=list_canvas.yview)
        list_inner = ttk.Frame(list_canvas)
        list_inner.bind(
            '<Configure>',
            lambda event: list_canvas.configure(scrollregion=list_canvas.bbox('all')))
        list_canvas.create_window((0, 0), window=list_inner, anchor='nw')
        list_canvas.configure(yscrollcommand=list_scrollbar.set)
        list_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        list_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.function_widgets = {}
        for category in CATEGORY_ORDER:
            row_frame = ttk.Frame(list_inner)
            row_frame.pack(fill=tk.X, pady=2)
            
            vis_var = tk.BooleanVar(value=self.category_visibility.get(category, True))
            ttk.Checkbutton(row_frame, text="Show", variable=vis_var,
                            command=lambda c=category, v=vis_var: self.toggle_category(c, v)).pack(side=tk.LEFT)
            
            text = category if len(category) <= 30 else category[:28] + '...'
            ttk.Label(row_frame, text=text, width=30, anchor=tk.W).pack(side=tk.LEFT, padx=(5, 5))
            
            color_btn = ttk.Button(row_frame, text="", width=10)
            color_btn.configure(command=lambda c=category, b=color_btn: self.choose_category_color(c, b))
            color_btn.pack(side=tk.LEFT)
            self.update_color_button(color_btn, self.get_category_color(category))
            
            self.function_widgets[category] = {'vis_var': vis_var, 'color_btn': color_btn}
        
        # Apply and Reset buttons
        button_frame = ttk.Frame(style_window)
        button_frame.pack(fill=tk.X, padx=10, pady=10)

        apply_btn = ttk.Button(button_frame, text="Apply",
                              command=lambda: self.apply_settings({
                                  'title_size': title_size.get(),
                                  'xlabel_size': xlabel_size.get(),
                                  'ylabel_size': ylabel_size.get(),
                                  'xtick_size': xtick_size.get(),
                                  'ytick_size': ytick_size.get(),
                                  'curve_width': curve_width.get(),
                                  'curve_style': curve_style.get(),
                                  'scatter_size': scatter_size.get(),
                                  'scatter_alpha': scatter_alpha.get(),
                                  'scatter_edge': scatter_edge.get(),
                                  'legend_size': legend_size.get()
                              }, style_window))
        apply_btn.pack(side=tk.LEFT, padx=5)
        
        reset_btn = ttk.Button(button_frame, text="Reset to Default",
                              command=lambda: self.reset_settings(style_window))
        reset_btn.pack(side=tk.RIGHT, padx=5)
        
        # Store references
        self.style_refs = {
            'title_color_btn': title_color_btn,
            'xlabel_color_btn': xlabel_color_btn,
            'ylabel_color_btn': ylabel_color_btn,
            'curve_color_btn': curve_color_btn,
            'curve_style': curve_style,
            'scatter_fill_btn': scatter_fill_btn,
            'scatter_edge_btn': scatter_edge_btn,
            'legend_color_btn': legend_color_btn,
            'function_widgets': self.function_widgets
        }
    
    def choose_color(self, element, button):
        """Choose color for an element"""
        color = colorchooser.askcolor(title="Choose Color")[1]
        if color:
            if element == 'curve':
                self.style_settings['curve']['color'] = color
            elif element == 'scatter':
                self.style_settings['scatter']['color'] = color
            elif element == 'scatter_edge':
                self.style_settings['scatter']['edgecolor'] = color
            elif element == 'title':
                self.style_settings['title']['color'] = color
            elif element == 'xlabel':
                self.style_settings['xlabel']['color'] = color
            elif element == 'ylabel':
                self.style_settings['ylabel']['color'] = color
            elif element == 'xticklabels':
                self.style_settings['xticklabels']['color'] = color
            elif element == 'yticklabels':
                self.style_settings['yticklabels']['color'] = color
            elif element == 'legend':
                self.style_settings['legend']['color'] = color

            self.update_color_button(button, color)
    
    def update_color_button(self, button, color):
        """Update color button appearance"""
        button.configure(text=color)
    
    def toggle_visibility(self, element, var):
        """Toggle visibility of element"""
        self.style_settings[element]['visible'] = var.get()
        print(f"Toggle {element} visibility: {self.style_settings[element]['visible']}")
    
    def apply_settings(self, settings, window):
        """Apply style settings"""
        self.style_settings['title']['fontsize'] = settings['title_size']
        self.style_settings['xlabel']['fontsize'] = settings['xlabel_size']
        self.style_settings['ylabel']['fontsize'] = settings['ylabel_size']
        self.style_settings['xticklabels']['fontsize'] = settings['xtick_size']
        self.style_settings['yticklabels']['fontsize'] = settings['ytick_size']
        self.style_settings['curve']['linewidth'] = settings['curve_width']
        self.style_settings['curve']['linestyle'] = settings['curve_style']
        self.style_settings['scatter']['size'] = settings['scatter_size']
        self.style_settings['scatter']['alpha'] = settings['scatter_alpha']
        self.style_settings['scatter']['linewidth'] = settings['scatter_edge']
        self.style_settings['legend']['fontsize'] = settings['legend_size']

        # Re-render if data exists
        self.redraw()

        messagebox.showinfo("Success", "Style settings applied.")
        window.destroy()
    
    def reset_settings(self, window):
        """Reset to default settings"""
        self.style_settings = {
            'title': {
                'fontsize': 14,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xlabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'ylabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'yticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'curve': {
                'color': 'red',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'scatter': {
                'color': 'blue',
                'edgecolor': 'darkblue',
                'size': 50,
                'alpha': 0.6,
                'linewidth': 0.5,
                'visible': True
            },
            'legend': {
                'fontsize': 10,
                'color': 'black',
                'visible': True
            }
        }
        
        # Reset gene function coloring to defaults
        self.color_by_function.set(False)
        self.category_colors = dict(DEFAULT_CATEGORY_COLORS)
        self.category_visibility = {cat: True for cat in CATEGORY_ORDER}

        # Reset the outlier gene labels to defaults
        self.label_outliers.set(False)
        self.label_method.set(LABEL_METHOD_CENTER)
        self.label_count.set(5)
        self.label_fontsize.set(7)
        self.label_color = 'darkgreen'
        self.update_color_button(self.label_color_btn, self.label_color)

        # Repaint with the restored defaults
        self.redraw()

        messagebox.showinfo("Success", "Settings reset to default.")
        window.destroy()


def main():
    root = tk.Tk()
    app = GC3ENCPlotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
