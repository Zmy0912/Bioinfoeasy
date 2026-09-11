import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os

# Set Chinese fonts
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def plot_codon_usage_table_and_chart(excel_file_path, output_path='codon_usage_chart.png', colormap_name='tab20c', show_codon_labels=True, figsize=(18, 13), legend_block_height=0.12, legend_block_width=0.6, show_legend_blocks=True, codon_label_size=10, xlabel_size=14, ylabel_size=14, tick_label_size=12, ytick_label_size=12, legend_text_size=8, border_width=1.0, legend_bottom_offset=0.22):
    """
    Draw a combined plot of codon usage frequency: stacked bar chart below

    Parameters:
        excel_file_path: Excel file path
        output_path: Output image path
        colormap_name: Color scheme name, options: 'tab20c', 'tab20b', 'tab20', 'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'viridis', 'plasma', 'inferno'
        show_codon_labels: Whether to display codon names in the bar chart, True to display, False to hide
        figsize: Canvas size, format as (width, height), default is (18, 13)
        legend_block_height: Height of codon stacked blocks, default is 0.12
        legend_block_width: Width of codon stacked blocks, default is 0.6
        show_legend_blocks: Whether to display codon stacked block legend, True to display, False to hide
        codon_label_size: Font size of codon names in bar chart, default is 10
        xlabel_size: Font size of x-axis label, default is 14
        ylabel_size: Font size of y-axis label, default is 14
        tick_label_size: Font size of x-axis tick labels, default is 12
        ytick_label_size: Font size of y-axis tick labels, default is 12
        legend_text_size: Font size of stacked block legend text, default is 8
        border_width: Chart border and tick line thickness, default is 1.0
        legend_bottom_offset: Distance between stacked blocks and X-axis, default is 0.22
    """

    # Read Excel file, try different engines
    try:
        df = pd.read_excel(excel_file_path, engine='openpyxl')
    except:
        try:
            df = pd.read_excel(excel_file_path, engine='xlrd')
        except:
            df = pd.read_excel(excel_file_path)

    # Get column names (adapt to different column name formats)
    print("Excel column names:", df.columns.tolist())

    # Check if there is data
    if len(df) == 0:
        raise ValueError("Excel file is empty, no data")

    # Try to identify codon, frequency, RSCU columns
    codon_col = None
    freq_col = None
    rscu_col = None
    aa_col = None

    for col in df.columns:
        col_lower = str(col).lower()
        col_str = str(col)
        if 'codon' in col_lower or '密码子' in col_str:
            codon_col = col
        elif 'freq' in col_lower or '频数' in col_str or 'count' in col_lower or '数量' in col_str:
            freq_col = col
        elif 'rscu' in col_lower or 'rscu' in col_str:
            rscu_col = col
        elif 'amino' in col_lower or 'aa' in col_lower or '氨基酸' in col_str:
            aa_col = col

    if codon_col is None:
        codon_col = df.columns[0]
    if freq_col is None:
        freq_col = df.columns[1]
    if rscu_col is None:
        rscu_col = df.columns[2]

    # If no amino acid column, derive from codon
    if aa_col is None:
        # Ensure codon column data is string type
        df[codon_col] = df[codon_col].astype(str)
        df['AminoAcid'] = df[codon_col].apply(codon_to_aa)
        aa_col = 'AminoAcid'

    # Rename columns for use
    rename_dict = {
        codon_col: 'Codon',
        freq_col: 'Frequency',
        rscu_col: 'RSCU'
    }
    if aa_col is not None:
        rename_dict[aa_col] = 'AminoAcid'

    df = df.rename(columns=rename_dict)

    # Ensure required columns exist
    if 'Codon' not in df.columns:
        raise ValueError(f"Codon column not found, available columns: {df.columns.tolist()}")
    if 'RSCU' not in df.columns and 'Frequency' not in df.columns:
        raise ValueError(f"RSCU or Frequency column not found, available columns: {df.columns.tolist()}")

    # Group by amino acid
    amino_order = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Tyr',
                   'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly']

    # Ensure data is sorted correctly by amino acid
    df['AminoAcid'] = pd.Categorical(df['AminoAcid'], categories=amino_order, ordered=True)
    df = df.sort_values('AminoAcid')

    # Create figure
    fig = plt.figure(figsize=figsize)

    # ===== Stacked Bar Chart =====
    ax_chart = fig.add_axes([0.08, 0.22, 0.85, 0.68])

    # Prepare stacked data for each amino acid
    bottom_values = {aa: 0 for aa in amino_order}
    codon_colors = {}
    position_colors = {}  # Store fixed colors for each position

    # Assign fixed colors for each position (positions 1, 2, 3... use different colors)
    # Use user-specified color map
    try:
        colormap = plt.get_cmap(colormap_name)
    except:
        print(f"Warning: Color scheme '{colormap_name}' not found, using default 'tab20c'")
        colormap = plt.get_cmap('tab20c')

    # Group by amino acid, assign colors for each amino acid's codons by position
    for aa in amino_order:
        aa_data = df[df['AminoAcid'] == aa].reset_index(drop=True)
        for idx, row in aa_data.iterrows():
            codon = row['Codon']
            position = idx  # Position of codon in amino acid
            if position not in position_colors:
                position_colors[position] = colormap(position % 20)
            codon_colors[codon] = position_colors[position]

    # Draw stacked bar chart
    x_positions = np.arange(len(amino_order))
    bar_width = 0.8

    # Store codon information for each amino acid, for labeling
    aa_codons_info = {}

    for aa in amino_order:
        aa_data = df[df['AminoAcid'] == aa]
        if len(aa_data) == 0:
            continue

        aa_codons = aa_data['Codon'].tolist()
        aa_rscu = aa_data['RSCU'].tolist()

        aa_codons_info[aa] = {
            'codons': aa_codons,
            'rscu': aa_rscu
        }

        for codon, rscu in zip(aa_codons, aa_rscu):
            ax_chart.bar(aa, rscu, bar_width,
                        bottom=bottom_values[aa],
                        color=codon_colors[codon],
                        edgecolor='white',
                        linewidth=0.5)
            bottom_values[aa] += rscu

    # Set axis labels
    ax_chart.set_xlabel('Amino Acid', fontsize=xlabel_size)
    ax_chart.set_ylabel('RSCU', fontsize=ylabel_size)
    ax_chart.set_xticks(x_positions)
    ax_chart.set_xticklabels(amino_order, fontsize=tick_label_size)
    # Set Y-axis tick label size
    ax_chart.tick_params(axis='y', labelsize=ytick_label_size)

    # Set border and tick line thickness
    for spine in ax_chart.spines.values():
        spine.set_linewidth(border_width)
    ax_chart.tick_params(width=border_width)

    # Calculate RSCU sum for each amino acid (total height after stacking)
    aa_rscu_sums = df.groupby('AminoAcid')['RSCU'].sum()
    max_total_rscu = aa_rscu_sums.max()

    # Set Y-axis upper limit based on maximum total height after stacking
    if max_total_rscu > 3.0:
        # If maximum value after stacking is greater than 3.0, set a slightly larger value
        y_upper_limit = max_total_rscu * 1.05  # Add 5% margin
    else:
        # If maximum value is less than or equal to 3.0, default to 3.0
        y_upper_limit = 3.0
    ax_chart.set_ylim(0, y_upper_limit)

    # Label corresponding codons below each amino acid bar
    if show_codon_labels:
        for idx, aa in enumerate(amino_order):
            if aa in aa_codons_info:
                codons = aa_codons_info[aa]['codons']
                rscu_values = aa_codons_info[aa]['rscu']

                # Calculate center position of each color block
                bottom = 0
                for codon, rscu in zip(codons, rscu_values):
                    # Center position of color block
                    center_y = bottom + rscu / 2

                    # Label codon name
                    ax_chart.text(idx, center_y, codon,
                                ha='center', va='center',
                                fontsize=codon_label_size,
                                color='black' if codon_colors[codon] > (0.5, 0.5, 0.5) else 'white')

                    bottom += rscu

    # Add grid lines
    ax_chart.grid(True, axis='y', alpha=0.3, linestyle='--')
    ax_chart.set_axisbelow(True)

    # Add title
    fig.text(0.5, 0.02, '密码子使用频率分析 / Codon Usage Analysis',
            ha='center', fontsize=16, fontweight='bold')

    # ===== Add codon stacked legend below X-axis =====
    if show_legend_blocks:
        ax_legend = fig.add_axes([0.08, 0.06, 0.85, 0.14], sharex=ax_chart)
        ax_legend.axis('off')
        ax_legend.set_xlim(ax_chart.get_xlim())

        # Set fixed legend height range
        ax_legend.set_ylim(0, 1.0)

        # Create vertical stacked block legend for each amino acid
        block_width = legend_block_width  # Use adjustable width
        block_height = legend_block_height  # Use adjustable height
        y_start = 0.02  # Start from bottom

        for idx, aa in enumerate(amino_order):
            if aa in aa_codons_info:
                codons = aa_codons_info[aa]['codons']

                x_center = idx
                x_left = x_center - block_width / 2

                for j, codon in enumerate(codons):
                    y_bottom = y_start + j * block_height
                    color = codon_colors[codon]

                    # Draw rectangle (vertical stack, same height)
                    rect = plt.Rectangle((x_left, y_bottom), block_width, block_height,
                                       facecolor=color, edgecolor='white', linewidth=0.8)
                    ax_legend.add_patch(rect)

                    # Add codon name
                    ax_legend.text(x_center, y_bottom + block_height / 2, codon,
                                 ha='center', va='center',
                                 fontsize=legend_text_size,
                                 color='black' if color > (0.5, 0.5, 0.5) else 'white')

            # Adjust bar chart X-axis label position, move up to avoid overlap with legend
            ax_chart.xaxis.set_label_coords(0.5, -legend_bottom_offset)
    else:
        # When not showing stacked blocks, adjust bar chart position
        ax_chart.xaxis.set_label_coords(0.5, -0.05)

    # Adjust layout
    fig.subplots_adjust(left=0.08, right=0.95, top=0.98, bottom=0.08)
    return fig, df, output_path

def codon_to_aa(codon):
    """Codon to amino acid mapping (using RNA base U)"""
    genetic_code = {
        # Phe
        'UUU': 'Phe', 'UUC': 'Phe',
        # Leu
        'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        # Ile
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        # Met
        'AUG': 'Met',
        # Val
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        # Ser
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
        # Pro
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        # Thr
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        # Ala
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        # Tyr
        'UAU': 'Tyr', 'UAC': 'Tyr',
        # Stop
        'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
        # His
        'CAU': 'His', 'CAC': 'His',
        # Gln
        'CAA': 'Gln', 'CAG': 'Gln',
        # Asn
        'AAU': 'Asn', 'AAC': 'Asn',
        # Lys
        'AAA': 'Lys', 'AAG': 'Lys',
        # Asp
        'GAU': 'Asp', 'GAC': 'Asp',
        # Glu
        'GAA': 'Glu', 'GAG': 'Glu',
        # Cys
        'UGU': 'Cys', 'UGC': 'Cys',
        # Trp
        'UGG': 'Trp',
        # Arg
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
        # Gly
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }
    # Also support DNA base T (convert to U)
    codon = codon.upper().replace('T', 'U')
    return genetic_code.get(codon, 'Unknown')


class CodonUsageApp:
    """Interactive codon usage frequency analysis application"""

    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Codon Usage Analysis Tool / 密码子使用频率分析工具")
        self.root.geometry("1200x800")
        self.root.minsize(1000, 700)

        self.excel_file_path = None
        self.fig = None
        self.df = None
        self.canvas = None
        self.lbl_file = None
        self.lbl_hint = None
        self.colormap_var = tk.StringVar(value='tab20c')
        self.show_codon_labels_var = tk.BooleanVar(value=True)
        self.show_legend_blocks_var = tk.BooleanVar(value=True)
        self.fig_width_var = tk.StringVar(value='18')
        self.fig_height_var = tk.StringVar(value='13')
        self.legend_block_height_var = tk.StringVar(value='0.12')
        self.legend_block_width_var = tk.StringVar(value='0.6')

        # Text size settings (bold feature removed)
        self.codon_label_size_var = tk.StringVar(value='10')
        self.xlabel_size_var = tk.StringVar(value='14')
        self.ylabel_size_var = tk.StringVar(value='14')
        self.tick_label_size_var = tk.StringVar(value='12')
        self.ytick_label_size_var = tk.StringVar(value='12')
        self.legend_text_size_var = tk.StringVar(value='8')
        self.border_width_var = tk.StringVar(value='1.0')
        self.legend_bottom_offset_var = tk.StringVar(value='0.22')  # Distance between stacked blocks and X-axis

    def create_widgets(self):
        """Create interface components"""
        # Top toolbar - First row
        toolbar_row1 = tk.Frame(self.root, bd=1, relief=tk.RAISED)
        toolbar_row1.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        # Select Excel file button
        btn_load = tk.Button(toolbar_row1, text="Select Excel File", command=self.load_excel_file,
                            bg='#4CAF50', fg='white', font=('Arial', 10, 'bold'),
                            padx=15, pady=5)
        btn_load.pack(side=tk.LEFT, padx=5)

        # File path display
        if self.lbl_file is None:
            self.lbl_file = tk.Label(toolbar_row1, text="No file selected", fg='gray')
        self.lbl_file.pack(side=tk.LEFT, padx=10)

        # Separator
        tk.Frame(toolbar_row1, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Save buttons
        btn_pdf = tk.Button(toolbar_row1, text="Save as PDF", command=self.save_as_pdf,
                           bg='#2196F3', fg='white', font=('Arial', 10, 'bold'),
                           padx=15, pady=5)
        btn_pdf.pack(side=tk.LEFT, padx=5)

        btn_svg = tk.Button(toolbar_row1, text="Save as SVG", command=self.save_as_svg,
                           bg='#FF9800', fg='white', font=('Arial', 10, 'bold'),
                           padx=15, pady=5)
        btn_svg.pack(side=tk.LEFT, padx=5)

        btn_png = tk.Button(toolbar_row1, text="Save as PNG", command=self.save_as_png,
                           bg='#9C27B0', fg='white', font=('Arial', 10, 'bold'),
                           padx=15, pady=5)
        btn_png.pack(side=tk.LEFT, padx=5)

        # Top toolbar - Second row
        toolbar_row2 = tk.Frame(self.root, bd=1, relief=tk.RAISED)
        toolbar_row2.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        # Color scheme selection
        tk.Label(toolbar_row2, text="Color Scheme:", font=('Arial', 10)).pack(side=tk.LEFT, padx=5)
        colormap_options = ['tab20c', 'tab20b', 'tab20', 'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'viridis', 'plasma', 'inferno']
        self.colormap_combo = ttk.Combobox(toolbar_row2, textvariable=self.colormap_var,
                                           values=colormap_options, state='readonly', width=10)
        self.colormap_combo.pack(side=tk.LEFT, padx=5)

        # Codon label display option
        self.chk_codon_labels = tk.Checkbutton(toolbar_row2, text="Show Codon Names",
                                               variable=self.show_codon_labels_var,
                                               font=('Arial', 10),
                                               command=self.on_codon_label_toggle)
        self.chk_codon_labels.pack(side=tk.LEFT, padx=10)

        # Stacked block display option
        self.chk_legend_blocks = tk.Checkbutton(toolbar_row2, text="Show Stacked Blocks",
                                               variable=self.show_legend_blocks_var,
                                               font=('Arial', 10),
                                               command=self.on_legend_blocks_toggle)
        self.chk_legend_blocks.pack(side=tk.LEFT, padx=10)

        # Separator
        tk.Frame(toolbar_row2, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Canvas size adjustment
        tk.Label(toolbar_row2, text="Canvas Size:", font=('Arial', 10)).pack(side=tk.LEFT, padx=2)
        tk.Label(toolbar_row2, text="Width:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_width = tk.Entry(toolbar_row2, textvariable=self.fig_width_var, width=5)
        self.entry_width.pack(side=tk.LEFT, padx=2)
        tk.Label(toolbar_row2, text="Height:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_height = tk.Entry(toolbar_row2, textvariable=self.fig_height_var, width=5)
        self.entry_height.pack(side=tk.LEFT, padx=2)

        # Codon stacked block height adjustment
        tk.Frame(toolbar_row2, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        tk.Label(toolbar_row2, text="Block Height:", font=('Arial', 10)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_height = tk.Entry(toolbar_row2, textvariable=self.legend_block_height_var, width=6)
        self.entry_legend_height.pack(side=tk.LEFT, padx=2)

        # Codon stacked block width adjustment
        tk.Label(toolbar_row2, text="Width:", font=('Arial', 10)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_width = tk.Entry(toolbar_row2, textvariable=self.legend_block_width_var, width=6)
        self.entry_legend_width.pack(side=tk.LEFT, padx=2)

        # Separator
        tk.Frame(toolbar_row2, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Redraw button
        btn_redraw = tk.Button(toolbar_row2, text="Redraw", command=self.redraw_chart,
                               bg='#607D8B', fg='white', font=('Arial', 10, 'bold'),
                               padx=15, pady=5)
        btn_redraw.pack(side=tk.LEFT, padx=5)

        # Top toolbar - Third row (text settings)
        toolbar_row3 = tk.Frame(self.root, bd=1, relief=tk.RAISED)
        toolbar_row3.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        # Codon name settings in bar chart
        tk.Label(toolbar_row3, text="Bar Text:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="Size:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_codon_label_size = tk.Entry(toolbar_row3, textvariable=self.codon_label_size_var, width=4)
        self.entry_codon_label_size.pack(side=tk.LEFT, padx=2)

        # Separator
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Axis label settings
        tk.Label(toolbar_row3, text="Axis:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="X Size:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_xlabel_size = tk.Entry(toolbar_row3, textvariable=self.xlabel_size_var, width=4)
        self.entry_xlabel_size.pack(side=tk.LEFT, padx=2)

        tk.Label(toolbar_row3, text="Y Size:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_ylabel_size = tk.Entry(toolbar_row3, textvariable=self.ylabel_size_var, width=4)
        self.entry_ylabel_size.pack(side=tk.LEFT, padx=2)

        # Separator
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Tick label settings
        tk.Label(toolbar_row3, text="Ticks:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="X Size:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_tick_label_size = tk.Entry(toolbar_row3, textvariable=self.tick_label_size_var, width=4)
        self.entry_tick_label_size.pack(side=tk.LEFT, padx=2)

        tk.Label(toolbar_row3, text="Y Size:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_ytick_label_size = tk.Entry(toolbar_row3, textvariable=self.ytick_label_size_var, width=4)
        self.entry_ytick_label_size.pack(side=tk.LEFT, padx=2)

        # Separator
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Border thickness settings
        tk.Label(toolbar_row3, text="Border:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="Width:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_border_width = tk.Entry(toolbar_row3, textvariable=self.border_width_var, width=5)
        self.entry_border_width.pack(side=tk.LEFT, padx=2)

        # Separator
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Stacked block position settings
        tk.Label(toolbar_row3, text="Block Pos:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="From X:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_bottom_offset = tk.Entry(toolbar_row3, textvariable=self.legend_bottom_offset_var, width=5)
        self.entry_legend_bottom_offset.pack(side=tk.LEFT, padx=2)

        # Separator
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # Stacked block text settings
        tk.Label(toolbar_row3, text="Block Text:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="Size:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_text_size = tk.Entry(toolbar_row3, textvariable=self.legend_text_size_var, width=4)
        self.entry_legend_text_size.pack(side=tk.LEFT, padx=2)

        # Main plotting area container (with scrollbars)
        self.plot_container = tk.Frame(self.root)
        self.plot_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create vertical scrollbar
        v_scrollbar = tk.Scrollbar(self.plot_container, orient=tk.VERTICAL, command=self._on_vscroll)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Create horizontal scrollbar
        h_scrollbar = tk.Scrollbar(self.plot_container, orient=tk.HORIZONTAL, command=self._on_hscroll)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)

        # Create scrollable canvas container
        self.scroll_canvas = tk.Canvas(self.plot_container, bg='white')
        self.scroll_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Configure scrollable canvas scroll region
        self.scroll_canvas.config(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)

        # Create canvas container (fixed size, for placing matplotlib canvas)
        self.canvas_frame = tk.Frame(self.scroll_canvas, bg='white')
        self.canvas_window = self.scroll_canvas.create_window((0, 0), window=self.canvas_frame, anchor='nw')

        # Bind canvas size change event
        self.canvas_frame.bind('<Configure>', self._on_frame_configure)

        # Hint message
        if self.lbl_hint is None:
            self.lbl_hint = tk.Label(self.root,
                                   text="Please select an Excel file containing codon data, the program will automatically generate a table and stacked bar chart",
                                   fg='blue', font=('Arial', 10))
        self.lbl_hint.pack(side=tk.BOTTOM, pady=5)

    def load_excel_file(self):
        """Load Excel file"""
        file_path = filedialog.askopenfilename(
            title="Select Excel File",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )

        if file_path:
            self.excel_file_path = file_path
            self.lbl_file.config(text=os.path.basename(file_path), fg='black')
            self.plot_chart()

    def plot_chart(self):
        """Draw chart"""
        try:
            # Clear old chart
            if self.canvas:
                self.canvas.get_tk_widget().destroy()
            if hasattr(self, 'toolbar'):
                self.toolbar.destroy()

            # Generate chart
            colormap_name = self.colormap_var.get()
            show_codon_labels = self.show_codon_labels_var.get()

            # Get canvas size
            try:
                fig_width = float(self.fig_width_var.get())
                fig_height = float(self.fig_height_var.get())
                if fig_width < 10 or fig_height < 8:
                    raise ValueError("Canvas size too small")
                if fig_width > 30 or fig_height > 30:
                    raise ValueError("Canvas size too large")
                figsize = (fig_width, fig_height)
            except ValueError as e:
                messagebox.showwarning("Warning", f"Canvas size format error:\n{str(e)}\n\nUsing default size (18, 13)")
                figsize = (18, 13)
                self.fig_width_var.set('18')
                self.fig_height_var.set('13')

            # Get codon stacked block height
            try:
                legend_block_height = float(self.legend_block_height_var.get())
                if legend_block_height < 0.05 or legend_block_height > 0.3:
                    raise ValueError("Block height should be between 0.05 and 0.3")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Block height format error:\n{str(e)}\n\nUsing default value 0.12")
                legend_block_height = 0.12
                self.legend_block_height_var.set('0.12')

            # Get codon stacked block width
            try:
                legend_block_width = float(self.legend_block_width_var.get())
                if legend_block_width < 0.1 or legend_block_width > 1.0:
                    raise ValueError("Block width should be between 0.1 and 1.0")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Block width format error:\n{str(e)}\n\nUsing default value 0.6")
                legend_block_width = 0.6
                self.legend_block_width_var.set('0.6')

            # Get text settings
            try:
                codon_label_size = float(self.codon_label_size_var.get())
                if codon_label_size < 6 or codon_label_size > 24:
                    raise ValueError("Bar chart text size should be between 6 and 24")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Bar chart text size format error:\n{str(e)}\n\nUsing default value 10")
                codon_label_size = 10
                self.codon_label_size_var.set('10')

            try:
                xlabel_size = float(self.xlabel_size_var.get())
                ylabel_size = float(self.ylabel_size_var.get())
                if xlabel_size < 8 or xlabel_size > 30 or ylabel_size < 8 or ylabel_size > 30:
                    raise ValueError("Axis label size should be between 8 and 30")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Axis label size format error:\n{str(e)}\n\nUsing default value 14")
                xlabel_size = 14
                ylabel_size = 14
                self.xlabel_size_var.set('14')
                self.ylabel_size_var.set('14')

            try:
                tick_label_size = float(self.tick_label_size_var.get())
                ytick_label_size = float(self.ytick_label_size_var.get())
                if tick_label_size < 6 or tick_label_size > 20 or ytick_label_size < 6 or ytick_label_size > 20:
                    raise ValueError("Tick label size should be between 6 and 20")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Tick label size format error:\n{str(e)}\n\nUsing default value 12")
                tick_label_size = 12
                ytick_label_size = 12
                self.tick_label_size_var.set('12')
                self.ytick_label_size_var.set('12')

            try:
                border_width = float(self.border_width_var.get())
                if border_width < 0.5 or border_width > 5.0:
                    raise ValueError("Border thickness should be between 0.5 and 5.0")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Border thickness format error:\n{str(e)}\n\nUsing default value 1.0")
                border_width = 1.0
                self.border_width_var.set('1.0')

            try:
                legend_text_size = float(self.legend_text_size_var.get())
                if legend_text_size < 4 or legend_text_size > 16:
                    raise ValueError("Block text size should be between 4 and 16")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Block text size format error:\n{str(e)}\n\nUsing default value 8")
                legend_text_size = 8
                self.legend_text_size_var.set('8')

            # Get distance between stacked blocks and X-axis
            try:
                legend_bottom_offset = float(self.legend_bottom_offset_var.get())
                if legend_bottom_offset < 0.05 or legend_bottom_offset > 0.5:
                    raise ValueError("Block distance from X-axis should be between 0.05 and 0.5")
            except ValueError as e:
                messagebox.showwarning("Warning", f"Block distance from X-axis format error:\n{str(e)}\n\nUsing default value 0.22")
                legend_bottom_offset = 0.22
                self.legend_bottom_offset_var.set('0.22')

            self.fig, self.df, _ = plot_codon_usage_table_and_chart(self.excel_file_path, colormap_name=colormap_name, show_codon_labels=show_codon_labels, figsize=figsize, legend_block_height=legend_block_height, legend_block_width=legend_block_width, show_legend_blocks=self.show_legend_blocks_var.get(), codon_label_size=codon_label_size, xlabel_size=xlabel_size, ylabel_size=ylabel_size, tick_label_size=tick_label_size, ytick_label_size=ytick_label_size, legend_text_size=legend_text_size, border_width=border_width, legend_bottom_offset=legend_bottom_offset)

            # Display in Tkinter (use fixed size, no window scaling)
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
            self.canvas.draw()

            # Add canvas to frame, use fixed width/height, no window scaling
            canvas_widget = self.canvas.get_tk_widget()
            canvas_widget.pack_propagate(False)  # Prevent frame auto-resizing
            canvas_widget.pack(fill=tk.BOTH, expand=False)

            # Set canvas fixed size (based on figsize and dpi)
            dpi = 100
            canvas_width = int(figsize[0] * dpi)
            canvas_height = int(figsize[1] * dpi)
            canvas_widget.config(width=canvas_width, height=canvas_height)

            # Update canvas_frame size to match canvas
            self.canvas_frame.config(width=canvas_width, height=canvas_height)
            self.scroll_canvas.config(scrollregion=(0, 0, canvas_width, canvas_height))

            # Add navigation toolbar (supports zoom, pan, etc.)
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.canvas_frame)
            self.toolbar.update()

            # Force update scroll region
            self.root.update_idletasks()
            self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox("all"))

            self.lbl_hint.config(text="Chart generated successfully! Use the toolbar below to zoom, pan the chart", fg='green')

        except FileNotFoundError:
            messagebox.showerror("Error", f"File does not exist:\n{self.excel_file_path}")
            self.lbl_hint.config(text="Error: File does not exist", fg='red')
        except PermissionError:
            messagebox.showerror("Error", f"File is in use or no permission:\n{self.excel_file_path}")
            self.lbl_hint.config(text="Error: File in use or no permission", fg='red')
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")
            self.lbl_hint.config(text=f"Error: {str(e)}", fg='red')

    def _on_frame_configure(self, event):
        """Update scroll region when canvas size changes"""
        self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox("all"))

    def _on_vscroll(self, *args):
        """Vertical scroll callback"""
        self.scroll_canvas.yview(*args)

    def _on_hscroll(self, *args):
        """Horizontal scroll callback"""
        self.scroll_canvas.xview(*args)

    def redraw_chart(self):
        """Redraw chart with new color scheme"""
        if self.excel_file_path is None:
            messagebox.showwarning("Warning", "Please select an Excel file first")
            return
        self.plot_chart()

    def on_codon_label_toggle(self):
        """Codon label show/hide toggle"""
        if self.excel_file_path is None:
            # If no file loaded yet, only update checkbox state
            return
        # Redraw chart
        self.redraw_chart()

    def on_legend_blocks_toggle(self):
        """Stacked block show/hide toggle"""
        if self.excel_file_path is None:
            # If no file loaded yet, only update checkbox state
            return
        # Redraw chart
        self.redraw_chart()

    def save_as_pdf(self):
        """Save as PDF"""
        if self.fig is None:
            messagebox.showwarning("Warning", "Please load Excel file and generate chart first")
            return

        file_path = filedialog.asksaveasfilename(
            title="Save as PDF",
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")]
        )

        if file_path:
            self.fig.savefig(file_path, format='pdf', dpi=300, bbox_inches='tight')
            messagebox.showinfo("Success", f"PDF file saved to:\n{file_path}")

    def save_as_svg(self):
        """Save as SVG"""
        if self.fig is None:
            messagebox.showwarning("Warning", "Please load Excel file and generate chart first")
            return

        file_path = filedialog.asksaveasfilename(
            title="Save as SVG",
            defaultextension=".svg",
            filetypes=[("SVG files", "*.svg"), ("All files", "*.*")]
        )

        if file_path:
            self.fig.savefig(file_path, format='svg', dpi=300, bbox_inches='tight')
            messagebox.showinfo("Success", f"SVG file saved to:\n{file_path}")

    def save_as_png(self):
        """Save as PNG"""
        if self.fig is None:
            messagebox.showwarning("Warning", "Please load Excel file and generate chart first")
            return

        file_path = filedialog.asksaveasfilename(
            title="Save as PNG",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
        )

        if file_path:
            self.fig.savefig(file_path, format='png', dpi=300, bbox_inches='tight')
            messagebox.showinfo("Success", f"PNG file saved to:\n{file_path}")

    def run(self):
        """Run application"""
        self.root.mainloop()


# Command line mode (backward compatible)
def run_cli():
    """Run in command line mode"""
    import sys
    if len(sys.argv) > 1:
        excel_path = sys.argv[1]
    else:
        excel_path = '新建文件夹/A1.xlsx'

    if len(sys.argv) > 2:
        output_path = sys.argv[2]
    else:
        output_path = 'codon_usage_chart.pdf'

    fig, _, saved_path = plot_codon_usage_table_and_chart(excel_path, output_path)
    print(f"Chart saved to: {saved_path}")
    plt.show()


# Debug function: view Excel file content
def debug_excel_file(file_path):
    """Debug function: view Excel file content"""
    try:
        print(f"\nReading file: {file_path}")
        df = pd.read_excel(file_path)
        print(f"\nFile column names: {df.columns.tolist()}")
        print(f"Data types: {df.dtypes}")
        print(f"\nFirst 10 rows:")
        print(df.head(10))
        print(f"\nData shape: {df.shape}")
        return df
    except Exception as e:
        print(f"Failed to read Excel file: {e}")
        import traceback
        traceback.print_exc()
        return None


# Usage example
if __name__ == '__main__':
    import sys

    # If command line argument is --debug, enter debug mode
    if len(sys.argv) > 1 and sys.argv[1] == '--debug':
        print("=== Debug Mode ===")
        debug_excel_file('新建文件夹/A1.xlsx')
        input("\nPress Enter to exit...")
    else:
        # Start interactive application
        try:
            app = CodonUsageApp()
            app.create_widgets()
            app.run()
        except Exception as e:
            print(f"Program error: {e}")
            import traceback
            traceback.print_exc()
            input("Press Enter to exit...")
