#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MISA SSR Analyzer
=================
A GUI tool for visualizing MISA (MIcroSAtellite) search results.

Features
--------
1. Heatmap (Panel B)  : SSR motif presence/absence across species, with motifs
                        grouped by repeat-unit length (1 - 10+ bp) and species
                        ordered by Ward hierarchical clustering.
2. Pie chart (Panel C): Percentage distribution of SSR categories (p1-p6, c, c*)
                        across all species.
3. Stacked bar (Panel D): Per-species distribution and abundance of SSR
                          categories, with species sorted by total SSR count
                          (descending).
4. Exports the underlying data tables (tab-separated, UTF-8-BOM) for use in
   GraphPad Prism, Excel, TBtools, R, Origin, etc.

Requirements
------------
Python >= 3.8, tkinter (bundled with most Python builds), matplotlib >= 3.5,
numpy.

Run
---
    python misa_ssr_viewer.py
"""

import os
import re
import sys
from collections import defaultdict, OrderedDict

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import numpy as np


# ---------------------------------------------------------------------------
# Data layer
# ---------------------------------------------------------------------------
MOTIF_RE = re.compile(r'\(([ACGTN]{1,60})\)(\d+)')

TYPE_INFO = OrderedDict([
    ('p1', ('Mono-nucleotide',  1)),
    ('p2', ('Di-nucleotide',   2)),
    ('p3', ('Tri-nucleotide',  3)),
    ('p4', ('Tetra-nucleotide',4)),
    ('p5', ('Penta-nucleotide',5)),
    ('p6', ('Hexa-nucleotide', 6)),
    ('c',  ('Compound (c)',    0)),
    ('c*', ('Compound (c*)',   0)),
])


def parse_misa(path):
    """Read a MISA .misa file and return a list of SSR records."""
    records = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        fh.readline()  # skip header
        for line in fh:
            line = line.rstrip('\r\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            records.append({'id': cols[0], 'type': cols[2], 'ssr': cols[3]})
    return records


def extract_motifs(ssr_string):
    """Extract every repeat-unit motif from a (possibly compound) SSR string."""
    return [m.upper() for m, _ in MOTIF_RE.findall(ssr_string)]


def ward_leaf_order(matrix):
    """Ward minimum-variance hierarchical clustering on a binary matrix.

    Distance = squared Euclidean (= number of differing features).  Returns
    the leaf order produced by Ward's linkage, suitable for plotting the
    dendrogram leaves in the same order as the resulting heatmap rows.
    """
    n = len(matrix)
    d = len(matrix[0]) if n else 0
    if n <= 1:
        return list(range(n))

    dist = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = sum((matrix[i][k] - matrix[j][k]) ** 2 for k in range(d))

    clusters = [{'id': i, 'mem': [i], 'tree': i, 'sz': 1} for i in range(n)]
    nid = n
    while len(clusters) > 1:
        ba, bb, bv = 0, 1, None
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                ca, cb = clusters[a]['id'], clusters[b]['id']
                v = dist[(min(ca, cb), max(ca, cb))]
                if bv is None or v < bv:
                    bv, ba, bb = v, a, b
        A, B = clusters[ba], clusters[bb]
        na, nb = A['sz'], B['sz']
        for C in clusters:
            if C is A or C is B:
                continue
            nc = C['sz']
            ka = (min(A['id'], C['id']), max(A['id'], C['id']))
            kb = (min(B['id'], C['id']), max(B['id'], C['id']))
            d_ab = dist[(min(A['id'], B['id']), max(A['id'], B['id']))]
            d_new = ((na + nc) / (na + nb + nc)) * dist[ka] \
                  + ((nb + nc) / (na + nb + nc)) * dist[kb] \
                  - (nc / (na + nb + nc)) * d_ab
            dist[(min(nid, C['id']), max(nid, C['id']))] = d_new
        merged = {
            'id':  nid,
            'mem': A['mem'] + B['mem'],
            'tree': (A['tree'], B['tree']),
            'sz':  na + nb,
        }
        clusters = [c for i, c in enumerate(clusters) if i not in (ba, bb)]
        clusters.append(merged)
        nid += 1

    def flatten(tree):
        if isinstance(tree, int):
            return [tree]
        out = []
        for sub in tree:
            out.extend(flatten(sub))
        return out

    return flatten(clusters[0]['tree'])


class Analyzer:
    """Holds processed MISA data and computes the matrices used by the GUI."""

    def __init__(self):
        self.records = []
        self.species = []
        self.species_type_count = defaultdict(lambda: defaultdict(int))
        self.species_motifs = defaultdict(set)
        self.type_total = defaultdict(int)
        self.total_ssr = 0
        self.all_motifs = []
        self.ward_species = []
        self.mat = None          # ndarray, shape (species, motif), Ward order

    def process(self, path):
        self.records = parse_misa(path)
        self.species = sorted({r['id'] for r in self.records})
        self.species_type_count = defaultdict(lambda: defaultdict(int))
        self.species_motifs = defaultdict(set)
        self.type_total = defaultdict(int)
        for r in self.records:
            self.species_type_count[r['id']][r['type']] += 1
            self.type_total[r['type']] += 1
            for mo in extract_motifs(r['ssr']):
                self.species_motifs[r['id']].add(mo)

        self.total_ssr = len(self.records)
        self.all_motifs = sorted(
            {m for s in self.species for m in self.species_motifs[s]},
            key=lambda x: (len(x), x),
        )

        mat_orig = [[1 if m in self.species_motifs[s] else 0 for m in self.all_motifs]
                    for s in self.species]
        order = ward_leaf_order(mat_orig)
        self.ward_species = [self.species[i] for i in order]
        self.mat = np.array(
            [[mat_orig[i][j] for j in range(len(self.all_motifs))] for i in order],
            dtype=int,
        )

    # ----- data-table export -----
    def export_tables(self, out_dir):
        os.makedirs(out_dir, exist_ok=True)

        def write(name, header, rows):
            p = os.path.join(out_dir, name)
            with open(p, 'w', encoding='utf-8-sig', newline='') as fh:
                fh.write('\t'.join(header) + '\n')
                for r in rows:
                    fh.write('\t'.join(str(x) for x in r) + '\n')

        motifs = self.all_motifs
        species = self.species
        ward_sp = self.ward_species
        sm = self.species_motifs
        stc = self.species_type_count
        tt = self.type_total
        total = self.total_ssr

        # B1
        write('B1_heatmap_presence_absence_matrix.txt',
              ['Species\\Motif'] + motifs,
              [[s] + [int(self.mat[i, j]) for j in range(len(motifs))]
               for i, s in enumerate(ward_sp)])
        write('B1b_heatmap_matrix_original_species_order.txt',
              ['Species\\Motif'] + motifs,
              [[s] + [1 if m in sm[s] else 0 for m in motifs] for s in species])
        # B2
        write('B2_motif_info.txt',
              ['MotifOrder', 'Motif', 'RepeatUnitLength_bp', 'SpeciesCount', 'Frequency_pct'],
              [[i + 1, mo, len(mo),
                sum(1 for s in species if mo in sm[s]),
                f"{100.0 * sum(1 for s in species if mo in sm[s]) / len(species):.2f}%"]
               for i, mo in enumerate(motifs)])
        # B3
        groups = defaultdict(int)
        for mo in motifs:
            key = str(len(mo)) if len(mo) <= 10 else '10+'
            groups[key] += 1
        order_keys = [str(i) for i in range(1, 11)] + ['10+']
        write('B3_motif_length_group_count.txt',
              ['RepeatUnitLength_bp', 'MotifCount'],
              [[k, groups[k]] for k in order_keys if k in groups])
        # B4
        write('B4_species_ward_cluster_order.txt',
              ['ClusterOrder', 'Species', 'Total_SSR', 'MotifCount'],
              [[i + 1, s, sum(stc[s].values()), len(sm[s])]
               for i, s in enumerate(ward_sp)])

        # C
        rows_c = [[t, TYPE_INFO[t][0], tt.get(t, 0),
                   f"{100.0 * tt.get(t, 0) / total:.2f}%"] for t in TYPE_INFO]
        write('C_pie_SSR_category_percentage.txt',
              ['MISA_Type', 'Category', 'SSR_Count', 'Percentage'], rows_c)
        merged = [[TYPE_INFO[t][0], tt.get(t, 0),
                   f"{100.0 * tt.get(t, 0) / total:.2f}%"]
                  for t in ['p1', 'p2', 'p3', 'p4', 'p5', 'p6']]
        comp = tt.get('c', 0) + tt.get('c*', 0)
        merged.append(['Compound SSR (c + c*)', comp,
                       f"{100.0 * comp / total:.2f}%"])
        write('C_pie_SSR_category_percentage_merged.txt',
              ['Category', 'SSR_Count', 'Percentage'], merged)

        # D
        sp_by_total = sorted(species, key=lambda s: -sum(stc[s].values()))
        cats = list(TYPE_INFO.keys())
        write('D1_stackedbar_species_category_counts.txt',
              ['Species'] + [f"{t}({TYPE_INFO[t][0]})" for t in cats] + ['Total_SSR'],
              [[s] + [stc[s].get(t, 0) for t in cats] + [sum(stc[s].values())]
               for s in sp_by_total])
        long_rows = []
        for s in sp_by_total:
            tot = sum(stc[s].values())
            for t in cats:
                long_rows.append([s, t, TYPE_INFO[t][0], stc[s].get(t, 0), tot])
        write('D2_stackedbar_long_format.txt',
              ['Species', 'MISA_Type', 'Category', 'Count', 'SpeciesTotal'],
              long_rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def get_cmap(name):
    """Return a colormap in a way that works across matplotlib versions."""
    try:
        return matplotlib.colormaps[name]
    except Exception:
        return cm.get_cmap(name)


def plot_heatmap(analyzer, fig):
    fig.clear()
    mat = analyzer.mat
    motifs = analyzer.all_motifs
    species = analyzer.ward_species
    n_sp, n_mo = mat.shape

    # Group motifs by repeat-unit length
    groups = OrderedDict()
    for i, m in enumerate(motifs):
        L = len(m)
        key = L if L <= 10 else '10+'
        groups.setdefault(key, []).append(i)

    # Layout: top subplot for the length grouping, bottom for the heatmap
    ax_top = fig.add_subplot(2, 1, 1)
    ax = fig.add_subplot(2, 1, 2)
    main_colors = get_cmap('tab10')

    for g_key, idxs in groups.items():
        if not idxs:
            continue
        x0, x1 = min(idxs) - 0.5, max(idxs) + 0.5
        centre = (x0 + x1) / 2
        width = x1 - x0
        color_idx = (g_key if isinstance(g_key, int) else 10) % 10
        ax_top.bar([centre], [len(idxs)], width=width * 0.92, align='center',
                   color=main_colors(color_idx), edgecolor='black', linewidth=0.4)
        ax_top.text(centre, len(idxs) + 0.12, f"{len(idxs)}",
                    ha='center', va='bottom', fontsize=8)
        ax_top.text(centre, -0.3, f"{g_key} bp", ha='center', va='top',
                    rotation=90, fontsize=8, color=main_colors(color_idx),
                    fontweight='bold')
    max_cnt = max(len(idxs) for idxs in groups.values())
    ax_top.set_xlim(-0.5, n_mo - 0.5)
    ax_top.set_ylim(-2.0, max_cnt * 1.25 + 0.4)
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_top.set_title('Motifs grouped by repeat-unit length (1 - 10+ bp)',
                     fontsize=10)

    # Main heatmap
    cmap = mcolors.ListedColormap(['white', '#1f77b4'])
    ax.imshow(mat, aspect='auto', cmap=cmap, vmin=0, vmax=1)
    for g_key, idxs in groups.items():
        if not idxs:
            continue
        x0, x1 = min(idxs) - 0.5, max(idxs) + 0.5
        color_idx = (g_key if isinstance(g_key, int) else 10) % 10
        ax.axvspan(x0, x1, alpha=0.10, color=main_colors(color_idx), zorder=0)
    ax.set_yticks(range(n_sp))
    ax.set_yticklabels(species, fontsize=9)
    ax.set_xticks(range(n_mo))
    ax.set_xticklabels(motifs, rotation=90, fontsize=8)
    ax.set_xlabel('SSR motif')
    ax.set_ylabel('Species')
    ax.set_title('B. SSR motif presence (blue = 1) / absence (white = 0)\n'
                 'Species ordered by Ward hierarchical clustering',
                 fontsize=11)
    fig.tight_layout()


def plot_pie(analyzer, fig):
    fig.clear()
    ax = fig.add_subplot(111)
    palette = get_cmap('tab10')

    # Collect (label, count) and skip zero categories
    entries = []
    for t, (name, _) in TYPE_INFO.items():
        c = analyzer.type_total.get(t, 0)
        if c > 0:
            pct = 100.0 * c / analyzer.total_ssr
            entries.append((name, c, pct))
    sizes = [e[1] for e in entries]
    colors = [palette(i % 10) for i in range(len(entries))]

    wedges, _ = ax.pie(
        sizes, colors=colors, startangle=90,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.0},
    )

    # Place labels outside with leader lines, avoiding overlap on each side.
    right, left = [], []
    for w, (name, cnt, pct) in zip(wedges, entries):
        ang = np.deg2rad((w.theta1 + w.theta2) / 2.0)
        x, y = float(np.cos(ang)), float(np.sin(ang))
        (right if x >= 0 else left).append([y, ang, x, name, cnt, pct])

    min_dy = 0.20  # minimum vertical spacing between labels
    for group in (right, left):
        group.sort(key=lambda t: t[0])
        for i in range(1, len(group)):
            if group[i][0] - group[i - 1][0] < min_dy:
                group[i][0] = group[i - 1][0] + min_dy
        for y, ang, x, name, cnt, pct in group:
            side = 1 if x >= 0 else -1
            x0, y0 = np.cos(ang) * 1.02, np.sin(ang) * 1.02
            x1, y1 = side * 1.30, y
            ax.annotate(
                f"{name}\n{cnt} ({pct:.2f}%)",
                xy=(x0, y0), xytext=(x1, y1),
                ha='left' if side > 0 else 'right', va='center',
                fontsize=9,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.7,
                                shrinkA=0, shrinkB=4),
            )

    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect('equal')
    ax.set_title('C. Percentage distribution of SSR categories\n'
                 'across all analysed species', fontsize=12)
    fig.tight_layout()


def plot_stacked_bar(analyzer, fig):
    fig.clear()
    ax = fig.add_subplot(111)
    species_by_total = sorted(
        analyzer.species,
        key=lambda s: -sum(analyzer.species_type_count[s].values()),
    )
    cats = list(TYPE_INFO.keys())
    counts = np.array(
        [[analyzer.species_type_count[s].get(t, 0) for t in cats]
         for s in species_by_total],
        dtype=int,
    )
    palette = get_cmap('tab10')
    bottoms = np.zeros(len(species_by_total), dtype=int)
    for j, t in enumerate(cats):
        ax.bar(species_by_total, counts[:, j], bottom=bottoms,
               label=TYPE_INFO[t][0], color=palette(j % 10),
               edgecolor='white', linewidth=0.4)
        bottoms += counts[:, j]
    ax.set_ylabel('SSR count')
    ax.set_xlabel('Species (sorted by total SSR count, descending)')
    ax.set_title('D. SSR category distribution and abundance per species',
                 fontsize=12)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    max_len = max((len(s) for s in species_by_total), default=0)
    if max_len > 3 or len(species_by_total) > 12:
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(45)
            lbl.set_ha('right')
    fig.tight_layout()


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MISA SSR Analyzer")
        self.geometry("1180x820")
        self.minsize(900, 650)

        self.analyzer = Analyzer()
        self.misa_path = tk.StringVar()
        self.out_dir = tk.StringVar(value=os.getcwd())
        self.status = tk.StringVar(value="Ready. Open a MISA .misa file to begin.")

        self._build_style()
        self._build_menu()
        self._build_toolbar()
        self._build_notebook()
        self._build_statusbar()

    # ----- style / layout -----
    def _build_style(self):
        try:
            style = ttk.Style(self)
            if 'clam' in style.theme_names():
                style.theme_use('clam')
            style.configure('TButton', padding=4)
        except tk.TclError:
            pass

    def _build_menu(self):
        menubar = tk.Menu(self)
        m_file = tk.Menu(menubar, tearoff=0)
        m_file.add_command(label="Open MISA file...", command=self.choose_misa,
                           accelerator="Ctrl+O")
        m_file.add_command(label="Set output directory...",
                           command=self.choose_outdir)
        m_file.add_separator()
        m_file.add_command(label="Export data tables", command=self.export_tables)
        m_file.add_separator()
        m_file.add_command(label="Quit", command=self.destroy,
                           accelerator="Ctrl+Q")
        menubar.add_cascade(label="File", menu=m_file)

        m_plot = tk.Menu(menubar, tearoff=0)
        m_plot.add_command(label="Save current figure as PNG...",
                           command=self.save_current_figure)
        m_plot.add_command(label="Save all figures as PNG...",
                           command=self.save_all_figures)
        menubar.add_cascade(label="Plots", menu=m_plot)

        m_help = tk.Menu(menubar, tearoff=0)
        m_help.add_command(label="About", command=self.show_about)
        menubar.add_cascade(label="Help", menu=m_help)

        self.config(menu=menubar)
        self.bind_all('<Control-o>', lambda e: self.choose_misa())
        self.bind_all('<Control-q>', lambda e: self.destroy())

    def _build_toolbar(self):
        bar = ttk.Frame(self, padding=(8, 6))
        bar.pack(fill='x')

        ttk.Label(bar, text="MISA file:").pack(side='left')
        ttk.Entry(bar, textvariable=self.misa_path, width=55)\
            .pack(side='left', padx=(4, 4))
        ttk.Button(bar, text="Browse...", command=self.choose_misa)\
            .pack(side='left')

        ttk.Label(bar, text="  Output dir:").pack(side='left')
        ttk.Entry(bar, textvariable=self.out_dir, width=30)\
            .pack(side='left', padx=(4, 4))
        ttk.Button(bar, text="Browse...", command=self.choose_outdir)\
            .pack(side='left')

        ttk.Button(bar, text="Run Analysis", command=self.run_analysis)\
            .pack(side='left', padx=(10, 4))
        ttk.Button(bar, text="Export Tables", command=self.export_tables)\
            .pack(side='left', padx=4)

    def _build_notebook(self):
        self.nb = ttk.Notebook(self)
        self.nb.pack(fill='both', expand=True, padx=8, pady=(2, 2))

        self.tabs = {}
        self.figs = {}
        self.canvases = {}
        for key, title in [('B', 'B. Heatmap'),
                           ('C', 'C. Pie chart'),
                           ('D', 'D. Stacked bar')]:
            frame = ttk.Frame(self.nb)
            self.nb.add(frame, text=title)
            fig = Figure(figsize=(10, 7), dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.get_tk_widget().pack(fill='both', expand=True)
            toolbar = NavigationToolbar2Tk(canvas, frame)
            toolbar.update()
            self.tabs[key] = frame
            self.figs[key] = fig
            self.canvases[key] = canvas

    def _build_statusbar(self):
        ttk.Separator(self, orient='horizontal').pack(fill='x')
        ttk.Label(self, textvariable=self.status, anchor='w',
                  relief='sunken', padding=(6, 2)).pack(fill='x')

    # ----- actions -----
    def choose_misa(self):
        path = filedialog.askopenfilename(
            title="Select a MISA result file",
            filetypes=[("MISA output", "*.misa"), ("All files", "*.*")],
        )
        if path:
            self.misa_path.set(path)
            self.status.set(f"Selected: {os.path.basename(path)}")

    def choose_outdir(self):
        d = filedialog.askdirectory(title="Select output directory",
                                    initialdir=self.out_dir.get() or os.getcwd())
        if d:
            self.out_dir.set(d)

    def run_analysis(self):
        path = self.misa_path.get().strip()
        if not path or not os.path.isfile(path):
            messagebox.showerror("Error", "Please choose a valid MISA .misa file.")
            return
        try:
            self.status.set("Processing...")
            self.update_idletasks()
            self.analyzer.process(path)
        except Exception as exc:
            messagebox.showerror("Processing failed", str(exc))
            self.status.set("Failed.")
            return

        plot_heatmap(self.analyzer, self.figs['B'])
        plot_pie(self.analyzer, self.figs['C'])
        plot_stacked_bar(self.analyzer, self.figs['D'])
        for c in self.canvases.values():
            c.draw_idle()

        self.status.set(
            f"Loaded {self.analyzer.total_ssr} SSRs from "
            f"{len(self.analyzer.species)} species "
            f"({len(self.analyzer.all_motifs)} distinct motifs)."
        )

    def export_tables(self):
        if not self.analyzer.records:
            messagebox.showinfo("No data",
                                "Run an analysis before exporting tables.")
            return
        out = self.out_dir.get().strip() or os.getcwd()
        try:
            self.analyzer.export_tables(out)
        except Exception as exc:
            messagebox.showerror("Export failed", str(exc))
            return
        self.status.set(f"Data tables exported to: {out}")
        messagebox.showinfo("Export complete",
                            f"Data tables written to:\n{out}")

    def _save_figure(self, key):
        if not self.analyzer.records:
            messagebox.showinfo("No data", "Run an analysis first.")
            return None
        default = os.path.join(self.out_dir.get() or os.getcwd(),
                               f"figure_{key}.png")
        path = filedialog.asksaveasfilename(
            title=f"Save figure {key}",
            defaultextension=".png",
            initialfile=os.path.basename(default),
            initialdir=os.path.dirname(default),
            filetypes=[("PNG image", "*.png"),
                       ("PDF document", "*.pdf"),
                       ("SVG vector", "*.svg")],
        )
        if not path:
            return None
        self.figs[key].savefig(path, dpi=300, bbox_inches='tight')
        return path

    def save_current_figure(self):
        idx = self.nb.index(self.nb.select())
        key = list(self.tabs.keys())[idx]
        path = self._save_figure(key)
        if path:
            self.status.set(f"Saved figure {key}: {path}")

    def save_all_figures(self):
        if not self.analyzer.records:
            messagebox.showinfo("No data", "Run an analysis first.")
            return
        out = self.out_dir.get() or os.getcwd()
        os.makedirs(out, exist_ok=True)
        saved = []
        for key in self.tabs:
            p = os.path.join(out, f"figure_{key}.png")
            self.figs[key].savefig(p, dpi=300, bbox_inches='tight')
            saved.append(p)
        self.status.set("Saved figures: " + "; ".join(saved))
        messagebox.showinfo("Figures saved",
                            "Saved:\n" + "\n".join(saved))

    def show_about(self):
        messagebox.showinfo(
            "About MISA SSR Analyzer",
            "MISA SSR Analyzer\n"
            "-----------------\n"
            "A visualization tool for MISA SSR search results.\n\n"
            "Panels:\n"
            "  B. Heatmap of SSR motif presence/absence\n"
            "  C. Pie chart of SSR category percentages\n"
            "  D. Stacked bar of category abundance per species\n\n"
            f"Python {sys.version.split()[0]}, "
            f"matplotlib {matplotlib.__version__}",
        )


def main():
    app = App()
    app.mainloop()


if __name__ == '__main__':
    main()