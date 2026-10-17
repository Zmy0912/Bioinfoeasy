# -*- coding: utf-8 -*-
"""GC12-GC3 Regression Table

A GUI tool that follows the calculation method of the v3.0.0
gc12_gc3_plotter.py:

  - Parse .out files (header line skipped); column 9 is GC3s, column 10 is GC
  - GC12 = (GC1 + GC2) / 2 = (3 * GC - GC3) / 2
  - Ordinary least squares regression of GC12 on GC3:
        GC12 = slope * GC3 + intercept
  - Report, in table form, the linear equation, slope, R², etc. of every
    species (file)

Input can be a single .out file or a folder containing several .out files.
The results are shown in a table (Treeview) and can be exported to CSV.
"""

import os
import glob
import csv
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np


# ---------------------------------------------------------------------------
# Data parsing and regression (identical to v3.0.0)
# ---------------------------------------------------------------------------

def parse_out_file(filepath):
    """Parse a single .out file.

    Returns a dict with the GC3 / GC12 data and the regression statistics,
    or ``None`` when no valid data could be extracted.
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Skip the header line
    data_lines = lines[1:]

    gc3_data, gc12_data = [], []

    for line in data_lines:
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        # Need at least 11 columns (indices 0-10): name ... GC3s (9), GC (10)
        if len(parts) < 11:
            continue

        try:
            gc3_str = parts[9]   # GC3s
            gc_str = parts[10]   # overall GC

            # Skip invalid data
            if gc3_str == '*****' or gc_str == '*****':
                continue

            gc3 = float(gc3_str)
            gc = float(gc_str)

            # Accept both fractions (0-1) and percentages (0-100)
            if gc3 > 1.0:
                gc3 /= 100.0
            if gc > 1.0:
                gc /= 100.0

            if not (0.0 <= gc3 <= 1.0 and 0.0 <= gc <= 1.0):
                continue

            # GC overall = (GC1 + GC2 + GC3) / 3
            # => GC12 = (GC1 + GC2) / 2 = (3 * GC - GC3) / 2
            gc12 = (3.0 * gc - gc3) / 2.0

            if not (0.0 <= gc12 <= 1.0):
                continue

            gc3_data.append(gc3 * 100.0)
            gc12_data.append(gc12 * 100.0)

        except (ValueError, IndexError):
            continue

    if not gc3_data:
        return None

    dataset = {
        'name': os.path.basename(filepath),
        'gc3': gc3_data,
        'gc12': gc12_data,
    }
    dataset['stats'] = compute_regression(gc3_data, gc12_data)
    return dataset


def compute_regression(gc3_data, gc12_data):
    """Linear regression of GC12 on GC3 (GC12 = slope * GC3 + intercept).

    Returns a dict with ``slope``, ``intercept``, ``r``, ``r_squared``,
    ``n``, the mean values and the standard error of the slope.
    """
    x = np.array(gc3_data, dtype=float)
    y = np.array(gc12_data, dtype=float)
    n = len(x)

    if n < 2 or np.std(x) == 0:
        slope = 0.0
        intercept = float(y[0]) if n else 0.0
        r = 0.0
        slope_se = 0.0
    else:
        slope, intercept = np.polyfit(x, y, 1)
        if np.std(y) == 0:
            r = 0.0
        else:
            r = float(np.corrcoef(x, y)[0, 1])
        # Standard error of the slope: SE = sqrt( SSE/(n-2) / Sxx )
        y_fit = slope * x + intercept
        sse = float(np.sum((y - y_fit) ** 2))
        sxx = float(np.sum((x - x.mean()) ** 2))
        if n > 2 and sxx > 0:
            slope_se = float(np.sqrt(sse / (n - 2) / sxx))
        else:
            slope_se = 0.0

    return {
        'slope': float(slope),
        'intercept': float(intercept),
        'r': r,
        'r_squared': r ** 2,
        'n': n,
        'gc3_mean': float(x.mean()) if n else 0.0,
        'gc12_mean': float(y.mean()) if n else 0.0,
        'slope_se': slope_se,
    }


def format_equation(slope, intercept):
    """Format the regression equation, e.g. y = 0.1234x + 45.6789."""
    if intercept >= 0:
        return f"y = {slope:.4f}x + {intercept:.4f}"
    return f"y = {slope:.4f}x - {abs(intercept):.4f}"


# ---------------------------------------------------------------------------
# Graphical user interface
# ---------------------------------------------------------------------------

class GC12GC3RegressionTable:
    # Table columns: (column id, heading, width, is numeric)
    TABLE_COLUMNS = [
        ('species',   'Species / File',        220, False),
        ('equation',  'Linear Equation',       220, False),
        ('slope',     'Slope',                 90, True),
        ('intercept', 'Intercept',            100, True),
        ('r',         'r',                      80, True),
        ('r2',        'R²',                     80, True),
        ('n',         'Genes (n)',              80, True),
        ('gc3_mean',  'GC3 mean (%)',           95, True),
        ('gc12_mean', 'GC12 mean (%)',          95, True),
        ('slope_se',  'Slope SE',               90, True),
    ]

    def __init__(self, root):
        self.root = root
        self.root.title("GC12-GC3 Regression Table")
        self.root.geometry("1150x640")

        # Datasets of the last analysis (kept for export)
        self.datasets = []

        self.create_widgets()

    # ------------------------- Widget construction -------------------------

    def create_widgets(self):
        control = ttk.LabelFrame(self.root, text="Input Selection", padding=10)
        control.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(control, text="Input mode:").grid(row=0, column=0, sticky=tk.W, pady=4)
        self.input_mode = tk.StringVar(value='folder')
        mode_frame = ttk.Frame(control)
        mode_frame.grid(row=0, column=1, columnspan=2, sticky=tk.W, padx=5)
        ttk.Radiobutton(mode_frame, text="Single file",
                        variable=self.input_mode, value='file',
                        command=self.on_mode_change).pack(side=tk.LEFT)
        ttk.Radiobutton(mode_frame, text="Folder (all .out files)",
                        variable=self.input_mode, value='folder',
                        command=self.on_mode_change).pack(side=tk.LEFT, padx=(15, 0))

        ttk.Label(control, text="Input file (.out):").grid(row=1, column=0, sticky=tk.W, pady=4)
        self.file_path = tk.StringVar()
        ttk.Entry(control, textvariable=self.file_path).grid(
            row=1, column=1, padx=5, sticky=tk.EW)
        ttk.Button(control, text="Browse File...", command=self.browse_file).grid(
            row=1, column=2, padx=5)

        ttk.Label(control, text="Input folder:").grid(row=2, column=0, sticky=tk.W, pady=4)
        self.folder_path = tk.StringVar(value=self._default_out_folder())
        ttk.Entry(control, textvariable=self.folder_path).grid(
            row=2, column=1, padx=5, sticky=tk.EW)
        ttk.Button(control, text="Browse Folder...", command=self.browse_folder).grid(
            row=2, column=2, padx=5)

        ttk.Label(control, text="Name prefix:").grid(row=3, column=0, sticky=tk.W, pady=4)
        self.name_prefix = tk.StringVar()
        ttk.Entry(control, textvariable=self.name_prefix, width=25).grid(
            row=3, column=1, sticky=tk.W, padx=5)

        ttk.Button(control, text="Analyze", command=self.run_analysis).grid(
            row=3, column=2, padx=5, sticky=tk.E)

        control.columnconfigure(1, weight=1)

        # ---------------- Results table ----------------
        table_frame = ttk.LabelFrame(
            self.root, text="Regression Statistics (click a column header to sort)",
            padding=5)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        columns = [c[0] for c in self.TABLE_COLUMNS]
        self.tree = ttk.Treeview(table_frame, columns=columns, show='headings')

        for col_id, heading, width, _ in self.TABLE_COLUMNS:
            self.tree.heading(col_id, text=heading,
                              command=lambda c=col_id: self.sort_by_column(c, False))
            self.tree.column(col_id, width=width,
                             anchor=tk.CENTER if col_id != 'species' else tk.W)

        vsb = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.tree.yview)
        hsb = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.grid(row=0, column=0, sticky='nsew')
        vsb.grid(row=0, column=1, sticky='ns')
        hsb.grid(row=1, column=0, sticky='ew')
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)

        # Status bar and export button
        bottom = ttk.Frame(self.root)
        bottom.pack(fill=tk.X, padx=10, pady=5)

        self.status_bar = ttk.Label(bottom, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.LEFT, fill=tk.X, expand=True)

        ttk.Button(bottom, text="Export CSV", command=self.export_csv).pack(
            side=tk.RIGHT, padx=(5, 0))

        self.on_mode_change()

    def _default_out_folder(self):
        """Locate the ``out`` folder next to this script (if it exists)."""
        try:
            base = os.path.dirname(os.path.abspath(__file__))
            candidate = os.path.join(os.path.dirname(base), 'out')
            if os.path.isdir(candidate):
                return candidate
            candidate = os.path.join(base, 'out')
            if os.path.isdir(candidate):
                return candidate
        except Exception:
            pass
        return ''

    # ------------------------- Event handlers -------------------------

    def on_mode_change(self):
        if self.input_mode.get() == 'folder':
            self.status_bar.config(
                text="Folder mode: the regression statistics of every .out "
                     "file in the folder will be computed")
        else:
            self.status_bar.config(text="Single-file mode")

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Select .out file",
            filetypes=[("OUT files", "*.out"), ("All files", "*.*")]
        )
        if filename:
            self.input_mode.set('file')
            self.file_path.set(filename)
            self.status_bar.config(text=f"Selected file: {os.path.basename(filename)}")

    def browse_folder(self):
        folder = filedialog.askdirectory(
            title="Select folder containing .out files")
        if folder:
            self.input_mode.set('folder')
            self.folder_path.set(folder)
            n = len(glob.glob(os.path.join(folder, '*.out')))
            self.status_bar.config(
                text=f"Selected folder: {folder}  ({n} .out file(s) found)")

    # ------------------------- Main workflow -------------------------

    def collect_files(self):
        """Return the .out files to analyze according to the input mode."""
        if self.input_mode.get() == 'file':
            path = self.file_path.get().strip()
            if not path or not os.path.isfile(path):
                messagebox.showerror("Error", "Please select a valid input file!")
                return []
            return [path]

        folder = self.folder_path.get().strip()
        if not folder or not os.path.isdir(folder):
            messagebox.showerror("Error", "Please select a valid input folder!")
            return []
        files = sorted(glob.glob(os.path.join(folder, '*.out')))
        if not files:
            messagebox.showerror("Error", "No .out files found in the selected folder!")
            return []
        return files

    def run_analysis(self):
        files = self.collect_files()
        if not files:
            return

        prefix = self.name_prefix.get().strip()
        datasets = []
        failed = []

        for filepath in files:
            try:
                dataset = parse_out_file(filepath)
            except Exception as e:
                print(f"Failed to read {filepath}: {e}")
                dataset = None
            if dataset is None:
                failed.append(os.path.basename(filepath))
                continue
            # Result name: optional prefix + file name without extension
            base = os.path.splitext(os.path.basename(filepath))[0]
            dataset['name'] = f"{prefix}{base}" if prefix else base
            datasets.append(dataset)

        if not datasets:
            messagebox.showerror(
                "Error", "No valid data could be extracted from any file!")
            return

        self.datasets = datasets
        self.fill_table()

        msg = f"Analysis complete: {len(datasets)} species / file(s)"
        if failed:
            msg += f" ({len(failed)} file(s) without valid data were skipped)"
        self.status_bar.config(text=msg)

    def fill_table(self, sort_col=None, reverse=False):
        """Fill the table with ``self.datasets``."""
        self.tree.delete(*self.tree.get_children())

        rows = []
        for ds in self.datasets:
            s = ds['stats']
            rows.append({
                'species': ds['name'],
                'equation': format_equation(s['slope'], s['intercept']),
                'slope': f"{s['slope']:.4f}",
                'intercept': f"{s['intercept']:.4f}",
                'r': f"{s['r']:.4f}",
                'r2': f"{s['r_squared']:.4f}",
                'n': str(s['n']),
                'gc3_mean': f"{s['gc3_mean']:.2f}",
                'gc12_mean': f"{s['gc12_mean']:.2f}",
                'slope_se': f"{s['slope_se']:.4f}",
            })

        if sort_col:
            numeric = {c[0] for c in self.TABLE_COLUMNS if c[3]}
            rows.sort(
                key=lambda r: float(r[sort_col]) if sort_col in numeric else r[sort_col],
                reverse=reverse)

        for row in rows:
            self.tree.insert('', tk.END, values=[row[c[0]] for c in self.TABLE_COLUMNS])

    def sort_by_column(self, col, reverse):
        """Sort the table by a column (click on the header)."""
        items = [(self.tree.set(k, col), k) for k in self.tree.get_children('')]
        numeric = {c[0] for c in self.TABLE_COLUMNS if c[3]}
        if col in numeric:
            items.sort(key=lambda t: float(t[0]), reverse=reverse)
        else:
            items.sort(key=lambda t: t[0], reverse=reverse)

        for index, (_, k) in enumerate(items):
            self.tree.move(k, '', index)

        # Clicking the same header again reverses the order
        self.tree.heading(col, command=lambda: self.sort_by_column(col, not reverse))

    # ------------------------- Export -------------------------

    def export_csv(self):
        if not self.datasets:
            messagebox.showwarning("Warning", "Please run the analysis first!")
            return

        filename = filedialog.asksaveasfilename(
            title="Export CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not filename:
            return

        try:
            with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
                writer = csv.writer(f)
                writer.writerow(['Species/File', 'Linear Equation', 'Slope',
                                 'Intercept', 'r', 'R²', 'Genes (n)',
                                 'GC3 mean (%)', 'GC12 mean (%)',
                                 'Slope SE'])
                for ds in self.datasets:
                    s = ds['stats']
                    writer.writerow([
                        ds['name'],
                        format_equation(s['slope'], s['intercept']),
                        f"{s['slope']:.4f}",
                        f"{s['intercept']:.4f}",
                        f"{s['r']:.4f}",
                        f"{s['r_squared']:.4f}",
                        s['n'],
                        f"{s['gc3_mean']:.2f}",
                        f"{s['gc12_mean']:.2f}",
                        f"{s['slope_se']:.4f}",
                    ])
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")
            return

        self.status_bar.config(text=f"Exported: {os.path.basename(filename)}")
        messagebox.showinfo("Success", f"Table saved to:\n{filename}")


def main():
    root = tk.Tk()
    GC12GC3RegressionTable(root)
    root.mainloop()


if __name__ == "__main__":
    main()
