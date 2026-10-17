# -*- coding: utf-8 -*-
"""
Direct / Palindromic Repeat Analyzer
=====================================
Parses Vmatch `repfind` output files (one per species, one file = one species).
Each data line has the format:
    length1  pos1  P/D  length2  pos2  distance  E-value
where P = Palindromic (reverse) repeat and D = Direct (forward) repeat.
Statistics are displayed in a table and can be exported to CSV (Excel-ready).
"""

import csv
import os
import statistics
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.constants import W


def parse_repeat_file(path):
    """Parse a single repfind output file.

    Returns a list of records: [(length, type 'P'/'D', evalue), ...]
    """
    records = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                length = int(parts[0])
                rtype = parts[2].upper()
                evalue = float(parts[6])
            except ValueError:
                continue
            if rtype not in ("P", "D"):
                continue
            records.append((length, rtype, evalue))
    return records


def species_name_from_file(filename):
    """Derive the species name from the file name:
    strip the extension and any trailing 'repeats' marker."""
    name = os.path.splitext(os.path.basename(filename))[0]
    lower = name.lower()
    if lower.endswith("repeats"):
        name = name[: -len("repeats")]
    return name.strip()


class RepeatAnalyzerApp:
    HEADERS = [
        "Species",
        "Direct (D) Repeats",
        "Palindromic (P) Repeats",
        "Total",
        "Total Length (bp)",
        "Mean Length (bp)",
        "Min Length (bp)",
        "Max Length (bp)",
    ]
    COL_WIDTHS = [120, 120, 140, 70, 120, 110, 110, 110]

    def __init__(self, root):
        self.root = root
        self.root.title("Direct / Palindromic Repeat Analyzer")
        self.root.geometry("1020x620")
        self.records = {}  # species name -> record list
        self.table_data = []
        self.total_row = None

        self._build_toolbar()
        self._build_table()
        self._build_statusbar()

    # ---------- UI ----------
    def _build_toolbar(self):
        top = ttk.Frame(self.root, padding=8)
        top.pack(fill="x")

        ttk.Button(top, text="Select Folder", command=self.choose_folder).pack(side="left")
        self.folder_var = tk.StringVar()
        ttk.Label(top, textvariable=self.folder_var, foreground="gray").pack(
            side="left", padx=8
        )

        ttk.Label(top, text="Min length:").pack(side="left", padx=(16, 4))
        self.min_len_var = tk.StringVar(value="10")
        ttk.Spinbox(
            top, from_=1, to=1000, textvariable=self.min_len_var, width=6
        ).pack(side="left")

        ttk.Label(top, text="Max E-value (blank = no limit):").pack(side="left", padx=(16, 4))
        self.max_e_var = tk.StringVar(value="")
        ttk.Entry(top, textvariable=self.max_e_var, width=8).pack(side="left")

        ttk.Button(top, text="Analyze", command=self.analyze).pack(side="left", padx=16)
        ttk.Button(top, text="Export CSV", command=self.export_csv).pack(side="left")

    def _build_table(self):
        frame = ttk.Frame(self.root)
        frame.pack(fill="both", expand=True, padx=8, pady=(0, 4))

        columns = list(range(len(self.HEADERS)))
        self.tree = ttk.Treeview(frame, columns=columns, show="headings")
        for i, (h, w) in enumerate(zip(self.HEADERS, self.COL_WIDTHS)):
            self.tree.heading(i, text=h, command=lambda c=i: self.sort_by(c))
            self.tree.column(i, width=w, anchor="center")
        self.tree.column(0, anchor="w")

        vsb = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

    def _build_statusbar(self):
        self.status_var = tk.StringVar(
            value="Select a folder containing repeat result files to begin"
        )
        ttk.Label(
            self.root, textvariable=self.status_var, relief="sunken", anchor=W
        ).pack(fill="x", side="bottom")

    # ---------- Logic ----------
    def choose_folder(self):
        folder = filedialog.askdirectory(title="Select folder with repeat result files")
        if folder:
            self.folder_var.set(folder)
            self.records = {}
            self._clear_table()
            self.status_var.set(f"Selected: {folder}  —  click 'Analyze'")

    def _clear_table(self):
        for item in self.tree.get_children():
            self.tree.delete(item)

    def _get_filters(self):
        try:
            min_len = int(self.min_len_var.get())
        except ValueError:
            messagebox.showerror("Invalid Parameter", "Minimum length must be an integer")
            return None
        max_e = None
        raw = self.max_e_var.get().strip()
        if raw:
            try:
                max_e = float(raw)
            except ValueError:
                messagebox.showerror(
                    "Invalid Parameter", "Max E-value must be a number (e.g. 1e-5)"
                )
                return None
        return min_len, max_e

    def analyze(self):
        folder = self.folder_var.get()
        if not folder or not os.path.isdir(folder):
            messagebox.showerror("Error", "Please select a valid folder first")
            return
        filters = self._get_filters()
        if filters is None:
            return
        min_len, max_e = filters

        self.records = {}
        files = sorted(
            f for f in os.listdir(folder)
            if f.lower().endswith(".txt") and os.path.isfile(os.path.join(folder, f))
        )
        if not files:
            messagebox.showwarning("Warning", "No .txt files found in the selected folder")
            return

        for f in files:
            recs = parse_repeat_file(os.path.join(folder, f))
            recs = [r for r in recs if r[0] >= min_len]
            if max_e is not None:
                recs = [r for r in recs if r[2] <= max_e]
            self.records[species_name_from_file(f)] = recs

        self._fill_table()

    def _fill_table(self):
        self._clear_table()
        self.table_data = []
        sum_d = sum_p = sum_len = 0
        all_lengths = []
        for sp in sorted(self.records):
            recs = self.records[sp]
            n_d = sum(1 for r in recs if r[1] == "D")
            n_p = sum(1 for r in recs if r[1] == "P")
            lengths = [r[0] for r in recs]
            total_bp = sum(lengths)
            row = [
                sp,
                n_d,
                n_p,
                len(recs),
                total_bp,
                round(statistics.mean(lengths), 1) if lengths else 0,
                min(lengths) if lengths else 0,
                max(lengths) if lengths else 0,
            ]
            self.table_data.append(row)
            self.tree.insert("", "end", values=row)
            sum_d += n_d
            sum_p += n_p
            sum_len += total_bp
            all_lengths.extend(lengths)

        if self.table_data:
            total_row = [
                f"TOTAL ({len(self.table_data)} species)",
                sum_d,
                sum_p,
                sum_d + sum_p,
                sum_len,
                round(statistics.mean(all_lengths), 1),
                min(all_lengths),
                max(all_lengths),
            ]
            self.total_row = total_row
            self.tree.insert("", "end", values=total_row, tags=("total",))
            self.tree.tag_configure("total", background="#e8f0fe")
            self.status_var.set(
                f"Done: {len(self.table_data)} species, "
                f"{sum_d} direct repeats, {sum_p} palindromic repeats"
            )
        else:
            self.total_row = None
            self.status_var.set(
                "No repeat records parsed — check filter settings or file contents"
            )

    def sort_by(self, col):
        """Sort by clicking a column header (species alphabetically, others numerically)."""
        if not self.table_data:
            return
        reverse = False
        if self.tree.heading(col)["text"].endswith(" ^"):
            reverse = True
        self.table_data.sort(
            key=lambda row: row[col] if col == 0 else float(row[col]),
            reverse=reverse,
        )
        for i, h in enumerate(self.HEADERS):
            arrow = ""
            if i == col:
                arrow = " v" if not reverse else " ^"
            self.tree.heading(i, text=h + arrow)
        self._clear_table()
        for row in self.table_data:
            self.tree.insert("", "end", values=row)
        if self.total_row:
            self.tree.insert("", "end", values=self.total_row, tags=("total",))

    def export_csv(self):
        if not self.table_data:
            messagebox.showwarning("Warning", "No data yet — run 'Analyze' first")
            return
        default_name = os.path.join(
            self.folder_var.get() or os.getcwd(), "repeat_statistics.csv"
        )
        path = filedialog.asksaveasfilename(
            title="Save Statistics",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile=default_name,
        )
        if not path:
            return
        try:
            with open(path, "w", newline="", encoding="utf-8-sig") as fh:
                writer = csv.writer(fh)
                writer.writerow(self.HEADERS)
                writer.writerows(self.table_data)
                # Summary row
                n_rows = len(self.table_data)
                sum_d = sum(r[1] for r in self.table_data)
                sum_p = sum(r[2] for r in self.table_data)
                sum_bp = sum(r[4] for r in self.table_data)
                writer.writerow(
                    [f"TOTAL ({n_rows} species)", sum_d, sum_p, sum_d + sum_p, sum_bp, "", "", ""]
                )
        except OSError as exc:
            messagebox.showerror("Save Failed", str(exc))
            return
        messagebox.showinfo("Success", f"Statistics saved to:\n{path}")


def main():
    root = tk.Tk()
    try:  # Sharper UI under high DPI (Windows)
        from ctypes import windll

        windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass
    style = ttk.Style()
    if "vista" in style.theme_names():
        style.theme_use("vista")
    RepeatAnalyzerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
