# -*- coding: utf-8 -*-
"""
CodonW2RSCUPlot -- Aggregate codonW RSCU results into the input format
required by the RSCU-Plot web app.

Workflow: codonW (compute per-species RSCU) -> this tool (merge multiple
species into one wide-table CSV) -> upload to
https://pcg-lab.shinyapps.io/RSCU-Plot/ for visualisation.

Input: one RSCU table file per species (MEGA-style .blk wide table or
       codonW-style .out text table). The file name (without extension)
       is used as the species column name.
Output: wide-table CSV
        col 1 = Amino Acid (three-letter abbreviation)
        col 2 = Codon      (RNA letters U/T both accepted)
        col 3+ = RSCU value of each species
"""

import os
import re
import csv
import glob
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


# ---------- Parser ----------

# Codon entry pattern: optional amino-acid abbreviation + triplet codon +
# count + RSCU. Used for both MEGA .blk (4 entries per line) and codonW
# .out (1 entry per line) styles.
CODON_RE = re.compile(r"((?:[A-Za-z]{3})\s+)?([ATUCG]{3})\s+(\d+)\s+([\d.]+)")

# Standard codon ordering (matches the RSCU-Plot template: grouped by amino acid)
STANDARD_ORDER = [
    ("Phe", ["UUU", "UUC"]),
    ("Leu", ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"]),
    ("Ile", ["AUU", "AUC", "AUA"]),
    ("Met", ["AUG"]),
    ("Val", ["GUU", "GUC", "GUA", "GUG"]),
    ("Ser", ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]),
    ("Pro", ["CCU", "CCC", "CCA", "CCG"]),
    ("Thr", ["ACU", "ACC", "ACA", "ACG"]),
    ("Ala", ["GCU", "GCC", "GCA", "GCG"]),
    ("Tyr", ["UAU", "UAC"]),
    ("TER", ["UAA", "UAG", "UGA"]),
    ("His", ["CAU", "CAC"]),
    ("Gln", ["CAA", "CAG"]),
    ("Asn", ["AAU", "AAC"]),
    ("Lys", ["AAA", "AAG"]),
    ("Asp", ["GAU", "GAC"]),
    ("Glu", ["GAA", "GAG"]),
    ("Cys", ["UGU", "UGC"]),
    ("Trp", ["UGG"]),
    ("Arg", ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]),
    ("Gly", ["GGU", "GGC", "GGA", "GGG"]),
]


def _normalize_codon(c):
    """Normalise DNA letter T to RNA letter U so files from different
    sources can be merged uniformly."""
    return c.upper().replace("T", "U")


def _parse_entries_from_text(text):
    """Extract (amino-acid, codon, count, RSCU) entries from text.

    Handles MEGA .blk style (4 entries per line) and codonW .out style
    (1 entry per line: AminoAcid  Codon  Count  RSCU).
    """
    entries = []
    for line in text.splitlines():
        if " codons in " in line:        # skip MEGA summary footer
            continue
        for m in CODON_RE.finditer(line):
            aa = m.group(1).strip() if m.group(1) else None
            codon = _normalize_codon(m.group(2))
            count = int(m.group(3))
            rscu = float(m.group(4))
            entries.append((aa, codon, count, rscu))
    return entries


def parse_input_file(filepath):
    """Read a single RSCU file; return list of (aa, codon, count, RSCU)."""
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    return _parse_entries_from_text(text)


def collect_input_files(input_dir):
    """Collect all RSCU result files in a directory (.blk/.out/.txt/.csv)."""
    patterns = ("*.blk", "*.out", "*.txt", "*.csv")
    seen = set()
    files = []
    for pat in patterns:
        for fp in sorted(glob.glob(os.path.join(input_dir, pat))):
            key = os.path.normcase(os.path.abspath(fp))
            if key not in seen:
                seen.add(key)
                files.append(fp)
    return files


def merge_rscu_files(input_dir, output_csv, include_ter=True,
                      status_callback=None):
    """Aggregate all RSCU result files in a directory into a single CSV
    suitable as RSCU-Plot input.

    Returns (n_species, n_codon_rows, list_of_empty_or_failed_files).
    """
    files = collect_input_files(input_dir)
    if not files:
        raise ValueError("No RSCU result files (.blk/.out/.txt/.csv) found "
                         "in the selected folder!")

    species_data = {}   # species name -> {codon: rscu}
    codon_aa = {}       # codon -> amino acid
    empty_files = []    # files from which no data could be parsed

    for i, fp in enumerate(files, 1):
        name = os.path.splitext(os.path.basename(fp))[0]
        if status_callback:
            status_callback("Parsing {}/{}: {}".format(i, len(files), name))
        try:
            entries = parse_input_file(fp)
        except Exception as e:
            empty_files.append("{} (read error: {})".format(name, e))
            continue
        if not entries:
            empty_files.append(name)
            if status_callback:
                status_callback("Warning: no data parsed from {}".format(name))
            continue
        data = {}
        for aa, codon, count, rscu in entries:
            if not include_ter and (aa == "TER" or codon_aa.get(codon) == "TER"):
                continue
            if aa:
                codon_aa[codon] = aa
            data[codon] = rscu
        species_data[name] = data

    if not species_data:
        raise ValueError("No RSCU data could be parsed from any file; "
                         "please check the input format.")

    # Row order: standard codon table first, unknown codons appended at the end
    known = {c for _, codons in STANDARD_ORDER for c in codons}
    ordered = [(aa, c) for aa, codons in STANDARD_ORDER for c in codons
               if any(c in d for d in species_data.values())]
    for d in species_data.values():
        for c in d:
            if c not in known and c not in [x[1] for x in ordered]:
                ordered.append((codon_aa.get(c, ""), c))

    if status_callback:
        status_callback("Writing aggregated table...")

    with open(output_csv, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(["Amino Acid", "Codon"] + list(species_data.keys()))
        for aa, codon in ordered:
            row = [aa, codon]
            for name in species_data:
                row.append(species_data[name].get(codon, ""))
            writer.writerow(row)

    return len(species_data), len(ordered), empty_files


# ---------- GUI ----------

APP_NAME = "CodonW2RSCUPlot"


class CodonW2RSCUPlotApp:
    def __init__(self, root):
        self.root = root
        self.root.title(APP_NAME + " -- codonW to RSCU-Plot converter")
        self.root.geometry("680x440")
        self.root.minsize(600, 400)

        main = ttk.Frame(self.root, padding=15)
        main.pack(fill=tk.BOTH, expand=True)

        # ---- Input folder ----
        frm_in = ttk.LabelFrame(
            main,
            text="(1) Folder containing codonW / MEGA RSCU files "
                 "(one file per species)",
            padding=10)
        frm_in.pack(fill=tk.X, pady=(0, 10))
        self.var_input = tk.StringVar()
        ttk.Entry(frm_in, textvariable=self.var_input).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 8))
        ttk.Button(frm_in, text="Browse...", command=self.choose_input).pack(
            side=tk.LEFT)

        # ---- Output CSV ----
        frm_out = ttk.LabelFrame(
            main,
            text="(2) Output CSV path (will be uploaded to RSCU-Plot)",
            padding=10)
        frm_out.pack(fill=tk.X, pady=(0, 10))
        self.var_output = tk.StringVar()
        ttk.Entry(frm_out, textvariable=self.var_output).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 8))
        ttk.Button(frm_out, text="Save as...", command=self.choose_output).pack(
            side=tk.LEFT)

        # ---- Options ----
        frm_opt = ttk.LabelFrame(main, text="Options", padding=10)
        frm_opt.pack(fill=tk.X, pady=(0, 10))
        self.var_include_ter = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            frm_opt, text="Include stop codons (TER)",
            variable=self.var_include_ter).pack(anchor=tk.W)

        # ---- Run button ----
        self.btn_run = ttk.Button(main, text="(3) Convert", command=self.run)
        self.btn_run.pack(pady=(0, 10))

        # ---- Log ----
        frm_log = ttk.LabelFrame(main, text="Log", padding=5)
        frm_log.pack(fill=tk.BOTH, expand=True)
        self.txt_log = tk.Text(frm_log, height=10, state=tk.DISABLED,
                               wrap=tk.WORD, font=("Consolas", 9))
        scroll = ttk.Scrollbar(frm_log, command=self.txt_log.yview)
        self.txt_log.configure(yscrollcommand=scroll.set)
        self.txt_log.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

    # ----- UI helpers -----
    def log(self, msg):
        self.txt_log.configure(state=tk.NORMAL)
        self.txt_log.insert(tk.END, msg + "\n")
        self.txt_log.see(tk.END)
        self.txt_log.configure(state=tk.DISABLED)

    def choose_input(self):
        d = filedialog.askdirectory(title="Select RSCU results folder",
                                    initialdir=self.var_input.get() or os.getcwd())
        if d:
            self.var_input.set(d)
            # auto-suggest an output path
            if not self.var_output.get():
                self.var_output.set(os.path.join(d, "rscu_input_for_plot.csv"))

    def choose_output(self):
        f = filedialog.asksaveasfilename(
            title="Save aggregated CSV", defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile="rscu_input_for_plot.csv",
            initialdir=(self.var_input.get() or os.getcwd()))
        if f:
            self.var_output.set(f)

    # ----- Run -----
    def run(self):
        input_dir = self.var_input.get().strip()
        output_csv = self.var_output.get().strip()

        if not input_dir or not os.path.isdir(input_dir):
            messagebox.showerror("Error", "Please select a valid RSCU folder!")
            return
        if not output_csv:
            messagebox.showerror("Error", "Please specify an output CSV path!")
            return
        if os.path.splitext(output_csv)[1].lower() != ".csv":
            output_csv += ".csv"
            self.var_output.set(output_csv)

        self.btn_run.configure(state=tk.DISABLED)
        self.txt_log.configure(state=tk.NORMAL)
        self.txt_log.delete("1.0", tk.END)
        self.txt_log.configure(state=tk.DISABLED)

        def status(msg):
            self.root.after(0, lambda: self.log(msg))

        def worker():
            try:
                result = merge_rscu_files(
                    input_dir, output_csv,
                    include_ter=self.var_include_ter.get(),
                    status_callback=status)
                self.root.after(0, lambda: self.on_done(result, output_csv))
            except Exception as e:
                self.root.after(0, lambda: self.on_error(str(e)))

        threading.Thread(target=worker, daemon=True).start()

    def on_done(self, result, output_csv):
        n_species, n_codons, empty_files = result
        self.log("Done! {} species, {} codon rows.".format(n_species, n_codons))
        if empty_files:
            self.log("Warnings (column kept empty): {}".format(
                ", ".join(empty_files)))
        self.log("Saved to: {}".format(output_csv))
        self.log("Next step: upload this CSV to "
                 "https://pcg-lab.shinyapps.io/RSCU-Plot/")
        self.btn_run.configure(state=tk.NORMAL)
        msg = ("Conversion complete!\n\n"
               "Species: {}\n"
               "Codon rows: {}\n\n"
               "File: {}\n\n"
               "You can now upload it to the RSCU-Plot web app for "
               "visualisation.".format(n_species, n_codons, output_csv))
        if empty_files:
            msg += "\n\nWarning: the following files were empty or could not be parsed:\n{}".format(
                "\n".join(empty_files))
        messagebox.showinfo("Done", msg)

    def on_error(self, msg):
        self.log("Error: " + msg)
        self.btn_run.configure(state=tk.NORMAL)
        messagebox.showerror("Error", msg)


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("vista")
    except tk.TclError:
        pass
    CodonW2RSCUPlotApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()