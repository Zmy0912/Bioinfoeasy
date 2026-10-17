# -*- coding: utf-8 -*-
"""
ENC-GC3s Neutral Expectation Curve Analyzer
===========================================
Addresses reviewer comment:
3- ENC-GC3s analysis does not specify how expected ENC values under neutrality
   were calculated; authors must provide Wright's formula for expected ENC
   calculation in Methods, quantify the proportion of genes falling below/on/above
   the expected curve, and use this to strengthen their conclusion about
   selection dominance.

Features:
1. Select a folder; the program automatically detects every species' .out file
   (codonw output) inside it.
2. Expected ENC under neutrality is calculated with Wright's (1990) formula:
       ENC_exp = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2)
3. Each gene's observed Nc is compared with its expected ENC and classified as:
       below  (Nc_obs < ENC_exp - tol): stronger bias than the neutral
               expectation, suggesting natural selection;
       on     (|Nc_obs - ENC_exp| <= tol): consistent with neutrality;
       above  (Nc_obs > ENC_exp + tol): weaker bias than the neutral
               expectation, consistent with mutation pressure dominance.
4. Results are reported as a per-species table (counts and percentages),
   a per-gene detail table, an ENC-GC3s plot with the expected curve,
   and can be exported to CSV.

Requires: Python 3.7+, matplotlib  (pip install matplotlib)
Run:      python ENC_GC3s_Analyzer.py
"""

import os
import csv
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

matplotlib.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


# ----------------------------------------------------------------------
# Core computation
# ----------------------------------------------------------------------
def expected_enc(gc3s):
    """Wright's (1990) formula for the expected ENC under neutrality."""
    gc3s = min(max(gc3s, 1e-9), 1.0 - 1e-9)
    return 2.0 + gc3s + 29.0 / (gc3s ** 2 + (1.0 - gc3s) ** 2)


def parse_out_file(path):
    """Parse a codonw .out file and return [(gene, nc, gc3s), ...]."""
    records = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()
    # Locate the header line (must contain both Nc and GC3s)
    header_idx, cols = None, []
    for i, line in enumerate(lines):
        if "Nc" in line and "GC3s" in line:
            header_idx = i
            cols = line.rstrip("\r\n").rstrip().split("\t")
            break
    if header_idx is None:
        return records
    try:
        nc_col = next(i for i, c in enumerate(cols) if c.strip().lower() in ("nc", "enc"))
        gc_col = next(i for i, c in enumerate(cols) if c.strip().lower() == "gc3s")
    except StopIteration:
        return records
    for line in lines[header_idx + 1:]:
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) <= max(nc_col, gc_col):
            continue
        nc_raw = parts[nc_col].strip()
        gc_raw = parts[gc_col].strip()
        try:
            nc = float(nc_raw)
            gc3s = float(gc_raw)
        except ValueError:
            continue  # skip invalid values such as "*****"
        gene = parts[0].split("_")[0].strip()
        if not gene:
            gene = parts[0].strip()[:20]
        records.append((gene, nc, gc3s))
    return records


def classify_genes(records, tol_mode="relative", tol=0.05):
    """Classify genes into below / on / above the expected curve."""
    below, on, above = [], [], []
    for gene, nc, gc3s in records:
        exp = expected_enc(gc3s)
        delta = tol * exp if tol_mode == "relative" else tol
        if nc < exp - delta:
            below.append((gene, nc, gc3s, exp))
        elif nc > exp + delta:
            above.append((gene, nc, gc3s, exp))
        else:
            on.append((gene, nc, gc3s, exp))
    return below, on, above


def analyze_folder(folder, tol_mode="relative", tol=0.05):
    """Scan a folder and return {species: {summary..., records...}}."""
    results = {}
    for fname in sorted(os.listdir(folder)):
        if not fname.lower().endswith(".out"):
            continue
        species = os.path.splitext(fname)[0]
        records = parse_out_file(os.path.join(folder, fname))
        below, on, above = classify_genes(records, tol_mode, tol)
        n = len(records)
        results[species] = {
            "file": fname,
            "total": n,
            "below": below, "on": on, "above": above,
            "n_below": len(below), "n_on": len(on), "n_above": len(above),
            "p_below": len(below) / n * 100 if n else 0.0,
            "p_on": len(on) / n * 100 if n else 0.0,
            "p_above": len(above) / n * 100 if n else 0.0,
            "records": records,
        }
    return results


# ----------------------------------------------------------------------
# GUI
# ----------------------------------------------------------------------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("ENC-GC3s Neutral Expectation Curve Analyzer (Wright, 1990)")
        self.geometry("1180x760")
        self.results = {}
        self.folder = tk.StringVar()

        self._build_ui()

    # ---------------- UI ----------------
    def _build_ui(self):
        top = ttk.Frame(self, padding=8)
        top.pack(fill="x")

        ttk.Label(top, text="Input folder:").grid(row=0, column=0, sticky="w")
        ttk.Entry(top, textvariable=self.folder, width=70).grid(row=0, column=1, padx=5)
        ttk.Button(top, text="Browse...", command=self.choose_folder).grid(row=0, column=2)

        ttk.Label(top, text="Tolerance:").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.tol_mode = tk.StringVar(value="relative")
        cb = ttk.Combobox(top, textvariable=self.tol_mode, width=14,
                          values=("relative", "absolute"), state="readonly")
        cb.grid(row=1, column=1, sticky="w", pady=(6, 0))
        self.tol = tk.DoubleVar(value=0.05)
        ttk.Spinbox(top, from_=0.001, to=1.0, increment=0.01, width=8,
                    textvariable=self.tol).grid(row=1, column=1, sticky="e", padx=(0, 90), pady=(6, 0))
        ttk.Label(top, text="(relative = % of expected ENC; absolute = ENC units)").grid(
            row=1, column=2, sticky="w", pady=(6, 0))

        ttk.Button(top, text="Run analysis", command=self.run_analysis).grid(
            row=2, column=0, columnspan=2, sticky="w", pady=(8, 0))
        ttk.Button(top, text="Export summary (CSV)", command=self.export_csv).grid(
            row=2, column=1, sticky="e", pady=(8, 0))
        ttk.Button(top, text="Export gene detail (CSV)", command=self.export_detail).grid(
            row=2, column=2, sticky="w", pady=(8, 0))

        self.status = tk.StringVar(value="Select an input folder and click 'Run analysis'.")
        ttk.Label(self, textvariable=self.status, relief="sunken",
                  anchor="w", padding=4).pack(fill="x", side="bottom")

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=8)

        # --- Summary table ---
        f1 = ttk.Frame(nb)
        nb.add(f1, text="Species summary")
        cols = ("species", "total", "n_below", "p_below", "n_on", "p_on",
                "n_above", "p_above")
        heads = ("Species / file", "Genes", "Below (n)", "Below (%)",
                 "On (n)", "On (%)", "Above (n)", "Above (%)")
        self.tree = ttk.Treeview(f1, columns=cols, show="headings", height=18)
        widths = (110, 70, 90, 90, 80, 90, 90, 90)
        for c, h, w in zip(cols, heads, widths):
            self.tree.heading(c, text=h)
            self.tree.column(c, width=w, anchor="center")
        vsb = ttk.Scrollbar(f1, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.tree.bind("<<TreeviewSelect>>", lambda e: self.draw_plot())

        # --- Gene detail table ---
        f2 = ttk.Frame(nb)
        nb.add(f2, text="Gene detail")
        dcols = ("species", "gene", "gc3s", "nc_obs", "nc_exp", "diff", "class")
        dheads = ("Species", "Gene", "GC3s", "Observed Nc", "Expected ENC",
                  "Nc - ENCexp", "Class")
        self.detail = ttk.Treeview(f2, columns=dcols, show="headings", height=24)
        dwidths = (90, 150, 80, 90, 100, 90, 80)
        for c, h, w in zip(dcols, dheads, dwidths):
            self.detail.heading(c, text=h)
            self.detail.column(c, width=w, anchor="center")
        dvsb = ttk.Scrollbar(f2, orient="vertical", command=self.detail.yview)
        self.detail.configure(yscrollcommand=dvsb.set)
        self.detail.pack(side="left", fill="both", expand=True)
        dvsb.pack(side="right", fill="y")

        # --- Plot ---
        f3 = ttk.Frame(nb)
        nb.add(f3, text="ENC-GC3s plot")
        self.fig = Figure(figsize=(8, 5.5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=f3)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # --- Methods note ---
        f4 = ttk.Frame(nb)
        nb.add(f4, text="Methods")
        txt = tk.Text(f4, wrap="word", height=30)
        txt.pack(fill="both", expand=True)
        txt.insert("1.0",
            "Methods (ready to paste into the manuscript, in reply to comment 3)\n"
            "=================================================================\n\n"
            "The expected ENC under neutrality (no selection) was calculated with "
            "Wright's (1990) formula:\n"
            "    ENC_exp = 2 + GC3s + 29 / [GC3s^2 + (1 - GC3s)^2]\n"
            "where GC3s is the GC content at synonymous codon third positions "
            "(computed with codonw).\n\n"
            "Classification criteria:\n"
            "  - Below the curve (below): Nc_obs < ENC_exp - tol. Codon usage bias is "
            "stronger than expected from mutation pressure alone, indicating that "
            "natural selection contributes to the codon usage bias of these genes.\n"
            "  - On the curve (on): |Nc_obs - ENC_exp| <= tol. Codon usage is "
            "consistent with the neutral expectation.\n"
            "  - Above the curve (above): Nc_obs > ENC_exp + tol. Codon usage bias is "
            "weaker than the neutral expectation, also consistent with mutation "
            "pressure dominance.\n"
            "The tolerance tol can be set in the GUI: in 'relative' mode it is a "
            "percentage of the expected ENC (default 5%); in 'absolute' mode it is an "
            "absolute ENC difference.\n\n"
            "Suggested conclusion wording: if the majority of genes lie on or above "
            "the expected curve and only a minority fall significantly below it, the "
            "results support the conclusion that mutation pressure is the dominant "
            "force shaping codon usage, with natural selection acting on a subset of "
            "genes; conversely, a majority of genes below the curve would indicate "
            "that selection dominates.\n\n"
            "Reference:\n"
            "Wright F. 1990. The 'effective number of codons' used in a gene. "
            "Gene 87:23-29.\n")
        txt.config(state="disabled")

    # ---------------- actions ----------------
    def choose_folder(self):
        d = filedialog.askdirectory(initialdir=self.folder.get() or "/")
        if d:
            self.folder.set(d)

    def run_analysis(self):
        folder = self.folder.get()
        if not os.path.isdir(folder):
            messagebox.showerror("Error", "Folder not found: %s" % folder)
            return
        tol_mode = self.tol_mode.get()
        tol = self.tol.get()
        self.status.set("Analyzing...")

        def work():
            res = analyze_folder(folder, tol_mode, tol)
            self.after(0, lambda: self._show(res, tol_mode, tol))
        threading.Thread(target=work, daemon=True).start()

    def _show(self, res, tol_mode, tol):
        self.results = res
        self.tol_used = (tol_mode, tol)
        self.tree.delete(*self.tree.get_children())
        self.detail.delete(*self.detail.get_children())
        if not res:
            self.status.set("No .out files found in the selected folder.")
            return
        tot_b = tot_o = tot_a = tot_n = 0
        for sp, r in sorted(res.items()):
            self.tree.insert("", "end", iid=sp, values=(
                sp, r["total"],
                r["n_below"], "%.2f" % r["p_below"],
                r["n_on"], "%.2f" % r["p_on"],
                r["n_above"], "%.2f" % r["p_above"]))
            tot_b += r["n_below"]; tot_o += r["n_on"]
            tot_a += r["n_above"]; tot_n += r["total"]
            for grp, cls in ((r["below"], "below"),
                             (r["on"], "on"),
                             (r["above"], "above")):
                for gene, nc, gc3s, exp in grp:
                    self.detail.insert("", "end", values=(
                        sp, gene, "%.3f" % gc3s, "%.2f" % nc, "%.2f" % exp,
                        "%+.2f" % (nc - exp), cls))
        if tot_n:
            self.tree.insert("", "end", values=(
                "TOTAL", tot_n, tot_b, "%.2f" % (tot_b / tot_n * 100),
                tot_o, "%.2f" % (tot_o / tot_n * 100),
                tot_a, "%.2f" % (tot_a / tot_n * 100)))
        self.status.set("Done: %d species, %d valid genes (tolerance: %s=%s)."
                        % (len(res), tot_n, tol_mode, tol))
        self.draw_plot()

    # ---------------- plot ----------------
    def draw_plot(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        xs = [i / 200.0 for i in range(201)]
        ax.plot(xs, [expected_enc(x) for x in xs], "k-", lw=1.5,
                label="Expected ENC (Wright, 1990)")
        sel = self.tree.selection()
        if sel and sel[0] in self.results:
            sp = sel[0]
            r = self.results[sp]
            for grp, color, lab in ((r["below"], "#d62728", "Below (selection)"),
                                    (r["on"], "#7f7f7f", "On curve (neutral)"),
                                    (r["above"], "#1f77b4", "Above (mutation)")):
                if grp:
                    ax.scatter([g[2] for g in grp], [g[1] for g in grp],
                               c=color, s=28, alpha=0.8, label=lab)
            ax.set_title("%s  (below %d / on %d / above %d)"
                         % (sp, r["n_below"], r["n_on"], r["n_above"]))
        else:
            for r in self.results.values():
                ax.scatter([g[2] for g in r["records"]],
                           [g[1] for g in r["records"]],
                           s=14, alpha=0.5)
            ax.set_title("All species (select a species row in the summary table "
                         "to view it individually)")
        ax.set_xlabel("GC3s")
        ax.set_ylabel("ENC (Nc)")
        ax.set_xlim(0, 1)
        ax.set_ylim(20, 62)
        ax.legend(loc="lower center", fontsize=8)
        ax.grid(alpha=0.3)
        self.canvas.draw()

    # ---------------- export ----------------
    def export_csv(self):
        if not self.results:
            messagebox.showinfo("Info", "Please run the analysis first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            initialfile="ENC_GC3s_summary.csv",
            filetypes=[("CSV files", "*.csv")])
        if not path:
            return
        tol_mode, tol = getattr(self, "tol_used", ("relative", 0.05))
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(["Wright(1990) expected ENC: ENC_exp = 2 + GC3s + 29/(GC3s^2+(1-GC3s)^2)"])
            w.writerow(["Tolerance: %s = %s" % (tol_mode, tol)])
            w.writerow([])
            w.writerow(["Species", "Total genes", "Below (n)", "Below (%)",
                        "On (n)", "On (%)", "Above (n)", "Above (%)"])
            tot = [0] * 4
            for sp, r in sorted(self.results.items()):
                w.writerow([sp, r["total"], r["n_below"], "%.2f" % r["p_below"],
                            r["n_on"], "%.2f" % r["p_on"],
                            r["n_above"], "%.2f" % r["p_above"]])
                tot[0] += r["total"]; tot[1] += r["n_below"]
                tot[2] += r["n_on"]; tot[3] += r["n_above"]
            if tot[0]:
                w.writerow(["TOTAL", tot[0], tot[1], "%.2f" % (tot[1] / tot[0] * 100),
                            tot[2], "%.2f" % (tot[2] / tot[0] * 100),
                            tot[3], "%.2f" % (tot[3] / tot[0] * 100)])
        messagebox.showinfo("Done", "Summary saved to:\n%s" % path)

    def export_detail(self):
        if not self.results:
            messagebox.showinfo("Info", "Please run the analysis first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            initialfile="ENC_GC3s_gene_detail.csv",
            filetypes=[("CSV files", "*.csv")])
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(["Species", "Gene", "GC3s", "Observed_Nc",
                        "Expected_ENC", "Nc_minus_expected", "Class"])
            for sp, r in sorted(self.results.items()):
                for grp, cls in ((r["below"], "below"),
                                 (r["on"], "on"),
                                 (r["above"], "above")):
                    for gene, nc, gc3s, exp in grp:
                        w.writerow([sp, gene, "%.4f" % gc3s, "%.2f" % nc,
                                    "%.2f" % exp, "%+.2f" % (nc - exp), cls])
        messagebox.showinfo("Done", "Gene detail saved to:\n%s" % path)


if __name__ == "__main__":
    app = App()
    app.mainloop()
