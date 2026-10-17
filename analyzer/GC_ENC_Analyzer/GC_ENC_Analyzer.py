# -*- coding: utf-8 -*-
"""
GC / ENC Codon Usage Bias Analyzer
====================================
Input a folder; each FASTA file inside is treated as one species
(species name = file name), and each sequence record within a file is
one gene (gene names are automatically extracted from the header line,
with coordinates, strand direction, "join" and other non-gene-name parts
removed).

For each CDS gene of each species, the tool computes:
    GC1 (%) / GC2 (%) / GC3 (%)   -- GC content at codon positions 1/2/3
    GC3s (%)                      -- GC content at the 3rd position of
                                     synonymous codons (CodonW convention)
    GC12 (%)                      -- average GC content of positions 1 and 2
    GC total (%)                  -- overall GC content of the sequence
    ENC                           -- effective number of codons, computed
                                     exactly like CodonW 1.4.4

A terminal stop codon is excluded from all statistics (CodonW
convention).

Three summary tables are produced:
    Table 1  All species x all genes (one row per species-gene combination,
             raw values)
    Table 2  All species (one row per species; values are the mean ± SD
             across all genes of that species; GC12-GC3 slope and R² come
             from the regression of GC12 on GC3 across genes)
    Table 3  All genes (one row per gene; values are the mean ± SD across
             all species; GC12-GC3 slope and R² come from the regression of
             GC12 on GC3 across species)

Outputs are written as CSV (UTF-8 with BOM, opens directly in Excel);
if openpyxl is installed, an xlsx version is also produced. Results are
saved in the "results" subfolder of the selected folder.

Dependencies: Python standard library only (tkinter); openpyxl optional.
"""

import os
import csv
import re
import queue
import threading
import traceback
from statistics import mean, stdev

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# ---------------------------------------------------------------------------
# Genetic code (standard codon table)
# ---------------------------------------------------------------------------
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

CODON2AA = {}
_i = 0
for _b1 in _BASES:
    for _b2 in _BASES:
        for _b3 in _BASES:
            CODON2AA[_b1 + _b2 + _b3] = _AAS[_i]
            _i += 1

# Group codons into amino-acid families (excluding Met, Trp, stop codons)
AA_FAMILIES = {}  # aa -> [codon, ...]
for _c, _a in CODON2AA.items():
    if _a in ("*", "M", "W"):
        continue
    AA_FAMILIES.setdefault(_a, []).append(_c)

DEGENERACY = {aa: len(cods) for aa, cods in AA_FAMILIES.items()}
# Number of amino-acid families per degeneracy class in the genetic code
# (Wright's coefficients):  ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
ENC_COEF = {2: 9.0, 3: 1.0, 4: 5.0, 6: 3.0}

# Codons belonging to degenerate (synonymous) families, i.e. all codons
# except Met, Trp and stop codons. CodonW's GC3s and L_sym use exactly
# this codon set.
SYN_CODONS = set()
for _aa, _cods in AA_FAMILIES.items():
    SYN_CODONS.update(_cods)

STOPS = {"TAA", "TAG", "TGA"}

VALID_BASES = set("ATCG")


# ---------------------------------------------------------------------------
# Gene name extraction
# ---------------------------------------------------------------------------
_GENE_SEP_RE = re.compile(r"(_\[|\[|_join\{|_complement|\||\s|\()")


def extract_gene_name(header: str) -> str:
    """Extract the gene name from a FASTA header, removing coordinates,
    strand direction, "join" and other non-gene-name parts.

    Examples:
        >psbD_[33312:34373](+)|NC01.fasta|NC01.fasta  ->  psbD
        >rps16_join{[5754:5793](-),[4686:4894](-)}|.. ->  rps16
        >rbcL_[56315:57757](+)                        ->  rbcL
    """
    name = header.strip()
    m = _GENE_SEP_RE.search(name)
    if m:
        name = name[: m.start()]
    name = name.strip().strip("_-")
    # Keep only a part that starts with a letter and consists of
    # alphanumeric characters and underscores
    m2 = re.match(r"^[A-Za-z][A-Za-z0-9_]*", name)
    if not m2:
        name = re.sub(r"[^A-Za-z0-9_]", "", name)
    return name


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------
def parse_fasta(path: str):
    """Minimal FASTA parser; returns [(header, seq), ...]"""
    records = []
    header, buf = None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(buf)))
                header, buf = line[1:], []
            else:
                buf.append(line)
    if header is not None:
        records.append((header, "".join(buf)))
    return records


# ---------------------------------------------------------------------------
# Sequence-level computations
# ---------------------------------------------------------------------------
def seq_to_codons(seq: str):
    """Clean the sequence and split it into codons (keep ACGT only;
    drop any trailing incomplete codon)."""
    seq = "".join(c for c in seq.upper() if c in VALID_BASES)
    n = len(seq) - len(seq) % 3
    return [seq[i:i + 3] for i in range(0, n, 3)]


def calc_gc_positions(codons):
    """Return (gc1, gc2, gc3): GC percentage at each codon position."""
    if not codons:
        return None, None, None
    c1 = c2 = c3 = 0
    for codon in codons:
        if codon[0] in "GC":
            c1 += 1
        if codon[1] in "GC":
            c2 += 1
        if codon[2] in "GC":
            c3 += 1
    n = len(codons)
    return c1 / n * 100, c2 / n * 100, c3 / n * 100


def calc_gc3s(codons):
    """GC percentage at the 3rd position of synonymous (degenerate)
    codons -- identical to CodonW's GC3s (Met, Trp and stop codons are
    excluded from both numerator and denominator)."""
    gc = tot = 0
    for codon in codons:
        if codon in SYN_CODONS:
            tot += 1
            if codon[2] in "GC":
                gc += 1
    return (gc / tot * 100) if tot else None


def calc_enc(codons):
    """Effective number of codons, ENC -- faithful re-implementation of
    CodonW 1.4.4 (enc() in codon_idx.c), which differs from a plain
    Wright (1990) implementation:

    * families observed only once (n <= 1) are excluded from the
      per-class homozygosity average;
    * families whose homozygosity bb is ~0 (perfectly even codon use)
      are also excluded (bb must be > 1e-7);
    * if the 3-fold class (Ile) has no valid family, the average of the
      2-fold and 4-fold class averages is used instead;
    * if any other class has no valid family, ENC is undefined (None,
      printed as '*****' by CodonW);
    * the result is only clamped at the upper bound 61 (no lower clamp).
    """
    counts = {}
    for c in codons:
        counts[c] = counts.get(c, 0) + 1

    totb = {2: 0.0, 3: 0.0, 4: 0.0, 6: 0.0}   # sum of bb per class
    numaa = {2: 0, 3: 0, 4: 0, 6: 0}          # valid families per class
    for aa, codon_list in AA_FAMILIES.items():
        n = sum(counts.get(c, 0) for c in codon_list)
        if n <= 1:
            bb = 0.0
        else:
            p2 = sum((counts.get(c, 0) / n) ** 2 for c in codon_list)
            bb = (n * p2 - 1.0) / (n - 1.0)
        if bb > 1e-7:
            k = DEGENERACY[aa]
            totb[k] += bb
            numaa[k] += 1

    enc = 2.0  # fold[1]: Met + Trp
    for z in (2, 3, 4, 6):
        if numaa[z] and totb[z] > 0:
            averb = totb[z] / numaa[z]
        elif z == 3 and numaa[2] and numaa[4]:
            # CodonW special case: no valid Ile family -> use the mean
            # of the 2-fold and 4-fold class averages
            averb = (totb[2] / numaa[2] + totb[4] / numaa[4]) * 0.5
        else:
            return None  # CodonW prints "*****" in this case
        enc += ENC_COEF[z] / averb
        if enc > 61.0:
            enc = 61.0
    return enc


def gene_metrics(seq: str):
    """Compute all metrics for one gene; returns a dict (values may be None).

    Following the CodonW convention, a terminal stop codon
    (TAA/TAG/TGA) is excluded from all statistics (internal stop
    codons, if any, are kept)."""
    codons = seq_to_codons(seq)
    if not codons:
        return None
    if codons[-1] in STOPS:
        codons = codons[:-1]
        if not codons:
            return None
    gc1, gc2, gc3 = calc_gc_positions(codons)
    gc3s = calc_gc3s(codons)
    gc12 = (gc1 + gc2) / 2.0
    allbases = "".join(codons)
    gc_all = sum(1 for b in allbases if b in "GC") / len(allbases) * 100
    enc = calc_enc(codons)
    return {
        "length": len(codons) * 3,
        "ncodons": len(codons),
        "GC1": gc1,
        "GC2": gc2,
        "GC3": gc3,
        "GC3s": gc3s,
        "GC12": gc12,
        "GCall": gc_all,
        "ENC": enc,
    }


# ---------------------------------------------------------------------------
# Regression
# ---------------------------------------------------------------------------
def linreg_slope_r2(xs, ys):
    """Least-squares linear regression of y on x; returns (slope, R²).
    Returns (None, None) when fewer than 3 points are available."""
    n = len(xs)
    if n < 3:
        return None, None
    mx, my = mean(xs), mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None, None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    syy = sum((y - my) ** 2 for y in ys)
    if syy == 0:
        r2 = 1.0 if slope != 0 else None
    else:
        r2 = (sxy * sxy) / (sxx * syy)
    return slope, r2


# ---------------------------------------------------------------------------
# Number formatting
# ---------------------------------------------------------------------------
def fmt_val(v, nd=1):
    return "" if v is None else f"{v:.{nd}f}"


def fmt_mean_std(values, nd=1):
    """Format as 'mean ± std'; show only the mean when there is a single value."""
    vs = [v for v in values if v is not None]
    if not vs:
        return ""
    if len(vs) == 1:
        return f"{vs[0]:.{nd}f}"
    return f"{mean(vs):.{nd}f} ± {stdev(vs):.{nd}f}"


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------
def process_folder(folder: str, log=print, progress=None):
    """Process a folder; returns (detail rows, per-species rows, per-gene rows,
    number of files, number of species-gene combinations)."""
    fasta_files = sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith((".fasta", ".fa", ".faa", ".fna", ".fst"))
    )
    if not fasta_files:
        raise RuntimeError(
            "No FASTA files (.fasta/.fa/.fna, etc.) were found in the "
            "selected folder.")

    detail_rows = []          # each row: (species, gene, metrics)
    gene_order = []           # order of first appearance of each gene

    for idx, fname in enumerate(fasta_files, 1):
        species = os.path.splitext(fname)[0]
        path = os.path.join(folder, fname)
        try:
            records = parse_fasta(path)
        except Exception as e:
            log(f"[Skipped] Cannot read {fname}: {e}")
            continue

        ngene = 0
        for header, seq in records:
            gene = extract_gene_name(header)
            if not gene:
                log(f"[Warning] {fname}: cannot recognize a gene name "
                    f"(header: {header[:60]}); record skipped")
                continue
            m = gene_metrics(seq)
            if m is None or m["ncodons"] < 3:
                log(f"[Warning] {fname} / {gene}: sequence too short or "
                    f"unparseable; record skipped")
                continue
            detail_rows.append((species, gene, m))
            if gene not in gene_order:
                gene_order.append(gene)
            ngene += 1

        log(f"({idx}/{len(fasta_files)}) {fname}: {species} -- {ngene} genes")
        if progress:
            progress(idx, len(fasta_files))

    # -------- Table 1: all species x all genes (detail) --------
    table1 = []
    for species, gene, m in detail_rows:
        table1.append([
            species, gene, m["length"], m["ncodons"],
            fmt_val(m["GC1"]), fmt_val(m["GC2"]), fmt_val(m["GC3"]),
            fmt_val(m["GC3s"]), fmt_val(m["GC12"]), fmt_val(m["GCall"]),
            fmt_val(m["ENC"]),
        ])

    # -------- Table 2: one row per species (averaged over genes) --------
    by_species = {}
    for species, gene, m in detail_rows:
        by_species.setdefault(species, []).append(m)

    species_names = sorted(by_species)
    table2 = []
    for species in species_names:
        ms = by_species[species]
        xs = [m["GC3"] for m in ms if m["GC3"] is not None]
        ys = [m["GC12"] for m in ms
              if m["GC3"] is not None and m["GC12"] is not None]
        slope, r2 = linreg_slope_r2(xs, ys)
        table2.append([
            species, len(ms),
            fmt_mean_std([m["GC1"] for m in ms]),
            fmt_mean_std([m["GC2"] for m in ms]),
            fmt_mean_std([m["GC3"] for m in ms]),
            fmt_mean_std([m["GC3s"] for m in ms]),
            fmt_mean_std([m["GC12"] for m in ms]),
            fmt_mean_std([m["ENC"] for m in ms]),
            fmt_val(slope, 3), fmt_val(r2, 3),
        ])

    # -------- Table 3: one row per gene (averaged over species) --------
    by_gene = {}
    for species, gene, m in detail_rows:
        by_gene.setdefault(gene, []).append(m)

    table3 = []
    for gene in gene_order:
        ms = by_gene[gene]
        xs = [m["GC3"] for m in ms if m["GC3"] is not None]
        ys = [m["GC12"] for m in ms
              if m["GC3"] is not None and m["GC12"] is not None]
        slope, r2 = linreg_slope_r2(xs, ys)
        table3.append([
            gene, len(ms),
            fmt_mean_std([m["GC1"] for m in ms]),
            fmt_mean_std([m["GC2"] for m in ms]),
            fmt_mean_std([m["GC3"] for m in ms]),
            fmt_mean_std([m["GC3s"] for m in ms]),
            fmt_mean_std([m["GC12"] for m in ms]),
            fmt_mean_std([m["ENC"] for m in ms]),
            fmt_val(slope, 3), fmt_val(r2, 3),
        ])

    return table1, table2, table3, len(fasta_files), len(detail_rows)


HEADER1 = ["Species", "Gene", "CDS length (bp)", "Codons",
           "GC1 (%)", "GC2 (%)", "GC3 (%)", "GC3s (%)", "GC12 (%)",
           "GC total (%)", "ENC"]
HEADER2 = ["Species", "Genes",
           "GC1 (%)", "GC2 (%)", "GC3 (%)", "GC3s (%)", "GC12 (%)", "ENC",
           "GC12-GC3 slope", "R2"]
HEADER3 = ["Gene", "Species",
           "GC1 (%)", "GC2 (%)", "GC3 (%)", "GC3s (%)", "GC12 (%)", "ENC",
           "GC12-GC3 slope", "R2"]


def save_results(out_dir, table1, table2, table3):
    """Save the three CSV tables (plus an xlsx file if openpyxl is installed)."""
    os.makedirs(out_dir, exist_ok=True)

    p1 = os.path.join(out_dir, "1_all_species_all_genes_detail.csv")
    p2 = os.path.join(out_dir, "2_species_averages.csv")
    p3 = os.path.join(out_dir, "3_gene_averages.csv")

    for path, header, rows in ((p1, HEADER1, table1),
                               (p2, HEADER2, table2),
                               (p3, HEADER3, table3)):
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)

    # Optional xlsx output
    xlsx_path = None
    try:
        from openpyxl import Workbook
        xlsx_path = os.path.join(out_dir, "GC_ENC_summary.xlsx")
        wb = Workbook()
        for sheet, header, rows in (
            ("All species all genes", HEADER1, table1),
            ("Species averages", HEADER2, table2),
            ("Gene averages", HEADER3, table3),
        ):
            ws = wb.active if sheet == "All species all genes" \
                else wb.create_sheet()
            ws.title = sheet
            ws.append(header)
            for r in rows:
                ws.append(r)
        wb.save(xlsx_path)
    except ImportError:
        xlsx_path = None
    except Exception:
        xlsx_path = None

    return [p1, p2, p3] + ([xlsx_path] if xlsx_path else [])


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("GC / ENC Codon Usage Bias Analyzer")
        self.geometry("780x560")
        self.minsize(680, 460)
        self._queue = queue.Queue()
        self.after(100, self._poll_queue)

        pad = {"padx": 10, "pady": 6}

        top = ttk.Frame(self)
        top.pack(fill="x", **pad)
        ttk.Label(top, text="Sequence folder:").pack(side="left")
        self.var_folder = tk.StringVar()
        ttk.Entry(top, textvariable=self.var_folder).pack(
            side="left", fill="x", expand=True, padx=6)
        ttk.Button(top, text="Browse...",
                   command=self.choose_folder).pack(side="left")
        ttk.Button(top, text="Start analysis",
                   command=self.start).pack(side="left", padx=6)

        self.progress = ttk.Progressbar(self, mode="determinate")
        self.progress.pack(fill="x", padx=10)

        ttk.Label(
            self,
            text="Notes: each FASTA file = one species (species name is the file "
                 "name); each record inside a file = one gene, whose name is "
                 "automatically extracted from the header (coordinates, strand "
                 "direction, 'join', etc. are removed).\n"
                 "GC12-GC3 slope / R2: slope and coefficient of determination "
                 "of the linear regression of GC12 on GC3; in the summary "
                 "tables each metric is reported as mean ± standard deviation.",
            foreground="#555", justify="left", wraplength=740,
        ).pack(fill="x", padx=10)

        self.log_text = tk.Text(self, height=18, state="disabled",
                                font=("Consolas", 9))
        self.log_text.pack(fill="both", expand=True, padx=10, pady=8)

        bottom = ttk.Frame(self)
        bottom.pack(fill="x", padx=10, pady=(0, 8))
        self.var_status = tk.StringVar(value="Ready")
        ttk.Label(bottom, textvariable=self.var_status,
                  foreground="#0a7").pack(side="left")
        ttk.Button(bottom, text="Open results folder",
                   command=self.open_outdir).pack(side="right")

        self.out_dir = None
        self._running = False

    # ---- Logging (thread-safe, forwarded to the main thread via a queue) ----
    def log(self, msg):
        self._queue.put(("log", msg))

    def _poll_queue(self):
        try:
            while True:
                item = self._queue.get_nowait()
                if item[0] == "log":
                    self.log_text.configure(state="normal")
                    self.log_text.insert("end", item[1] + "\n")
                    self.log_text.see("end")
                    self.log_text.configure(state="disabled")
                elif item[0] == "prog":
                    cur, total = item[1], item[2]
                    self.progress.configure(maximum=total, value=cur)
                elif item[0] == "done":
                    self.var_status.set(item[1])
                    self._running = False
        except queue.Empty:
            pass
        self.after(100, self._poll_queue)

    # ---- Interaction ----
    def choose_folder(self):
        d = filedialog.askdirectory(
            title="Select the folder containing species FASTA files")
        if d:
            self.var_folder.set(os.path.normpath(d))

    def open_outdir(self):
        if self.out_dir and os.path.isdir(self.out_dir):
            os.startfile(self.out_dir)

    def start(self):
        if self._running:
            return
        folder = self.var_folder.get().strip()
        if not folder or not os.path.isdir(folder):
            messagebox.showwarning(
                "Notice", "Please select a valid sequence folder first!")
            return
        self._running = True
        self.var_status.set("Analyzing...")
        self.progress.configure(maximum=100, value=0)
        threading.Thread(target=self._work, args=(folder,), daemon=True).start()

    def _work(self, folder):
        try:
            self.log("=" * 60)
            self.log(f"Input folder: {folder}")
            t1, t2, t3, nfile, npair = process_folder(
                folder, log=self.log, progress=self._progress)
            out_dir = os.path.join(folder, "results")
            paths = save_results(out_dir, t1, t2, t3)
            self.out_dir = out_dir
            self.log("-" * 60)
            self.log(f"Done! {nfile} species files, {npair} species-gene "
                     f"combinations in total.")
            for p in paths:
                self.log(f"  Output: {p}")
            self._queue.put(
                ("done", f"Done: {nfile} species / {npair} species-gene "
                         f"combinations. Results saved to {out_dir}"))
        except Exception as e:
            err = "".join(traceback.format_exception_only(type(e), e)).strip()
            self.log(f"[Error] {err}")
            self._queue.put(("done", f"Analysis failed: {err}"))

    def _progress(self, cur, total):
        self._queue.put(("prog", cur, total))


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
