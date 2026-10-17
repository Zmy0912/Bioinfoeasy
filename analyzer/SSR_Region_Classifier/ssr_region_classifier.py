#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SSR Region Classifier
=====================

Overview
--------
This program reads a MISA SSR result file together with the GenBank (GB)
annotation files of the corresponding species, and classifies every SSR
according to its genomic coordinate into one of three region types:

  - CDS    : coding region (CDS exons)
  - Intron : introns inferred from the gene span minus the CDS exons
  - IGS    : intergenic spacer (neither CDS nor Intron)

The output files can be loaded directly by Excel, Origin, R or Python
(matplotlib/seaborn) for plotting or further statistics.

Input Requirements
------------------
1. MISA file:
   - Tab-separated. The header must contain at least:
     ID, SSR nr., SSR type, SSR, size, start, end
   - The first column (ID, e.g. NC01, NC02, ...) must match the basename
     of the GB file (NC01.gb, NC02.gb, ...).

2. GB folder:
   - A folder (usually named the same as the MISA file) that contains
     one GB file per species, e.g. NC01.gb, NC02.gb, ...
   - Each GB file must contain ``gene`` and ``CDS`` features in its
     FEATURES table.

Outputs
-------
1. <out>/SSR_with_region.txt        : MISA contents + Region / RegionDetail columns
2. <out>/SSR_region_summary.tsv     : per-species region counts and percentages (long, ggplot-friendly)
3. <out>/SSR_region_summary_wide.tsv: per-species region counts (wide)
4. <out>/SSR_type_by_region.tsv     : region x SSR-type counts (good for stacked bar charts)
5. <out>/log.txt                    : run log

Usage
-----
Double-click the launcher, or run from the command line:

    python ssr_region_classifier.py
"""

import os
import re
import csv
import sys
import time
import threading
from collections import defaultdict, OrderedDict

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext


# ---------------------------------------------------------------------------
# 1. GenBank file parsing
# ---------------------------------------------------------------------------

_LOCATION_RE = re.compile(r"(\d+)\.\.(\d+)")


def parse_location(loc: str):
    """
    Parse a GenBank Location string into a list of ``(start, end)`` intervals.

    Handles single-segment, ``join(...)`` and combinations, including the
    ``complement(...)`` wrapper.

    The returned intervals are always sorted in ascending order; the
    complement information is returned separately via a boolean flag.
    """
    raw = loc.strip()
    is_complement = "complement" in raw

    # Strip complement(...) / join(...) wrappers (possibly nested).
    s = raw
    while True:
        changed = False
        if s.startswith("complement(") and s.endswith(")"):
            s = s[len("complement("):-1].strip()
            changed = True
        if s.startswith("join(") and s.endswith(")"):
            s = s[len("join("):-1].strip()
            changed = True
        if not changed:
            break

    intervals = []
    for m in _LOCATION_RE.finditer(s):
        a, b = int(m.group(1)), int(m.group(2))
        intervals.append((min(a, b), max(a, b)))

    intervals.sort()
    return intervals, is_complement


def _read_features_block(gb_text: str):
    """
    Read the FEATURES section of a GenBank file and merge continuation
    lines into single records.

    Returns a list of ``(feature_key, location_str)`` tuples in file order.
    """
    lines = gb_text.splitlines()
    features_start = None
    for i, line in enumerate(lines):
        if line.startswith("FEATURES"):
            features_start = i + 1
            break
    if features_start is None:
        return []

    records = []
    i = features_start
    while i < len(lines):
        line = lines[i]
        # Stop when we hit the next section (ORIGIN, //).
        if line.startswith("ORIGIN") or line.startswith("//"):
            break
        m = re.match(r"^\s{5}(\S+)\s+(.*)$", line)
        if m:
            key = m.group(1)
            loc = m.group(2).strip()
            # Collect continuation lines (21-space indented).
            j = i + 1
            while j < len(lines):
                cont = lines[j]
                if re.match(r"^\s{5}\S+\s", cont):
                    break  # new feature key
                if re.match(r"^\s{21}\S", cont) or cont.startswith(" " * 21):
                    loc += " " + cont[21:].strip()
                    j += 1
                else:
                    break
            records.append((key, loc))
            i = j
        else:
            i += 1

    return records


def parse_genome_length(gb_text: str) -> int:
    m = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", gb_text, re.MULTILINE)
    if m:
        return int(m.group(1))
    return 0


def parse_gb_file(gb_path: str):
    """
    Parse a single GenBank file and return:

        {
            'length':       genome length in bp,
            'cds_exons':    [(start, end), ...],  # all CDS exons
            'gene_ranges':  [(start, end), ...],  # all gene spans
            'cds_by_gene':  OrderedDict{ gene_range -> [exon_ranges] }
                                              # used for intron inference
        }
    """
    with open(gb_path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    length = parse_genome_length(text)
    records = _read_features_block(text)

    cds_exons = []
    gene_ranges = []
    cds_per_gene = OrderedDict()  # gene_range -> list[exon_range]

    # First pass: collect all gene and CDS features.
    for key, loc in records:
        if key == "gene":
            intervals, _ = parse_location(loc)
            for seg in intervals:
                if seg not in gene_ranges:
                    gene_ranges.append(seg)
                cds_per_gene.setdefault(seg, [])
        elif key == "CDS":
            intervals, _ = parse_location(loc)
            for seg in intervals:
                cds_exons.append(seg)

    # Associate CDS exons with the most-overlapping gene span.
    for gene_seg in cds_per_gene:
        best = []
        for ex in cds_exons:
            ov = max(0, min(gene_seg[1], ex[1]) - max(gene_seg[0], ex[0]) + 1)
            if ov > 0:
                best.append((ov, ex))
        best.sort(key=lambda x: (-x[0], x[1][0]))
        cds_per_gene[gene_seg] = [b[1] for b in best]

    return {
        "length": length,
        "cds_exons": cds_exons,
        "gene_ranges": gene_ranges,
        "cds_by_gene": cds_per_gene,
    }


def build_intron_intervals(gene_ranges, cds_by_gene):
    """
    Infer intron intervals from gene spans and CDS exons.

    Introns are strictly defined as the gaps between adjacent CDS exons
    within the same gene (5'/3' UTRs are not counted as introns).
    """
    introns = []
    for gene_seg, exons in cds_by_gene.items():
        if not exons:
            continue
        # Only keep exons that lie strictly inside the gene span.
        in_gene = [e for e in exons
                   if e[0] >= gene_seg[0] and e[1] <= gene_seg[1]]
        if len(in_gene) < 2:
            continue
        in_gene.sort()
        for i in range(len(in_gene) - 1):
            a_end = in_gene[i][1]
            b_start = in_gene[i + 1][0]
            if b_start > a_end + 1:
                # Half-open interval (a_end, b_start); an SSR is considered
                # intronic when start >= a_end and end <= b_start.
                introns.append((a_end, b_start))
    return introns


# ---------------------------------------------------------------------------
# 2. MISA parsing and region classification
# ---------------------------------------------------------------------------

def parse_misa(misa_path: str):
    """Read a MISA file. Returns ``(header, rows)`` where each row is a dict."""
    rows = []
    header = None
    with open(misa_path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            parts = line.split("\t")
            if header is None:
                header = [c.strip() for c in parts]
                continue
            row = {}
            for i, col in enumerate(parts):
                if i < len(header):
                    row[header[i]] = col.strip()
            rows.append(row)
    return header, rows


def classify_ssr(start: int, end: int, cds_exons, intron_intervals):
    """
    Decide which region an SSR belongs to.

    Priority:
      1. CDS if fully contained (or center falls inside a CDS exon).
      2. Intron if fully contained (or center falls inside an intron).
      3. Otherwise IGS.
    """
    # 1) Strict containment.
    for s, e in cds_exons:
        if s <= start and end <= e:
            return "CDS", f"CDS:{s}..{e}"
    for s, e in intron_intervals:
        if s <= start and end <= e:
            return "Intron", f"Intron:{s}..{e}"

    # 2) Fallback: classify by the midpoint of the SSR.
    mid = (start + end) / 2.0
    for s, e in cds_exons:
        if s <= mid <= e:
            return "CDS", f"CDS:{s}..{e}(center)"
    for s, e in intron_intervals:
        if s <= mid <= e:
            return "Intron", f"Intron:{s}..{e}(center)"

    # 3) No hit -> IGS.
    return "IGS", "IGS"


# ---------------------------------------------------------------------------
# 3. Main analysis pipeline
# ---------------------------------------------------------------------------

def analyze(misa_path: str, gb_dir: str, out_dir: str,
            log_func=print, cancel_flag=None):
    """
    Run the full analysis pipeline.

    log_func(str)   : callback used to emit log lines.
    cancel_flag     : threading.Event(); set it to abort the run.
    """
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, "log.txt")
    log_fp = open(log_path, "w", encoding="utf-8")

    def log(msg):
        ts = time.strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        log_func(line)
        log_fp.write(line + "\n")
        log_fp.flush()

    log("Analysis started.")
    log(f"MISA file: {misa_path}")
    log(f"GB folder: {gb_dir}")
    log(f"Output dir: {out_dir}")

    header, rows = parse_misa(misa_path)
    if not rows:
        raise RuntimeError("MISA file is empty or has an invalid format.")
    log(f"Loaded {len(rows)} SSR records from the MISA file.")

    # Collect unique species in file order.
    species_list = []
    seen = set()
    for r in rows:
        sp = r.get("ID", "").strip()
        if sp and sp not in seen:
            seen.add(sp)
            species_list.append(sp)
    log(f"Detected {len(species_list)} species: {', '.join(species_list)}")

    # Map species -> GB file path.
    gb_cache = {}
    for sp in species_list:
        candidates = [
            os.path.join(gb_dir, sp + ".gb"),
            os.path.join(gb_dir, sp + ".GB"),
            os.path.join(gb_dir, sp + ".genbank"),
        ]
        gb_path = None
        for c in candidates:
            if os.path.isfile(c):
                gb_path = c
                break
        if gb_path is None:
            # Recurse once into sub-directories.
            for root_dir, _, files in os.walk(gb_dir):
                for f in files:
                    name, ext = os.path.splitext(f)
                    if name == sp and ext.lower() in (".gb", ".genbank"):
                        gb_path = os.path.join(root_dir, f)
                        break
                if gb_path:
                    break
        if gb_path is None:
            log(f"  [Warning] No GB file found for species {sp}; skipped.")
            continue
        gb_cache[sp] = gb_path
        log(f"  - {sp} -> {gb_path}")

    if not gb_cache:
        raise RuntimeError("No GB file could be matched. "
                           "Please check the folder path and file names.")

    # Classify SSRs and accumulate statistics.
    detail_rows = []
    stats_species_region = defaultdict(lambda: defaultdict(int))  # species -> region -> count
    stats_region_type = defaultdict(lambda: defaultdict(int))      # region -> type -> count
    stats_species_total = defaultdict(int)
    stats_species_unmapped = defaultdict(int)

    total = len(rows)
    for idx, r in enumerate(rows, 1):
        if cancel_flag is not None and cancel_flag.is_set():
            log("Analysis cancelled by user.")
            log_fp.close()
            return False

        sp = r.get("ID", "").strip()
        try:
            start = int(r.get("start", "0"))
            end = int(r.get("end", "0"))
        except Exception:
            start, end = 0, 0

        r2 = dict(r)
        if sp not in gb_cache:
            r2["Region"] = "Unmapped"
            r2["RegionDetail"] = "missing_gb"
            stats_species_unmapped[sp] += 1
        else:
            gb_info = parse_gb_file(gb_cache[sp])
            introns = build_intron_intervals(gb_info["gene_ranges"],
                                             gb_info["cds_by_gene"])
            region, detail = classify_ssr(start, end,
                                          gb_info["cds_exons"], introns)
            r2["Region"] = region
            r2["RegionDetail"] = detail
            stats_species_region[sp][region] += 1
            stats_region_type[region][r.get("SSR type", "").strip()] += 1
            stats_species_total[sp] += 1

        detail_rows.append(r2)
        if idx % 200 == 0 or idx == total:
            log(f"  Progress: {idx}/{total} ({idx*100//total}%)")

    # ---- Output 1: detailed per-SSR file with region info ----
    out_detail = os.path.join(out_dir, "SSR_with_region.txt")
    cols = list(header) + ["Region", "RegionDetail"]
    with open(out_detail, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in detail_rows:
            w.writerow([r.get(c, "") for c in cols])
    log(f"Wrote detailed file: {out_detail}")

    # ---- Output 2: long-form summary (species, region, count, percent) ----
    out_long = os.path.join(out_dir, "SSR_region_summary.tsv")
    with open(out_long, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Species", "Region", "Count", "Percent"])
        species_order = [s for s in species_list if s in stats_species_total]
        for sp in species_order:
            total_sp = stats_species_total[sp]
            for region in ("CDS", "Intron", "IGS"):
                c = stats_species_region[sp].get(region, 0)
                pct = (c / total_sp * 100.0) if total_sp else 0.0
                w.writerow([sp, region, c, f"{pct:.2f}"])
    log(f"Wrote long-form summary: {out_long}")

    # ---- Output 3: wide-form summary (species x region) ----
    out_wide = os.path.join(out_dir, "SSR_region_summary_wide.tsv")
    with open(out_wide, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Species", "CDS", "Intron", "IGS", "Total"])
        for sp in species_order:
            cds = stats_species_region[sp].get("CDS", 0)
            intr = stats_species_region[sp].get("Intron", 0)
            igs = stats_species_region[sp].get("IGS", 0)
            w.writerow([sp, cds, intr, igs, cds + intr + igs])
    log(f"Wrote wide-form summary: {out_wide}")

    # ---- Output 4: region x SSR-type matrix (stacked-bar friendly) ----
    out_type = os.path.join(out_dir, "SSR_type_by_region.tsv")
    all_types = sorted({t for reg in stats_region_type
                        for t in stats_region_type[reg]})
    with open(out_type, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Region"] + all_types + ["Total"])
        for region in ("CDS", "Intron", "IGS"):
            row = [region]
            tot = 0
            for t in all_types:
                c = stats_region_type[region].get(t, 0)
                row.append(c)
                tot += c
            row.append(tot)
            w.writerow(row)
    log(f"Wrote region x SSR-type matrix: {out_type}")

    log("Analysis finished.")
    log_fp.close()
    return True


# ---------------------------------------------------------------------------
# 4. Tkinter GUI
# ---------------------------------------------------------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SSR Region Classifier (CDS / Intron / IGS)")
        self.geometry("820x620")
        self.minsize(760, 560)

        # State
        self.misa_path = tk.StringVar()
        self.gb_dir = tk.StringVar()
        self.out_dir = tk.StringVar()
        self.status = tk.StringVar(value="Ready")
        self.cancel_event = threading.Event()
        self.worker = None

        self._build_ui()

    # ---- UI construction ----
    def _build_ui(self):
        pad = {"padx": 8, "pady": 4}

        frm_top = ttk.LabelFrame(self, text="Input and Output")
        frm_top.pack(fill="x", padx=10, pady=8)

        # MISA file
        ttk.Label(frm_top, text="MISA file:").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(frm_top, textvariable=self.misa_path, width=70)\
            .grid(row=0, column=1, sticky="ew", **pad)
        ttk.Button(frm_top, text="Browse...", command=self.choose_misa)\
            .grid(row=0, column=2, **pad)

        # GB folder
        ttk.Label(frm_top, text="GB folder:").grid(row=1, column=0, sticky="w", **pad)
        ttk.Entry(frm_top, textvariable=self.gb_dir, width=70)\
            .grid(row=1, column=1, sticky="ew", **pad)
        ttk.Button(frm_top, text="Browse...", command=self.choose_gb)\
            .grid(row=1, column=2, **pad)

        # Output folder
        ttk.Label(frm_top, text="Output dir:").grid(row=2, column=0, sticky="w", **pad)
        ttk.Entry(frm_top, textvariable=self.out_dir, width=70)\
            .grid(row=2, column=1, sticky="ew", **pad)
        ttk.Button(frm_top, text="Browse...", command=self.choose_out)\
            .grid(row=2, column=2, **pad)

        frm_top.columnconfigure(1, weight=1)

        # Buttons
        frm_btn = ttk.Frame(self)
        frm_btn.pack(fill="x", padx=10)
        self.btn_run = ttk.Button(frm_btn, text="Run analysis",
                                  command=self.on_run)
        self.btn_run.pack(side="left", padx=4, pady=4)
        self.btn_cancel = ttk.Button(frm_btn, text="Cancel",
                                     command=self.on_cancel,
                                     state="disabled")
        self.btn_cancel.pack(side="left", padx=4, pady=4)
        ttk.Button(frm_btn, text="Clear log", command=self.clear_log)\
            .pack(side="left", padx=4, pady=4)
        ttk.Button(frm_btn, text="Open output folder",
                   command=self.open_out)\
            .pack(side="left", padx=4, pady=4)
        ttk.Button(frm_btn, text="Help", command=self.show_help)\
            .pack(side="right", padx=4, pady=4)

        # Progress bar
        self.progress = ttk.Progressbar(self, mode="determinate")
        self.progress.pack(fill="x", padx=10, pady=2)

        # Log
        frm_log = ttk.LabelFrame(self, text="Run log")
        frm_log.pack(fill="both", expand=True, padx=10, pady=8)
        self.txt_log = scrolledtext.ScrolledText(frm_log, height=20,
                                                 font=("Consolas", 10))
        self.txt_log.pack(fill="both", expand=True, padx=4, pady=4)

        # Status bar
        ttk.Label(self, textvariable=self.status, anchor="w",
                  relief="sunken").pack(fill="x", side="bottom")

    # ---- Path choosers ----
    def choose_misa(self):
        path = filedialog.askopenfilename(
            title="Select the MISA result file",
            filetypes=[("MISA files", "*.misa *.MISA *.txt *.tsv"),
                       ("All files", "*.*")])
        if path:
            self.misa_path.set(path)
            # Auto-suggest the GB folder and output directory.
            misa_dir = os.path.dirname(path)
            base = os.path.splitext(os.path.basename(path))[0]
            gb_guess = os.path.join(misa_dir, base)
            if os.path.isdir(gb_guess) and not self.gb_dir.get():
                self.gb_dir.set(gb_guess)
            if not self.out_dir.get():
                self.out_dir.set(os.path.join(misa_dir, base + "_region_result"))

    def choose_gb(self):
        d = filedialog.askdirectory(title="Select the GB folder")
        if d:
            self.gb_dir.set(d)

    def choose_out(self):
        d = filedialog.askdirectory(title="Select the output directory")
        if d:
            self.out_dir.set(d)

    def open_out(self):
        d = self.out_dir.get().strip()
        if d and os.path.isdir(d):
            try:
                if sys.platform.startswith("win"):
                    os.startfile(d)  # noqa
                elif sys.platform == "darwin":
                    os.system(f'open "{d}"')
                else:
                    os.system(f'xdg-open "{d}"')
            except Exception as e:
                messagebox.showerror("Open failed", str(e))
        else:
            messagebox.showinfo("Notice", "The output directory does not exist "
                                          "or has not been set.")

    # ---- Logging ----
    def log(self, msg: str):
        self.txt_log.insert("end", msg + "\n")
        self.txt_log.see("end")
        self.update_idletasks()

    def clear_log(self):
        self.txt_log.delete("1.0", "end")

    # ---- Run / cancel ----
    def on_run(self):
        misa = self.misa_path.get().strip()
        gb = self.gb_dir.get().strip()
        out = self.out_dir.get().strip()

        if not misa or not os.path.isfile(misa):
            messagebox.showerror("Error", "Please select a valid MISA file.")
            return
        if not gb or not os.path.isdir(gb):
            messagebox.showerror("Error", "Please select a valid GB folder.")
            return
        if not out:
            messagebox.showerror("Error", "Please select an output directory.")
            return
        os.makedirs(out, exist_ok=True)

        self.cancel_event.clear()
        self.btn_run.config(state="disabled")
        self.btn_cancel.config(state="normal")
        self.status.set("Analyzing...")

        self.worker = threading.Thread(
            target=self._worker, args=(misa, gb, out), daemon=True)
        self.worker.start()

    def _worker(self, misa, gb, out):
        try:
            ok = analyze(misa, gb, out,
                         log_func=lambda m: self.after(0, self.log, m),
                         cancel_flag=self.cancel_event)
            self.after(0, self._on_done, ok, None)
        except Exception as e:
            self.after(0, self._on_done, False, e)

    def _on_done(self, ok, err):
        self.btn_run.config(state="normal")
        self.btn_cancel.config(state="disabled")
        if err is not None:
            self.status.set("Error")
            messagebox.showerror("Runtime error", str(err))
        elif ok:
            self.status.set("Done")
            messagebox.showinfo("Done",
                                f"Analysis finished!\nResults saved to:\n"
                                f"{self.out_dir.get()}")
        else:
            self.status.set("Cancelled")

    def on_cancel(self):
        self.cancel_event.set()
        self.status.set("Cancelling...")

    # ---- Help ----
    def show_help(self):
        help_text = (
            "User Guide\n"
            "==========\n"
            "1. Select the MISA result file (*.misa / *.txt).\n"
            "2. Select the GB folder named after the MISA file "
            "(containing NC01.gb, NC02.gb, ...).\n"
            "3. Select an output directory.\n"
            "4. Click 'Run analysis'.\n\n"
            "Output files:\n"
            "  - SSR_with_region.txt          adds a Region column to every SSR\n"
            "  - SSR_region_summary.tsv       long form: species x region x count x %\n"
            "  - SSR_region_summary_wide.tsv  wide form: species x {CDS, Intron, IGS}\n"
            "  - SSR_type_by_region.tsv       region x SSR type (stacked-bar friendly)\n"
            "  - log.txt                      run log\n\n"
            "Classification rules:\n"
            "  CDS    : the SSR lies (or its midpoint lies) inside a CDS exon\n"
            "  Intron : the SSR lies inside an intron inferred from gene - CDS\n"
            "  IGS    : everything else is treated as intergenic spacer\n"
        )
        messagebox.showinfo("User Guide", help_text)


# ---------------------------------------------------------------------------
# 5. Entry point
# ---------------------------------------------------------------------------

def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()