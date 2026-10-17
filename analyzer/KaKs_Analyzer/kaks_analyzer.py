# -*- coding: utf-8 -*-
"""
KaKs Result Analysis Tool
=========================
Features:
    1. Auto-load either an organized .xlsx file or the raw .kaks text file.
    2. Group rows by gene, then drop genes whose missing-value ratio
       exceeds a user-selected threshold (<20% / 30% / 50% / 70%).
    3. Impute missing values using the per-gene median of valid values
       from other species.
       (Ka/Ks == 0 is treated as missing because, biologically, it is
       a degenerate output produced when Ka fails to compute.)
    4. Simultaneously export an analysis/filtering report
       (xlsx + txt).

Author:  MiniMax-M3
Requires: Python 3.8+, pandas, openpyxl, tkinter (standard library)
"""
import os
import re
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
APP_DIR = Path(__file__).resolve().parent     # Folder containing this script
EXAMPLE_FILE = APP_DIR / "example_kaks.xlsx"  # Bundled example file path


# ---------------------------------------------------------------------------
# Data I/O and core processing
# ---------------------------------------------------------------------------
class KaksProcessor:
    """Encapsulates KaKs file loading, filtering, imputation and reporting."""

    COL_KAKS = "Ka/Ks"
    COL_KA = "Ka"
    COL_KS = "Ks"

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------
    @staticmethod
    def load_file(file_path: str) -> pd.DataFrame:
        """Auto-detect file type (.xlsx organized, or .kaks raw)."""
        ext = os.path.splitext(file_path)[1].lower()
        if ext in (".xlsx", ".xls"):
            return KaksProcessor._load_xlsx(file_path)
        if ext in (".kaks", ".txt", ".tsv", ".axt.kaks", ""):
            return KaksProcessor._load_kaks(file_path)
        raise ValueError(f"Unsupported file type: {ext}")

    @staticmethod
    def _load_xlsx(path: str) -> pd.DataFrame:
        raw = pd.read_excel(path)
        if raw.shape[1] < 6:
            raise ValueError("xlsx must contain at least 6 columns "
                             "(Sp1, Sp2, Gene, Ka, Ks, Ka/Ks)")
        raw.columns = ["Species1", "Species2", "Gene", "Ka", "Ks", "Ka/Ks"]
        # Drop rows whose Species1 contains no alphabetic characters
        # (typically summary / blank rows in the original spreadsheet).
        raw = raw[raw["Species1"].astype(str).str.contains(
            r"[A-Za-z]", na=False)].copy()
        for col in ["Ka", "Ks", "Ka/Ks"]:
            raw[col] = pd.to_numeric(raw[col], errors="coerce")
        raw = KaksProcessor._deduplicate(raw)
        return raw.reset_index(drop=True)

    @staticmethod
    def _load_kaks(path: str) -> pd.DataFrame:
        raw = pd.read_csv(path, sep="\t")
        if "Sequence" not in raw.columns:
            raise ValueError(".kaks file is missing the 'Sequence' column")

        def parse_gene(s):
            # Forms like "accD_mafft" / "ndhB_copy2" / "ndhB_copy2_mafft"
            name = str(s).split("|")[-1].strip()
            # Strip trailing suffixes such as _mafft or _copy2 iteratively.
            for _ in range(5):
                new = re.sub(r"(_mafft|_copy\d+)$", "", name)
                if new == name:
                    break
                name = new
            return name

        def parse_species(s):
            head = str(s).split("|")[0]
            parts = head.split("&")
            if len(parts) >= 2:
                return parts[0].strip(), parts[1].strip()
            return parts[0].strip(), ""

        parsed = raw["Sequence"].apply(parse_species)
        raw["Species1"] = [p[0] for p in parsed]
        raw["Species2"] = [p[1] for p in parsed]
        raw["Gene"] = raw["Sequence"].apply(parse_gene)
        for col in ["Ka", "Ks", "Ka/Ks"]:
            raw[col] = pd.to_numeric(raw[col], errors="coerce")
        raw = KaksProcessor._deduplicate(
            raw[["Species1", "Species2", "Gene", "Ka", "Ks", "Ka/Ks"]])
        return raw.reset_index(drop=True)

    @staticmethod
    def _deduplicate(df: pd.DataFrame) -> pd.DataFrame:
        """Deduplicate rows by (Species1, Species2, Gene).

        Prefer rows whose Ka/Ks is a valid value (>0); if all are invalid,
        keep an arbitrary one.
        """
        if df.duplicated(subset=["Species1", "Species2", "Gene"]).any():
            df = df.copy()
            df["__score__"] = df[KaksProcessor.COL_KAKS].apply(
                lambda x: 2 if (pd.notna(x) and x > 0)
                else (1 if pd.notna(x) else 0))
            df = df.sort_values(
                ["Species1", "Species2", "Gene", "__score__"],
                ascending=[True, True, True, False])
            df = df.drop_duplicates(
                subset=["Species1", "Species2", "Gene"], keep="first")
            df = df.drop(columns=["__score__"])
        return df

    # ------------------------------------------------------------------
    # Missing-value detection (core)
    # ------------------------------------------------------------------
    @staticmethod
    def _is_valid(s: pd.Series) -> pd.Series:
        """A value is "valid" when it is not NaN and is greater than 0.

        Note: Ka/Ks == 0 is treated as invalid because, in KaKs software,
        it represents a degenerate output emitted when Ka cannot be
        computed; it carries no biological meaning.
        """
        return s.notna() & (s > 0)

    @staticmethod
    def _is_missing(s: pd.Series) -> pd.Series:
        """A value is "missing" when it is NaN or equal to 0."""
        return s.isna() | (s == 0)

    # ------------------------------------------------------------------
    # Core: filtering + imputation
    # ------------------------------------------------------------------
    @staticmethod
    def analyze(df: pd.DataFrame, threshold: float, fill_columns: list):
        """Compute the per-gene missing ratio, drop genes over the
        threshold, then impute remaining missing values with the median.

        Missing-value definitions:
            Ka/Ks column : NaN or 0 (Ka/Ks == 0 is the degenerate value
                           emitted when Ka fails to compute).
            Ka column    : NaN
            Ks column    : NaN

        Median computation:
            Uses only "valid values" (Ka/Ks > 0, Ka non-NaN, Ks non-NaN)
            so that degenerate zeros do not skew the median.

        Parameters
        ----------
        df : DataFrame
            Must contain the columns
            Species1, Species2, Gene, Ka, Ks, Ka/Ks.
        threshold : float
            Missing-ratio cutoff in [0, 1]. Genes whose missing ratio
            is greater than or equal to this value are removed.
        fill_columns : list
            Columns to impute, e.g. ["Ka/Ks"].

        Returns
        -------
        dict:
            total_genes   : total number of genes
            kept_genes    : number of genes kept
            removed_genes : list of removed gene names
            kept_df       : DataFrame containing only the kept genes
                            with imputation applied
            report_rows   : per-gene statistics
        """
        result = {
            "total_genes": 0,
            "kept_genes": 0,
            "removed_genes": [],
            "kept_df": None,
            "report_rows": [],
        }

        work = df.copy()
        work["__missing_kaks__"] = KaksProcessor._is_missing(work[KaksProcessor.COL_KAKS])
        work["__missing_ka__"] = work[KaksProcessor.COL_KA].isna()
        work["__missing_ks__"] = work[KaksProcessor.COL_KS].isna()

        grouped = work.groupby("Gene", sort=False)

        report_rows = []
        kept_dfs = []
        removed = []
        for gene, sub in grouped:
            total = len(sub)
            missing = int(sub["__missing_kaks__"].sum())
            missing_ka = int(sub["__missing_ka__"].sum())
            missing_ks = int(sub["__missing_ks__"].sum())
            rate = missing / total if total else 0.0

            # Always compute the median of all three columns (over valid
            # values) so that it can be shown in the report.
            medians = {}
            for col in [KaksProcessor.COL_KA,
                        KaksProcessor.COL_KS,
                        KaksProcessor.COL_KAKS]:
                if col == KaksProcessor.COL_KAKS:
                    valid = sub[col][sub[col] > 0]    # Exclude 0 and NaN
                else:
                    valid = sub[col].dropna()
                medians[col] = float(valid.median()) if len(valid) else np.nan

            status = "Kept" if rate < threshold else "Removed"
            if status == "Removed":
                removed.append(gene)
            else:
                # Median imputation
                sub_filled = sub.copy()
                for col in fill_columns:
                    val = medians[col]
                    if pd.isna(val):
                        continue
                    mask = sub_filled[col].isna() | (sub_filled[col] == 0)
                    sub_filled.loc[mask, col] = val
                kept_dfs.append(sub_filled)

            report_rows.append({
                "Gene": gene,
                "TotalRows": total,
                "MissingCount": missing,
                "MissingRate": round(rate, 4),
                "MissingKa": missing_ka,
                "MissingKs": missing_ks,
                "FilledCount": total - missing,
                "MedianKa": medians.get("Ka"),
                "MedianKs": medians.get("Ks"),
                "MedianKaKs": medians.get("Ka/Ks"),
                "Status": status,
            })

        result["total_genes"] = len(report_rows)
        result["removed_genes"] = removed
        result["kept_genes"] = result["total_genes"] - len(removed)
        result["report_rows"] = report_rows
        if kept_dfs:
            kept = pd.concat(kept_dfs, ignore_index=True)
            for c in ["__missing_kaks__", "__missing_ka__", "__missing_ks__"]:
                if c in kept.columns:
                    kept = kept.drop(columns=[c])
        else:
            kept = pd.DataFrame(columns=df.columns)
        result["kept_df"] = kept
        return result

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------
    @staticmethod
    def export(result: dict, output_dir: str, threshold: float,
               source_file: str):
        """Export the processed data plus the analysis report.

        Files written:
            <name>_filtered.xlsx            Filtered + imputed data
            <name>_report.xlsx              Per-gene statistics report
            <name>_summary.txt              Concise text report
        """
        os.makedirs(output_dir, exist_ok=True)
        base = os.path.splitext(os.path.basename(source_file))[0]
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        prefix = f"{base}_filtered_{int(threshold*100)}pct_{ts}"

        # 1. Filtered data
        xlsx_data = os.path.join(output_dir, f"{prefix}.xlsx")
        if result["kept_df"] is not None and not result["kept_df"].empty:
            result["kept_df"].to_excel(xlsx_data, index=False)
        else:
            pd.DataFrame().to_excel(xlsx_data, index=False)

        # 2. Report (one row per gene)
        report_df = pd.DataFrame(result["report_rows"])
        xlsx_report = os.path.join(output_dir, f"{prefix}_report.xlsx")
        report_df.to_excel(xlsx_report, index=False)

        # 3. Concise text report
        txt_report = os.path.join(output_dir, f"{prefix}_summary.txt")
        with open(txt_report, "w", encoding="utf-8") as fh:
            fh.write("=" * 80 + "\n")
            fh.write("KaKs Analysis & Filtering Report\n")
            fh.write("=" * 80 + "\n")
            fh.write(f"Source file       : {source_file}\n")
            fh.write(f"Generated at      : "
                     f"{datetime.now():%Y-%m-%d %H:%M:%S}\n")
            fh.write(f"Missing ratio     : < {threshold * 100:.0f}%\n")
            fh.write("Missing rule      : Ka/Ks NaN or == 0 "
                     "(Ka/Ks==0 is treated as a degenerate value "
                     "produced when Ka fails to compute)\n")
            fh.write("Imputation        : Median of valid values from "
                     "the same gene's other species\n")
            fh.write("-" * 80 + "\n")
            fh.write(f"Total genes       : {result['total_genes']}\n")
            fh.write(f"Kept genes        : {result['kept_genes']}\n")
            fh.write(f"Removed genes     : {len(result['removed_genes'])}\n")
            fh.write("-" * 80 + "\n")
            fh.write(f"Removed genes ({len(result['removed_genes'])}):\n")
            if result["removed_genes"]:
                for g in result["removed_genes"]:
                    fh.write(f"   - {g}\n")
            else:
                fh.write("   (none)\n")
            fh.write("-" * 80 + "\n")
            fh.write("Per-gene statistics:\n")
            fh.write(f"{'Gene':<10}{'Total':>6}{'MissKa/Ks':>10}"
                     f"{'MissKa':>7}{'MissKs':>7}{'Rate':>8}"
                     f"{'+Filled':>8}{'MedKa':>11}{'MedKs':>11}"
                     f"{'MedKa/Ks':>11}{'Status':>9}\n")
            for row in result["report_rows"]:

                def fmt(v, w):
                    if v is None or (isinstance(v, float) and v != v):
                        return "NA".rjust(w)
                    return f"{v:>{w}.4g}"

                fh.write(
                    f"{row['Gene']:<10}"
                    f"{row['TotalRows']:>6}"
                    f"{row['MissingCount']:>10}"
                    f"{row['MissingKa']:>7}"
                    f"{row['MissingKs']:>7}"
                    f"{row['MissingRate']*100:>7.2f}%"
                    f"{row['FilledCount']:>8}"
                    f"{fmt(row['MedianKa'], 11)}"
                    f"{fmt(row['MedianKs'], 11)}"
                    f"{fmt(row['MedianKaKs'], 11)}"
                    f"{row['Status']:>9}\n"
                )
            fh.write("=" * 80 + "\n")
            fh.write("\nField descriptions:\n")
            fh.write("  Total       : Number of rows for the gene "
                     "(i.e. species pairs).\n")
            fh.write("  MissKa/Ks   : Missing count for the Ka/Ks column "
                     "(NaN or == 0).\n")
            fh.write("  MissKa      : Missing count for the Ka column "
                     "(NaN).\n")
            fh.write("  MissKs      : Missing count for the Ks column "
                     "(NaN).\n")
            fh.write("  Rate        : Ka/Ks missing ratio = "
                     "MissKa/Ks / Total.\n")
            fh.write("  +Filled     : Number of valid rows after "
                     "imputation.\n")
            fh.write("  MedKa/MedKs/MedKa/Ks : Median used for imputation "
                     "(computed over valid values).\n")
            fh.write("  Status      : Kept or Removed.\n")

        return {
            "data": xlsx_data,
            "report": xlsx_report,
            "summary": txt_report,
        }


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------
class App(tk.Tk):
    THRESHOLDS = [
        ("20% (drop genes with missing ratio >= 20%)", 0.20),
        ("30% (drop genes with missing ratio >= 30%)", 0.30),
        ("50% (drop genes with missing ratio >= 50%)", 0.50),
        ("70% (drop genes with missing ratio >= 70%)", 0.70),
    ]

    def __init__(self):
        super().__init__()
        self.title("KaKs Result Analysis & Filtering Tool v2.0")
        self.geometry("950x720")
        self.minsize(860, 660)

        self._df = None
        self._file_path = tk.StringVar()
        self._out_dir = tk.StringVar(value=str(APP_DIR))
        self._threshold_var = tk.IntVar(value=2)   # Default 50%
        self._fill_ka = tk.BooleanVar(value=False)
        self._fill_ks = tk.BooleanVar(value=False)
        self._fill_kaks = tk.BooleanVar(value=True)

        self._build_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        try:
            ttk.Style().theme_use("vista")
        except Exception:
            pass

        pad = {"padx": 8, "pady": 4}

        # ===== Top: File selection =====
        frm_file = ttk.LabelFrame(self, text="(1) Select Input File")
        frm_file.pack(fill="x", **pad)
        ttk.Entry(frm_file, textvariable=self._file_path, width=80).pack(
            side="left", fill="x", expand=True, padx=6, pady=6)
        ttk.Button(frm_file, text="Browse...", command=self._choose_file).pack(
            side="left", padx=4, pady=6)
        ttk.Button(frm_file, text="Load",
                   command=self._load_file).pack(side="left", padx=4, pady=6)
        ttk.Button(frm_file, text="Use Example",
                   command=self._use_example).pack(side="left", padx=4, pady=6)

        # ===== Parameters =====
        frm_param = ttk.LabelFrame(self, text="(2) Filter Parameters")
        frm_param.pack(fill="x", **pad)

        ttk.Label(frm_param,
                  text="Missing-ratio threshold "
                       "(genes with Ka/Ks missing ratio >= threshold "
                       "are dropped):"
                  ).grid(row=0, column=0, sticky="w", padx=6, pady=4)
        for i, (label, val) in enumerate(self.THRESHOLDS):
            ttk.Radiobutton(frm_param, text=label, value=i,
                            variable=self._threshold_var
                            ).grid(row=i + 1, column=0, sticky="w", padx=20)

        ttk.Label(frm_param, text="Columns to impute with median:"
                  ).grid(row=0, column=1, sticky="w", padx=20, pady=4)
        ttk.Checkbutton(frm_param, text="Impute Ka (NaN -> median)",
                        variable=self._fill_ka
                        ).grid(row=1, column=1, sticky="w", padx=24)
        ttk.Checkbutton(frm_param, text="Impute Ks (NaN -> median)",
                        variable=self._fill_ks
                        ).grid(row=2, column=1, sticky="w", padx=24)
        ttk.Checkbutton(frm_param, text="Impute Ka/Ks (NaN/0 -> median)",
                        variable=self._fill_kaks
                        ).grid(row=3, column=1, sticky="w", padx=24)
        ttk.Label(frm_param,
                  text="* Ka/Ks == 0 is treated as missing "
                       "(KaKs software emits 0 when Ka computation fails).",
                  foreground="gray"
                  ).grid(row=4, column=1, sticky="w", padx=24)

        # ===== Output directory =====
        frm_out = ttk.LabelFrame(self, text="(3) Output Directory")
        frm_out.pack(fill="x", **pad)
        ttk.Entry(frm_out, textvariable=self._out_dir, width=80).pack(
            side="left", fill="x", expand=True, padx=6, pady=6)
        ttk.Button(frm_out, text="Browse...",
                   command=self._choose_outdir).pack(side="left", padx=4, pady=6)

        # ===== Action buttons =====
        frm_act = ttk.Frame(self)
        frm_act.pack(fill="x", **pad)
        self.btn_run = ttk.Button(frm_act, text="Start Analysis",
                                  command=self._run)
        self.btn_run.pack(side="left", padx=6, pady=4)
        ttk.Button(frm_act, text="Clear Log",
                   command=lambda: self._log.delete("1.0", "end")
                   ).pack(side="left", padx=4)
        ttk.Button(frm_act, text="Open Output Folder",
                   command=self._open_outdir).pack(side="left", padx=4)
        ttk.Button(frm_act, text="About",
                   command=self._show_about).pack(side="right", padx=4)

        # ===== Progress bar =====
        self.progress = ttk.Progressbar(self, mode="determinate", maximum=100)
        self.progress.pack(fill="x", padx=10, pady=2)

        # ===== Notebook: Log + Report preview + Data preview =====
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=6)

        # Log
        self._log = scrolledtext.ScrolledText(nb, height=10, wrap="word")
        nb.add(self._log, text="Run Log")

        # Report preview
        self.tree, frm_tree = self._make_treeview(nb)
        nb.add(frm_tree, text="Report Preview")

        # Data preview
        self.tree_data, frm_data = self._make_treeview(nb)
        nb.add(frm_data, text="Data Preview")

    def _make_treeview(self, parent):
        frm = ttk.Frame(parent)
        tv = ttk.Treeview(frm, show="headings", height=14)
        vsb = ttk.Scrollbar(frm, orient="vertical", command=tv.yview)
        hsb = ttk.Scrollbar(frm, orient="horizontal", command=tv.xview)
        tv.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tv.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frm.rowconfigure(0, weight=1)
        frm.columnconfigure(0, weight=1)
        return tv, frm

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------
    def _choose_file(self):
        p = filedialog.askopenfilename(
            title="Select a KaKs result file",
            filetypes=[("KaKs results", "*.xlsx *.xls *.kaks *.txt *.tsv"),
                       ("All files", "*.*")])
        if p:
            self._file_path.set(p)

    def _choose_outdir(self):
        p = filedialog.askdirectory(title="Select output directory",
                                    initialdir=self._out_dir.get())
        if p:
            self._out_dir.set(p)

    def _open_outdir(self):
        d = self._out_dir.get()
        if os.path.isdir(d):
            try:
                os.startfile(d)   # Windows
            except AttributeError:
                import subprocess
                subprocess.Popen(["xdg-open", d])

    def _use_example(self):
        """Load the bundled example file."""
        if EXAMPLE_FILE.exists():
            self._file_path.set(str(EXAMPLE_FILE))
            self._log_write(f"Example file selected: {EXAMPLE_FILE}")
        else:
            messagebox.showwarning(
                "Example file not found",
                f"Could not locate example file: {EXAMPLE_FILE}\n"
                "Please use 'Browse...' to select a .xlsx or .kaks file.")

    def _load_file(self):
        path = self._file_path.get().strip()
        if not path:
            messagebox.showwarning("Notice", "Please choose an input file first")
            return
        if not os.path.exists(path):
            messagebox.showerror("Error", f"File does not exist:\n{path}")
            return
        try:
            self._log_write(f"Loading file: {path}")
            self._df = KaksProcessor.load_file(path)
            self._log_write(
                f"Loaded successfully: {len(self._df)} rows, "
                f"{self._df['Gene'].nunique()} genes, "
                f"{self._df['Species2'].nunique()} compared species.")
            self._populate_data_preview()
            self._log_write("Data preview updated.")
        except Exception as exc:
            self._log_write(f"[Error] Failed to load: {exc}")
            messagebox.showerror("Load failed", str(exc))

    def _populate_data_preview(self):
        self.tree_data.delete(*self.tree_data.get_children())
        if self._df is None:
            return
        cols = list(self._df.columns)
        self.tree_data["columns"] = cols
        for c in cols:
            self.tree_data.heading(c, text=c)
            self.tree_data.column(c, width=110, anchor="w")
        for _, row in self._df.head(500).iterrows():
            self.tree_data.insert("", "end",
                                  values=[str(v) for v in row.tolist()])

    def _get_fill_columns(self):
        cols = []
        if self._fill_ka.get():
            cols.append("Ka")
        if self._fill_ks.get():
            cols.append("Ks")
        if self._fill_kaks.get():
            cols.append("Ka/Ks")
        if not cols:
            messagebox.showwarning("Notice",
                                   "Please select at least one column to impute")
            return None
        return cols

    def _run(self):
        if self._df is None:
            messagebox.showwarning("Notice", "Please load a file first")
            return
        cols = self._get_fill_columns()
        if cols is None:
            return
        threshold = self.THRESHOLDS[self._threshold_var.get()][1]

        self.btn_run.config(state="disabled")
        self.progress["value"] = 0
        self._log_write(
            f"Starting | threshold={threshold*100:.0f}% "
            f"| impute columns={cols}")

        def task():
            try:
                self._post(20, lambda: self.progress.configure(value=30))
                result = KaksProcessor.analyze(
                    self._df, threshold, cols)
                self._post(40, lambda: self.progress.configure(value=60))
                files = KaksProcessor.export(
                    result, self._out_dir.get(), threshold,
                    self._file_path.get())
                self._post(80, lambda: self.progress.configure(value=90))
                self._last_result = result
                self._last_files = files
                self._post(0, self._populate_report_preview)
                self._post(0, lambda: self._log_write(
                    f"Done!\nTotal genes: {result['total_genes']}\n"
                    f"Kept: {result['kept_genes']}\n"
                    f"Removed: {len(result['removed_genes'])}\n"
                    f"Generated files:\n  - {files['data']}\n"
                    f"  - {files['report']}\n"
                    f"  - {files['summary']}\n"))
                self._post(0, lambda: self.progress.configure(value=100))
            except Exception as exc:
                self._post(0, lambda: self._log_write(f"[Error] {exc}"))
                self._post(0, lambda: messagebox.showerror("Processing failed", str(exc)))
            finally:
                self._post(0, lambda: self.btn_run.config(state="normal"))

        threading.Thread(target=task, daemon=True).start()

    def _post(self, delay_ms, fn):
        """Thread-safe GUI update; dispatches the callback to Tk's main loop."""
        if delay_ms > 0:
            self.after(delay_ms, fn)
        else:
            self.after_idle(fn)

    def _populate_report_preview(self):
        if not hasattr(self, "_last_result"):
            return
        report_rows = self._last_result["report_rows"]
        cols = ["Gene", "TotalRows", "MissingCount", "MissingKa",
                "MissingKs", "MissingRate", "FilledCount",
                "MedianKa", "MedianKs", "MedianKaKs", "Status"]
        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = cols
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=90, anchor="w")
        for r in report_rows:
            vals = []
            for c in cols:
                v = r[c]
                if c == "MissingRate":
                    vals.append(f"{v*100:.2f}%")
                elif c in ("MedianKa", "MedianKs", "MedianKaKs"):
                    vals.append("NA" if (v is None or
                                          (isinstance(v, float) and v != v))
                                else f"{v:.5g}")
                else:
                    vals.append(str(v))
            self.tree.insert("", "end", values=vals)

    def _show_about(self):
        messagebox.showinfo(
            "About",
            "KaKs Result Analysis & Filtering Tool v2.0\n\n"
            "Features:\n"
            " - Loads .xlsx organized sheets or raw .kaks text files\n"
            " - Drops genes whose missing ratio meets/exceeds the selected "
            "threshold (20% / 30% / 50% / 70%)\n"
            " - Imputes missing data with the per-gene median of valid "
            "values from other species\n"
            " - Exports filtered data + per-gene report + text summary\n\n"
            "Missing-value rule:\n"
            " - Ka/Ks column: NaN or 0 (Ka/Ks==0 is the degenerate value "
            "produced when Ka computation fails)\n"
            " - Ka, Ks columns: NaN\n\n"
            "Median imputation:\n"
            " - Computed over valid values only; missing/degenerate values "
            "do not influence the median\n\n"
            "Dependencies: pandas / openpyxl "
            "(tkinter is part of the standard library)")

    def _log_write(self, msg):
        self._log.insert("end", f"[{datetime.now():%H:%M:%S}] {msg}\n")
        self._log.see("end")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    app = App()
    app.mainloop()