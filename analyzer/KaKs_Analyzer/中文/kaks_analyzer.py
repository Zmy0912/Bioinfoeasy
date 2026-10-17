# -*- coding: utf-8 -*-
"""
KaKs 结果分析工具
=================
功能:
    1. 自动加载整理后的 .xlsx 文件 或 原始 .kaks 文本文件
    2. 按基因分组, 按缺失比例阈值(<20%/30%/50%/70%)筛选并删除缺失过多的基因
    3. 用每个基因其他物种有效值的中位数填充缺失值
       (Ka/Ks=0 被视为缺失 — 该值在生物学上代表 Ka 无法计算的退化值)
    4. 同时输出分析/筛选报告 (xlsx + txt)

作者:  MiniMax-M3
依赖:  Python 3.8+, pandas, openpyxl, tkinter(标准库)
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
# 常量
# ---------------------------------------------------------------------------
APP_DIR = Path(__file__).resolve().parent     # 程序所在目录
EXAMPLE_FILE = APP_DIR / "example_kaks.xlsx"  # 示例文件路径 (与程序同目录)


# ---------------------------------------------------------------------------
# 数据读写与核心处理
# ---------------------------------------------------------------------------
class KaksProcessor:
    """封装 KaKs 文件加载、筛选、填充、报告生成的业务逻辑."""

    COL_KAKS = "Ka/Ks"
    COL_KA = "Ka"
    COL_KS = "Ks"

    # ------------------------------------------------------------------
    # 加载
    # ------------------------------------------------------------------
    @staticmethod
    def load_file(file_path: str) -> pd.DataFrame:
        """根据扩展名自动识别 xlsx(整理后) 或 .kaks(原始) 格式."""
        ext = os.path.splitext(file_path)[1].lower()
        if ext in (".xlsx", ".xls"):
            return KaksProcessor._load_xlsx(file_path)
        if ext in (".kaks", ".txt", ".tsv", ".axt.kaks", ""):
            return KaksProcessor._load_kaks(file_path)
        raise ValueError(f"不支持的文件类型: {ext}")

    @staticmethod
    def _load_xlsx(path: str) -> pd.DataFrame:
        raw = pd.read_excel(path)
        if raw.shape[1] < 6:
            raise ValueError("xlsx 文件至少需要 6 列 (Sp1,Sp2,Gene,Ka,Ks,Ka/Ks)")
        raw.columns = ["Species1", "Species2", "Gene", "Ka", "Ks", "Ka/Ks"]
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
            raise ValueError(".kaks 文件缺少 'Sequence' 列")

        def parse_gene(s):
            # 形如 "accD_mafft" / "ndhB_copy2" / "ndhB_copy2_mafft"
            name = str(s).split("|")[-1].strip()
            # 反复去掉尾部后缀 (如 _mafft 或 _copy2)
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
        """对 (Species1, Species2, Gene) 组合去重.
        优先保留 Ka/Ks 有效值(>0)的行; 若都缺失则任选其一.
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
    # 缺失识别 (核心)
    # ------------------------------------------------------------------
    @staticmethod
    def _is_valid(s: pd.Series) -> pd.Series:
        """判断一列是否为'有效值': 非 NaN 且 > 0.
        注意: Ka/Ks 列的 0 值视为缺失, 因为该值在 KaKs 软件中代表
              Ka 计算失败的退化值, 没有生物学意义.
        """
        return s.notna() & (s > 0)

    @staticmethod
    def _is_missing(s: pd.Series) -> pd.Series:
        """判断一列是否为'缺失': NaN 或 == 0."""
        return s.isna() | (s == 0)

    # ------------------------------------------------------------------
    # 核心: 筛选 + 填充
    # ------------------------------------------------------------------
    @staticmethod
    def analyze(df: pd.DataFrame, threshold: float, fill_columns: list):
        """按基因统计缺失比例, 删除超过阈值的基因, 然后用中位数填充.

        缺失定义:
            Ka/Ks 列的 NaN 或 0 值 (Ka/Ks=0 表示 Ka 计算失败的退化值)
            Ka   列的 NaN 值
            Ks   列的 NaN 值

        中位数计算:
            只使用"有效值"(Ka/Ks>0, Ka 非NaN, Ks 非NaN) 计算中位数
            排除缺失/退化值的影响, 使填充值更符合生物学意义.

        Parameters
        ----------
        df : DataFrame
            包含 Species1, Species2, Gene, Ka, Ks, Ka/Ks 列
        threshold : float
            缺失比例阈值 (0~1). 缺失率 >= threshold 的基因被删除
        fill_columns : list
            需要进行填充的列名列表, 例如 ["Ka/Ks"]

        Returns
        -------
        dict:
            total_genes   : 总基因数
            kept_genes    : 保留基因数
            removed_genes : 被删除的基因名列表
            kept_df       : 保留基因 + 已填充的 DataFrame
            report_rows   : 每个基因的统计行
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

            # 始终计算三个列的中位数 (基于有效值), 便于在报告中展示
            medians = {}
            for col in [KaksProcessor.COL_KA, KaksProcessor.COL_KS, KaksProcessor.COL_KAKS]:
                if col == KaksProcessor.COL_KAKS:
                    valid = sub[col][sub[col] > 0]    # 排除 0 和 NaN
                else:
                    valid = sub[col].dropna()
                medians[col] = float(valid.median()) if len(valid) else np.nan

            status = "保留" if rate < threshold else "删除"
            if status == "删除":
                removed.append(gene)
            else:
                # 中位数填充
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
    # 导出
    # ------------------------------------------------------------------
    @staticmethod
    def export(result: dict, output_dir: str, threshold: float,
               source_file: str):
        """导出处理后的数据 + 报告.

        生成文件:
            <name>_filtered.xlsx         筛选+填充后的数据
            <name>_report.xlsx           各基因统计报告
            <name>_summary.txt           简要文本报告
        """
        os.makedirs(output_dir, exist_ok=True)
        base = os.path.splitext(os.path.basename(source_file))[0]
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        prefix = f"{base}_filtered_{int(threshold*100)}pct_{ts}"

        # 1. 保留数据
        xlsx_data = os.path.join(output_dir, f"{prefix}.xlsx")
        if result["kept_df"] is not None and not result["kept_df"].empty:
            result["kept_df"].to_excel(xlsx_data, index=False)
        else:
            pd.DataFrame().to_excel(xlsx_data, index=False)

        # 2. 报告 (每个基因)
        report_df = pd.DataFrame(result["report_rows"])
        xlsx_report = os.path.join(output_dir, f"{prefix}_report.xlsx")
        report_df.to_excel(xlsx_report, index=False)

        # 3. 简要文本报告
        txt_report = os.path.join(output_dir, f"{prefix}_summary.txt")
        with open(txt_report, "w", encoding="utf-8") as fh:
            fh.write("=" * 80 + "\n")
            fh.write("KaKs 分析与筛选报告\n")
            fh.write("=" * 80 + "\n")
            fh.write(f"源文件         : {source_file}\n")
            fh.write(f"生成时间       : {datetime.now():%Y-%m-%d %H:%M:%S}\n")
            fh.write(f"缺失比例阈值   : < {threshold * 100:.0f}%\n")
            fh.write("缺失识别       : Ka/Ks 列的 NaN 或 0 "
                     "(Ka/Ks=0 视为 Ka 计算失败的退化值)\n")
            fh.write("填充策略       : 用本基因其他物种有效值的中位数填充\n")
            fh.write("-" * 80 + "\n")
            fh.write(f"总基因数       : {result['total_genes']}\n")
            fh.write(f"保留基因数     : {result['kept_genes']}\n")
            fh.write(f"删除基因数     : {len(result['removed_genes'])}\n")
            fh.write("-" * 80 + "\n")
            fh.write(f"删除的基因 ({len(result['removed_genes'])}):\n")
            if result["removed_genes"]:
                for g in result["removed_genes"]:
                    fh.write(f"   - {g}\n")
            else:
                fh.write("   (无)\n")
            fh.write("-" * 80 + "\n")
            fh.write("各基因统计明细:\n")
            fh.write(f"{'Gene':<10}{'Total':>6}{'MissKa/Ks':>10}{'MissKa':>7}"
                     f"{'MissKs':>7}{'Rate':>8}"
                     f"{'+Filled':>8}{'MedKa':>11}{'MedKs':>11}"
                     f"{'MedKa/Ks':>11}{'Status':>8}\n")
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
                    f"{row['Status']:>8}\n"
                )
            fh.write("=" * 80 + "\n")
            fh.write("\n字段说明:\n")
            fh.write("  Total      : 该基因在数据中的总行数 (物种对数)\n")
            fh.write("  MissKa/Ks  : Ka/Ks 列缺失数 (NaN 或 =0)\n")
            fh.write("  MissKa     : Ka 列缺失数 (NaN)\n")
            fh.write("  MissKs     : Ks 列缺失数 (NaN)\n")
            fh.write("  Rate       : Ka/Ks 缺失比例 = MissKa/Ks / Total\n")
            fh.write("  +Filled    : 填充后有效数据行数\n")
            fh.write("  MedKa/MedKs/MedKa/Ks : 填充所用的中位数 (基于有效值)\n")
            fh.write("  Status     : 保留 或 删除\n")

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
        ("20% (缺失率 >= 20% 删除)", 0.20),
        ("30% (缺失率 >= 30% 删除)", 0.30),
        ("50% (缺失率 >= 50% 删除)", 0.50),
        ("70% (缺失率 >= 70% 删除)", 0.70),
    ]

    def __init__(self):
        super().__init__()
        self.title("KaKs 结果分析与筛选工具 v2.0")
        self.geometry("950x720")
        self.minsize(860, 660)

        self._df = None
        self._file_path = tk.StringVar()
        self._out_dir = tk.StringVar(value=str(APP_DIR))
        self._threshold_var = tk.IntVar(value=2)   # 默认 50%
        self._fill_ka = tk.BooleanVar(value=False)
        self._fill_ks = tk.BooleanVar(value=False)
        self._fill_kaks = tk.BooleanVar(value=True)

        self._build_ui()

    # ------------------------------------------------------------------
    # 界面构建
    # ------------------------------------------------------------------
    def _build_ui(self):
        try:
            ttk.Style().theme_use("vista")
        except Exception:
            pass

        pad = {"padx": 8, "pady": 4}

        # ===== 顶部: 文件选择 =====
        frm_file = ttk.LabelFrame(self, text="① 选择输入文件")
        frm_file.pack(fill="x", **pad)
        ttk.Entry(frm_file, textvariable=self._file_path, width=80).pack(
            side="left", fill="x", expand=True, padx=6, pady=6)
        ttk.Button(frm_file, text="浏览…", command=self._choose_file).pack(
            side="left", padx=4, pady=6)
        ttk.Button(frm_file, text="加载", command=self._load_file).pack(
            side="left", padx=4, pady=6)
        ttk.Button(frm_file, text="使用示例", command=self._use_example).pack(
            side="left", padx=4, pady=6)

        # ===== 参数 =====
        frm_param = ttk.LabelFrame(self, text="② 筛选参数")
        frm_param.pack(fill="x", **pad)

        ttk.Label(frm_param, text="缺失比例阈值 (删除 Ka/Ks 缺失比例 ≥ 该阈值的基因):"
                 ).grid(row=0, column=0, sticky="w", padx=6, pady=4)
        for i, (label, val) in enumerate(self.THRESHOLDS):
            ttk.Radiobutton(frm_param, text=label, value=i,
                            variable=self._threshold_var
                            ).grid(row=i + 1, column=0, sticky="w", padx=20)

        ttk.Label(frm_param, text="中位数填充列:").grid(
            row=0, column=1, sticky="w", padx=20, pady=4)
        ttk.Checkbutton(frm_param, text="填充 Ka (NaN → 中位数)",
                        variable=self._fill_ka
                        ).grid(row=1, column=1, sticky="w", padx=24)
        ttk.Checkbutton(frm_param, text="填充 Ks (NaN → 中位数)",
                        variable=self._fill_ks
                        ).grid(row=2, column=1, sticky="w", padx=24)
        ttk.Checkbutton(frm_param, text="填充 Ka/Ks (NaN/0 → 中位数)",
                        variable=self._fill_kaks
                        ).grid(row=3, column=1, sticky="w", padx=24)
        ttk.Label(frm_param,
                  text="* Ka/Ks=0 视为缺失(因 Ka 计算失败时软件输出 0)",
                  foreground="gray"
                  ).grid(row=4, column=1, sticky="w", padx=24)

        # ===== 输出目录 =====
        frm_out = ttk.LabelFrame(self, text="③ 输出目录")
        frm_out.pack(fill="x", **pad)
        ttk.Entry(frm_out, textvariable=self._out_dir, width=80).pack(
            side="left", fill="x", expand=True, padx=6, pady=6)
        ttk.Button(frm_out, text="浏览…", command=self._choose_outdir).pack(
            side="left", padx=4, pady=6)

        # ===== 操作按钮 =====
        frm_act = ttk.Frame(self)
        frm_act.pack(fill="x", **pad)
        self.btn_run = ttk.Button(frm_act, text="开始分析", command=self._run)
        self.btn_run.pack(side="left", padx=6, pady=4)
        ttk.Button(frm_act, text="清空日志",
                   command=lambda: self._log.delete("1.0", "end")
                   ).pack(side="left", padx=4)
        ttk.Button(frm_act, text="打开输出目录",
                   command=self._open_outdir).pack(side="left", padx=4)
        ttk.Button(frm_act, text="关于",
                   command=self._show_about).pack(side="right", padx=4)

        # ===== 进度条 =====
        self.progress = ttk.Progressbar(self, mode="determinate", maximum=100)
        self.progress.pack(fill="x", padx=10, pady=2)

        # ===== Notebook: 日志 + 报告预览 + 数据预览 =====
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=6)

        # 日志
        self._log = scrolledtext.ScrolledText(nb, height=10, wrap="word")
        nb.add(self._log, text="运行日志")

        # 报告预览
        self.tree, frm_tree = self._make_treeview(nb)
        nb.add(frm_tree, text="报告预览")

        # 数据预览
        self.tree_data, frm_data = self._make_treeview(nb)
        nb.add(frm_data, text="数据预览")

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
    # 事件
    # ------------------------------------------------------------------
    def _choose_file(self):
        p = filedialog.askopenfilename(
            title="选择 KaKs 结果文件",
            filetypes=[("KaKs 结果", "*.xlsx *.xls *.kaks *.txt *.tsv"),
                       ("所有文件", "*.*")])
        if p:
            self._file_path.set(p)

    def _choose_outdir(self):
        p = filedialog.askdirectory(title="选择输出目录",
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
        """使用程序自带示例文件."""
        if EXAMPLE_FILE.exists():
            self._file_path.set(str(EXAMPLE_FILE))
            self._log_write(f"已选择示例文件: {EXAMPLE_FILE}")
        else:
            messagebox.showwarning(
                "示例文件不存在",
                f"未找到示例文件: {EXAMPLE_FILE}\n"
                f"请使用「浏览」手动选择 .xlsx 或 .kaks 文件.")

    def _load_file(self):
        path = self._file_path.get().strip()
        if not path:
            messagebox.showwarning("提示", "请先选择输入文件")
            return
        if not os.path.exists(path):
            messagebox.showerror("错误", f"文件不存在:\n{path}")
            return
        try:
            self._log_write(f"正在加载文件: {path}")
            self._df = KaksProcessor.load_file(path)
            self._log_write(
                f"加载成功, 共 {len(self._df)} 行, "
                f"{self._df['Gene'].nunique()} 个基因, "
                f"{self._df['Species2'].nunique()} 个对比物种.")
            self._populate_data_preview()
            self._log_write("数据预览已更新.")
        except Exception as exc:
            self._log_write(f"[错误] 加载失败: {exc}")
            messagebox.showerror("加载失败", str(exc))

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
            messagebox.showwarning("提示", "请至少选择一列进行填充")
            return None
        return cols

    def _run(self):
        if self._df is None:
            messagebox.showwarning("提示", "请先加载文件")
            return
        cols = self._get_fill_columns()
        if cols is None:
            return
        threshold = self.THRESHOLDS[self._threshold_var.get()][1]

        self.btn_run.config(state="disabled")
        self.progress["value"] = 0
        self._log_write(
            f"开始处理 | 阈值={threshold*100:.0f}% | 填充列={cols}")

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
                    f"完成!\n总基因数: {result['total_genes']}\n"
                    f"保留: {result['kept_genes']}\n"
                    f"删除: {len(result['removed_genes'])}\n"
                    f"已生成文件:\n  - {files['data']}\n"
                    f"  - {files['report']}\n"
                    f"  - {files['summary']}\n"))
                self._post(0, lambda: self.progress.configure(value=100))
            except Exception as exc:
                self._post(0, lambda: self._log_write(f"[错误] {exc}"))
                self._post(0, lambda: messagebox.showerror("处理失败", str(exc)))
            finally:
                self._post(0, lambda: self.btn_run.config(state="normal"))

        threading.Thread(target=task, daemon=True).start()

    def _post(self, delay_ms, fn):
        """线程安全的 GUI 更新, 把回调投递到 Tk 主线程."""
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
            "关于",
            "KaKs 结果分析与筛选工具 v2.0\n\n"
            "功能:\n"
            " • 加载 .xlsx 整理表 或 原始 .kaks 文本\n"
            " • 按缺失比例阈值(20%/30%/50%/70%)删除基因\n"
            " • 用本基因其他物种有效值的中位数填充缺失\n"
            " • 导出: 筛选后数据 + 各基因报告 + 文本报告\n\n"
            "缺失识别:\n"
            " • Ka/Ks 列的 NaN 或 0 (Ka/Ks=0 是 Ka 计算失败的退化值)\n"
            " • Ka, Ks 列的 NaN\n\n"
            "中位数填充:\n"
            " • 仅基于有效值计算, 不受缺失/退化值影响\n\n"
            "依赖: pandas / openpyxl (标准库: tkinter)")

    def _log_write(self, msg):
        self._log.insert("end", f"[{datetime.now():%H:%M:%S}] {msg}\n")
        self._log.see("end")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    app = App()
    app.mainloop()