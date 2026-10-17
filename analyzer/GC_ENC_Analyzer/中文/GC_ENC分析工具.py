# -*- coding: utf-8 -*-
"""
GC / ENC 密码子使用偏好分析工具
====================================
输入一个文件夹，文件夹内每个 FASTA 文件视为一个物种（物种名 = 文件名），
文件内每条序列记录为一个基因（基因名自动从标题行提取，清除坐标、链方向、
join 等非基因名部分）。

对每个物种的每个 CDS 基因计算：
    GC1 (%) / GC2 (%) / GC3 (%)   —— 密码子第 1/2/3 位 GC 含量
    GC3s (%)                      —— 四倍简并密码子第 3 位 GC 含量
    GC12 (%)                      —— 第 1、2 位平均 GC 含量
    GC总含量 (%)                   —— 全序列 GC
    ENC                           —— 有效密码子数 (Wright, 1990)

并输出三个汇总表格：
    表1  所有物种 × 所有基因（每行一个 物种-基因 组合，具体数值）
    表2  所有物种（每行一个物种，值为该物种所有基因的 平均值 ± 标准差，
         GC12–GC3 slope 与 R² 为该物种内各基因 GC12 对 GC3 的回归）
    表3  所有基因（每行一个基因，值为该基因所有物种的 平均值 ± 标准差，
         GC12–GC3 slope 与 R² 为该基因在各物种间 GC12 对 GC3 的回归）

输出为 CSV（UTF-8-BOM，可直接用 Excel 打开）；若安装了 openpyxl，
则同时输出 xlsx 版本。结果保存在所选文件夹下的 "分析结果" 子文件夹中。

依赖：仅 Python 标准库（tkinter）；openpyxl 为可选。
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
# 遗传密码表（标准密码子表）
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

# 按氨基酸归并密码子家族（排除 Met、Trp、终止密码子）
AA_FAMILIES = {}  # aa -> [codon, ...]
for _c, _a in CODON2AA.items():
    if _a in ("*", "M", "W"):
        continue
    AA_FAMILIES.setdefault(_a, []).append(_c)

DEGENERACY = {aa: len(cods) for aa, cods in AA_FAMILIES.items()}
# ENC 公式中各简并度类的系数:  ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
ENC_COEF = {2: 9.0, 3: 1.0, 4: 5.0, 6: 3.0}

# 四倍简并密码子家族（前两位碱基）：CT,CC,CG,GC,AC,TC,GG,GT
FOURFOLD_PREFIXES = {"CT", "CC", "CG", "GC", "AC", "TC", "GG", "GT"}

VALID_BASES = set("ATCG")


# ---------------------------------------------------------------------------
# 基因名提取
# ---------------------------------------------------------------------------
_GENE_SEP_RE = re.compile(r"(_\[|\[|_join\{|_complement|\||\s|\()")


def extract_gene_name(header: str) -> str:
    """从 FASTA 标题行提取基因名，去除坐标、链方向、join 等非基因名部分。

    示例:
        >psbD_[33312:34373](+)|NC01.fasta|NC01.fasta  ->  psbD
        >rps16_join{[5754:5793](-),[4686:4894](-)}|.. ->  rps16
        >rbcL_[56315:57757](+)                        ->  rbcL
    """
    name = header.strip()
    m = _GENE_SEP_RE.search(name)
    if m:
        name = name[: m.start()]
    name = name.strip().strip("_-")
    # 只保留以字母开头、由字母数字下划线组成的部分
    m2 = re.match(r"^[A-Za-z][A-Za-z0-9_]*", name)
    if not m2:
        name = re.sub(r"[^A-Za-z0-9_]", "", name)
    return name


# ---------------------------------------------------------------------------
# FASTA 解析
# ---------------------------------------------------------------------------
def parse_fasta(path: str):
    """简易 FASTA 解析器，返回 [(header, seq), ...]"""
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
# 序列水平计算
# ---------------------------------------------------------------------------
def seq_to_codons(seq: str):
    """清理序列并切分为密码子列表（只保留 ACGT；截掉末尾不完整密码子）。"""
    seq = "".join(c for c in seq.upper() if c in VALID_BASES)
    n = len(seq) - len(seq) % 3
    return [seq[i:i + 3] for i in range(0, n, 3)]


def calc_gc_positions(codons):
    """返回 (gc1, gc2, gc3)：各密码子位置的 GC 百分比。"""
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
    """四倍简并密码子第 3 位的 GC 百分比。"""
    gc = tot = 0
    for codon in codons:
        if codon[:2] in FOURFOLD_PREFIXES:
            tot += 1
            if codon[2] in "GC":
                gc += 1
    return (gc / tot * 100) if tot else None


def calc_enc(codons):
    """有效密码子数 ENC (Wright 1990)。"""
    counts = {}
    for c in codons:
        counts[c] = counts.get(c, 0) + 1

    sum_by_deg = {2: [], 3: [], 4: [], 6: []}
    for aa, codon_list in AA_FAMILIES.items():
        n = sum(counts.get(c, 0) for c in codon_list)
        if n <= 0:
            continue
        k = DEGENERACY[aa]
        if n == 1:
            F = 1.0
        else:
            p2 = sum((counts.get(c, 0) / n) ** 2 for c in codon_list)
            F = (n * p2 - 1.0) / (n - 1.0)
        sum_by_deg[k].append(F)

    terms = 0.0
    n_classes = 0
    for k, flist in sum_by_deg.items():
        if flist:
            avg_f = sum(flist) / len(flist)
            # F=0（密码子完全均一使用）意味着 ENC 趋向 61，用极小值保护再截断
            terms += ENC_COEF[k] / (avg_f if avg_f > 1e-12 else 1e-12)
            n_classes += 1
    if n_classes == 0:
        return None
    enc = 2.0 + terms
    # ENC 理论范围 20 ~ 61
    if enc < 20:
        enc = 20.0
    elif enc > 61:
        enc = 61.0
    return enc


def gene_metrics(seq: str):
    """计算单个基因的全部指标，返回 dict（值可能为 None）。"""
    codons = seq_to_codons(seq)
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
# 回归
# ---------------------------------------------------------------------------
def linreg_slope_r2(xs, ys):
    """y 对 x 的最小二乘线性回归，返回 (slope, R²)。点数 < 3 时返回 (None, None)。"""
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
# 数值格式化
# ---------------------------------------------------------------------------
def fmt_val(v, nd=1):
    return "" if v is None else f"{v:.{nd}f}"


def fmt_mean_std(values, nd=1):
    """格式化为 'mean ± std'；样本不足时只显示平均值。"""
    vs = [v for v in values if v is not None]
    if not vs:
        return ""
    if len(vs) == 1:
        return f"{vs[0]:.{nd}f}"
    return f"{mean(vs):.{nd}f} ± {stdev(vs):.{nd}f}"


# ---------------------------------------------------------------------------
# 核心分析
# ---------------------------------------------------------------------------
def process_folder(folder: str, log=print, progress=None):
    """处理文件夹，返回 (明细行, 物种汇总行, 基因汇总行, 文件数, 基因组合数)。"""
    fasta_files = sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith((".fasta", ".fa", ".faa", ".fna", ".fst"))
    )
    if not fasta_files:
        raise RuntimeError("所选文件夹中没有找到 FASTA 文件（.fasta/.fa/.fna 等）。")

    detail_rows = []          # 每行: (species, gene, metrics)
    gene_order = []           # 基因出现顺序

    for idx, fname in enumerate(fasta_files, 1):
        species = os.path.splitext(fname)[0]
        path = os.path.join(folder, fname)
        try:
            records = parse_fasta(path)
        except Exception as e:
            log(f"[跳过] 无法读取 {fname}: {e}")
            continue

        ngene = 0
        for header, seq in records:
            gene = extract_gene_name(header)
            if not gene:
                log(f"[警告] {fname}: 无法识别基因名（标题行: {header[:60]}），已跳过")
                continue
            m = gene_metrics(seq)
            if m is None or m["ncodons"] < 3:
                log(f"[警告] {fname} / {gene}: 序列过短或无法解析，已跳过")
                continue
            detail_rows.append((species, gene, m))
            if gene not in gene_order:
                gene_order.append(gene)
            ngene += 1

        log(f"({idx}/{len(fasta_files)}) {fname}: {species} —— {ngene} 个基因")
        if progress:
            progress(idx, len(fasta_files))

    # ---------------- 表1：所有物种 × 所有基因（明细） ----------------
    table1 = []
    for species, gene, m in detail_rows:
        table1.append([
            species, gene, m["length"], m["ncodons"],
            fmt_val(m["GC1"]), fmt_val(m["GC2"]), fmt_val(m["GC3"]),
            fmt_val(m["GC3s"]), fmt_val(m["GC12"]), fmt_val(m["GCall"]),
            fmt_val(m["ENC"]),
        ])

    # ---------------- 表2：每个物种一行（对基因求均值） ----------------
    by_species = {}
    for species, gene, m in detail_rows:
        by_species.setdefault(species, []).append(m)

    species_names = sorted(by_species)
    table2 = []
    for species in species_names:
        ms = by_species[species]
        xs = [m["GC3"] for m in ms if m["GC3"] is not None]
        ys = [m["GC12"] for m in ms if m["GC3"] is not None and m["GC12"] is not None]
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

    # ---------------- 表3：每个基因一行（对物种求均值） ----------------
    by_gene = {}
    for species, gene, m in detail_rows:
        by_gene.setdefault(gene, []).append(m)

    table3 = []
    for gene in gene_order:
        ms = by_gene[gene]
        xs = [m["GC3"] for m in ms if m["GC3"] is not None]
        ys = [m["GC12"] for m in ms if m["GC3"] is not None and m["GC12"] is not None]
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


HEADER1 = ["Species", "Gene", "CDS长度", "密码子数",
           "GC1 (%)", "GC2 (%)", "GC3 (%)", "GC3s (%)", "GC12 (%)",
           "GC总含量 (%)", "ENC"]
HEADER2 = ["Species", "基因数",
           "GC1 (%)", "GC2 (%)", "GC3 (%)", "GC3s (%)", "GC12 (%)", "ENC",
           "GC12–GC3 slope", "R²"]
HEADER3 = ["Gene", "物种数",
           "GC1 (%)", "GC2 (%)", "GC3 (%)", "GC3s (%)", "GC12 (%)", "ENC",
           "GC12–GC3 slope", "R²"]


def save_results(out_dir, table1, table2, table3):
    """保存三个 CSV（若安装 openpyxl，则同时保存 xlsx）。"""
    os.makedirs(out_dir, exist_ok=True)

    p1 = os.path.join(out_dir, "1_所有物种所有基因_明细.csv")
    p2 = os.path.join(out_dir, "2_所有物种平均值.csv")
    p3 = os.path.join(out_dir, "3_所有基因平均值.csv")

    for path, header, rows in ((p1, HEADER1, table1),
                               (p2, HEADER2, table2),
                               (p3, HEADER3, table3)):
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)

    # 可选的 xlsx 输出
    xlsx_path = None
    try:
        from openpyxl import Workbook
        xlsx_path = os.path.join(out_dir, "GC_ENC分析汇总.xlsx")
        wb = Workbook()
        for sheet, header, rows in (
            ("所有物种所有基因", HEADER1, table1),
            ("物种平均值", HEADER2, table2),
            ("基因平均值", HEADER3, table3),
        ):
            ws = wb.active if sheet == "所有物种所有基因" else wb.create_sheet()
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
        self.title("GC / ENC 密码子使用偏好分析工具")
        self.geometry("780x560")
        self.minsize(680, 460)
        self._queue = queue.Queue()
        self.after(100, self._poll_queue)

        pad = {"padx": 10, "pady": 6}

        top = ttk.Frame(self)
        top.pack(fill="x", **pad)
        ttk.Label(top, text="序列文件夹：").pack(side="left")
        self.var_folder = tk.StringVar()
        ttk.Entry(top, textvariable=self.var_folder).pack(
            side="left", fill="x", expand=True, padx=6)
        ttk.Button(top, text="浏览…", command=self.choose_folder).pack(side="left")
        ttk.Button(top, text="开始分析", command=self.start).pack(side="left", padx=6)

        self.progress = ttk.Progressbar(self, mode="determinate")
        self.progress.pack(fill="x", padx=10)

        ttk.Label(
            self,
            text="说明：每个 FASTA 文件 = 一个物种（物种名取文件名）；文件内每条记录 = 一个基因，"
                 "基因名自动从标题行提取（清除坐标、链方向、join 等非基因名部分）。\n"
                 "GC12–GC3 slope / R²：GC12 对 GC3 的线性回归斜率与决定系数；"
                 "汇总表中各指标为 平均值 ± 标准差。",
            foreground="#555", justify="left", wraplength=740,
        ).pack(fill="x", padx=10)

        self.log_text = tk.Text(self, height=18, state="disabled",
                                font=("Consolas", 9))
        self.log_text.pack(fill="both", expand=True, padx=10, pady=8)

        bottom = ttk.Frame(self)
        bottom.pack(fill="x", padx=10, pady=(0, 8))
        self.var_status = tk.StringVar(value="就绪")
        ttk.Label(bottom, textvariable=self.var_status,
                  foreground="#0a7").pack(side="left")
        ttk.Button(bottom, text="打开结果文件夹",
                   command=self.open_outdir).pack(side="right")

        self.out_dir = None
        self._running = False

    # ---- 日志（线程安全，经队列转发到主线程） ----
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

    # ---- 交互 ----
    def choose_folder(self):
        d = filedialog.askdirectory(title="选择包含物种 FASTA 文件的文件夹")
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
            messagebox.showwarning("提示", "请先选择有效的序列文件夹！")
            return
        self._running = True
        self.var_status.set("正在分析…")
        self.progress.configure(maximum=100, value=0)
        threading.Thread(target=self._work, args=(folder,), daemon=True).start()

    def _work(self, folder):
        try:
            self.log("=" * 60)
            self.log(f"输入文件夹: {folder}")
            t1, t2, t3, nfile, npair = process_folder(
                folder, log=self.log, progress=self._progress)
            out_dir = os.path.join(folder, "分析结果")
            paths = save_results(out_dir, t1, t2, t3)
            self.out_dir = out_dir
            self.log("-" * 60)
            self.log(f"完成！共 {nfile} 个物种文件，{npair} 个 物种-基因 组合。")
            for p in paths:
                self.log(f"  输出: {p}")
            self._queue.put(("done", f"完成：{nfile} 个物种 / {npair} 个物种-基因组合，"
                                     f"结果已保存到 {out_dir}"))
        except Exception as e:
            err = "".join(traceback.format_exception_only(type(e), e)).strip()
            self.log(f"[错误] {err}")
            self._queue.put(("done", f"分析失败：{err}"))

    def _progress(self, cur, total):
        self._queue.put(("prog", cur, total))


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
