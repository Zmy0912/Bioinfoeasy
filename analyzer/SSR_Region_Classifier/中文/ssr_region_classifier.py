#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SSR 区域分类分析工具
====================

功能概述
--------
本程序读取 MISA 输出的 SSR 结果文件，结合同一物种对应的 GenBank (GB) 注释文件，
将每条 SSR 按其基因组坐标归类到三类区域之一：

  - CDS    : 编码区（CDS 外显子）
  - Intron : 基因内的内含子（由 gene 范围与 CDS 外显子差集推断）
  - IGS    : 基因间间隔区（既不在 CDS 也不在 Intron 的区域）

输出文件可直接使用 Excel、Origin、R、Python(matplotlib/seaborn) 等软件进行绘图。

输入要求
--------
1. MISA 文件：
   - 制表符分隔，至少包含列：ID、SSR nr.、SSR type、SSR、size、start、end
   - 第 1 列 ID（如 NC01、NC02 ...）与 GB 文件主名一一对应。

2. GB 文件夹：
   - 与 MISA 文件同名（同编号）的子文件夹，
   - 子文件夹内含 ID.gb（如 NC01.gb, NC02.gb, ...），
   - GB 文件需包含 FEATURES 表的 gene 与 CDS 特征。

输出
----
1. <out>/SSR_with_region.txt       ：MISA 原内容 + Region / RegionDetail 列
2. <out>/SSR_region_summary.tsv    ：每个物种各区域 SSR 数量与占比（长表，便于 ggplot）
3. <out>/SSR_region_summary_wide.tsv：每个物种各区域 SSR 数量（宽表）
4. <out>/SSR_type_by_region.tsv    ：每区域 × SSR 类型的数量（便于堆叠柱状图）
5. <out>/log.txt                   ：运行日志

使用
---
直接双击或在命令行中运行：
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
# 1. GB 文件解析
# ---------------------------------------------------------------------------

_LOCATION_RE = re.compile(r"(\d+)\.\.(\d+)")


def parse_location(loc: str):
    """
    解析 GenBank Location 字符串为 [(start, end), ...] 区间列表。
    支持单段、join 多段及其组合，含 complement 包装。

    返回的区间始终按 start..end 升序（不反转方向），方向信息单独返回。
    """
    raw = loc.strip()
    is_complement = "complement" in raw

    # 去除 complement(...) / join(...) 的外层包装
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
    读取 GB 文件的 FEATURES 区段，将可能跨行的特征合并为单行记录。
    返回 [(feature_key, location_str), ...] 的列表，保留文件中出现的顺序。
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
        # 遇到下一段起始（如 ORIGIN、//）则结束
        if line.startswith("ORIGIN") or line.startswith("//"):
            break
        m = re.match(r"^\s{5}(\S+)\s+(.*)$", line)
        if m:
            key = m.group(1)
            loc = m.group(2).strip()
            # 收集续行（以 21 个空格起头的为特征表的续行）
            j = i + 1
            while j < len(lines):
                cont = lines[j]
                if re.match(r"^\s{5}\S+\s", cont):
                    break  # 新的特征键
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
    解析单个 GB 文件，返回：
      {
        'length': 基因组长度,
        'cds_exons':   [(start, end), ...],   # 所有 CDS 外显子（合并所有 CDS 特征）
        'gene_ranges': [(start, end), ...],   # 所有 gene 总范围
        'cds_by_gene': OrderedDict{ gene_range -> [exon_ranges] }  # 用于推断内含子
      }
    """
    with open(gb_path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    length = parse_genome_length(text)
    records = _read_features_block(text)

    cds_exons = []
    gene_ranges = []
    cds_per_gene = OrderedDict()  # gene_range -> list[exon_range]

    # 先收集所有 gene 和 CDS
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

    # 关联 CDS 外显子到对应的 gene（用最大重叠来匹配）
    used = set()
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
    由 gene 与 CDS 外显子推断内含子区间。
    内含子严格定义为同一基因内、相邻 CDS 外显子之间的间隔（不含 5'/3'UTR）。
    """
    introns = []
    for gene_seg, exons in cds_by_gene.items():
        if not exons:
            continue
        # 只保留严格落在 gene 区间内的外显子（保险起见）
        in_gene = [e for e in exons if e[0] >= gene_seg[0] and e[1] <= gene_seg[1]]
        if len(in_gene) < 2:
            continue
        in_gene.sort()
        for i in range(len(in_gene) - 1):
            a_start, a_end = in_gene[i]
            b_start, b_end = in_gene[i + 1]
            if b_start > a_end + 1:
                # 使用 (a_end, b_start) 半开半闭区间；分类时按 start/end 落入判断
                introns.append((a_end, b_start))
    return introns


# ---------------------------------------------------------------------------
# 2. MISA 解析与区域分类
# ---------------------------------------------------------------------------

def parse_misa(misa_path: str):
    """读取 MISA 文件，返回 list[dict] 与表头。"""
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
    判断 SSR 属于哪个区域。优先匹配 CDS，再匹配 Intron，否则 IGS。
    完全包含（start >= seg_start 且 end <= seg_end）才算严格命中；
    否则按中心点判断。
    """
    # 1) 严格包含匹配
    for s, e in cds_exons:
        if s <= start and end <= e:
            return "CDS", f"CDS:{s}..{e}"
    for s, e in intron_intervals:
        if s <= start and end <= e:
            return "Intron", f"Intron:{s}..{e}"

    # 2) 按 SSR 中心点判断
    mid = (start + end) / 2.0
    best_label, best_key, best_dist = None, None, None
    for s, e in cds_exons:
        if s <= mid <= e:
            return "CDS", f"CDS:{s}..{e}(center)"
    for s, e in intron_intervals:
        if s <= mid <= e:
            return "Intron", f"Intron:{s}..{e}(center)"

    # 3) 都没命中则 IGS，附带最接近的基因注释
    return "IGS", "IGS"


# ---------------------------------------------------------------------------
# 3. 主分析流程
# ---------------------------------------------------------------------------

def analyze(misa_path: str, gb_dir: str, out_dir: str,
            log_func=print, cancel_flag=None):
    """
    执行完整分析流程。
    log_func(str): 日志回调
    cancel_flag: threading.Event()，用于中途取消
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

    log(f"开始分析")
    log(f"MISA 文件: {misa_path}")
    log(f"GB 文件夹: {gb_dir}")
    log(f"输出目录 : {out_dir}")

    header, rows = parse_misa(misa_path)
    if not rows:
        raise RuntimeError("MISA 文件为空或格式错误。")
    log(f"已加载 MISA 数据，共 {len(rows)} 条 SSR 记录。")

    # 提取所有唯一物种
    species_list = []
    seen = set()
    for r in rows:
        sp = r.get("ID", "").strip()
        if sp and sp not in seen:
            seen.add(sp)
            species_list.append(sp)
    log(f"检测到 {len(species_list)} 个物种: {', '.join(species_list)}")

    # 物种 -> GB 文件路径
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
            # 在子目录中再找一次
            for root_dir, _, files in os.walk(gb_dir):
                for f in files:
                    name, ext = os.path.splitext(f)
                    if name == sp and ext.lower() in (".gb", ".genbank"):
                        gb_path = os.path.join(root_dir, f)
                        break
                if gb_path:
                    break
        if gb_path is None:
            log(f"  [警告] 未找到物种 {sp} 对应的 GB 文件，已跳过该物种。")
            continue
        gb_cache[sp] = gb_path
        log(f"  - {sp} -> {gb_path}")

    if not gb_cache:
        raise RuntimeError("未匹配到任何 GB 文件，请检查文件夹路径与命名。")

    # 区域分类
    detail_rows = []
    stats_species_region = defaultdict(lambda: defaultdict(int))  # species -> region -> count
    stats_region_type = defaultdict(lambda: defaultdict(int))      # region -> type -> count
    stats_species_total = defaultdict(int)
    stats_species_unmapped = defaultdict(int)

    total = len(rows)
    for idx, r in enumerate(rows, 1):
        if cancel_flag is not None and cancel_flag.is_set():
            log("用户取消分析。")
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
            log(f"  进度: {idx}/{total} ({idx*100//total}%)")

    # ---- 输出 1：带区域信息的明细 ----
    out_detail = os.path.join(out_dir, "SSR_with_region.txt")
    cols = list(header) + ["Region", "RegionDetail"]
    with open(out_detail, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in detail_rows:
            w.writerow([r.get(c, "") for c in cols])
    log(f"已写出明细文件: {out_detail}")

    # ---- 输出 2：长表 (species, region, count, percent) ----
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
    log(f"已写出长表统计: {out_long}")

    # ---- 输出 3：宽表 (species × region) ----
    out_wide = os.path.join(out_dir, "SSR_region_summary_wide.tsv")
    with open(out_wide, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Species", "CDS", "Intron", "IGS", "Total"])
        for sp in species_order:
            cds = stats_species_region[sp].get("CDS", 0)
            intr = stats_species_region[sp].get("Intron", 0)
            igs = stats_species_region[sp].get("IGS", 0)
            w.writerow([sp, cds, intr, igs, cds + intr + igs])
    log(f"已写出宽表统计: {out_wide}")

    # ---- 输出 4：区域 × SSR 类型（堆叠柱状图友好） ----
    out_type = os.path.join(out_dir, "SSR_type_by_region.tsv")
    all_types = sorted({t for reg in stats_region_type for t in stats_region_type[reg]})
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
    log(f"已写出区域×SSR类型统计: {out_type}")

    log("分析完成！")
    log_fp.close()
    return True


# ---------------------------------------------------------------------------
# 4. Tkinter GUI
# ---------------------------------------------------------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SSR 区域分类分析工具 (CDS / Intron / IGS)")
        self.geometry("820x620")
        self.minsize(760, 560)

        # 状态
        self.misa_path = tk.StringVar()
        self.gb_dir = tk.StringVar()
        self.out_dir = tk.StringVar()
        self.status = tk.StringVar(value="就绪")
        self.cancel_event = threading.Event()
        self.worker = None

        self._build_ui()

    # ---- UI 构造 ----
    def _build_ui(self):
        pad = {"padx": 8, "pady": 4}

        frm_top = ttk.LabelFrame(self, text="输入与输出")
        frm_top.pack(fill="x", padx=10, pady=8)

        # MISA
        ttk.Label(frm_top, text="MISA 文件：").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(frm_top, textvariable=self.misa_path, width=70)\
            .grid(row=0, column=1, sticky="ew", **pad)
        ttk.Button(frm_top, text="选择…", command=self.choose_misa)\
            .grid(row=0, column=2, **pad)

        # GB 文件夹
        ttk.Label(frm_top, text="GB 文件夹：").grid(row=1, column=0, sticky="w", **pad)
        ttk.Entry(frm_top, textvariable=self.gb_dir, width=70)\
            .grid(row=1, column=1, sticky="ew", **pad)
        ttk.Button(frm_top, text="选择…", command=self.choose_gb)\
            .grid(row=1, column=2, **pad)

        # 输出目录
        ttk.Label(frm_top, text="输出目录：").grid(row=2, column=0, sticky="w", **pad)
        ttk.Entry(frm_top, textvariable=self.out_dir, width=70)\
            .grid(row=2, column=1, sticky="ew", **pad)
        ttk.Button(frm_top, text="选择…", command=self.choose_out)\
            .grid(row=2, column=2, **pad)

        frm_top.columnconfigure(1, weight=1)

        # 按钮区
        frm_btn = ttk.Frame(self)
        frm_btn.pack(fill="x", padx=10)
        self.btn_run = ttk.Button(frm_btn, text="开始分析", command=self.on_run)
        self.btn_run.pack(side="left", padx=4, pady=4)
        self.btn_cancel = ttk.Button(frm_btn, text="取消", command=self.on_cancel,
                                     state="disabled")
        self.btn_cancel.pack(side="left", padx=4, pady=4)
        ttk.Button(frm_btn, text="清空日志", command=self.clear_log)\
            .pack(side="left", padx=4, pady=4)
        ttk.Button(frm_btn, text="打开输出目录", command=self.open_out)\
            .pack(side="left", padx=4, pady=4)
        ttk.Button(frm_btn, text="使用说明", command=self.show_help)\
            .pack(side="right", padx=4, pady=4)

        # 进度条
        self.progress = ttk.Progressbar(self, mode="determinate")
        self.progress.pack(fill="x", padx=10, pady=2)

        # 日志
        frm_log = ttk.LabelFrame(self, text="运行日志")
        frm_log.pack(fill="both", expand=True, padx=10, pady=8)
        self.txt_log = scrolledtext.ScrolledText(frm_log, height=20,
                                                 font=("Consolas", 10))
        self.txt_log.pack(fill="both", expand=True, padx=4, pady=4)

        # 状态栏
        ttk.Label(self, textvariable=self.status, anchor="w",
                  relief="sunken").pack(fill="x", side="bottom")

    # ---- 路径选择 ----
    def choose_misa(self):
        path = filedialog.askopenfilename(
            title="选择 MISA 结果文件",
            filetypes=[("MISA files", "*.misa *.MISA *.txt *.tsv"),
                       ("All files", "*.*")])
        if path:
            self.misa_path.set(path)
            # 自动推荐 GB 文件夹与输出目录
            misa_dir = os.path.dirname(path)
            base = os.path.splitext(os.path.basename(path))[0]
            gb_guess = os.path.join(misa_dir, base)
            if os.path.isdir(gb_guess) and not self.gb_dir.get():
                self.gb_dir.set(gb_guess)
            if not self.out_dir.get():
                self.out_dir.set(os.path.join(misa_dir, base + "_region_result"))

    def choose_gb(self):
        d = filedialog.askdirectory(title="选择 GB 文件所在文件夹")
        if d:
            self.gb_dir.set(d)

    def choose_out(self):
        d = filedialog.askdirectory(title="选择结果输出目录")
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
                messagebox.showerror("打开失败", str(e))
        else:
            messagebox.showinfo("提示", "输出目录不存在或未设置。")

    # ---- 日志 ----
    def log(self, msg: str):
        self.txt_log.insert("end", msg + "\n")
        self.txt_log.see("end")
        self.update_idletasks()

    def clear_log(self):
        self.txt_log.delete("1.0", "end")

    # ---- 运行 ----
    def on_run(self):
        misa = self.misa_path.get().strip()
        gb = self.gb_dir.get().strip()
        out = self.out_dir.get().strip()

        if not misa or not os.path.isfile(misa):
            messagebox.showerror("错误", "请选择有效的 MISA 文件。")
            return
        if not gb or not os.path.isdir(gb):
            messagebox.showerror("错误", "请选择有效的 GB 文件夹。")
            return
        if not out:
            messagebox.showerror("错误", "请选择输出目录。")
            return
        os.makedirs(out, exist_ok=True)

        self.cancel_event.clear()
        self.btn_run.config(state="disabled")
        self.btn_cancel.config(state="normal")
        self.status.set("正在分析…")

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
            self.status.set("出错")
            messagebox.showerror("运行错误", str(err))
        elif ok:
            self.status.set("完成")
            messagebox.showinfo("完成",
                                f"分析完成！结果保存在：\n{self.out_dir.get()}")
        else:
            self.status.set("已取消")

    def on_cancel(self):
        self.cancel_event.set()
        self.status.set("正在取消…")

    # ---- 帮助 ----
    def show_help(self):
        help_text = (
            "使用说明\n"
            "========\n"
            "1. 选择 MISA 结果文件（*.misa / *.txt）。\n"
            "2. 选择与 MISA 同名的 GB 文件夹（内含 NC01.gb、NC02.gb 等）。\n"
            "3. 选择结果输出目录。\n"
            "4. 点击「开始分析」。\n\n"
            "输出文件：\n"
            "  - SSR_with_region.txt         每条 SSR 增加 Region 列\n"
            "  - SSR_region_summary.tsv      长表：物种 × 区域 × 数量 × 占比\n"
            "  - SSR_region_summary_wide.tsv 宽表：物种 × {CDS, Intron, IGS}\n"
            "  - SSR_type_by_region.tsv      区域 × SSR 类型（堆叠图）\n"
            "  - log.txt                     运行日志\n\n"
            "区域判断规则：\n"
            "  CDS    : SSR 完全落在 CDS 外显子区间内（或中心点在外显子内）\n"
            "  Intron : SSR 落在由 gene − CDS 推断的内含子中\n"
            "  IGS    : 其余均为基因间间隔区\n"
        )
        messagebox.showinfo("使用说明", help_text)


# ---------------------------------------------------------------------------
# 5. 入口
# ---------------------------------------------------------------------------

def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()