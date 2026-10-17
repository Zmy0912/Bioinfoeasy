# -*- coding: utf-8 -*-
"""
CodonW2RSCUPlot —— 将 codonW 的 RSCU 结果汇总为 RSCU-Plot 应用所需输入格式。

工作流：codonW（按物种计算 RSCU）→  本工具（多物种合并为宽表 CSV）
       → 上传至 https://pcg-lab.shinyapps.io/RSCU-Plot/ 进行可视化

输入：每个物种一个 RSCU 表文件（MEGA 风格的 .blk 宽表 / codonW 风格的
       .out 文本表），文件名（去扩展名）将作为物种列名。
输出：CSV，宽表布局
        列 1 = Amino Acid   （氨基酸，三字母缩写）
        列 2 = Codon        （密码子，RNA 字母 U/T 均可）
        列 3+ = 各物种的 RSCU 数值
"""

import os
import re
import csv
import glob
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


# ---------- 解析器 ----------

# MEGA 风格 / codonW 风格的密码子条目：可选氨基酸缩写 + 三联体密码子 + 计数 + RSCU
CODON_RE = re.compile(r"((?:[A-Za-z]{3})\s+)?([ATUCG]{3})\s+(\d+)\s+([\d.]+)")

# 标准密码子表顺序（与 RSCU-Plot 模板一致：按氨基酸分组）
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
    """将 DNA 字母 T 统一为 RNA 字母 U，便于跨格式汇总。"""
    return c.upper().replace("T", "U")


def _parse_entries_from_text(text):
    """从文本中提取 (氨基酸, 密码子, 计数, RSCU) 条目列表。

    支持 MEGA .blk 风格（每行 4 列密码子）和 codonW .out 风格
    （每行 1 组：AminoAcid  Codon  Count  RSCU）。
    """
    entries = []
    for line in text.splitlines():
        if " codons in " in line:        # 跳过 MEGA 的汇总说明行
            continue
        for m in CODON_RE.finditer(line):
            aa = m.group(1).strip() if m.group(1) else None
            codon = _normalize_codon(m.group(2))
            count = int(m.group(3))
            rscu = float(m.group(4))
            entries.append((aa, codon, count, rscu))
    return entries


def parse_input_file(filepath):
    """读取单个 RSCU 文件，返回 (氨基酸, 密码子, 计数, RSCU) 列表。"""
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    return _parse_entries_from_text(text)


def collect_input_files(input_dir):
    """收集文件夹中所有支持的 RSCU 结果文件。"""
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
    """汇总文件夹中所有 RSCU 结果文件到 RSCU-Plot 输入格式的 CSV。

    返回 (物种数, 密码子行数, 空文件列表)。
    """
    files = collect_input_files(input_dir)
    if not files:
        raise ValueError("所选文件夹中没有找到任何 RSCU 结果文件（.blk/.out/.txt/.csv）！")

    species_data = {}   # 物种名 -> {codon: rscu}
    codon_aa = {}       # 密码子 -> 氨基酸
    empty_files = []    # 未解析到数据的文件

    for i, fp in enumerate(files, 1):
        name = os.path.splitext(os.path.basename(fp))[0]
        if status_callback:
            status_callback("正在解析 {}/{}：{}".format(i, len(files), name))
        try:
            entries = parse_input_file(fp)
        except Exception as e:
            empty_files.append("{}（读取失败：{}）".format(name, e))
            continue
        if not entries:
            empty_files.append(name)
            if status_callback:
                status_callback("警告：{} 中未解析到任何数据".format(name))
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
        raise ValueError("所有文件都未能解析出 RSCU 数据，请检查输入格式。")

    # 行顺序：按标准密码子表；未知密码子追加到末尾
    known = {c for _, codons in STANDARD_ORDER for c in codons}
    ordered = [(aa, c) for aa, codons in STANDARD_ORDER for c in codons
               if any(c in d for d in species_data.values())]
    for d in species_data.values():
        for c in d:
            if c not in known and c not in [x[1] for x in ordered]:
                ordered.append((codon_aa.get(c, ""), c))

    if status_callback:
        status_callback("正在写入汇总表格...")

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
        self.root.title(APP_NAME + " —— codonW → RSCU-Plot 转换工具")
        self.root.geometry("680x440")
        self.root.minsize(600, 400)

        main = ttk.Frame(self.root, padding=15)
        main.pack(fill=tk.BOTH, expand=True)

        # ---- 输入文件夹 ----
        frm_in = ttk.LabelFrame(
            main,
            text="① 选择 codonW / MEGA RSCU 结果文件夹（每个物种一个文件）",
            padding=10)
        frm_in.pack(fill=tk.X, pady=(0, 10))
        self.var_input = tk.StringVar()
        ttk.Entry(frm_in, textvariable=self.var_input).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 8))
        ttk.Button(frm_in, text="浏览...", command=self.choose_input).pack(
            side=tk.LEFT)

        # ---- 输出 CSV ----
        frm_out = ttk.LabelFrame(
            main,
            text="② 设置输出 CSV 路径（将作为 RSCU-Plot 的输入）",
            padding=10)
        frm_out.pack(fill=tk.X, pady=(0, 10))
        self.var_output = tk.StringVar()
        ttk.Entry(frm_out, textvariable=self.var_output).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 8))
        ttk.Button(frm_out, text="另存为...", command=self.choose_output).pack(
            side=tk.LEFT)

        # ---- 选项 ----
        frm_opt = ttk.LabelFrame(main, text="选项", padding=10)
        frm_opt.pack(fill=tk.X, pady=(0, 10))
        self.var_include_ter = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            frm_opt, text="包含终止密码子（TER / Stop）",
            variable=self.var_include_ter).pack(anchor=tk.W)

        # ---- 运行按钮 ----
        self.btn_run = ttk.Button(main, text="③ 开始转换", command=self.run)
        self.btn_run.pack(pady=(0, 10))

        # ---- 日志 ----
        frm_log = ttk.LabelFrame(main, text="运行日志", padding=5)
        frm_log.pack(fill=tk.BOTH, expand=True)
        self.txt_log = tk.Text(frm_log, height=10, state=tk.DISABLED,
                               wrap=tk.WORD, font=("Consolas", 9))
        scroll = ttk.Scrollbar(frm_log, command=self.txt_log.yview)
        self.txt_log.configure(yscrollcommand=scroll.set)
        self.txt_log.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

    # ----- UI 辅助 -----
    def log(self, msg):
        self.txt_log.configure(state=tk.NORMAL)
        self.txt_log.insert(tk.END, msg + "\n")
        self.txt_log.see(tk.END)
        self.txt_log.configure(state=tk.DISABLED)

    def choose_input(self):
        d = filedialog.askdirectory(title="选择 RSCU 结果文件夹",
                                    initialdir=self.var_input.get() or os.getcwd())
        if d:
            self.var_input.set(d)
            # 自动建议输出路径
            if not self.var_output.get():
                self.var_output.set(os.path.join(d, "rscu_input_for_plot.csv"))

    def choose_output(self):
        f = filedialog.asksaveasfilename(
            title="保存汇总 CSV", defaultextension=".csv",
            filetypes=[("CSV 文件", "*.csv"), ("所有文件", "*.*")],
            initialfile="rscu_input_for_plot.csv",
            initialdir=(self.var_input.get() or os.getcwd()))
        if f:
            self.var_output.set(f)

    # ----- 执行 -----
    def run(self):
        input_dir = self.var_input.get().strip()
        output_csv = self.var_output.get().strip()

        if not input_dir or not os.path.isdir(input_dir):
            messagebox.showerror("错误", "请选择有效的 RSCU 结果文件夹！")
            return
        if not output_csv:
            messagebox.showerror("错误", "请设置输出 CSV 文件路径！")
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
        self.log("完成！共 {} 个物种、{} 个密码子行。".format(n_species, n_codons))
        if empty_files:
            self.log("警告（对应列保留为空）：{}".format(", ".join(empty_files)))
        self.log("结果已保存至：{}".format(output_csv))
        self.log("下一步：将该 CSV 上传至 https://pcg-lab.shinyapps.io/RSCU-Plot/")
        self.btn_run.configure(state=tk.NORMAL)
        msg = ("转换完成！\n\n物种数：{}\n密码子行数：{}\n\n文件：{}\n\n"
               "可直接上传至 RSCU-Plot 应用进行可视化。".format(
                   n_species, n_codons, output_csv))
        if empty_files:
            msg += "\n\n警告：以下文件为空或未解析到数据：\n{}".format(
                "\n".join(empty_files))
        messagebox.showinfo("完成", msg)

    def on_error(self, msg):
        self.log("出错：" + msg)
        self.btn_run.configure(state=tk.NORMAL)
        messagebox.showerror("错误", msg)


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