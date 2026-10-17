# -*- coding: utf-8 -*-
"""
FASTA 标识头重命名工具 (Fasta Header Renamer)
------------------------------------------------
功能:
1. 选择文件夹, 自动识别 *.fasta / *.fa / *.fas / *.fna 等序列文件
2. 解析每个文件中的序列条数, 并展示首条序列的标识头
3. 在主界面表格中可单独编辑"新标识头"
4. 一键批量应用所有新标识头; 一键将文件名同步为首条序列的新标识头
5. 对于多条序列的文件, 可进入"多序列编辑器"逐条修改, 并支持:
   - 仅修改所选序列的头部
   - 一键将文件名同步为首条序列的新标识头
6. 写盘采用临时文件 -> 原子替换, 避免中途崩溃导致原文件损坏
7. 主流程在后台线程执行, UI 不卡顿

依赖: 仅使用 Python 3.8+ 标准库 (tkinter)
"""

import os
import re
import sys
import shutil
import threading
import tempfile
import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import OrderedDict


# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------
FASTA_EXTS = (".fasta", ".fa", ".fas", ".fna", ".ffn", ".faa", ".frn")
HEADER_PREFIX = ">"
WINDOW_TITLE = "FASTA 标识头重命名工具"
APP_VERSION = "v1.0"


# ---------------------------------------------------------------------------
# 解析与写盘相关函数
# ---------------------------------------------------------------------------
def list_fasta_files(folder):
    """递归列出文件夹下所有 FASTA 文件, 返回相对于 folder 的相对路径列表"""
    results = []
    for root, _, files in os.walk(folder):
        for f in files:
            if f.lower().endswith(FASTA_EXTS):
                full = os.path.join(root, f)
                rel = os.path.relpath(full, folder)
                results.append(rel)
    results.sort(key=lambda x: x.lower())
    return results


def parse_fasta(path):
    """
    解析一个 FASTA 文件.
    返回 OrderedDict: {header_id: {"header": 完整标识头(不含>和首尾空白),
                                  "seq_lines": [行, ...],
                                  "seq_text": 序列拼接字符串}}
    注: 同一标识头出现多次时会自动加 _2 _3 后缀以保证 key 唯一.
    """
    records = OrderedDict()
    cur_key = None
    cur_header = None
    cur_lines = None

    def _commit():
        """把当前累积的序列行提交到 records"""
        if cur_key is None:
            return
        records[cur_key] = {
            "header": cur_header or "",
            "seq_lines": list(cur_lines or []),
            "seq_text": "".join(cur_lines or []),
        }

    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        for line in fp:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(HEADER_PREFIX):
                # 先提交上一条
                _commit()
                header = line[len(HEADER_PREFIX):].strip()
                base_key = (header.split()[0] if header else "seq") or "seq"
                key = base_key
                n = 2
                while key in records:
                    key = f"{base_key}_{n}"
                    n += 1
                cur_key = key
                cur_header = header
                cur_lines = []
            else:
                if cur_key is not None:
                    cur_lines.append(line)
        # 文件末尾再 commit 一次
        _commit()

    return records


def write_fasta_atomic(path, records, header_order):
    """将 records 按 header_order 顺序写到 path, 原子替换."""
    tmp_fd, tmp_path = tempfile.mkstemp(
        prefix=".rename_", suffix=".tmp", dir=os.path.dirname(path) or None
    )
    try:
        with os.fdopen(tmp_fd, "w", encoding="utf-8", newline="\n") as fp:
            for key in header_order:
                rec = records[key]
                header = rec["header"].strip()
                fp.write(f">{header}\n")
                seq = rec.get("seq_text", "".join(rec.get("seq_lines", [])))
                # 按 70 字符一行换行 (FASTA 习惯)
                line_width = 70
                for i in range(0, len(seq), line_width):
                    fp.write(seq[i:i + line_width] + "\n")
        # 原子替换
        if os.name == "nt":
            # Windows 下 os.replace 默认是原子的
            if os.path.exists(path):
                os.remove(path)
            os.rename(tmp_path, path)
        else:
            os.replace(tmp_path, path)
    except Exception:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except Exception:
            pass
        raise


def sanitize_for_filename(name):
    """去掉/替换不适合做文件名的字符, 但保留中文/字母/数字/常见符号"""
    # Windows 禁用字符:  \ / : * ? " < > | 以及控制字符
    if not name:
        return "unnamed"
    # 先截断过长的
    name = name.strip()
    bad = '<>:"/\\|?*\n\r\t'
    for ch in bad:
        name = name.replace(ch, "_")
    # 控制字符
    name = "".join(ch for ch in name if ord(ch) >= 32)
    # 不以 . 或空格 结尾
    name = name.rstrip(" .")
    if not name:
        return "unnamed"
    return name


def pick_filename_stem(header, accession_only=True):
    """从一条标识头中提取建议的文件名主体 (不含扩展名).

    accession_only=True:  取首段空格之前的部分 (例如 "OQ851342.1")
    accession_only=False: 保留完整标识头 (最长 120 字符避免过长)
    """
    h = (header or "").strip()
    if not h:
        return "unnamed"
    if accession_only:
        stem = h.split()[0] if h.split() else "unnamed"
    else:
        stem = h[:120]
    return sanitize_for_filename(stem) or "unnamed"


def unique_target_path(target_full):
    """若 target_full 已存在, 自动追加 _1 _2 … 返回尚未被占用的新完整路径."""
    if not os.path.exists(target_full):
        return target_full
    base, ext = os.path.splitext(target_full)
    i = 1
    while os.path.exists(f"{base}_{i}{ext}"):
        i += 1
    return f"{base}_{i}{ext}"


def header_from_filename(filepath, keep_extension=False):
    """从文件路径提取建议的标识头.
    默认去掉扩展名 (如 NC01.fasta -> NC01);
    keep_extension=True 时包含扩展名 (NC01.fasta -> NC01.fasta).
    会经过 sanitize_for_filename 清洗非法字符.
    """
    base = os.path.basename(filepath or "")
    if not base:
        return "unnamed"
    if not keep_extension:
        base = os.path.splitext(base)[0]
    return sanitize_for_filename(base) or "unnamed"


# ---------------------------------------------------------------------------
# GUI 主程序
# ---------------------------------------------------------------------------
class FastaRenamerApp:
    def __init__(self, root):
        self.root = root
        self.root.title(f"{WINDOW_TITLE} {APP_VERSION}")
        self.root.geometry("1180x760")
        self.root.minsize(960, 620)

        self.folder = tk.StringVar()
        self.status = tk.StringVar(value="请选择一个包含 FASTA 文件的文件夹。")
        self.files_info = {}  # rel_path -> {path, full_path, header_count, first_key, first_header, records}
        # "改名为标识头"时, 是否仅使用标识首段 (即 accession / 第一个空格之前的部分)
        # 默认 True: 大多数场景下用 accession 作为文件名最实用
        self.use_accession_only = tk.BooleanVar(value=True)
        # "改标识头=文件名"时, 是否对多序列文件中的所有序列都改名 (默认仅首条)
        self.apply_to_all_seqs = tk.BooleanVar(value=False)
        # "改标识头=文件名"时, 是否保留扩展名 (默认去掉, 即 NC01.fasta -> NC01)
        self.keep_extension = tk.BooleanVar(value=False)

        self._build_ui()

    # ---------------------- UI 构建 ----------------------
    def _build_ui(self):
        style = ttk.Style()
        # 使用更现代的主题 (Windows 上 'vista' 一般可用)
        for theme in ("vista", "clam", "default"):
            try:
                style.theme_use(theme)
                break
            except tk.TclError:
                continue
        style.configure("Treeview", rowheight=24, font=("Microsoft YaHei", 10))
        style.configure("Treeview.Heading", font=("Microsoft YaHei", 10, "bold"))
        style.configure("TButton", padding=(8, 4))
        style.configure("Accent.TButton", padding=(10, 5))

        # 顶部: 文件夹选择
        top = ttk.Frame(self.root, padding=(10, 10, 10, 5))
        top.pack(fill=tk.X)
        ttk.Label(top, text="文件夹:").pack(side=tk.LEFT)
        entry = ttk.Entry(top, textvariable=self.folder)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(6, 6))
        ttk.Button(top, text="浏览…", command=self.browse_folder).pack(side=tk.LEFT)
        ttk.Button(top, text="重新扫描", command=self.reload_folder).pack(
            side=tk.LEFT, padx=(6, 0)
        )

        # 中部: 文件列表 + 操作按钮
        body = ttk.Frame(self.root, padding=(10, 0, 10, 0))
        body.pack(fill=tk.BOTH, expand=True)
        body.columnconfigure(0, weight=1)
        body.rowconfigure(0, weight=1)

        # 左: Treeview
        list_frame = ttk.LabelFrame(body, text="FASTA 文件列表 (双击进入多序列编辑器)")
        list_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 6))
        list_frame.rowconfigure(0, weight=1)
        list_frame.columnconfigure(0, weight=1)

        cols = ("filename", "seqcount", "firstheader", "newheader")
        self.tree = ttk.Treeview(
            list_frame, columns=cols, show="headings", selectmode="extended"
        )
        self.tree.heading("filename", text="文件名")
        self.tree.heading("seqcount", text="序列数")
        self.tree.heading("firstheader", text="当前首条标识头")
        self.tree.heading("newheader", text="新标识头 (留空 = 不改)")
        self.tree.column("filename", width=240, anchor="w")
        self.tree.column("seqcount", width=70, anchor="center")
        self.tree.column("firstheader", width=330, anchor="w")
        self.tree.column("newheader", width=360, anchor="w")

        vsb = ttk.Scrollbar(list_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")

        self.tree.bind("<Double-1>", self.on_double_click)
        # 右键菜单 - 也打开编辑器
        self.tree.bind("<Button-3>", self.on_right_click)

        # 右: 操作按钮 (可滚动,避免按钮过多时被遮挡)
        right_outer = ttk.Frame(body)
        right_outer.grid(row=0, column=1, sticky="nsew", padx=(6, 0))
        right_outer.columnconfigure(0, weight=1)
        right_outer.rowconfigure(0, weight=1)

        right_canvas = tk.Canvas(right_outer, borderwidth=0, highlightthickness=0)
        right_canvas.grid(row=0, column=0, sticky="nsew")
        right_vsb = ttk.Scrollbar(
            right_outer, orient="vertical", command=right_canvas.yview
        )
        right_vsb.grid(row=0, column=1, sticky="ns")
        right_canvas.configure(yscrollcommand=right_vsb.set)

        # 内层 LabelFrame (保留"操作"标题)
        btn_frame = ttk.LabelFrame(
            right_canvas,
            text="操作 (有选中时仅作用于选中行)",
            padding=(8, 4),
        )
        btn_window = right_canvas.create_window(
            (0, 0), window=btn_frame, anchor="nw"
        )

        def _sync_canvas_width(event):
            right_canvas.itemconfigure(btn_window, width=event.width)

        def _sync_scroll_region(event):
            right_canvas.configure(
                scrollregion=right_canvas.bbox("all")
            )

        right_canvas.bind("<Configure>", _sync_canvas_width)
        btn_frame.bind("<Configure>", _sync_scroll_region)

        # 鼠标滚轮支持 (仅当鼠标在右侧面板时)
        def _on_wheel(event):
            right_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        right_canvas.bind(
            "<Enter>",
            lambda e: right_canvas.bind_all("<MouseWheel>", _on_wheel),
        )
        right_canvas.bind(
            "<Leave>",
            lambda e: right_canvas.unbind_all("<MouseWheel>"),
        )

        btn_frame.columnconfigure(0, weight=1)

        # ====== 区块 1: 头 → 文件名 (已存在) ======
        ttk.Button(
            btn_frame, text="📝 应用新标识头", style="Accent.TButton",
            command=self.apply_new_headers,
        ).grid(row=0, column=0, sticky="ew", padx=8, pady=(8, 4))
        ttk.Button(
            btn_frame, text="📁 用标识头改名为文件名",
            command=self.rename_files_by_headers,
        ).grid(row=1, column=0, sticky="ew", padx=8, pady=4)
        ttk.Button(
            btn_frame, text="🧬 编辑多序列…",
            command=self.open_multi_seq_editor,
        ).grid(row=2, column=0, sticky="ew", padx=8, pady=4)
        ttk.Separator(btn_frame, orient="horizontal").grid(
            row=3, column=0, sticky="ew", padx=8, pady=8
        )
        ttk.Button(
            btn_frame, text="↻ 全部清空新标识头",
            command=self.clear_all_new_headers,
        ).grid(row=4, column=0, sticky="ew", padx=8, pady=4)
        ttk.Button(
            btn_frame, text="从首条标识头填充",
            command=self.fill_new_from_first,
        ).grid(row=5, column=0, sticky="ew", padx=8, pady=4)
        ttk.Separator(btn_frame, orient="horizontal").grid(
            row=6, column=0, sticky="ew", padx=8, pady=8
        )
        ttk.Button(
            btn_frame, text="✅ 选中行: 改名为标识头",
            command=self.rename_selected_by_headers,
        ).grid(row=7, column=0, sticky="ew", padx=8, pady=4)
        ttk.Button(
            btn_frame, text="✅ 选中行: 应用头并改文件名",
            command=self.apply_and_rename_selected,
        ).grid(row=8, column=0, sticky="ew", padx=8, pady=4)

        # ====== 区块 2: 文件名 → 头 (新) ======
        ttk.Separator(btn_frame, orient="horizontal").grid(
            row=9, column=0, sticky="ew", padx=8, pady=8
        )
        ttk.Label(
            btn_frame, text="↓ 反向: 文件名 → 标识头 ↓",
            font=("Microsoft YaHei", 9, "bold"),
        ).grid(row=10, column=0, sticky="w", padx=8)
        ttk.Button(
            btn_frame, text="✅ 选中行: 改标识头=文件名",
            style="Accent.TButton",
            command=self.rename_selected_headers_to_filenames,
        ).grid(row=11, column=0, sticky="ew", padx=8, pady=4)
        ttk.Button(
            btn_frame, text="📥 选中行: 用文件名填充新标识头",
            command=self.fill_new_headers_from_filename,
        ).grid(row=12, column=0, sticky="ew", padx=8, pady=4)
        # 反向模式选项
        reverse_frame = ttk.LabelFrame(
            btn_frame, text="反向改名选项", padding=(6, 4)
        )
        reverse_frame.grid(row=13, column=0, sticky="ew", padx=8, pady=(6, 4))
        reverse_frame.columnconfigure(0, weight=1)
        ttk.Checkbutton(
            reverse_frame,
            text="多序列文件中所有序列都改",
            variable=self.apply_to_all_seqs,
        ).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(
            reverse_frame,
            text="标识头保留扩展名 (如 NC01.fasta)",
            variable=self.keep_extension,
        ).grid(row=1, column=0, sticky="w")
        ttk.Button(
            reverse_frame, text="模式说明…",
            command=self.show_header_change_mode_hint,
        ).grid(row=2, column=0, sticky="w", pady=(2, 0))

        # ====== 区块 3: 头→文件名 命名模式 ======
        ttk.Separator(btn_frame, orient="horizontal").grid(
            row=14, column=0, sticky="ew", padx=8, pady=8
        )
        mode_frame = ttk.Frame(btn_frame)
        mode_frame.grid(row=15, column=0, sticky="ew", padx=8, pady=(0, 4))
        ttk.Checkbutton(
            mode_frame,
            text="仅用 accession 首段作文件名",
            variable=self.use_accession_only,
        ).pack(anchor="w")
        ttk.Button(
            mode_frame, text="命名模式说明…",
            command=self.show_naming_mode_hint,
        ).pack(anchor="w", pady=(2, 0))

        # ====== 区块 4: 一键 + 导出 ======
        ttk.Separator(btn_frame, orient="horizontal").grid(
            row=16, column=0, sticky="ew", padx=8, pady=8
        )
        ttk.Button(
            btn_frame, text="💾 应用并重命名 (一键)",
            command=self.apply_headers_and_rename,
        ).grid(row=17, column=0, sticky="ew", padx=8, pady=4)
        ttk.Button(
            btn_frame, text="导出全部到新文件夹…",
            command=self.export_all_to_new_folder,
        ).grid(row=18, column=0, sticky="ew", padx=8, pady=(4, 8))

        # 下部: 状态栏
        bar = ttk.Frame(self.root, padding=(10, 6, 10, 8))
        bar.pack(fill=tk.X)
        self.progress = ttk.Progressbar(bar, mode="determinate", length=200)
        self.progress.pack(side=tk.RIGHT)
        ttk.Label(bar, textvariable=self.status).pack(side=tk.LEFT, fill=tk.X, expand=True)

        # 编辑态 - 用户在 Treeview 单元格上原地编辑
        self.tree.bind("<Double-1>", self.on_double_click)

    # ---------------------- 事件 ----------------------
    def browse_folder(self):
        d = filedialog.askdirectory(title="选择包含 FASTA 文件的文件夹")
        if d:
            self.folder.set(d)
            self.reload_folder()

    def reload_folder(self):
        folder = self.folder.get().strip()
        if not folder or not os.path.isdir(folder):
            messagebox.showerror("错误", "请先选择有效的文件夹")
            return
        # 在后台线程扫描, 避免大文件夹卡死 UI
        self.set_status(f"正在扫描 {folder} …")
        self.progress.config(mode="indeterminate")
        self.progress.start(10)
        threading.Thread(
            target=self._scan_worker, args=(folder,), daemon=True
        ).start()

    def _scan_worker(self, folder):
        try:
            files = list_fasta_files(folder)
            infos = OrderedDict()
            for rel in files:
                full = os.path.join(folder, rel)
                try:
                    records = parse_fasta(full)
                    keys = list(records.keys())
                    first_header = records[keys[0]]["header"] if keys else ""
                    infos[rel] = {
                        "full_path": full,
                        "records": records,
                        "header_count": len(keys),
                        "first_key": keys[0] if keys else None,
                        "first_header": first_header,
                    }
                except Exception as e:
                    infos[rel] = {
                        "full_path": full,
                        "records": OrderedDict(),
                        "header_count": 0,
                        "first_key": None,
                        "first_header": f"<解析失败: {e}>",
                        "error": True,
                    }
            self._schedule_ui(lambda: self._apply_scan_result(folder, infos))
        except Exception:
            err = traceback.format_exc()
            self._schedule_ui(lambda: messagebox.showerror("扫描失败", err))

    def _schedule_ui(self, func):
        """线程安全地调度 UI 更新 (主线程 mainloop 中有效; mainloop 未运行时直接调用)"""
        try:
            self.root.after(0, func)
        except RuntimeError:
            # mainloop 未启动 (如测试环境) — 直接调用
            try:
                func()
            except Exception:
                pass

    def _apply_scan_result(self, folder, infos):
        self.files_info = infos
        # 清空 tree
        for iid in self.tree.get_children():
            self.tree.delete(iid)
        # 填入
        for idx, (rel, info) in enumerate(infos.items(), 1):
            seqcnt = info.get("header_count", 0)
            self.tree.insert(
                "", "end", iid=rel,
                values=(rel, seqcnt, info.get("first_header", ""), ""),
            )
        self.progress.stop()
        self.progress.config(mode="determinate", value=0)
        n = len(infos)
        total_seqs = sum(i.get("header_count", 0) for i in infos.values())
        self.set_status(f"已扫描 {folder}: 共 {n} 个 FASTA 文件, 合计 {total_seqs} 条序列。")

    def fill_new_from_first(self):
        """将每行的新标识头默认填充为当前首条标识头 (便于快速批量小改)"""
        count = 0
        for iid in self.tree.get_children():
            vals = list(self.tree.item(iid, "values"))
            if vals:
                cur = vals[2]
                vals[3] = cur
                self.tree.item(iid, values=vals)
                count += 1
        self.set_status(f"已把 {count} 个文件的“新标识头”填充为当前首条标识头。")

    def clear_all_new_headers(self):
        for iid in self.tree.get_children():
            vals = list(self.tree.item(iid, "values"))
            vals[3] = ""
            self.tree.item(iid, values=vals)
        self.set_status("已清空所有“新标识头”。")

    def on_double_click(self, event=None):
        """双击进入多序列编辑器"""
        sel = self.tree.selection()
        if not sel:
            return
        rel = sel[0]
        info = self.files_info.get(rel)
        if not info or info.get("error"):
            messagebox.showerror("错误", "此文件解析失败, 无法编辑。")
            return
        MultiSeqEditor(self, rel, info)

    def on_right_click(self, event):
        iid = self.tree.identify_row(event.y)
        if iid:
            self.tree.selection_set(iid)
            self.on_double_click()

    # ---------------------- 写盘 ----------------------
    def _write_one_file(self, rel, header_updates, header_order=None):
        """对单个文件应用 header 更新 (dict: key -> new_header).
        header_order: 期望的写入顺序 (OrderedDict 的 keys 顺序即可). 默认为 records 当前顺序.
        返回 (success, message)
        """
        info = self.files_info[rel]
        records = info["records"]
        full_path = info["full_path"]

        updated_keys = []
        for key, new_h in header_updates.items():
            if key in records:
                old = records[key]["header"]
                new = new_h.strip()
                if new != old:
                    records[key]["header"] = new
                    updated_keys.append(key)

        if not updated_keys:
            return True, "无改动"

        order = header_order or list(records.keys())
        write_fasta_atomic(full_path, records, order)

        # 刷新缓存
        keys = list(records.keys())
        info["first_header"] = records[keys[0]]["header"] if keys else ""

        return True, f"已更新 {len(updated_keys)} 条序列头"

    def apply_new_headers(self):
        """应用主界面“新标识头” - 首条序列改名 (智能作用域: 有选中则仅作用于选中)."""
        targets = []
        scope_label = "当前选中行"
        for iid, info, vals in self._iter_target_files(force_selected=False):
            new_h = vals[3].strip()
            if new_h:
                targets.append((iid, vals[0], new_h))
        if not targets:
            messagebox.showinfo("提示", "没有需要应用的“新标识头”。")
            return
        sel = self.tree.selection()
        if sel:
            scope_label = f"选中的 {len(sel)} 行"
        else:
            scope_label = f"全部 {len(self.tree.get_children())} 行"
        if not messagebox.askyesno(
            "确认", f"将对 {scope_label} 中 {len(targets)} 个文件应用新标识头, 是否继续?"
        ):
            return
        # 执行
        ok, fail = 0, 0
        for rel, _, new_h in targets:
            info = self.files_info[rel]
            if not info.get("records") or info.get("error"):
                fail += 1
                continue
            first_key = info["first_key"]
            if not first_key:
                fail += 1
                continue
            try:
                self._write_one_file(rel, {first_key: new_h})
                ok += 1
            except Exception as e:
                fail += 1
                print("写入失败:", rel, e)
        # 刷新显示
        self._refresh_view()
        messagebox.showinfo("完成", f"成功: {ok}, 失败: {fail}")

    def _refresh_view(self):
        for iid in self.tree.get_children():
            info = self.files_info.get(iid)
            if not info:
                continue
            vals = list(self.tree.item(iid, "values"))
            keys = list(info["records"].keys())
            vals[1] = info.get("header_count", 0)
            vals[2] = info.get("first_header", "")
            self.tree.item(iid, values=vals)

    # ---------------------- 重命名逻辑 ----------------------
    def _iter_target_files(self, force_selected=False):
        """返回一个迭代器, 产出 (iid, info, vals).
        若有选中且 force_selected=False: 仅产出选中的行.
        若无选中: 产出全部.
        若 force_selected=True: 即便无选中也不返回任何东西.
        """
        iids = list(self.tree.selection())
        if not iids and force_selected:
            return
        if not iids:
            iids = list(self.tree.get_children())
        for iid in iids:
            info = self.files_info.get(iid)
            if not info:
                continue
            vals = list(self.tree.item(iid, "values"))
            yield iid, info, vals

    def _collect_rename_pairs(self, force_selected=False):
        """构造待重命名文件列表 [(iid, old_full, new_full, header), ...].
        使用 self.use_accession_only 决定文件名主体.
        跳过无法生成有效名、或新名==旧名的情况.
        """
        accession_only = bool(self.use_accession_only.get())
        pairs = []
        for iid, info, vals in self._iter_target_files(force_selected=force_selected):
            target_header = (vals[3].strip() or info.get("first_header", "")).strip()
            if not target_header:
                continue
            stem = pick_filename_stem(target_header, accession_only=accession_only)
            if not stem or stem == "unnamed":
                # 退化: 跳过不能成为合法文件名的
                continue
            ext = os.path.splitext(info["full_path"])[1]
            old_full = info["full_path"]
            new_full = os.path.join(os.path.dirname(old_full), stem + ext)
            pairs.append((iid, old_full, new_full, target_header, stem))
        return pairs

    def _perform_renames(self, pairs, *, after_msg=None):
        """执行 os.rename, 自动处理重名冲突 (加 _1 _2…)."""
        if not pairs:
            messagebox.showinfo("提示", "没有可重命名的文件。")
            return
        msg = f"将重命名 {len(pairs)} 个文件:\n"
        # 仅展示最多 8 个示例
        for _, _, new_full, header, _ in pairs[:8]:
            msg += f"  · {header[:40]} → {os.path.basename(new_full)}\n"
        if len(pairs) > 8:
            msg += f"  …还有 {len(pairs) - 8} 个\n"
        msg += "\n文件名重命名是不可逆操作, 确认继续?"
        if not messagebox.askyesno("确认", msg):
            return
        ok = skip = fail = 0
        for _, old_full, new_full, _, _ in pairs:
            if new_full == old_full:
                skip += 1
                continue
            target = unique_target_path(new_full)
            try:
                os.rename(old_full, target)
                ok += 1
            except Exception:
                fail += 1
        # 文件已被改名, 重新扫描保证 Treeview 与磁盘一致
        body = f"重命名完成: 成功 {ok}, 跳过(同名) {skip}, 失败 {fail}"
        if after_msg:
            body += "\n" + after_msg
        body += "\n已自动重新扫描文件夹。"
        messagebox.showinfo("完成", body)
        self.reload_folder()

    # ---------------------- 重命名按钮 ----------------------
    def rename_files_by_headers(self):
        """智能: 有选中则仅作用于选中; 否则作用于全部 (无新标识头也能用当前首条头)."""
        pairs = self._collect_rename_pairs(force_selected=False)
        if not pairs:
            messagebox.showinfo("提示", "没有可重命名的目标。\n可能原因: 文件没有首条标识头、或重命名后无变化。")
            return
        self._perform_renames(pairs)

    def rename_selected_by_headers(self):
        """强制只作用于当前选中的行; 若无选中给出提示."""
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("提示", "请先在左侧列表中选择 1 个或多个文件 (按住 Ctrl/Shift 多选)。")
            return
        pairs = self._collect_rename_pairs(force_selected=True)
        if not pairs:
            messagebox.showinfo("提示", "所选文件中没有可重命名的目标。")
            return
        self._perform_renames(pairs)

    def apply_and_rename_selected(self):
        """仅对选中的行: 应用新标识头, 然后把文件改名为该标识头."""
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("提示", "请先选择要处理的行 (Ctrl/Shift 多选)。")
            return
        # 1) 仅对选中行应用 新标识头
        applied = 0
        for iid, info, vals in self._iter_target_files(force_selected=True):
            new_h = vals[3].strip()
            if new_h:
                first_key = info.get("first_key")
                if not first_key or info.get("error"):
                    continue
                try:
                    self._write_one_file(iid, {first_key: new_h})
                    applied += 1
                except Exception:
                    pass
        # 2) 仅对选中行做重命名
        pairs = self._collect_rename_pairs(force_selected=True)
        if not pairs and applied == 0:
            messagebox.showinfo("提示", "选中行没有可应用的“新标识头”或可重命名的目标。")
            return
        if pairs:
            self._perform_renames(pairs, after_msg=f"已先将 {applied} 个文件的首条标识头写入。")
        else:
            self._refresh_view()
            messagebox.showinfo("完成", f"已写入 {applied} 个文件的首条标识头, 但无可重命名目标。")

    def apply_headers_and_rename(self):
        """一键: 先应用所有新标识头, 再批量重命名 — 智能作用于选中/全部."""
        # 1) 应用 (智能作用域)
        targets = []
        for iid, info, vals in self._iter_target_files(force_selected=False):
            new_h = vals[3].strip()
            if new_h and not info.get("error"):
                targets.append((iid, new_h))
        if targets:
            for rel, new_h in targets:
                info = self.files_info[rel]
                first_key = info.get("first_key")
                if first_key:
                    try:
                        self._write_one_file(rel, {first_key: new_h})
                    except Exception:
                        pass
        # 2) 重命名 (智能作用域)
        pairs = self._collect_rename_pairs(force_selected=False)
        if not pairs:
            self._refresh_view()
            messagebox.showinfo("完成", f"已应用 {len(targets)} 个文件的标识头, 但没有可重命名的目标。")
            return
        self._perform_renames(pairs, after_msg=f"已应用 {len(targets)} 个文件的首条标识头。")

    def show_naming_mode_hint(self):
        """弹窗解释当前命名模式与切换方式."""
        is_acc = bool(self.use_accession_only.get())
        msg = (
            "当前命名模式: "
            + ("仅 accession 首段 (例如 OQ851342.1.fasta)"
               if is_acc else
               "完整标识头 (例如 OQ851342.1 Asparagus officinalis ... .fasta)")
            + "\n\n"
            + "勾选/取消勾选上面那个复选框即可切换。"
            + "\n推荐使用 accession 模式: 文件名更短、更易于命令行引用。"
        )
        messagebox.showinfo("文件名命名模式", msg)

    # ---------------------- 改标识头 = 文件名 (反向) ----------------------
    def _new_headers_from_filename(self, info, apply_all=None, keep_ext=None):
        """根据一个文件的 records, 返回 OrderedDict 形式的 {key: new_header}.
        仅首条序列: 返回 {first_key: filename}.
        所有序列 (apply_all=True):
            - 每个序列的 header -> {filename}, {filename}_2, {filename}_3 ...
            - 第一条 仍为 {filename} (不带后缀)
        若首条序列原本就是这个文件名, 自动视为 "无变化" (返回空 dict).
        """
        if apply_all is None:
            apply_all = bool(self.apply_to_all_seqs.get())
        if keep_ext is None:
            keep_ext = bool(self.keep_extension.get())

        records = info.get("records") or {}
        keys = list(records.keys())
        if not keys:
            return {}

        full = info.get("full_path", "")
        base = header_from_filename(full, keep_extension=keep_ext)
        # 默认仅改第一条
        first_key = keys[0]
        first_old = records[first_key]["header"]
        if base == first_old:
            # 首条无变化
            result = {}
        else:
            result = {first_key: base}

        if apply_all and len(keys) > 1:
            idx = 2
            for k in keys[1:]:
                candidate = f"{base}_{idx}"
                # 避免与已生成的命名冲突 (极端少见)
                while any(v == candidate for v in result.values()):
                    candidate = f"{base}_{idx}"
                    idx += 1
                if candidate != records[k]["header"]:
                    result[k] = candidate
                idx += 1
        return result

    def fill_new_headers_from_filename(self):
        """用文件名 (可去掉扩展名) 填充选中行的“新标识头”列, 便于预览/编辑后再点“应用”."""
        keep_ext = bool(self.keep_extension.get())
        apply_all = bool(self.apply_to_all_seqs.get())
        scope_label = "当前选中行"
        n_filled = 0
        files_info = {}
        for iid, info, vals in self._iter_target_files(force_selected=False):
            full = info.get("full_path", "")
            if not full or info.get("error") or not info.get("records"):
                continue
            base = header_from_filename(full, keep_extension=keep_ext)
            vals[3] = base
            self.tree.item(iid, values=vals)
            files_info[iid] = (info, base)
            n_filled += 1
        if n_filled == 0:
            messagebox.showinfo("提示", "没有可填充的行, 请先选择文件。")
            return
        sel = self.tree.selection()
        if sel:
            scope_label = f"选中的 {len(sel)} 行"
        else:
            scope_label = f"全部 {len(self.tree.get_children())} 行"
        if apply_all:
            self.set_status(
                f"已把 {scope_label} 的“新标识头”填充为文件名。"
                f"\n  多序列文件中其余序列的标识头可在多序列编辑器中查看/修改。"
            )
        else:
            self.set_status(
                f"已把 {scope_label} 的“新标识头”填充为文件名 (仅首条序列生效)。"
            )

    def rename_headers_to_filenames(self):
        """智能: 有选中仅作用于选中; 否则作用于全部 (把标识头改为文件名)."""
        keep_ext = bool(self.keep_extension.get())
        apply_all = bool(self.apply_to_all_seqs.get())
        # 收集每个文件的 updates
        updates_per_file = {}  # iid -> {key: new_h}
        for iid, info, vals in self._iter_target_files(force_selected=False):
            if info.get("error") or not info.get("records"):
                continue
            up = self._new_headers_from_filename(
                info, apply_all=apply_all, keep_ext=keep_ext
            )
            if up:
                updates_per_file[iid] = up
        if not updates_per_file:
            messagebox.showinfo("提示", "所选文件的首条标识头已经与文件名一致, 无需修改。")
            return
        if not self._confirm_header_change(updates_per_file):
            return
        ok = fail = 0
        total_changed = 0
        for iid, updates in updates_per_file.items():
            try:
                self._write_one_file(iid, updates)
                ok += 1
                total_changed += len(updates)
            except Exception:
                fail += 1
        self._refresh_view()
        messagebox.showinfo(
            "完成",
            f"成功: {ok} 个文件 (共 {total_changed} 条序列头), 失败: {fail}"
        )

    def rename_selected_headers_to_filenames(self):
        """强制只作用于当前选中行; 把标识头改为文件名."""
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo(
                "提示",
                "请先在左侧列表中选择 1 个或多个文件 (按住 Ctrl/Shift 多选)。"
            )
            return
        keep_ext = bool(self.keep_extension.get())
        apply_all = bool(self.apply_to_all_seqs.get())
        updates_per_file = {}
        for iid in sel:
            info = self.files_info.get(iid)
            if not info or info.get("error") or not info.get("records"):
                continue
            up = self._new_headers_from_filename(
                info, apply_all=apply_all, keep_ext=keep_ext
            )
            if up:
                updates_per_file[iid] = up
        if not updates_per_file:
            messagebox.showinfo(
                "提示",
                "所选文件的首条标识头已经与文件名一致, 无需修改。"
            )
            return
        if not self._confirm_header_change(updates_per_file):
            return
        ok = fail = 0
        total_changed = 0
        for iid, updates in updates_per_file.items():
            try:
                self._write_one_file(iid, updates)
                ok += 1
                total_changed += len(updates)
            except Exception:
                fail += 1
        self._refresh_view()
        messagebox.showinfo(
            "完成",
            f"成功: {ok} 个文件 (共 {total_changed} 条序列头), 失败: {fail}"
        )

    def _confirm_header_change(self, updates_per_file):
        """对即将做的“标识头 = 文件名”给出预览确认. updates_per_file: iid -> {key: new_h}."""
        msg = f"将把以下 {len(updates_per_file)} 个文件的首条标识头 (或多条) 改为文件名:\n"
        preview_items = list(updates_per_file.items())[:8]
        for iid, updates in preview_items:
            info = self.files_info.get(iid)
            if not info:
                continue
            vals = self.tree.item(iid, "values")
            fname = os.path.basename(info["full_path"])
            keys = list((info.get("records") or {}).keys())
            for k, new_h in updates.items():
                old = (info["records"][k]["header"]) if k in info["records"] else "?"
                msg += f"  · {fname}:  {old[:40]}  →  {new_h[:40]}\n"
            if len(updates) > 1:
                msg += f"      (本文件中其它 {len(updates) - 1} 条也已列出/将一并修改)\n"
        if len(updates_per_file) > 8:
            msg += f"  …还有 {len(updates_per_file) - 8} 个文件\n"
        msg += "\n写盘操作不可逆, 确认继续?"
        return messagebox.askyesno("确认", msg)

    def show_header_change_mode_hint(self):
        """弹窗提示如何切换 “改标识头 = 文件名” 的两种选项."""
        apply_all = bool(self.apply_to_all_seqs.get())
        keep_ext = bool(self.keep_extension.get())
        msg = (
            "改“标识头 = 文件名”模式:\n\n"
            "· 多序列文件: "
            + ("全部序列都改 (第一条 = 文件名, 后续 = 文件名_2, _3, …)"
               if apply_all else
               "仅改第一条序列 (推荐, 与文件名惯例一致)")
            + "\n"
            "· 文件名包含扩展名: "
            + ("是 (如 NC01.fasta)"
               if keep_ext else
               "否 (如 NC01, 推荐)")
            + "\n\n"
            "在右侧复选框勾选/取消勾选即可切换。"
        )
        messagebox.showinfo("改标识头 = 文件名 模式", msg)

    def export_all_to_new_folder(self):
        """把所有 (改后) 的 FASTA 导出到一个新文件夹, 文件名用首条标识头"""
        if not self.files_info:
            messagebox.showinfo("提示", "请先加载文件夹")
            return
        out = filedialog.askdirectory(title="选择导出目录")
        if not out:
            return
        ok = 0
        fail = 0
        for rel, info in self.files_info.items():
            records = info.get("records")
            if not records:
                continue
            try:
                keys = list(records.keys())
                first_h = records[keys[0]]["header"].strip()
                new_name = sanitize_for_filename(first_h or os.path.splitext(rel)[0])
                ext = os.path.splitext(info["full_path"])[1]
                dst = os.path.join(out, new_name + ext)
                base = os.path.splitext(dst)[0]
                i = 1
                while os.path.exists(dst):
                    dst = f"{base}_{i}{ext}"
                    i += 1
                write_fasta_atomic(dst, records, keys)
                ok += 1
            except Exception:
                fail += 1
        messagebox.showinfo("导出完成", f"成功: {ok}, 失败: {fail}")

    def open_multi_seq_editor(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("提示", "请在左侧列表中选择 (可多选) 一个文件后再点击此按钮。")
            return
        rel = sel[0]
        info = self.files_info.get(rel)
        if not info or info.get("error"):
            messagebox.showerror("错误", "该文件无法编辑。")
            return
        if info.get("header_count", 0) == 0:
            messagebox.showinfo("提示", "该文件中没有可识别的序列。")
            return
        MultiSeqEditor(self, rel, info)

    # ---------------------- 状态 ----------------------
    def set_status(self, msg):
        self.status.set(msg)


# ---------------------------------------------------------------------------
# 多序列编辑器
# ---------------------------------------------------------------------------
class MultiSeqEditor(tk.Toplevel):
    """弹窗: 对单个文件内每条序列的标识头进行单独编辑"""

    def __init__(self, master, rel, info):
        super().__init__(master.root)
        self.master_app = master
        self.rel = rel
        self.info = info
        self.records = info["records"]
        self.title(f"多序列编辑器 - {rel}  ({len(self.records)} 条序列)")
        self.geometry("900x560")
        self.minsize(700, 400)

        # 顶部: 文件信息
        top = ttk.Frame(self, padding=(10, 10, 10, 4))
        top.pack(fill=tk.X)
        ttk.Label(
            top,
            text=f"文件: {rel}    序列数: {len(self.records)}",
            font=("Microsoft YaHei", 10, "bold"),
        ).pack(side=tk.LEFT)

        # 顶部: 把首条序列改名并同步文件名用的输入框
        top2 = ttk.LabelFrame(self, text="首条序列同步重命名该文件 (可选)", padding=(10, 6))
        top2.pack(fill=tk.X, padx=10, pady=(4, 6))
        self.first_new_var = tk.StringVar(value=self.records[list(self.records.keys())[0]]["header"])
        self.first_use_acc = tk.BooleanVar(value=bool(master.use_accession_only.get()))
        ttk.Label(top2, text="新的首条标识头:").pack(side=tk.LEFT)
        ttk.Entry(top2, textvariable=self.first_new_var).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=(6, 6)
        )
        ttk.Checkbutton(
            top2, text="仅 accession 首段", variable=self.first_use_acc
        ).pack(side=tk.LEFT, padx=(4, 4))
        ttk.Button(top2, text="改名文件", command=self.rename_file).pack(side=tk.LEFT)

        # 中部: 可滚动序列列表
        body = ttk.Frame(self, padding=(10, 0, 10, 6))
        body.pack(fill=tk.BOTH, expand=True)
        body.rowconfigure(0, weight=1)
        body.columnconfigure(0, weight=1)

        canvas = tk.Canvas(body, highlightthickness=0)
        vsb = ttk.Scrollbar(body, orient="vertical", command=canvas.yview)
        self.inner = ttk.Frame(canvas)
        self.inner.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=self.inner, anchor="nw")
        canvas.configure(yscrollcommand=vsb.set)
        canvas.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")

        # 列标题
        ttk.Label(self.inner, text="#", width=4, font=("Microsoft YaHei", 9, "bold")).grid(
            row=0, column=0, padx=4
        )
        ttk.Label(self.inner, text="当前标识头", width=24, font=("Microsoft YaHei", 9, "bold")).grid(
            row=0, column=1, padx=4
        )
        ttk.Label(self.inner, text="新标识头", width=46, font=("Microsoft YaHei", 9, "bold")).grid(
            row=0, column=2, padx=4
        )
        ttk.Label(self.inner, text="序列长度", width=10, font=("Microsoft YaHei", 9, "bold")).grid(
            row=0, column=3, padx=4
        )

        self.entry_vars = {}  # key -> StringVar
        self.key_order = list(self.records.keys())
        for i, key in enumerate(self.key_order, 1):
            rec = self.records[key]
            ttk.Label(self.inner, text=str(i)).grid(row=i, column=0, padx=4, pady=2)
            ttk.Label(
                self.inner, text=rec["header"][:30] + ("…" if len(rec["header"]) > 30 else "")
            ).grid(row=i, column=1, sticky="w", padx=4, pady=2)
            var = tk.StringVar(value=rec["header"])
            ttk.Entry(self.inner, textvariable=var, width=60).grid(
                row=i, column=2, sticky="ew", padx=4, pady=2
            )
            self.entry_vars[key] = var
            ttk.Label(self.inner, text=str(len(rec["seq_text"]))).grid(
                row=i, column=3, padx=4, pady=2
            )

        self.inner.columnconfigure(2, weight=1)

        # 支持鼠标滚轮
        def _on_wheel(e):
            canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")
        self.bind_all("<MouseWheel>", _on_wheel)

        # 底部按钮
        bottom = ttk.Frame(self, padding=(10, 4, 10, 10))
        bottom.pack(fill=tk.X)
        ttk.Button(
            bottom, text="💾 保存修改到文件", style="Accent.TButton", command=self.save_all
        ).pack(side=tk.LEFT)
        ttk.Button(
            bottom, text="从“当前标识头”填充“新标识头”",
            command=self.fill_new_from_current,
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(
            bottom, text="清空所有“新标识头”", command=self.clear_new_headers
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(bottom, text="关闭", command=self.destroy).pack(side=tk.RIGHT)

    def fill_new_from_current(self):
        for key, var in self.entry_vars.items():
            cur = self.records[key]["header"]
            var.set(cur)

    def clear_new_headers(self):
        for var in self.entry_vars.values():
            var.set("")

    def save_all(self):
        updates = {}
        for key, var in self.entry_vars.items():
            new_h = var.get().strip()
            old_h = self.records[key]["header"]
            if new_h and new_h != old_h:
                updates[key] = new_h
        if not updates:
            messagebox.showinfo("提示", "未检测到任何修改。", parent=self)
            return
        if not messagebox.askyesno(
            "确认保存", f"将更新 {len(updates)} 条序列标识头, 是否继续?", parent=self
        ):
            return
        try:
            self.master_app._write_one_file(self.rel, updates, header_order=self.key_order)
        except Exception as e:
            messagebox.showerror("保存失败", str(e), parent=self)
            return
        # 刷新主界面
        self.master_app._refresh_view()
        messagebox.showinfo("保存成功", f"已写入 {self.rel}", parent=self)
        # 同步首条标识头显示
        first_key = self.key_order[0]
        self.first_new_var.set(self.records[first_key]["header"])

    def rename_file(self):
        new_header = self.first_new_var.get().strip()
        if not new_header:
            messagebox.showwarning("提示", "请先填写新的首条标识头。", parent=self)
            return
        old_full = self.info["full_path"]
        folder = os.path.dirname(old_full)
        ext = os.path.splitext(old_full)[1]
        # 与主窗口统一的命名模式
        stem = pick_filename_stem(new_header, accession_only=bool(self.first_use_acc.get()))
        new_full = os.path.join(folder, stem + ext)
        if new_full == old_full:
            messagebox.showinfo("提示", "新文件名与原文件名相同。", parent=self)
            return
        new_full = unique_target_path(new_full)
        try:
            os.rename(old_full, new_full)
        except Exception as e:
            messagebox.showerror("重命名失败", str(e), parent=self)
            return
        messagebox.showinfo(
            "成功", f"已重命名为:\n{os.path.basename(new_full)}", parent=self
        )
        self.master_app.reload_folder()
        self.destroy()


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------
def main():
    root = tk.Tk()
    app = FastaRenamerApp(root)
    if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
        app.folder.set(sys.argv[1])
        app.reload_folder()
    root.mainloop()


if __name__ == "__main__":
    main()
