import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import re
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class FASTACleanerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA 文件验证和清理工具")
        self.root.geometry("800x600")
        self.root.minsize(700, 500)

        self.input_files = []
        self.setup_ui()

    def setup_ui(self):
        """设置用户界面"""
        # 创建主框架
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # 配置网格权重
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(7, weight=1)

        # 标题
        title_label = ttk.Label(
            main_frame,
            text="FASTA 文件验证和清理工具",
            font=("Arial", 16, "bold")
        )
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))

        # 文件选择区域
        ttk.Label(main_frame, text="选择 FASTA 文件:").grid(row=1, column=0, sticky=tk.W, pady=5)

        file_frame = ttk.Frame(main_frame)
        file_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        self.file_listbox = tk.Listbox(file_frame, height=6, selectmode=tk.MULTIPLE)
        self.file_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        scrollbar = ttk.Scrollbar(file_frame, orient=tk.VERTICAL, command=self.file_listbox.yview)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.file_listbox.configure(yscrollcommand=scrollbar.set)

        file_frame.columnconfigure(0, weight=1)

        # 按钮区域
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=3, column=0, columnspan=3, pady=10)

        ttk.Button(button_frame, text="添加文件", command=self.add_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="添加文件夹", command=self.add_folder).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="移除选中", command=self.remove_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清空列表", command=self.clear_files).pack(side=tk.LEFT, padx=5)

        # 输出目录
        ttk.Label(main_frame, text="输出目录:").grid(row=4, column=0, sticky=tk.W, pady=5)

        output_frame = ttk.Frame(main_frame)
        output_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        output_frame.columnconfigure(0, weight=1)

        self.output_path_var = tk.StringVar(value="./cleaned_output")
        self.output_entry = ttk.Entry(output_frame, textvariable=self.output_path_var)
        self.output_entry.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        ttk.Button(output_frame, text="浏览", command=self.browse_output_dir).grid(row=0, column=1)

        # 选项
        options_frame = ttk.LabelFrame(main_frame, text="处理选项", padding="10")
        options_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        self.validate_var = tk.BooleanVar(value=True)
        self.clean_var = tk.BooleanVar(value=True)
        self.create_report_var = tk.BooleanVar(value=True)

        ttk.Checkbutton(options_frame, text="验证文件", variable=self.validate_var).grid(row=0, column=0, sticky=tk.W, padx=5)
        ttk.Checkbutton(options_frame, text="清理非法字符", variable=self.clean_var).grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Checkbutton(options_frame, text="生成报告", variable=self.create_report_var).grid(row=0, column=2, sticky=tk.W, padx=5)

        # 日志输出区域
        ttk.Label(main_frame, text="处理日志:").grid(row=7, column=0, sticky=tk.W, pady=(10, 5))

        self.log_text = scrolledtext.ScrolledText(main_frame, height=15, wrap=tk.WORD, font=("Consolas", 9))
        self.log_text.grid(row=8, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)

        # 底部按钮
        action_frame = ttk.Frame(main_frame)
        action_frame.grid(row=9, column=0, columnspan=3, pady=15)

        ttk.Button(action_frame, text="开始处理", command=self.process_files, style="Accent.TButton").pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="清空日志", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="退出", command=self.root.quit).pack(side=tk.LEFT, padx=5)

        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=10, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 0))

    def add_files(self):
        """添加文件"""
        files = filedialog.askopenfilenames(
            title="选择 FASTA 文件",
            filetypes=[
                ("FASTA Files", "*.fasta *.fa *.fna *.fas *.fsa"),
                ("All Files", "*.*")
            ]
        )
        for file in files:
            if file not in self.input_files:
                self.input_files.append(file)
                self.file_listbox.insert(tk.END, Path(file).name)
        self.status_var.set(f"已添加 {len(files)} 个文件")

    def add_folder(self):
        """添加文件夹"""
        folder = filedialog.askdirectory(title="选择包含 FASTA 文件的文件夹")
        if folder:
            extensions = {'.fasta', '.fa', '.fna', '.fas', '.fsa'}
            folder_path = Path(folder)
            count = 0
            for ext in extensions:
                for file in folder_path.glob(f"*{ext}"):
                    if str(file) not in self.input_files:
                        self.input_files.append(str(file))
                        self.file_listbox.insert(tk.END, file.name)
                        count += 1
            self.status_var.set(f"从文件夹添加了 {count} 个文件")

    def remove_files(self):
        """移除选中的文件"""
        selection = self.file_listbox.curselection()
        for index in reversed(selection):
            self.file_listbox.delete(index)
            del self.input_files[index]
        self.status_var.set("已移除选中的文件")

    def clear_files(self):
        """清空文件列表"""
        self.file_listbox.delete(0, tk.END)
        self.input_files.clear()
        self.status_var.set("文件列表已清空")

    def browse_output_dir(self):
        """浏览输出目录"""
        folder = filedialog.askdirectory(title="选择输出目录")
        if folder:
            self.output_path_var.set(folder)

    def log(self, message, level="INFO"):
        """添加日志"""
        prefix = {
            "INFO": "ℹ",
            "SUCCESS": "✓",
            "WARNING": "⚠",
            "ERROR": "✗",
            "HEADER": "═"
        }[level]
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(tk.END, f"{prefix} [{timestamp}] {message}\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def clear_log(self):
        """清空日志"""
        self.log_text.delete(1.0, tk.END)

    def validate_fasta(self, filepath):
        """验证 FASTA 文件"""
        invalid_chars = set()
        sequence_count = 0
        total_length = 0

        try:
            for record in SeqIO.parse(filepath, "fasta"):
                sequence_count += 1
                seq_str = str(record.seq).upper()
                total_length += len(seq_str)
                chars_in_seq = set(re.findall(r'[^-AGCT]', seq_str))
                if chars_in_seq:
                    invalid_chars.update(chars_in_seq)

            return {
                'valid': len(invalid_chars) == 0,
                'invalid_chars': sorted(invalid_chars),
                'sequence_count': sequence_count,
                'total_length': total_length
            }
        except Exception as e:
            return {'error': str(e), 'valid': False}

    def clean_fasta(self, input_path, output_path):
        """清理 FASTA 文件"""
        cleaned_records = []
        removed_chars_count = 0
        original_length = 0

        try:
            for record in SeqIO.parse(input_path, "fasta"):
                original_seq = str(record.seq)
                original_length += len(original_seq)
                clean_seq = re.sub(r'[^AGCT]', '', original_seq.upper())
                removed_chars_count += (len(original_seq) - len(clean_seq))

                cleaned_record = SeqRecord(
                    Seq(clean_seq),
                    id=record.id,
                    description=record.description
                )
                cleaned_records.append(cleaned_record)

            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(cleaned_records, output_path, "fasta")

            return {
                'success': True,
                'sequence_count': len(cleaned_records),
                'removed_chars': removed_chars_count,
                'original_length': original_length,
                'cleaned_length': sum(len(r.seq) for r in cleaned_records)
            }
        except Exception as e:
            return {'success': False, 'error': str(e)}

    def process_files(self):
        """处理文件"""
        if not self.input_files:
            messagebox.showwarning("警告", "请先添加 FASTA 文件！")
            return

        output_dir = Path(self.output_path_var.get())
        output_dir.mkdir(parents=True, exist_ok=True)

        report_data = []

        self.log("=" * 60, "HEADER")
        self.log("开始处理文件...", "INFO")
        self.log("=" * 60, "HEADER")

        for i, filepath in enumerate(self.input_files, 1):
            filename = Path(filepath).name
            self.log(f"\n[{i}/{len(self.input_files)}] 处理文件: {filename}", "INFO")

            # 验证
            if self.validate_var.get():
                result = self.validate_fasta(filepath)
                if 'error' in result:
                    self.log(f"  验证失败: {result['error']}", "ERROR")
                elif result['valid']:
                    self.log(f"  验证通过 - {result['sequence_count']} 个序列, {result['total_length']} bp", "SUCCESS")
                else:
                    self.log(f"  验证失败 - 发现非法字符: {', '.join(result['invalid_chars'])}", "WARNING")
                    self.log(f"  序列数: {result['sequence_count']}, 总长度: {result['total_length']} bp", "INFO")
                report_data.append({
                    'filename': filename,
                    'validation': result
                })

            # 清理
            if self.clean_var.get() and 'validation' in report_data[-1] and not report_data[-1]['validation'].get('valid', True):
                output_path = output_dir / f"{Path(filepath).stem}_cleaned{Path(filepath).suffix}"
                self.log(f"  开始清理...", "INFO")
                result = self.clean_fasta(filepath, output_path)

                if result['success']:
                    self.log(f"  清理完成 - 移除 {result['removed_chars']} 个非法字符", "SUCCESS")
                    self.log(f"  输出文件: {output_path.name}", "INFO")
                    self.log(f"  长度变化: {result['original_length']} → {result['cleaned_length']} bp", "INFO")
                    report_data[-1]['cleaning'] = result
                else:
                    self.log(f"  清理失败: {result.get('error', '未知错误')}", "ERROR")

        self.log("\n" + "=" * 60, "HEADER")
        self.log("处理完成！", "SUCCESS")
        self.log(f"共处理 {len(self.input_files)} 个文件", "INFO")
        self.log(f"输出目录: {output_dir}", "INFO")
        self.log("=" * 60, "HEADER")

        # 生成报告
        if self.create_report_var.get():
            self.generate_report(report_data, output_dir)

        self.status_var.set("处理完成")

    def generate_report(self, data, output_dir):
        """生成处理报告"""
        report_path = output_dir / "processing_report.txt"
        try:
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write("=" * 80 + "\n")
                f.write("FASTA 文件处理报告\n")
                f.write("=" * 80 + "\n\n")

                f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"处理文件数: {len(data)}\n\n")

                for item in data:
                    f.write("-" * 80 + "\n")
                    f.write(f"文件名: {item['filename']}\n\n")

                    # 验证结果
                    val = item['validation']
                    f.write("【验证结果】\n")
                    if 'error' in val:
                        f.write(f"  状态: 失败\n")
                        f.write(f"  错误: {val['error']}\n")
                    elif val['valid']:
                        f.write(f"  状态: 通过\n")
                        f.write(f"  序列数: {val['sequence_count']}\n")
                        f.write(f"  总长度: {val['total_length']} bp\n")
                    else:
                        f.write(f"  状态: 失败\n")
                        f.write(f"  非法字符: {', '.join(val['invalid_chars'])}\n")
                        f.write(f"  序列数: {val['sequence_count']}\n")
                        f.write(f"  总长度: {val['total_length']} bp\n")

                    # 清理结果
                    if 'cleaning' in item:
                        clean = item['cleaning']
                        f.write(f"\n【清理结果】\n")
                        if clean['success']:
                            f.write(f"  状态: 成功\n")
                            f.write(f"  处理序列数: {clean['sequence_count']}\n")
                            f.write(f"  移除字符数: {clean['removed_chars']}\n")
                            f.write(f"  原始长度: {clean['original_length']} bp\n")
                            f.write(f"  清理后长度: {clean['cleaned_length']} bp\n")
                        else:
                            f.write(f"  状态: 失败\n")
                            f.write(f"  错误: {clean.get('error', '未知错误')}\n")

                    f.write("\n")

                f.write("=" * 80 + "\n")
                f.write("报告结束\n")
                f.write("=" * 80 + "\n")

            self.log(f"\n报告已保存到: {report_path}", "SUCCESS")
        except Exception as e:
            self.log(f"生成报告失败: {e}", "ERROR")


from datetime import datetime


if __name__ == "__main__":
    root = tk.Tk()
    app = FASTACleanerGUI(root)
    root.mainloop()
