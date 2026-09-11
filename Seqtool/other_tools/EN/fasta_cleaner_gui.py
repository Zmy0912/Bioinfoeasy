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
        self.root.title("FASTA File Validator and Cleaner")
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

        # Title
        title_label = ttk.Label(
            main_frame,
            text="FASTA File Validator and Cleaner",
            font=("Arial", 16, "bold")
        )
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))

        # File selection area
        ttk.Label(main_frame, text="Select FASTA Files:").grid(row=1, column=0, sticky=tk.W, pady=5)

        file_frame = ttk.Frame(main_frame)
        file_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)

        self.file_listbox = tk.Listbox(file_frame, height=6, selectmode=tk.MULTIPLE)
        self.file_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        scrollbar = ttk.Scrollbar(file_frame, orient=tk.VERTICAL, command=self.file_listbox.yview)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.file_listbox.configure(yscrollcommand=scrollbar.set)

        file_frame.columnconfigure(0, weight=1)

        # Button area
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=3, column=0, columnspan=3, pady=10)

        ttk.Button(button_frame, text="Add Files", command=self.add_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Add Folder", command=self.add_folder).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Remove Selected", command=self.remove_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear List", command=self.clear_files).pack(side=tk.LEFT, padx=5)

        # Output directory
        ttk.Label(main_frame, text="Output Directory:").grid(row=4, column=0, sticky=tk.W, pady=5)

        output_frame = ttk.Frame(main_frame)
        output_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        output_frame.columnconfigure(0, weight=1)

        self.output_path_var = tk.StringVar(value="./cleaned_output")
        self.output_entry = ttk.Entry(output_frame, textvariable=self.output_path_var)
        self.output_entry.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        ttk.Button(output_frame, text="Browse", command=self.browse_output_dir).grid(row=0, column=1)

        # Options
        options_frame = ttk.LabelFrame(main_frame, text="Processing Options", padding="10")
        options_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=10)

        self.validate_var = tk.BooleanVar(value=True)
        self.clean_var = tk.BooleanVar(value=True)
        self.create_report_var = tk.BooleanVar(value=True)

        ttk.Checkbutton(options_frame, text="Validate Files", variable=self.validate_var).grid(row=0, column=0, sticky=tk.W, padx=5)
        ttk.Checkbutton(options_frame, text="Clean Invalid Characters", variable=self.clean_var).grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Checkbutton(options_frame, text="Generate Report", variable=self.create_report_var).grid(row=0, column=2, sticky=tk.W, padx=5)

        # Log output area
        ttk.Label(main_frame, text="Processing Log:").grid(row=7, column=0, sticky=tk.W, pady=(10, 5))

        self.log_text = scrolledtext.ScrolledText(main_frame, height=15, wrap=tk.WORD, font=("Consolas", 9))
        self.log_text.grid(row=8, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)

        # Bottom buttons
        action_frame = ttk.Frame(main_frame)
        action_frame.grid(row=9, column=0, columnspan=3, pady=15)

        ttk.Button(action_frame, text="Start Processing", command=self.process_files, style="Accent.TButton").pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Clear Log", command=self.clear_log).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Exit", command=self.root.quit).pack(side=tk.LEFT, padx=5)

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=10, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 0))

    def add_files(self):
        """Add files"""
        files = filedialog.askopenfilenames(
            title="Select FASTA Files",
            filetypes=[
                ("FASTA Files", "*.fasta *.fa *.fna *.fas *.fsa"),
                ("All Files", "*.*")
            ]
        )
        for file in files:
            if file not in self.input_files:
                self.input_files.append(file)
                self.file_listbox.insert(tk.END, Path(file).name)
        self.status_var.set(f"Added {len(files)} file(s)")

    def add_folder(self):
        """Add folder"""
        folder = filedialog.askdirectory(title="Select Folder with FASTA Files")
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
            self.status_var.set(f"Added {count} file(s) from folder")

    def remove_files(self):
        """Remove selected files"""
        selection = self.file_listbox.curselection()
        for index in reversed(selection):
            self.file_listbox.delete(index)
            del self.input_files[index]
        self.status_var.set("Removed selected files")

    def clear_files(self):
        """Clear file list"""
        self.file_listbox.delete(0, tk.END)
        self.input_files.clear()
        self.status_var.set("File list cleared")

    def browse_output_dir(self):
        """Browse output directory"""
        folder = filedialog.askdirectory(title="Select Output Directory")
        if folder:
            self.output_path_var.set(folder)

    def log(self, message, level="INFO"):
        """Add log"""
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
        """Clear log"""
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
        """Process files"""
        if not self.input_files:
            messagebox.showwarning("Warning", "Please add FASTA files first!")
            return

        output_dir = Path(self.output_path_var.get())
        output_dir.mkdir(parents=True, exist_ok=True)

        report_data = []

        self.log("=" * 60, "HEADER")
        self.log("Starting file processing...", "INFO")
        self.log("=" * 60, "HEADER")

        for i, filepath in enumerate(self.input_files, 1):
            filename = Path(filepath).name
            self.log(f"\n[{i}/{len(self.input_files)}] Processing file: {filename}", "INFO")

            # Validation
            if self.validate_var.get():
                result = self.validate_fasta(filepath)
                if 'error' in result:
                    self.log(f"  Validation failed: {result['error']}", "ERROR")
                elif result['valid']:
                    self.log(f"  Validation passed - {result['sequence_count']} sequence(s), {result['total_length']} bp", "SUCCESS")
                else:
                    self.log(f"  Validation failed - Found invalid characters: {', '.join(result['invalid_chars'])}", "WARNING")
                    self.log(f"  Sequence count: {result['sequence_count']}, Total length: {result['total_length']} bp", "INFO")
                report_data.append({
                    'filename': filename,
                    'validation': result
                })

            # Cleaning
            if self.clean_var.get() and 'validation' in report_data[-1] and not report_data[-1]['validation'].get('valid', True):
                output_path = output_dir / f"{Path(filepath).stem}_cleaned{Path(filepath).suffix}"
                self.log(f"  Starting cleanup...", "INFO")
                result = self.clean_fasta(filepath, output_path)

                if result['success']:
                    self.log(f"  Cleanup completed - Removed {result['removed_chars']} invalid character(s)", "SUCCESS")
                    self.log(f"  Output file: {output_path.name}", "INFO")
                    self.log(f"  Length change: {result['original_length']} → {result['cleaned_length']} bp", "INFO")
                    report_data[-1]['cleaning'] = result
                else:
                    self.log(f"  Cleanup failed: {result.get('error', 'Unknown error')}", "ERROR")

        self.log("\n" + "=" * 60, "HEADER")
        self.log("Processing completed!", "SUCCESS")
        self.log(f"Processed {len(self.input_files)} file(s) in total", "INFO")
        self.log(f"Output directory: {output_dir}", "INFO")
        self.log("=" * 60, "HEADER")

        # Generate report
        if self.create_report_var.get():
            self.generate_report(report_data, output_dir)

        self.status_var.set("Processing completed")

    def generate_report(self, data, output_dir):
        """Generate processing report"""
        report_path = output_dir / "processing_report.txt"
        try:
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write("=" * 80 + "\n")
                f.write("FASTA File Processing Report\n")
                f.write("=" * 80 + "\n\n")

                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Total files: {len(data)}\n\n")

                for item in data:
                    f.write("-" * 80 + "\n")
                    f.write(f"Filename: {item['filename']}\n\n")

                    # Validation results
                    val = item['validation']
                    f.write("[Validation Results]\n")
                    if 'error' in val:
                        f.write(f"  Status: Failed\n")
                        f.write(f"  Error: {val['error']}\n")
                    elif val['valid']:
                        f.write(f"  Status: Passed\n")
                        f.write(f"  Sequence count: {val['sequence_count']}\n")
                        f.write(f"  Total length: {val['total_length']} bp\n")
                    else:
                        f.write(f"  Status: Failed\n")
                        f.write(f"  Invalid characters: {', '.join(val['invalid_chars'])}\n")
                        f.write(f"  Sequence count: {val['sequence_count']}\n")
                        f.write(f"  Total length: {val['total_length']} bp\n")

                    # Cleaning results
                    if 'cleaning' in item:
                        clean = item['cleaning']
                        f.write(f"\n[Cleaning Results]\n")
                        if clean['success']:
                            f.write(f"  Status: Success\n")
                            f.write(f"  Sequences processed: {clean['sequence_count']}\n")
                            f.write(f"  Characters removed: {clean['removed_chars']}\n")
                            f.write(f"  Original length: {clean['original_length']} bp\n")
                            f.write(f"  Cleaned length: {clean['cleaned_length']} bp\n")
                        else:
                            f.write(f"  Status: Failed\n")
                            f.write(f"  Error: {clean.get('error', 'Unknown error')}\n")

                    f.write("\n")

                f.write("=" * 80 + "\n")
                f.write("End of Report\n")
                f.write("=" * 80 + "\n")

            self.log(f"\nReport saved to: {report_path}", "SUCCESS")
        except Exception as e:
            self.log(f"Failed to generate report: {e}", "ERROR")


from datetime import datetime


if __name__ == "__main__":
    root = tk.Tk()
    app = FASTACleanerGUI(root)
    root.mainloop()
