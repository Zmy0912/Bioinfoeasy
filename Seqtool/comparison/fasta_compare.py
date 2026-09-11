#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA Sequence Comparison Tool
For comparing sequences across multiple FASTA files to identify common and unique sequences
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
from pathlib import Path
from collections import defaultdict


class FastaParser:
    """FASTA file parser"""

    @staticmethod
    def parse_file(filepath):
        """Parse a single FASTA file, return sequence dictionary {header: sequence}"""
        sequences = {}
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                current_header = None
                current_seq = []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_header:
                            sequences[current_header] = ''.join(current_seq)
                        current_header = line[1:]  # Remove '>'
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
        except Exception as e:
            print(f"Error parsing {filepath}: {e}")
        return sequences

    @staticmethod
    def get_header_key(header):
        """Extract unique identifier from complete header (remove position info, etc.)"""
        # Extract gene name (e.g., extract "psbA" from "psbA_[76:1137](-)")
        if '_join{' in header:
            return header.split('_join{')[0]
        elif '[_' in header:
            return header.split('[_')[0]
        else:
            return header.split('[')[0] if '[' in header else header


class SequenceComparator:
    """Sequence comparator"""

    def __init__(self, fasta_files):
        self.fasta_files = fasta_files
        self.file_data = {}  # {filename: {header_key: {full_header: sequence}}}
        self.load_files()

    def load_files(self):
        """Load all FASTA files"""
        for filepath in self.fasta_files:
            filename = os.path.basename(filepath)
            sequences = FastaParser.parse_file(filepath)
            # Organize by header_key
            file_dict = defaultdict(dict)
            for header, seq in sequences.items():
                key = FastaParser.get_header_key(header)
                file_dict[key][header] = seq
            self.file_data[filename] = dict(file_dict)

    def compare(self):
        """Compare sequences across all files"""
        # Collect sequence keys from all files
        all_keys = set()
        file_keys = {}  # {filename: set(keys)}
        for filename, data in self.file_data.items():
            keys = set(data.keys())
            file_keys[filename] = keys
            all_keys.update(keys)

        # Count how many files each key appears in
        key_count = defaultdict(int)
        key_files = defaultdict(list)  # {key: [filenames]}
        for key in all_keys:
            for filename, keys in file_keys.items():
                if key in keys:
                    key_count[key] += 1
                    key_files[key].append(filename)

        # Classify
        total_files = len(self.fasta_files)
        common_sequences = {k: key_files[k] for k, v in key_count.items() if v == total_files}
        unique_sequences = {}
        for filename, keys in file_keys.items():
            unique_in_file = keys - set(k for k, v in key_count.items() if v > 1)
            unique_sequences[filename] = unique_in_file

        return {
            'all_keys': all_keys,
            'common': common_sequences,
            'unique': unique_sequences,
            'file_keys': file_keys
        }


class FastaCompareApp:
    """Graphical interface for FASTA sequence comparison tool"""

    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Sequence Comparison Tool")
        self.root.geometry("1200x800")

        self.fasta_files = []
        self.comparator = None
        self.comparison_result = None

        self.create_widgets()

    def create_widgets(self):
        """Create interface components"""
        # Top toolbar
        toolbar = ttk.Frame(self.root, padding=5)
        toolbar.pack(fill=tk.X)

        ttk.Button(toolbar, text="Select Folder", command=self.select_folder).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Start Comparison", command=self.start_compare).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Export Common", command=self.export_common).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Export Unique", command=self.export_unique).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=5)

        # File list
        files_frame = ttk.LabelFrame(self.root, text="Loaded Files", padding=5)
        files_frame.pack(fill=tk.X, padx=5, pady=5)

        self.files_listbox = tk.Listbox(files_frame, height=4)
        self.files_listbox.pack(fill=tk.X, expand=True)

        # Create Notebook for multiple tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # 1. Common sequences tab
        self.common_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.common_frame, text="Common")

        self.common_tree = ttk.Treeview(self.common_frame, columns=('Count', 'Files'), show='tree headings')
        self.common_tree.heading('#0', text='Sequence Name')
        self.common_tree.heading('Count', text='Count')
        self.common_tree.heading('Files', text='Files')
        self.common_tree.column('#0', width=300)
        self.common_tree.column('Count', width=100)
        self.common_tree.column('Files', width=500)
        self.common_tree.pack(fill=tk.BOTH, expand=True)
        self.common_tree.bind('<Double-1>', self.on_double_click)

        # 2. Unique sequences tab
        self.unique_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.unique_frame, text="Unique")

        self.unique_tree = ttk.Treeview(self.unique_frame, columns=('File'), show='tree headings')
        self.unique_tree.heading('#0', text='Sequence Name')
        self.unique_tree.heading('File', text='Source File')
        self.unique_tree.column('#0', width=400)
        self.unique_tree.column('File', width=600)
        self.unique_tree.pack(fill=tk.BOTH, expand=True)
        self.unique_tree.bind('<Double-1>', self.on_double_click)

        # 3. Detailed statistics tab
        self.stats_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.stats_frame, text="Statistics")

        self.stats_text = scrolledtext.ScrolledText(self.stats_frame, wrap=tk.WORD)
        self.stats_text.pack(fill=tk.BOTH, expand=True)

        # 4. Sequence viewer tab
        self.view_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.view_frame, text="View")

        view_controls = ttk.Frame(self.view_frame)
        view_controls.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(view_controls, text="Sequence Name:").pack(side=tk.LEFT)
        self.seq_name_var = tk.StringVar()
        self.seq_name_combo = ttk.Combobox(view_controls, textvariable=self.seq_name_var, width=50)
        self.seq_name_combo.pack(side=tk.LEFT, padx=5)

        ttk.Label(view_controls, text="Source File:").pack(side=tk.LEFT)
        self.source_file_var = tk.StringVar()
        self.source_file_combo = ttk.Combobox(view_controls, textvariable=self.source_file_var, width=20)
        self.source_file_combo.pack(side=tk.LEFT, padx=5)

        ttk.Button(view_controls, text="View", command=self.view_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Button(view_controls, text="Copy", command=self.copy_sequence).pack(side=tk.LEFT, padx=5)
        ttk.Button(view_controls, text="Export", command=self.export_current_sequence).pack(side=tk.LEFT, padx=5)

        self.seq_view_text = scrolledtext.ScrolledText(self.view_frame, wrap=tk.WORD)
        self.seq_view_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Please select folder containing FASTA files")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X)

    def select_folder(self):
        """Select folder"""
        folder = filedialog.askdirectory(title="Select folder containing FASTA files")
        if folder:
            self.load_fasta_files(folder)

    def load_fasta_files(self, folder):
        """Load FASTA files from folder"""
        self.fasta_files = []
        folder_path = Path(folder)
        fasta_files = list(folder_path.glob("*.fasta")) + list(folder_path.glob("*.fa"))

        if not fasta_files:
            messagebox.showwarning("Warning", "No FASTA files found in folder!")
            return

        self.fasta_files = [str(f) for f in fasta_files]

        # Update file list
        self.files_listbox.delete(0, tk.END)
        for f in self.fasta_files:
            self.files_listbox.insert(tk.END, os.path.basename(f))

        # Update source file dropdown
        self.source_file_combo['values'] = [os.path.basename(f) for f in self.fasta_files]
        if self.fasta_files:
            self.source_file_combo.current(0)

        self.status_var.set(f"Loaded {len(self.fasta_files)} FASTA files")

    def start_compare(self):
        """Start comparison"""
        if not self.fasta_files:
            messagebox.showwarning("Warning", "Please select folder containing FASTA files first!")
            return

        try:
            self.status_var.set("Comparing sequences...")
            self.root.update()

            self.comparator = SequenceComparator(self.fasta_files)
            self.comparison_result = self.comparator.compare()

            self.display_results()
            self.status_var.set("Comparison complete!")

        except Exception as e:
            messagebox.showerror("Error", f"Error during comparison: {e}")
            self.status_var.set("Comparison failed!")

    def display_results(self):
        """Display comparison results"""
        # Clear trees
        for item in self.common_tree.get_children():
            self.common_tree.delete(item)
        for item in self.unique_tree.get_children():
            self.unique_tree.delete(item)

        # Display common sequences
        for key, files in self.comparison_result['common'].items():
            self.common_tree.insert('', 'end', text=key, values=(len(files), ', '.join(files)))

        # Display unique sequences
        for filename, keys in self.comparison_result['unique'].items():
            for key in keys:
                self.unique_tree.insert('', 'end', text=key, values=(filename,))

        # Display statistics
        self.display_statistics()

        # Update sequence name dropdown
        all_keys = list(self.comparison_result['all_keys'])
        self.seq_name_combo['values'] = sorted(all_keys)

    def display_statistics(self):
        """Display statistics"""
        if not self.comparison_result:
            return

        stats = []
        stats.append("=" * 60)
        stats.append("FASTA Sequence Comparison Statistics")
        stats.append("=" * 60)
        stats.append(f"\nNumber of files compared: {len(self.fasta_files)}")

        for filename in self.fasta_files:
            fname = os.path.basename(filename)
            count = len(self.comparator.file_data[fname])
            stats.append(f"  - {fname}: {count} sequences")

        stats.append(f"\nTotal distinct sequences: {len(self.comparison_result['all_keys'])}")
        stats.append(f"Common sequences (in all files): {len(self.comparison_result['common'])}")

        stats.append("\nUnique sequences:")
        for filename, keys in self.comparison_result['unique'].items():
            stats.append(f"  - {filename}: {len(keys)} unique sequences")

        # Count sequence occurrence frequency
        stats.append("\nSequence frequency distribution:")
        key_count = {}
        for key, files in self.comparison_result['common'].items():
            key_count[key] = len(files)
        # Add non-common sequences
        for filename, keys in self.comparison_result['unique'].items():
            for key in keys:
                key_count[key] = 1

        count_distribution = defaultdict(int)
        for key, count in key_count.items():
            count_distribution[count] += 1

        for count in sorted(count_distribution.keys(), reverse=True):
            stats.append(f"  - Appearing in {count} file(s): {count_distribution[count]} sequences")

        self.stats_text.delete(1.0, tk.END)
        self.stats_text.insert(tk.END, '\n'.join(stats))

    def on_double_click(self, event):
        """Double-click to display sequence"""
        tree = event.widget
        item = tree.selection()[0] if tree.selection() else None
        if item:
            key = tree.item(item, 'text')
            self.seq_name_var.set(key)
            self.notebook.select(self.view_frame)
            self.view_sequence()

    def view_sequence(self):
        """View sequence"""
        key = self.seq_name_var.get().strip()
        source_file = self.source_file_var.get().strip()

        if not key:
            messagebox.showwarning("Warning", "Please enter sequence name!")
            return

        if not self.comparator:
            messagebox.showwarning("Warning", "Please perform comparison first!")
            return

        # Find sequence
        found = False
        result = []

        if source_file:
            # Search sequence in specified file
            files_to_search = [source_file]
        else:
            # Search sequence in all files
            files_to_search = self.comparator.file_data.keys()

        for fname in files_to_search:
            if key in self.comparator.file_data[fname]:
                found = True
                for full_header, seq in self.comparator.file_data[fname][key].items():
                    result.append(f"File: {fname}")
                    result.append(f"Sequence ID: {full_header}")
                    result.append(f"Sequence Length: {len(seq)} bp")
                    result.append("\nSequence:")
                    result.append(seq)
                    result.append("\n" + "=" * 60 + "\n")

        if not found:
            self.seq_view_text.delete(1.0, tk.END)
            self.seq_view_text.insert(tk.END, f"Sequence not found: {key}")
        else:
            self.seq_view_text.delete(1.0, tk.END)
            self.seq_view_text.insert(tk.END, '\n'.join(result))

    def copy_sequence(self):
        """Copy sequence to clipboard"""
        content = self.seq_view_text.get(1.0, tk.END)
        self.root.clipboard_clear()
        self.root.clipboard_append(content)
        messagebox.showinfo("Success", "Copied to clipboard!")

    def export_common(self):
        """Export common sequences"""
        if not self.comparison_result:
            messagebox.showwarning("Warning", "Please perform comparison first!")
            return

        output_file = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")],
            title="Save common sequences"
        )

        if not output_file:
            return

        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                for key in self.comparison_result['common'].keys():
                    for filename in self.fasta_files:
                        fname = os.path.basename(filename)
                        if key in self.comparator.file_data[fname]:
                            for full_header, seq in self.comparator.file_data[fname][key].items():
                                f.write(f">{full_header}|{fname}\n")
                                f.write(f"{seq}\n")
                                break  # Write only one

            messagebox.showinfo("Success", f"Exported {len(self.comparison_result['common'])} common sequences!")
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")

    def export_unique(self):
        """Export unique sequences"""
        if not self.comparison_result:
            messagebox.showwarning("Warning", "Please perform comparison first!")
            return

        output_dir = filedialog.askdirectory(title="Select folder to save unique sequences")
        if not output_dir:
            return

        try:
            total_exported = 0
            for filename, keys in self.comparison_result['unique'].items():
                if keys:  # If there are unique sequences
                    output_file = os.path.join(output_dir, f"{filename}_unique.fasta")
                    with open(output_file, 'w', encoding='utf-8') as f:
                        for key in keys:
                            for full_header, seq in self.comparator.file_data[filename][key].items():
                                f.write(f">{full_header}|{filename}_unique\n")
                                f.write(f"{seq}\n")
                                break
                    total_exported += len(keys)

            if total_exported > 0:
                messagebox.showinfo("Success", f"Exported {total_exported} unique sequences to {output_dir}!")
            else:
                messagebox.showinfo("Info", "No unique sequences to export!")
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")

    def export_current_sequence(self):
        """Export currently viewed sequence"""
        key = self.seq_name_var.get().strip()
        if not key:
            messagebox.showwarning("Warning", "Please view a sequence first!")
            return

        output_file = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")],
            title="Save sequence"
        )

        if not output_file:
            return

        try:
            content = self.seq_view_text.get(1.0, tk.END)
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(content)
            messagebox.showinfo("Success", f"Saved to {output_file}!")
        except Exception as e:
            messagebox.showerror("Error", f"Save failed: {e}")

    def clear_all(self):
        """Clear all data"""
        self.fasta_files = []
        self.comparator = None
        self.comparison_result = None
        self.files_listbox.delete(0, tk.END)
        self.common_tree.delete(*self.common_tree.get_children())
        self.unique_tree.delete(*self.unique_tree.get_children())
        self.stats_text.delete(1.0, tk.END)
        self.seq_view_text.delete(1.0, tk.END)
        self.seq_name_var.set("")
        self.source_file_var.set("")
        self.status_var.set("All data cleared")


def main():
    root = tk.Tk()
    app = FastaCompareApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
