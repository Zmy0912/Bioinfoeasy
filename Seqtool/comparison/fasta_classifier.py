#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA Sample Classification Tool
Classify sequences in mixed FASTA files by sample into separate files
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
from pathlib import Path
from collections import defaultdict


class FastaSampleClassifier:
    """FASTA sample sequence classifier"""

    def __init__(self, filepath):
        self.filepath = filepath
        self.sequences = []  # [(header, sequence, sample_name)]
        self.sample_sequences = defaultdict(list)  # {sample_name: [(header, sequence)]}
        self.parse_file()

    def parse_file(self):
        """Parse FASTA file"""
        try:
            with open(self.filepath, 'r', encoding='utf-8') as f:
                current_header = None
                current_seq = []
                sample_name = None

                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_header and current_seq:
                            seq = ''.join(current_seq)
                            self.sequences.append((current_header, seq, sample_name))
                            if sample_name:
                                self.sample_sequences[sample_name].append((current_header, seq))

                        # Parse new header
                        current_header = line[1:]
                        current_seq = []

                        # Extract sample name (from end of header)
                        sample_name = self.extract_sample_name(current_header)
                    else:
                        current_seq.append(line)

                # Save last sequence
                if current_header and current_seq:
                    seq = ''.join(current_seq)
                    self.sequences.append((current_header, seq, sample_name))
                    if sample_name:
                        self.sample_sequences[sample_name].append((current_header, seq))

        except Exception as e:
            raise Exception(f"File parsing failed: {e}")

    @staticmethod
    def extract_sample_name(header):
        """Extract sample name from header"""
        # Try multiple extraction methods
        # Method 1: End of header is |sample.fasta format
        if '|' in header:
            parts = header.split('|')
            last_part = parts[-1].strip()
            if last_part.endswith('.fasta') or last_part.endswith('.fa'):
                return last_part.replace('.fasta', '').replace('.fa', '')
            return last_part

        # Method 2: End of header is sample.fasta format (no | separator)
        if header.endswith('.fasta') or header.endswith('.fa'):
            # Search after last | or space
            for sep in ['|', ' ']:
                if sep in header:
                    parts = header.split(sep)
                    for part in reversed(parts):
                        part = part.strip()
                        if part.endswith('.fasta') or part.endswith('.fa'):
                            return part.replace('.fasta', '').replace('.fa', '')
                    break
            else:
                # No separator, extract directly
                return header.rsplit('.', 1)[0]

        # Method 3: Try to extract last part that looks like a sample name
        for sep in ['|', ' ']:
            if sep in header:
                parts = header.split(sep)
                last_part = parts[-1].strip()
                # If last part contains alphanumeric combination, might be sample name
                if any(c.isalpha() and c.isdigit() for c in last_part):
                    return last_part
                return last_part

        return "unknown"

    def get_samples(self):
        """Get all sample names"""
        return sorted(self.sample_sequences.keys())

    def get_sample_info(self, sample_name):
        """Get sample information"""
        if sample_name not in self.sample_sequences:
            return None
        sequences = self.sample_sequences[sample_name]
        total_length = sum(len(seq) for _, seq in sequences)
        return {
            'count': len(sequences),
            'total_length': total_length,
            'sequences': sequences
        }


class FastaClassifierApp:
    """Graphical interface for FASTA sample classification tool"""

    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Sample Classification Tool")
        self.root.geometry("1200x800")

        self.classifier = None
        self.selected_samples = set()

        self.create_widgets()

    def create_widgets(self):
        """Create interface components"""
        # Top toolbar
        toolbar = ttk.Frame(self.root, padding=5)
        toolbar.pack(fill=tk.X)

        ttk.Button(toolbar, text="Select File", command=self.select_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Export All", command=self.export_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Export Selected", command=self.export_selected).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Select All", command=self.select_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Deselect All", command=self.deselect_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=5)

        # File information
        info_frame = ttk.LabelFrame(self.root, text="File Information", padding=5)
        info_frame.pack(fill=tk.X, padx=5, pady=5)

        self.file_info_label = ttk.Label(info_frame, text="No file selected")
        self.file_info_label.pack(anchor=tk.W)

        # Create PanedWindow split view
        paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left: Sample list
        left_frame = ttk.Frame(paned)
        paned.add(left_frame, weight=1)

        ttk.Label(left_frame, text="Sample List").pack(anchor=tk.W, pady=5)

        # Create sample tree with checkboxes
        self.sample_tree = ttk.Treeview(left_frame, columns=('Count', 'Length'), show='tree headings', selectmode='extended')
        self.sample_tree.heading('#0', text='Sample Name')
        self.sample_tree.heading('Count', text='Seq Count')
        self.sample_tree.heading('Length', text='Total Length')
        self.sample_tree.column('#0', width=200)
        self.sample_tree.column('Count', width=80)
        self.sample_tree.column('Length', width=100)
        self.sample_tree.pack(fill=tk.BOTH, expand=True)

        # Add checkbox column
        self.sample_tree.bind('<ButtonRelease-1>', self.on_tree_click)

        # Right: Sequence details
        right_frame = ttk.Frame(paned)
        paned.add(right_frame, weight=2)

        ttk.Label(right_frame, text="Sequence Details").pack(anchor=tk.W, pady=5)

        # Sequence list
        ttk.Label(right_frame, text="Sequences in sample:").pack(anchor=tk.W)
        self.seq_listbox = tk.Listbox(right_frame, height=8)
        self.seq_listbox.pack(fill=tk.X, padx=5, pady=5)
        self.seq_listbox.bind('<<ListboxSelect>>', self.on_sequence_select)

        # Sequence content viewer
        ttk.Label(right_frame, text="Sequence content:").pack(anchor=tk.W, padx=5, pady=(10, 0))
        self.seq_text = scrolledtext.ScrolledText(right_frame, wrap=tk.WORD)
        self.seq_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Bottom status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Please select FASTA file")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(fill=tk.X)

    def on_tree_click(self, event):
        """Handle tree click event, implement checkbox functionality"""
        item = self.sample_tree.identify_row(event.y)
        column = self.sample_tree.identify_column(event.x)

        if item and column == '#1':
            # Toggle selection state
            sample_name = self.sample_tree.item(item, 'text')
            values = self.sample_tree.item(item, 'values')
            if values[0] == '☑':
                self.sample_tree.item(item, values=('☐',) + values[1:])
                self.selected_samples.discard(sample_name)
            else:
                self.sample_tree.item(item, values=('☑',) + values[1:])
                self.selected_samples.add(sample_name)

    def select_file(self):
        """Select FASTA file"""
        filepath = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        )

        if filepath:
            try:
                self.classifier = FastaSampleClassifier(filepath)
                self.display_samples()
                self.file_info_label.config(text=f"File: {os.path.basename(filepath)}")
                self.status_var.set(f"Loaded {len(self.classifier.sequences)} sequences, {len(self.classifier.get_samples())} samples")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {e}")
                self.status_var.set("File loading failed")

    def display_samples(self):
        """Display sample list"""
        # Clear list
        for item in self.sample_tree.get_children():
            self.sample_tree.delete(item)

        # Display samples
        for sample in self.classifier.get_samples():
            info = self.classifier.get_sample_info(sample)
            self.sample_tree.insert('', 'end', text=sample, values=('☐', info['count'], info['total_length']))

        self.selected_samples.clear()

    def on_sequence_select(self, event):
        """Display selected sequence"""
        selection = self.seq_listbox.curselection()
        if selection:
            index = selection[0]
            if hasattr(self, 'current_sample_sequences') and index < len(self.current_sample_sequences):
                header, seq = self.current_sample_sequences[index]
                content = f"> {header}\n\n"
                content += f"Sequence length: {len(seq)} bp\n\n"
                content += seq
                self.seq_text.delete(1.0, tk.END)
                self.seq_text.insert(tk.END, content)

    def select_all(self):
        """Select all samples"""
        for item in self.sample_tree.get_children():
            values = self.sample_tree.item(item, 'values')
            sample_name = self.sample_tree.item(item, 'text')
            if values[0] == '☐':
                self.sample_tree.item(item, values=('☑',) + values[1:])
                self.selected_samples.add(sample_name)

    def deselect_all(self):
        """Deselect all"""
        for item in self.sample_tree.get_children():
            values = self.sample_tree.item(item, 'values')
            sample_name = self.sample_tree.item(item, 'text')
            if values[0] == '☑':
                self.sample_tree.item(item, values=('☐',) + values[1:])
                self.selected_samples.discard(sample_name)

    def export_all(self):
        """Export all samples"""
        if not self.classifier:
            messagebox.showwarning("Warning", "Please select FASTA file first!")
            return

        output_dir = filedialog.askdirectory(title="Select save folder")
        if not output_dir:
            return

        try:
            total_samples = 0
            total_sequences = 0

            for sample in self.classifier.get_samples():
                info = self.classifier.get_sample_info(sample)
                output_file = os.path.join(output_dir, f"{sample}.fasta")

                with open(output_file, 'w', encoding='utf-8') as f:
                    for header, seq in info['sequences']:
                        f.write(f">{header}\n")
                        f.write(f"{seq}\n")

                total_samples += 1
                total_sequences += info['count']

            messagebox.showinfo("Success", f"Exported {total_samples} samples, {total_sequences} sequences!")
            self.status_var.set(f"Export complete: {total_samples} samples, {total_sequences} sequences")

        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")

    def export_selected(self):
        """Export selected samples"""
        if not self.classifier:
            messagebox.showwarning("Warning", "Please select FASTA file first!")
            return

        if not self.selected_samples:
            messagebox.showwarning("Warning", "Please select samples to export!")
            return

        output_dir = filedialog.askdirectory(title="Select save folder")
        if not output_dir:
            return

        try:
            total_sequences = 0
            exported_samples = []

            for sample in self.selected_samples:
                info = self.classifier.get_sample_info(sample)
                output_file = os.path.join(output_dir, f"{sample}.fasta")

                with open(output_file, 'w', encoding='utf-8') as f:
                    for header, seq in info['sequences']:
                        f.write(f">{header}\n")
                        f.write(f"{seq}\n")

                total_sequences += info['count']
                exported_samples.append(sample)

            messagebox.showinfo("Success", f"Exported {len(exported_samples)} samples, {total_sequences} sequences!")
            self.status_var.set(f"Export complete: {len(exported_samples)} samples, {total_sequences} sequences")

        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")

    def clear_all(self):
        """Clear all data"""
        self.classifier = None
        self.selected_samples.clear()

        for item in self.sample_tree.get_children():
            self.sample_tree.delete(item)

        self.seq_listbox.delete(0, tk.END)
        self.seq_text.delete(1.0, tk.END)
        self.file_info_label.config(text="No file selected")
        self.status_var.set("All data cleared")


def main():
    root = tk.Tk()
    app = FastaClassifierApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
