#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GFF3 to BED Converter
A graphical tool for converting GFF3 format to BED format
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os


class GFF3ToBEDConverter:
    def __init__(self, root):
        self.root = root
        self.root.title("GFF3 to BED Converter")
        self.root.geometry("800x750")
        self.gff3_file = ""
        self.bed_file = ""
        self.parsed_features = []
        self.gene_mapping = {}
        self.bed_preview_content = []
        
        self.setup_ui()
    
    def setup_ui(self):
        # Title
        title_label = tk.Label(
            self.root,
            text="GFF3 to BED Converter",
            font=("Arial", 16, "bold"),
            fg="#2c3e50",
            pady=10
        )
        title_label.pack()
        
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # File selection area
        file_frame = ttk.LabelFrame(main_frame, text="File Selection", padding="10")
        file_frame.pack(fill=tk.X, pady=5)
        
        # GFF3 file selection
        gff3_label = ttk.Label(file_frame, text="GFF3 File:", font=("Arial", 10))
        gff3_label.grid(row=0, column=0, sticky=tk.W, pady=5)
        
        gff3_entry = ttk.Entry(file_frame, width=45)
        gff3_entry.grid(row=0, column=1, padx=5, pady=5)
        self.gff3_entry = gff3_entry
        
        gff3_btn = ttk.Button(
            file_frame,
            text="Browse",
            command=self.select_gff3_file
        )
        gff3_btn.grid(row=0, column=2, padx=5, pady=5)
        
        # BED file save
        bed_label = ttk.Label(file_frame, text="BED File:", font=("Arial", 10))
        bed_label.grid(row=1, column=0, sticky=tk.W, pady=5)
        
        bed_entry = ttk.Entry(file_frame, width=45)
        bed_entry.grid(row=1, column=1, padx=5, pady=5)
        self.bed_entry = bed_entry
        
        bed_btn = ttk.Button(
            file_frame,
            text="Browse",
            command=self.select_bed_file
        )
        bed_btn.grid(row=1, column=2, padx=5, pady=5)
        
        # Options area
        option_frame = ttk.LabelFrame(main_frame, text="Conversion Options", padding="10")
        option_frame.pack(fill=tk.X, pady=5)
        
        # Feature type selection
        self.feature_types = {
            "gene": tk.BooleanVar(value=True),
            "mRNA": tk.BooleanVar(value=True),
            "exon": tk.BooleanVar(value=True),
            "CDS": tk.BooleanVar(value=True),
        }
        
        type_frame = ttk.Frame(option_frame)
        type_frame.pack(fill=tk.X, pady=5)
        
        row, col = 0, 0
        for feature, var in self.feature_types.items():
            checkbox = ttk.Checkbutton(type_frame, text=f"Convert {feature} features", variable=var)
            checkbox.grid(row=row, column=col, padx=10, pady=2, sticky=tk.W)
            col += 1
            if col > 1:
                col = 0
                row += 1
        
        # Gene name options
        name_frame = ttk.Frame(option_frame)
        name_frame.pack(fill=tk.X, pady=5)
        
        self.use_custom_name = tk.BooleanVar(value=False)
        custom_name_check = ttk.Checkbutton(
            name_frame,
            text="Use custom gene names (do not use names from GFF3)",
            variable=self.use_custom_name
        )
        custom_name_check.grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        # Conversion buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=10)
        
        preview_btn = ttk.Button(
            button_frame,
            text="Preview Output",
            command=self.preview_output,
            width=15
        )
        preview_btn.pack(side=tk.LEFT, padx=5)
        
        convert_btn = ttk.Button(
            button_frame,
            text="Start Conversion",
            command=self.convert_gff3_to_bed,
            width=15
        )
        convert_btn.pack(side=tk.LEFT, padx=5)
        
        clear_btn = ttk.Button(
            button_frame,
            text="Clear",
            command=self.clear_all,
            width=10
        )
        clear_btn.pack(side=tk.LEFT, padx=5)
        
        # Log area
        log_frame = ttk.LabelFrame(main_frame, text="Conversion Log", padding="10")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        self.log_text = scrolledtext.ScrolledText(
            log_frame,
            height=8,
            width=70,
            wrap=tk.WORD,
            font=("Consolas", 9)
        )
        self.log_text.pack(fill=tk.BOTH, expand=True)
        
        # Configure log colors
        self.log_text.tag_config("info", foreground="blue")
        self.log_text.tag_config("success", foreground="green")
        self.log_text.tag_config("error", foreground="red")
        self.log_text.tag_config("warning", foreground="orange")
    
    def select_gff3_file(self):
        file_path = filedialog.askopenfilename(
            title="Select GFF3 File",
            filetypes=[("GFF3 files", "*.gff3 *.gff"), ("All files", "*.*")]
        )
        if file_path:
            self.gff3_file = file_path
            self.gff3_entry.delete(0, tk.END)
            self.gff3_entry.insert(0, file_path)
            # Automatically set BED file path
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            bed_path = os.path.join(os.path.dirname(file_path), f"{base_name}.bed")
            self.bed_entry.delete(0, tk.END)
            self.bed_entry.insert(0, bed_path)
            self.bed_file = bed_path
            self.log(f"Selected GFF3 file: {file_path}", "info")
            # Automatically parse file to get gene information
            self.parse_gff3_file()
    
    def select_bed_file(self):
        file_path = filedialog.asksaveasfilename(
            title="Save BED File",
            defaultextension=".bed",
            filetypes=[("BED files", "*.bed"), ("All files", "*.*")]
        )
        if file_path:
            self.bed_file = file_path
            self.bed_entry.delete(0, tk.END)
            self.bed_entry.insert(0, file_path)
            self.log(f"Set BED output file: {file_path}", "info")
    
    def parse_gff3_file(self):
        """Parse GFF3 file"""
        self.parsed_features = []
        self.gene_mapping = {}
        
        try:
            with open(self.gff3_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    # Skip comment lines and empty lines
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        seqid = parts[0]
                        source = parts[1]
                        feature_type = parts[2]
                        start = int(parts[3])
                        end = int(parts[4])
                        score = parts[5]
                        strand = parts[6]
                        phase = parts[7]
                        attributes = parts[8]
                        
                        feature = {
                            'line_num': line_num,
                            'seqid': seqid,
                            'source': source,
                            'type': feature_type,
                            'start': start,
                            'end': end,
                            'score': score,
                            'strand': strand,
                            'phase': phase,
                            'attributes': attributes
                        }
                        
                        self.parsed_features.append(feature)
                        
                        # Extract gene ID
                        gene_id = self.extract_gene_id(feature)
                        if gene_id and feature_type == 'gene':
                            self.gene_mapping[gene_id] = gene_id
            
            self.log(f"Successfully parsed {len(self.parsed_features)} features", "success")
            
            # Display gene statistics
            if self.gene_mapping:
                self.log(f"Detected {len(self.gene_mapping)} genes:", "info")
                gene_list = sorted(self.gene_mapping.keys())
                for i, gene in enumerate(gene_list):
                    if i < 5:  # Only show first 5
                        self.log(f"  - {gene}", "info")
                    elif i == 5:
                        self.log(f"  ... and {len(gene_list) - 5} more genes", "info")
                        break
            
        except Exception as e:
            self.log(f"Error parsing GFF3 file: {str(e)}", "error")
            self.parsed_features = []
            self.gene_mapping = {}
    
    def extract_gene_id(self, feature):
        """Extract gene ID from feature"""
        attributes = feature['attributes']
        gene_id = ""
        
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('ID='):
                gene_id = attr.split('=', 1)[1]
                break
            elif attr.startswith('Name='):
                gene_id = attr.split('=', 1)[1]
                break
            elif attr.startswith('gene='):
                gene_id = attr.split('=', 1)[1]
                break
        
        return gene_id
    
    def find_gene_for_feature(self, feature, all_features):
        """Find corresponding gene for non-gene features"""
        gene_id = ""
        gene_name = ""
        
        # First try to find from Parent attribute
        attributes = feature['attributes']
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('Parent=gene-'):
                gene_name = attr.split('=', 1)[1].replace('gene-', '')
                break
            elif attr.startswith('Parent='):
                parent_id = attr.split('=', 1)[1]
                # Check if this parent ID is a gene
                for f in all_features:
                    if f['type'] == 'gene':
                        f_gene_id = self.extract_gene_id(f)
                        if f_gene_id == parent_id or parent_id.startswith('gene-' + f_gene_id):
                            gene_name = f_gene_id
                            break
                if gene_name:
                    break
        
        # If not found, try to find from gene attribute
        if not gene_name:
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene='):
                    gene_name = attr.split('=', 1)[1]
                    break
        
        # If still not found, search in all genes
        if not gene_name:
            feature_start = feature['start']
            feature_end = feature['end']
            feature_seqid = feature['seqid']
            
            for f in all_features:
                if f['type'] == 'gene' and f['seqid'] == feature_seqid:
                    # Check if feature is within gene range
                    if (f['start'] <= feature_start <= f['end'] or
                        f['start'] <= feature_end <= f['end'] or
                        (feature_start <= f['start'] and feature_end >= f['end'])):
                        gene_name = self.extract_gene_id(f)
                        break
        
        return gene_name if gene_name else feature['seqid']
    
    def generate_bed_entries(self):
        """Generate BED entries"""
        bed_entries = []
        
        selected_types = [ft for ft, var in self.feature_types.items() if var.get()]
        
        if not selected_types:
            self.log("Error: Please select at least one feature type!", "error")
            return None
        
        # Determine gene name for each feature
        for feature in self.parsed_features:
            if feature['type'] not in selected_types:
                continue
            
            # Determine gene name
            if feature['type'] == 'gene':
                gene_name = self.extract_gene_id(feature)
            else:
                gene_name = self.find_gene_for_feature(feature, self.parsed_features)
            
            # Start and end positions
            start = feature['start']
            end = feature['end']
            
            # Feature type
            feature_type = feature['type']
            
            # phase (optional column)
            phase = feature['phase'] if feature['phase'] != '.' else ''
            
            # Build BED entry
            if phase:
                bed_entry = f"{gene_name}\t{start}\t{end}\t{feature_type}\t{phase}"
            else:
                bed_entry = f"{gene_name}\t{start}\t{end}\t{feature_type}"
            
            bed_entries.append(bed_entry)
        
        return bed_entries
    
    def show_preview_window(self):
        """Show preview window"""
        preview_window = tk.Toplevel(self.root)
        preview_window.title("BED Output Preview")
        preview_window.geometry("900x600")
        
        # Toolbar
        toolbar = ttk.Frame(preview_window, padding="5")
        toolbar.pack(fill=tk.X)
        
        save_btn = ttk.Button(
            toolbar,
            text="Save to File",
            command=lambda: self.save_from_preview(preview_window)
        )
        save_btn.pack(side=tk.LEFT, padx=5)
        
        close_btn = ttk.Button(
            toolbar,
            text="Close",
            command=preview_window.destroy
        )
        close_btn.pack(side=tk.LEFT, padx=5)
        
        # Info label
        info_label = ttk.Label(
            toolbar,
            text=f"Total {len(self.bed_preview_content)} lines (showing first 500 lines)",
            font=("Arial", 9)
        )
        info_label.pack(side=tk.LEFT, padx=20)
        
        # Preview text box
        text_frame = ttk.Frame(preview_window)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        preview_text = scrolledtext.ScrolledText(
            text_frame,
            wrap=tk.NONE,
            font=("Consolas", 10)
        )
        preview_text.pack(fill=tk.BOTH, expand=True)
        
        # Display content
        preview_text.insert(tk.END, "# BED format file converted from GFF3\n")
        preview_text.insert(tk.END, "# Columns: geneID/transcriptID  start  end  featureType  phase(optional)\n")
        preview_text.insert(tk.END, "\n")
        
        display_lines = self.bed_preview_content[:500]  # Show max 500 lines
        for line in display_lines:
            preview_text.insert(tk.END, line + "\n")
        
        if len(self.bed_preview_content) > 500:
            preview_text.insert(tk.END, f"\n... and {len(self.bed_preview_content) - 500} more lines not shown\n")
        
        preview_text.config(state=tk.DISABLED)
    
    def save_from_preview(self, preview_window):
        """Save from preview window"""
        if not self.bed_file:
            messagebox.showwarning("Warning", "Please set BED output file first!")
            return
        
        if self.save_bed_file(self.bed_preview_content, self.bed_file):
            messagebox.showinfo("Success", f"Conversion successful!\nBED file saved to:\n{self.bed_file}")
            preview_window.destroy()
    
    def preview_output(self):
        """Preview output content"""
        if not self.gff3_file:
            messagebox.showwarning("Warning", "Please select GFF3 file first!")
            return
        
        if not self.parsed_features:
            self.log("Parsing GFF3 file...", "info")
            self.parse_gff3_file()
        
        if not self.parsed_features:
            messagebox.showerror("Error", "Unable to parse GFF3 file!")
            return
        
        self.log("Generating BED preview...", "info")
        self.bed_preview_content = self.generate_bed_entries()
        
        if self.bed_preview_content is None:
            return
        
        self.log(f"Generated {len(self.bed_preview_content)} BED entries", "success")
        
        # Show preview window
        self.show_preview_window()
    
    def save_bed_file(self, bed_entries, file_path):
        """Save BED file"""
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                # Write BED file header
                f.write("# BED format file converted from GFF3\n")
                f.write("# Columns: geneID/transcriptID  start  end  featureType  phase(optional)\n")
                f.write("\n")
                
                for entry in bed_entries:
                    f.write(entry + '\n')
            
            return True
        except Exception as e:
            self.log(f"Error saving BED file: {str(e)}", "error")
            return False
    
    def convert_gff3_to_bed(self):
        """Execute conversion"""
        if not self.gff3_file:
            messagebox.showwarning("Warning", "Please select GFF3 file first!")
            return
        
        if not self.bed_file:
            messagebox.showwarning("Warning", "Please set BED output file!")
            return
        
        if not os.path.exists(self.gff3_file):
            messagebox.showerror("Error", f"GFF3 file does not exist: {self.gff3_file}")
            return
        
        if not self.parsed_features:
            self.log("Parsing GFF3 file...", "info")
            self.parse_gff3_file()
        
        if not self.parsed_features:
            messagebox.showerror("Error", "Unable to parse GFF3 file!")
            return
        
        self.log("=" * 60, "info")
        self.log("Starting conversion...", "info")
        self.log(f"GFF3 file: {self.gff3_file}", "info")
        self.log(f"BED file: {self.bed_file}", "info")
        
        # Feature type statistics
        type_count = {}
        for f in self.parsed_features:
            ftype = f['type']
            type_count[ftype] = type_count.get(ftype, 0) + 1
        
        self.log("Feature statistics:", "info")
        for ftype, count in sorted(type_count.items()):
            self.log(f"  {ftype}: {count}", "info")
        
        # Generate BED entries
        self.log("Generating BED entries...", "info")
        bed_entries = self.generate_bed_entries()
        
        if bed_entries is None:
            return
        
        self.log(f"Generated {len(bed_entries)} BED entries", "success")
        
        # Display gene name statistics
        gene_names = set()
        for entry in bed_entries:
            gene_name = entry.split('\t')[0]
            gene_names.add(gene_name)
        
        self.log(f"Total {len(gene_names)} unique gene names", "success")
        
        # Save BED file
        self.log("Saving BED file...", "info")
        if self.save_bed_file(bed_entries, self.bed_file):
            self.log("Conversion successful! BED file saved to:", "success")
            self.log(f"  {self.bed_file}", "success")
            messagebox.showinfo("Success", f"Conversion successful!\nBED file saved to:\n{self.bed_file}")
        else:
            messagebox.showerror("Error", "Failed to save BED file!")
        
        self.log("=" * 60, "info")
    
    def clear_all(self):
        """Clear all inputs"""
        self.gff3_file = ""
        self.bed_file = ""
        self.parsed_features = []
        self.gene_mapping = {}
        self.bed_preview_content = []
        self.gff3_entry.delete(0, tk.END)
        self.bed_entry.delete(0, tk.END)
        self.log_text.delete(1.0, tk.END)
        self.log("Cleared all content", "info")
    
    def log(self, message, level="info"):
        """Add log"""
        self.log_text.insert(tk.END, message + "\n", level)
        self.log_text.see(tk.END)
        self.root.update()


def main():
    root = tk.Tk()
    app = GFF3ToBEDConverter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
