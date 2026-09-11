import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
import re


class BEDConverter:
    def __init__(self, root):
        self.root = root
        self.root.title("BED/GFF3 Coordinate Converter")
        self.root.geometry("900x700")
        
        self.input_file_path = ""
        self.output_data = []
        self.file_type = "bed"  # bed or gff3
        
        self.create_widgets()
    
    def create_widgets(self):
        # Title
        title_label = tk.Label(
            self.root,
            text="BED/GFF3 Coordinate Converter - Normalize gene coordinates to start from 0",
            font=("Arial", 14, "bold")
        )
        title_label.pack(pady=10)
        
        # File type selection
        type_frame = ttk.LabelFrame(self.root, text="File Type")
        type_frame.pack(fill="x", padx=10, pady=5)
        
        self.file_type_var = tk.StringVar(value="bed")
        ttk.Radiobutton(type_frame, text="BED Format", variable=self.file_type_var, 
                        value="bed", command=self.on_file_type_change).pack(side="left", padx=20, pady=5)
        ttk.Radiobutton(type_frame, text="GFF3 Format", variable=self.file_type_var, 
                        value="gff3", command=self.on_file_type_change).pack(side="left", padx=20, pady=5)
        
        # File selection frame
        file_frame = ttk.LabelFrame(self.root, text="File Selection")
        file_frame.pack(fill="x", padx=10, pady=5)
        
        # Input file
        input_frame = ttk.Frame(file_frame)
        input_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(input_frame, text="Input File:").pack(side="left", padx=5)
        self.input_entry = ttk.Entry(input_frame)
        self.input_entry.pack(side="left", fill="x", expand=True, padx=5)
        
        ttk.Button(input_frame, text="Browse", command=self.browse_input).pack(side="left", padx=5)
        
        # Action buttons
        button_frame = ttk.Frame(self.root)
        button_frame.pack(pady=10)
        
        ttk.Button(button_frame, text="Load File", command=self.load_file).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Convert", command=self.convert).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Save File", command=self.save_file).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Clear", command=self.clear_all).pack(side="left", padx=5)
        
        # Statistics info
        self.stats_label = ttk.Label(self.root, text="")
        self.stats_label.pack(pady=5)
        
        # Preview area
        preview_frame = ttk.LabelFrame(self.root, text="Conversion Result Preview")
        preview_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.preview_text = scrolledtext.ScrolledText(preview_frame, wrap=tk.NONE)
        self.preview_text.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Add horizontal and vertical scrollbars
        h_scroll = ttk.Scrollbar(self.preview_text, orient="horizontal", command=self.preview_text.xview)
        h_scroll.pack(side="bottom", fill="x")
        self.preview_text.configure(xscrollcommand=h_scroll.set)
    
    def on_file_type_change(self):
        self.file_type = self.file_type_var.get()
        # Clear current data
        self.input_entry.delete(0, tk.END)
        self.output_data = []
        self.preview_text.delete(1.0, tk.END)
        self.stats_label.config(text="")
    
    def browse_input(self):
        if self.file_type == "bed":
            filetypes = [("BED Files", "*.bed *.BED"), ("Text Files", "*.txt"), ("All Files", "*.*")]
            title = "Select BED File"
        else:
            filetypes = [("GFF3 Files", "*.gff *.gff3"), ("Text Files", "*.txt"), ("All Files", "*.*")]
            title = "Select GFF3 File"
        
        file_path = filedialog.askopenfilename(
            title=title,
            filetypes=filetypes
        )
        if file_path:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, file_path)
            self.input_file_path = file_path
    
    def load_file(self):
        self.input_file_path = self.input_entry.get().strip()
        
        if not self.input_file_path:
            messagebox.showwarning("Warning", "Please select an input file first!")
            return
        
        if not os.path.exists(self.input_file_path):
            messagebox.showerror("Error", "File does not exist!")
            return
        
        try:
            with open(self.input_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(tk.END, content)
            file_type_name = "BED" if self.file_type == "bed" else "GFF3"
            self.stats_label.config(text=f"Original {file_type_name} file loaded")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {str(e)}")
    
    def convert(self):
        if self.file_type == "bed":
            self.convert_bed()
        else:
            self.convert_gff3()
    
    def convert_bed(self):
        self.input_file_path = self.input_entry.get().strip()
        
        if not self.input_file_path:
            messagebox.showwarning("Warning", "Please select an input file first!")
            return
        
        if not os.path.exists(self.input_file_path):
            messagebox.showerror("Error", "File does not exist!")
            return
        
        try:
            # Read file
            with open(self.input_file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Parse data
            gene_data = {}
            comments = []
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    comments.append(line)
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                gene_id = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                feature_type = parts[3]
                phase = parts[4] if len(parts) > 4 else ''
                
                if gene_id not in gene_data:
                    gene_data[gene_id] = {
                        'gene_start': None,
                        'gene_end': None,
                        'features': []
                    }
                
                # Update gene start and end positions
                if gene_data[gene_id]['gene_start'] is None:
                    gene_data[gene_id]['gene_start'] = start
                    gene_data[gene_id]['gene_end'] = end
                else:
                    gene_data[gene_id]['gene_start'] = min(gene_data[gene_id]['gene_start'], start)
                    gene_data[gene_id]['gene_end'] = max(gene_data[gene_id]['gene_end'], end)
                
                # Save feature information
                gene_data[gene_id]['features'].append({
                    'start': start,
                    'end': end,
                    'type': feature_type,
                    'phase': phase,
                    'original_line': parts
                })
            
            # Convert coordinates
            self.output_data = []
            gene_count = 0
            feature_count = 0
            
            for gene_id, data in gene_data.items():
                gene_start = data['gene_start']
                gene_count += 1
                
                # Process features in original order
                for feature in data['features']:
                    new_start = feature['start'] - gene_start
                    new_end = feature['end'] - gene_start
                    
                    new_line = f"{gene_id}\t{new_start}\t{new_end}\t{feature['type']}"
                    if feature['phase']:
                        new_line += f"\t{feature['phase']}"
                    
                    self.output_data.append(new_line)
                    feature_count += 1
            
            # Show preview
            self.preview_text.delete(1.0, tk.END)
            
            # Add comment lines
            for comment in comments:
                self.preview_text.insert(tk.END, comment + "\n")
            
            # Add converted data
            for line in self.output_data:
                self.preview_text.insert(tk.END, line + "\n")
            
            self.stats_label.config(text=f"Conversion complete: {gene_count} genes, {feature_count} features")
            messagebox.showinfo("Success", f"Conversion complete!\nGenes: {gene_count}\nFeatures: {feature_count}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Conversion failed: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def convert_gff3(self):
        self.input_file_path = self.input_entry.get().strip()
        
        if not self.input_file_path:
            messagebox.showwarning("Warning", "Please select an input file first!")
            return
        
        if not os.path.exists(self.input_file_path):
            messagebox.showerror("Error", "File does not exist!")
            return
        
        try:
            # Read file
            with open(self.input_file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Parse data
            gene_data = {}
            gene_id_map = {}  # Store gene rows for easy lookup
            transcript_map = {}  # Store transcript-gene mapping
            
            for line_idx, line in enumerate(lines):
                line = line.rstrip('\n')
                
                # Skip empty lines and comment lines (but save comment lines)
                if not line.strip():
                    continue
                
                if line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 8:
                    continue
                
                seqid = parts[0]
                source = parts[1]
                feature_type = parts[2]
                start = int(parts[3]) - 1  # GFF3 is 1-based, convert to 0-based
                end = int(parts[4])
                score = parts[5]
                strand = parts[6]
                phase = parts[7]
                attributes = parts[8] if len(parts) > 8 else ''
                
                # Parse attributes
                attrs_dict = self.parse_gff3_attributes(attributes)
                
                gene_id = attrs_dict.get('gene_id') or attrs_dict.get('gene', '')
                transcript_id = attrs_dict.get('transcript_id') or attrs_dict.get('Parent', '')
                
                # Store gene rows
                if feature_type == 'gene' and gene_id:
                    if gene_id not in gene_id_map:
                        gene_id_map[gene_id] = []
                    gene_id_map[gene_id].append({
                        'line_idx': line_idx,
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'original_line': line
                    })
                
                # Establish transcript-gene mapping
                if feature_type == 'mRNA' or feature_type == 'transcript':
                    if transcript_id and gene_id:
                        transcript_map[transcript_id] = gene_id
                
                # Process feature data
                if gene_id:
                    if gene_id not in gene_data:
                        gene_data[gene_id] = {
                            'gene_start': None,
                            'gene_end': None,
                            'seqid': seqid,
                            'strand': strand,
                            'features': []
                        }
                    
                    # Update gene start and end positions
                    if gene_data[gene_id]['gene_start'] is None:
                        gene_data[gene_id]['gene_start'] = start
                        gene_data[gene_id]['gene_end'] = end
                    else:
                        gene_data[gene_id]['gene_start'] = min(gene_data[gene_id]['gene_start'], start)
                        gene_data[gene_id]['gene_end'] = max(gene_data[gene_id]['gene_end'], end)
                    
                    # Save feature information
                    gene_data[gene_id]['features'].append({
                        'line_idx': line_idx,
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'type': feature_type,
                        'score': score,
                        'strand': strand,
                        'phase': phase,
                        'attributes': attributes,
                        'original_line': line
                    })
            
            # Convert coordinates
            self.output_data = []
            gene_count = 0
            feature_count = 0
            
            # Sort by original line number
            for gene_id, data in gene_data.items():
                gene_start = data['gene_start']
                gene_count += 1
                
                # Process features in original order
                data['features'].sort(key=lambda x: x['line_idx'])
                
                for feature in data['features']:
                    new_start = feature['start'] - gene_start
                    new_end = feature['end'] - gene_start
                    
                    # Convert to 1-based (GFF3 standard)
                    new_start_gff = new_start + 1
                    new_end_gff = new_end
                    
                    # Build new line
                    new_line = "\t".join([
                        feature['seqid'],
                        feature['score'],  # source column, use original score value or '.'
                        feature['type'],
                        str(new_start_gff),
                        str(new_end_gff),
                        feature['score'],
                        feature['strand'],
                        feature['phase'],
                        feature['attributes']
                    ])
                    
                    self.output_data.append({
                        'line_idx': feature['line_idx'],
                        'content': new_line
                    })
                    feature_count += 1
            
            # Sort output by original line number
            self.output_data.sort(key=lambda x: x['line_idx'])
            
            # Rebuild complete file content (including comment lines and empty lines)
            converted_lines = []
            output_dict = {item['line_idx']: item['content'] for item in self.output_data}
            
            for line_idx, line in enumerate(lines):
                line_stripped = line.strip()
                
                # Keep comment lines
                if line_stripped.startswith('#'):
                    converted_lines.append(line.rstrip('\n'))
                # Keep empty lines
                elif not line_stripped:
                    converted_lines.append('')
                # Convert data lines
                elif line_idx in output_dict:
                    converted_lines.append(output_dict[line_idx])
            
            self.output_content = converted_lines
            
            # Show preview
            self.preview_text.delete(1.0, tk.END)
            for line in converted_lines:
                self.preview_text.insert(tk.END, line + "\n")
            
            self.stats_label.config(text=f"GFF3 conversion complete: {gene_count} genes, {feature_count} features")
            messagebox.showinfo("Success", f"GFF3 conversion complete!\nGenes: {gene_count}\nFeatures: {feature_count}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Conversion failed: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def parse_gff3_attributes(self, attributes_str):
        """Parse GFF3 attributes field"""
        attrs_dict = {}
        if not attributes_str or attributes_str == '.':
            return attrs_dict
        
        parts = attributes_str.split(';')
        for part in parts:
            part = part.strip()
            if '=' in part:
                key, value = part.split('=', 1)
                attrs_dict[key] = value
        
        return attrs_dict
    
    def save_file(self):
        if self.file_type == "bed":
            if not self.output_data:
                messagebox.showwarning("Warning", "No data to save. Please perform conversion first!")
                return
        else:
            if not hasattr(self, 'output_content') or not self.output_content:
                messagebox.showwarning("Warning", "No data to save. Please perform conversion first!")
                return
        
        # Set default extension based on file type
        if self.file_type == "bed":
            default_extension = ".bed"
            initialfile = "converted.bed"
            filetypes = [("BED Files", "*.bed"), ("Text Files", "*.txt"), ("All Files", "*.*")]
            title = "Save Converted BED File"
        else:
            default_extension = ".gff3"
            initialfile = "converted.gff3"
            filetypes = [("GFF3 Files", "*.gff3 *.gff"), ("Text Files", "*.txt"), ("All Files", "*.*")]
            title = "Save Converted GFF3 File"
        
        file_path = filedialog.asksaveasfilename(
            title=title,
            defaultextension=default_extension,
            initialfile=initialfile,
            filetypes=filetypes
        )
        
        if not file_path:
            return
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                if self.file_type == "bed":
                    # BED format save
                    with open(self.input_file_path, 'r', encoding='utf-8') as orig_f:
                        original_lines = orig_f.readlines()
                    
                    # Extract comment lines
                    comments = []
                    for line in original_lines:
                        if line.strip().startswith('#'):
                            comments.append(line)
                    
                    # Write comment lines
                    for comment in comments:
                        f.write(comment)
                    
                    # Write converted data
                    for line in self.output_data:
                        f.write(line + "\n")
                else:
                    # GFF3 format save
                    for line in self.output_content:
                        if line == '':
                            f.write("\n")
                        else:
                            f.write(line + "\n")
            
            messagebox.showinfo("Success", f"File saved to:\n{file_path}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {str(e)}")
    
    def clear_all(self):
        self.input_entry.delete(0, tk.END)
        self.input_file_path = ""
        self.output_data = []
        if hasattr(self, 'output_content'):
            delattr(self, 'output_content')
        self.preview_text.delete(1.0, tk.END)
        self.stats_label.config(text="")


def main():
    root = tk.Tk()
    app = BEDConverter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
