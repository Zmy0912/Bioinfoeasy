import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import re


class PiGeneAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Gene Nucleotide Polymorphism (Pi) Analyzer")
        self.root.geometry("800x600")
        
        self.pi_file_path = None
        self.gff_file_path = None
        self.result_df = None
        
        self.create_widgets()
    
    def create_widgets(self):
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        
        title_label = ttk.Label(
            main_frame, 
            text="Gene Nucleotide Polymorphism (Pi) Analyzer", 
            font=('Arial', 16, 'bold')
        )
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        ttk.Label(main_frame, text="Pi File:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.pi_entry = ttk.Entry(main_frame, width=60)
        self.pi_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        ttk.Button(main_frame, text="Browse", command=self.browse_pi_file).grid(row=1, column=2, padx=5, pady=5)
        
        ttk.Label(main_frame, text="GFF3 File:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.gff_entry = ttk.Entry(main_frame, width=60)
        self.gff_entry.grid(row=2, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        ttk.Button(main_frame, text="Browse", command=self.browse_gff_file).grid(row=2, column=2, padx=5, pady=5)
        
        ttk.Separator(main_frame, orient='horizontal').grid(
            row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=20
        )
        
        ttk.Button(
            main_frame,
            text="Analyze",
            command=self.analyze,
            style='Accent.TButton'
        ).grid(row=4, column=0, columnspan=3, pady=10)
        
        ttk.Label(main_frame, text="Results Preview:").grid(row=5, column=0, sticky=tk.W, pady=5)
        
        tree_frame = ttk.Frame(main_frame)
        tree_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.rowconfigure(0, weight=1)
        main_frame.rowconfigure(6, weight=1)
        
        self.tree = ttk.Treeview(tree_frame, show='headings')
        self.tree['columns'] = ('Gene Name', 'Start', 'End', 'Midpoint', 'Pi at Start', 'Pi at End', 'Pi at Midpoint', 'Average Pi')
        
        for col in self.tree['columns']:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=100, anchor=tk.CENTER)
        
        scrollbar = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        
        self.tree.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        ttk.Separator(main_frame, orient='horizontal').grid(
            row=7, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=20
        )
        
        ttk.Button(main_frame, text="Export Results", command=self.export_results).grid(row=8, column=0, columnspan=3, pady=5)
        
        self.status_label = ttk.Label(main_frame, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_label.grid(row=9, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(20, 0))
        
        self.set_default_files()
    
    def set_default_files(self):
        default_pi = r"C:\Users\MYzhang\Desktop\测试\Pi.txt"
        default_gff = r"C:\Users\MYzhang\Desktop\测试\sequence.gff3"
        if os.path.exists(default_pi):
            self.pi_entry.insert(0, default_pi)
            self.pi_file_path = default_pi
        if os.path.exists(default_gff):
            self.gff_entry.insert(0, default_gff)
            self.gff_file_path = default_gff
    
    def browse_pi_file(self):
        file_path = filedialog.askopenfilename(
            title="Select Pi File",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if file_path:
            self.pi_entry.delete(0, tk.END)
            self.pi_entry.insert(0, file_path)
            self.pi_file_path = file_path
    
    def browse_gff_file(self):
        file_path = filedialog.askopenfilename(
            title="Select GFF3 File",
            filetypes=[("GFF3 Files", "*.gff3"), ("All Files", "*.*")]
        )
        if file_path:
            self.gff_entry.delete(0, tk.END)
            self.gff_entry.insert(0, file_path)
            self.gff_file_path = file_path
    
    def parse_pi_file(self, file_path):
        pi_data = []
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    if not line or line.startswith('Window'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 3:
                        window = parts[0]
                        midpoint = float(parts[1])
                        pi = float(parts[2])
                        
                        start_end = window.split('-')
                        if len(start_end) == 2:
                            start = int(start_end[0])
                            end = int(start_end[1])
                            pi_data.append({
                                'start': start,
                                'end': end,
                                'midpoint': midpoint,
                                'pi': pi
                            })
        except Exception as e:
            raise Exception(f"Error parsing Pi file: {str(e)}")
        
        return pi_data
    
    def parse_gff_file(self, file_path):
        genes = {}
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue
                    
                    feature_type = parts[2]
                    if feature_type != 'gene':
                        continue
                    
                    start = int(parts[3])
                    end = int(parts[4])
                    
                    attributes = parts[8]
                    gene_name = self.extract_gene_name(attributes)
                    
                    if gene_name:
                        if gene_name not in genes:
                            genes[gene_name] = {'start': start, 'end': end}
                        else:
                            if start < genes[gene_name]['start']:
                                genes[gene_name]['start'] = start
                            if end > genes[gene_name]['end']:
                                genes[gene_name]['end'] = end
        except Exception as e:
            raise Exception(f"Error parsing GFF3 file: {str(e)}")
        
        return genes
    
    def extract_gene_name(self, attributes):
        match = re.search(r'gene=([^;]+)', attributes)
        if match:
            return match.group(1)
        match = re.search(r'Name=([^;]+)', attributes)
        if match:
            return match.group(1)
        return None
    
    def find_pi_at_position(self, position, pi_data):
        for window in pi_data:
            if window['start'] <= position <= window['end']:
                return window['pi']
        return None
    
    def calculate_gene_pi(self, gene_start, gene_end, pi_data):
        relevant_windows = []
        for window in pi_data:
            if window['end'] < gene_start or window['start'] > gene_end:
                continue
            
            overlap_start = max(window['start'], gene_start)
            overlap_end = min(window['end'], gene_end)
            overlap_length = overlap_end - overlap_start + 1
            
            relevant_windows.append({
                'pi': window['pi'],
                'length': overlap_length
            })
        
        if not relevant_windows:
            return None, None, None, None, None
        
        gene_mid = (gene_start + gene_end) // 2
        
        pi_at_start = self.find_pi_at_position(gene_start, pi_data)
        pi_at_end = self.find_pi_at_position(gene_end, pi_data)
        pi_at_mid = self.find_pi_at_position(gene_mid, pi_data)
        
        total_length = sum(w['length'] for w in relevant_windows)
        weighted_pi = sum(w['pi'] * w['length'] for w in relevant_windows)
        avg_pi = weighted_pi / total_length if total_length > 0 else 0
        
        return avg_pi, total_length, pi_at_start, pi_at_end, pi_at_mid
    
    def analyze(self):
        if not self.pi_file_path or not self.gff_file_path:
            messagebox.showerror("Error", "Please select both Pi file and GFF3 file!")
            return

        self.status_label.config(text="Analyzing...")
        self.root.update()

        try:
            pi_data = self.parse_pi_file(self.pi_file_path)
            genes = self.parse_gff_file(self.gff_file_path)

            if not pi_data:
                messagebox.showerror("Error", "No valid data found in Pi file!")
                self.status_label.config(text="Analysis failed")
                return

            if not genes:
                messagebox.showerror("Error", "No genes found in GFF3 file!")
                self.status_label.config(text="Analysis failed")
                return

            results = []
            for gene_name, gene_info in sorted(genes.items()):
                gene_start = gene_info['start']
                gene_end = gene_info['end']
                gene_mid = (gene_start + gene_end) // 2

                avg_pi, coverage, pi_at_start, pi_at_end, pi_at_mid = self.calculate_gene_pi(gene_start, gene_end, pi_data)

                if avg_pi is not None:
                    results.append({
                        'Gene Name': gene_name,
                        'Start': gene_start,
                        'End': gene_end,
                        'Midpoint': gene_mid,
                        'Pi at Start': pi_at_start if pi_at_start is not None else 0,
                        'Pi at End': pi_at_end if pi_at_end is not None else 0,
                        'Pi at Midpoint': pi_at_mid if pi_at_mid is not None else 0,
                        'Average Pi': avg_pi,
                        'Coverage': coverage
                    })

            self.result_df = pd.DataFrame(results)

            self.update_tree(self.result_df)

            self.status_label.config(text=f"Analysis complete! Processed {len(results)} genes")

        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text="Analysis failed")
    
    def update_tree(self, df):
        for item in self.tree.get_children():
            self.tree.delete(item)

        for _, row in df.iterrows():
            self.tree.insert('', tk.END, values=(
                row['Gene Name'],
                row['Start'],
                row['End'],
                row['Midpoint'],
                f"{row['Pi at Start']:.5f}" if row['Pi at Start'] > 0 else 'N/A',
                f"{row['Pi at End']:.5f}" if row['Pi at End'] > 0 else 'N/A',
                f"{row['Pi at Midpoint']:.5f}" if row['Pi at Midpoint'] > 0 else 'N/A',
                f"{row['Average Pi']:.5f}"
            ))
    
    def export_results(self):
        if self.result_df is None or self.result_df.empty:
            messagebox.showwarning("Warning", "No results to export!")
            return

        file_path = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".xlsx",
            filetypes=[
                ("Excel Files", "*.xlsx"),
                ("CSV Files", "*.csv"),
                ("All Files", "*.*")
            ]
        )

        if file_path:
            try:
                export_columns = ['Gene Name', 'Start', 'End', 'Midpoint', 'Pi at Start', 'Pi at End', 'Pi at Midpoint', 'Average Pi']
                if file_path.endswith('.xlsx'):
                    self.result_df[export_columns].to_excel(
                        file_path, index=False, float_format='%.5f'
                    )
                elif file_path.endswith('.csv'):
                    self.result_df[export_columns].to_csv(
                        file_path, index=False, float_format='%.5f'
                    )
                messagebox.showinfo("Success", f"Results saved to:\n{file_path}")
                self.status_label.config(text=f"Results exported to: {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {str(e)}")


import os

if __name__ == "__main__":
    root = tk.Tk()
    app = PiGeneAnalyzer(root)
    root.mainloop()
