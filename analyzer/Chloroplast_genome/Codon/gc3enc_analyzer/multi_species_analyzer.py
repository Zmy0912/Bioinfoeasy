import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import numpy as np
from scipy import stats


class MultiSpeciesENCAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Multi-Species ENC-GC3 Analyzer")
        self.root.geometry("1000x700")
        
        self.results_data = []
        
        self.create_widgets()
    
    def create_widgets(self):
        # Control panel
        control_frame = ttk.LabelFrame(self.root, text="Folder Selection", padding=10)
        control_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Folder path
        ttk.Label(control_frame, text="Data folder:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.folder_path = tk.StringVar()
        folder_entry = ttk.Entry(control_frame, textvariable=self.folder_path, width=50)
        folder_entry.grid(row=0, column=1, padx=5, pady=5)
        
        browse_btn = ttk.Button(control_frame, text="Browse...", command=self.browse_folder)
        browse_btn.grid(row=0, column=2, padx=5, pady=5)
        
        # Analyze button
        analyze_btn = ttk.Button(control_frame, text="Analyze", command=self.analyze_folder)
        analyze_btn.grid(row=0, column=3, padx=5, pady=5)
        
        # Export button
        export_btn = ttk.Button(control_frame, text="Export CSV", command=self.export_results)
        export_btn.grid(row=0, column=4, padx=5, pady=5)
        
        # Results panel
        results_frame = ttk.LabelFrame(self.root, text="Results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Treeview with scrollbar
        tree_frame = ttk.Frame(results_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        scrollbar_y = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL)
        scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
        
        scrollbar_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL)
        scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)
        
        columns = ('Species', 'Gene Count', 'Mean ENC', 'Std ENC', 'Min ENC', 'Max ENC',
                   'Mean GC3', 'Std GC3', 'Min GC3', 'Max GC3')
        self.tree = ttk.Treeview(tree_frame, columns=columns, show='headings',
                                 yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)
        
        for col in columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=100, anchor=tk.CENTER)
        
        self.tree.column('Species', width=150, anchor=tk.W)
        
        scrollbar_y.config(command=self.tree.yview)
        scrollbar_x.config(command=self.tree.xview)
        
        self.tree.pack(fill=tk.BOTH, expand=True)
        
        # Summary frame
        summary_frame = ttk.LabelFrame(self.root, text="Summary Statistics", padding=10)
        summary_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.summary_text = tk.Text(summary_frame, height=6, width=120, font=('Courier', 9))
        self.summary_text.pack(fill=tk.X)
        
        # Status bar
        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(fill=tk.X, padx=5, pady=5)
    
    def browse_folder(self):
        folder = filedialog.askdirectory(title="Select folder containing .out files")
        if folder:
            self.folder_path.set(folder)
            self.status_bar.config(text=f"Selected: {folder}")
    
    def parse_out_file(self, filepath):
        """Parse a single .out file and extract ENC and GC3 data"""
        enc_values = []
        gc3_values = []
        gene_count = 0
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Skip header line
            for line in lines[1:]:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) < 10:
                    continue
                
                try:
                    enc_str = parts[8]
                    gc3_str = parts[9]
                    
                    # Skip invalid data
                    if enc_str == '*****' or gc3_str == '*****':
                        continue
                    
                    enc = float(enc_str)
                    gc3 = float(gc3_str)
                    
                    enc_values.append(enc)
                    gc3_values.append(gc3)
                    gene_count += 1
                    
                except (ValueError, IndexError):
                    continue
            
        except Exception as e:
            print(f"Error parsing {filepath}: {e}")
            return None
        
        if gene_count == 0:
            return None
        
        return {
            'enc': np.array(enc_values),
            'gc3': np.array(gc3_values),
            'count': gene_count
        }
    
    def analyze_folder(self):
        folder = self.folder_path.get()
        
        if not folder or not os.path.exists(folder):
            messagebox.showerror("Error", "Please select a valid folder!")
            return
        
        # Find all .out files
        out_files = [f for f in os.listdir(folder) if f.endswith('.out')]
        
        if not out_files:
            messagebox.showerror("Error", "No .out files found in the selected folder!")
            return
        
        self.status_bar.config(text=f"Analyzing {len(out_files)} files...")
        self.root.update()
        
        # Clear previous results
        for item in self.tree.get_children():
            self.tree.delete(item)
        
        self.results_data = []
        
        all_enc_values = []
        all_gc3_values = []
        
        for filename in sorted(out_files):
            filepath = os.path.join(folder, filename)
            data = self.parse_out_file(filepath)
            
            if data is None:
                continue
            
            # Get species name from filename (remove extension)
            species_name = os.path.splitext(filename)[0]
            
            # Calculate statistics
            mean_enc = np.mean(data['enc'])
            std_enc = np.std(data['enc'])
            min_enc = np.min(data['enc'])
            max_enc = np.max(data['enc'])
            
            mean_gc3 = np.mean(data['gc3'])
            std_gc3 = np.std(data['gc3'])
            min_gc3 = np.min(data['gc3'])
            max_gc3 = np.max(data['gc3'])
            
            # Store result
            result = {
                'species': species_name,
                'count': data['count'],
                'mean_enc': mean_enc,
                'std_enc': std_enc,
                'min_enc': min_enc,
                'max_enc': max_enc,
                'mean_gc3': mean_gc3,
                'std_gc3': std_gc3,
                'min_gc3': min_gc3,
                'max_gc3': max_gc3
            }
            self.results_data.append(result)
            
            # Collect all values for overall statistics
            all_enc_values.extend(data['enc'])
            all_gc3_values.extend(data['gc3'])
            
            # Add to treeview
            self.tree.insert('', tk.END, values=(
                species_name,
                data['count'],
                f"{mean_enc:.2f}",
                f"{std_enc:.2f}",
                f"{min_enc:.2f}",
                f"{max_enc:.2f}",
                f"{mean_gc3:.3f}",
                f"{std_gc3:.3f}",
                f"{min_gc3:.3f}",
                f"{max_gc3:.3f}"
            ))
        
        # Update summary
        self.update_summary(all_enc_values, all_gc3_values)
        
        self.status_bar.config(text=f"Analysis complete - {len(self.results_data)} species analyzed")
    
    def update_summary(self, all_enc, all_gc3):
        """Update summary statistics text"""
        if not all_enc:
            return
        
        all_enc = np.array(all_enc)
        all_gc3 = np.array(all_gc3)
        
        summary = []
        summary.append("=" * 80)
        summary.append("OVERALL STATISTICS ACROSS ALL SPECIES")
        summary.append("=" * 80)
        summary.append("")
        summary.append(f"Total species analyzed:    {len(self.results_data)}")
        summary.append(f"Total genes analyzed:       {len(all_enc)}")
        summary.append("")
        summary.append("ENC Statistics:")
        summary.append(f"  Overall Mean ENC:        {np.mean(all_enc):.2f}")
        summary.append(f"  Overall Std ENC:         {np.std(all_enc):.2f}")
        summary.append(f"  Overall Min ENC:         {np.min(all_enc):.2f}")
        summary.append(f"  Overall Max ENC:         {np.max(all_enc):.2f}")
        summary.append(f"  Median ENC:              {np.median(all_enc):.2f}")
        summary.append("")
        summary.append("GC3 Statistics:")
        summary.append(f"  Overall Mean GC3:        {np.mean(all_gc3):.3f}")
        summary.append(f"  Overall Std GC3:        {np.std(all_gc3):.3f}")
        summary.append(f"  Overall Min GC3:        {np.min(all_gc3):.3f}")
        summary.append(f"  Overall Max GC3:         {np.max(all_gc3):.3f}")
        summary.append(f"  Median GC3:              {np.median(all_gc3):.3f}")
        summary.append("")
        
        # Spearman correlation between mean ENC and mean GC3
        if len(self.results_data) > 2:
            mean_encs = [r['mean_enc'] for r in self.results_data]
            mean_gc3s = [r['mean_gc3'] for r in self.results_data]
            corr, p_value = stats.spearmanr(mean_encs, mean_gc3s)
            summary.append(f"Spearman correlation (mean ENC vs mean GC3): {corr:.4f} (p={p_value:.4e})")
        
        summary.append("")
        summary.append("=" * 80)
        
        self.summary_text.delete('1.0', tk.END)
        self.summary_text.insert('1.0', '\n'.join(summary))
    
    def export_results(self):
        """Export results to CSV file"""
        if not self.results_data:
            messagebox.showwarning("Warning", "No results to export. Please analyze data first!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export Results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if not filename:
            return
        
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                # Header
                f.write("Species,Gene Count,Mean ENC,Std ENC,Min ENC,Max ENC,")
                f.write("Mean GC3,Std GC3,Min GC3,Max GC3\n")
                
                # Data
                for r in self.results_data:
                    f.write(f"{r['species']},{r['count']},{r['mean_enc']:.2f},{r['std_enc']:.2f},")
                    f.write(f"{r['min_enc']:.2f},{r['max_enc']:.2f},")
                    f.write(f"{r['mean_gc3']:.3f},{r['std_gc3']:.3f},")
                    f.write(f"{r['min_gc3']:.3f},{r['max_gc3']:.3f}\n")
            
            messagebox.showinfo("Success", f"Results exported to:\n{filename}")
            self.status_bar.config(text=f"Exported: {os.path.basename(filename)}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {str(e)}")


def main():
    root = tk.Tk()
    app = MultiSpeciesENCAnalyzer(root)
    root.mainloop()


if __name__ == "__main__":
    main()
