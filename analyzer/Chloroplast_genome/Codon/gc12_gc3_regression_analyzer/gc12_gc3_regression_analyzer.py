import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import numpy as np
from scipy import stats
import csv


class GC12GC3Analyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("GC12-GC3 Regression Analyzer")
        self.root.geometry("1100x650")
        
        self.results_data = []
        self.input_folder = None
        
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
        
        # Save button
        save_btn = ttk.Button(control_frame, text="Save Results", command=self.save_results)
        save_btn.grid(row=0, column=4, padx=5, pady=5)
        
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
        
        columns = ('Species', 'Gene Count', 'Slope', 'Intercept', 'R²', 
                   'P-value', 'Std Error', 'Pearson r', 'Mean GC12', 'Mean GC3')
        self.tree = ttk.Treeview(tree_frame, columns=columns, show='headings',
                                 yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)
        
        column_widths = {
            'Species': 100,
            'Gene Count': 80,
            'Slope': 70,
            'Intercept': 75,
            'R²': 60,
            'P-value': 100,
            'Std Error': 80,
            'Pearson r': 80,
            'Mean GC12': 80,
            'Mean GC3': 80
        }
        
        for col in columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=column_widths[col], anchor=tk.CENTER)
        
        self.tree.column('Species', anchor=tk.W)
        
        scrollbar_y.config(command=self.tree.yview)
        scrollbar_x.config(command=self.tree.xview)
        
        self.tree.pack(fill=tk.BOTH, expand=True)
        
        # Status bar
        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(fill=tk.X, padx=5, pady=5)
    
    def browse_folder(self):
        folder = filedialog.askdirectory(title="Select folder containing GC analysis CSV files")
        if folder:
            self.folder_path.set(folder)
            self.input_folder = folder
            self.status_bar.config(text=f"Selected: {folder}")
    
    def parse_csv_file(self, filepath):
        """Parse a CSV file and extract GC12 and GC3 data"""
        gc12_values = []
        gc3_values = []
        gene_count = 0
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                header = next(reader)  # Skip header
                
                # Find column indices
                gc12_idx = None
                gc3_idx = None
                
                for i, col in enumerate(header):
                    col_clean = col.strip().lower()
                    # Match "gc3 (%)" column (not gc3s)
                    if col_clean == 'gc3 (%)':
                        gc3_idx = i
                    # Match "gc12 (%)" column
                    elif col_clean == 'gc12 (%)':
                        gc12_idx = i
                
                # Fallback: try pattern matching if exact match failed
                if gc12_idx is None or gc3_idx is None:
                    for i, col in enumerate(header):
                        col_lower = col.lower().strip()
                        if 'gc12' in col_lower and '(' in col_lower and 'gc3' not in col_lower:
                            gc12_idx = i
                        elif 'gc3' in col_lower and '(' in col_lower and 'gc3s' not in col_lower:
                            gc3_idx = i
                
                if gc12_idx is None or gc3_idx is None:
                    print(f"Could not find GC12 or GC3 columns in {filepath}")
                    print(f"Header: {header}")
                    return None
                
                for row in reader:
                    if len(row) <= max(gc12_idx, gc3_idx):
                        continue
                    
                    try:
                        gc12 = float(row[gc12_idx])
                        gc3 = float(row[gc3_idx])
                        
                        gc12_values.append(gc12)
                        gc3_values.append(gc3)
                        gene_count += 1
                        
                    except (ValueError, IndexError):
                        continue
            
        except Exception as e:
            print(f"Error parsing {filepath}: {e}")
            return None
        
        if gene_count < 2:
            return None
        
        return {
            'gc12': np.array(gc12_values),
            'gc3': np.array(gc3_values),
            'count': gene_count
        }
    
    def perform_regression(self, gc3, gc12):
        """
        Perform linear regression using scipy.stats.linregress
        Regression: GC12 = slope * GC3 + intercept (GC3 as X, GC12 as Y)
        This matches the reference program gc_analysis_plotter.py
        """
        # Use scipy.stats.linregress (GC3 as x, GC12 as y)
        slope, intercept, r_value, p_value, std_err = stats.linregress(gc3, gc12)
        
        # R-squared
        r_squared = r_value ** 2
        
        # Pearson correlation (same as r_value from linregress)
        pearson_r = r_value
        
        return {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'p_value': p_value,
            'std_err': std_err,
            'pearson_r': pearson_r
        }
    
    def analyze_folder(self):
        folder = self.folder_path.get()
        
        if not folder or not os.path.exists(folder):
            messagebox.showerror("Error", "Please select a valid folder!")
            return
        
        # Find all CSV files
        csv_files = [f for f in os.listdir(folder) if f.endswith('.csv')]
        
        if not csv_files:
            messagebox.showerror("Error", "No CSV files found in the selected folder!")
            return
        
        self.input_folder = folder
        self.status_bar.config(text=f"Analyzing {len(csv_files)} files...")
        self.root.update()
        
        # Clear previous results
        for item in self.tree.get_children():
            self.tree.delete(item)
        
        self.results_data = []
        
        for filename in sorted(csv_files):
            filepath = os.path.join(folder, filename)
            data = self.parse_csv_file(filepath)
            
            if data is None:
                continue
            
            # Get species name from filename (remove extension)
            species_name = os.path.splitext(filename)[0]
            
            # Perform regression: GC12 = slope * GC3 + intercept
            reg_results = self.perform_regression(data['gc3'], data['gc12'])
            
            # Calculate means
            mean_gc12 = np.mean(data['gc12'])
            mean_gc3 = np.mean(data['gc3'])
            
            # Store result
            result = {
                'species': species_name,
                'count': data['count'],
                'slope': reg_results['slope'],
                'intercept': reg_results['intercept'],
                'r_squared': reg_results['r_squared'],
                'p_value': reg_results['p_value'],
                'std_err': reg_results['std_err'],
                'pearson_r': reg_results['pearson_r'],
                'mean_gc12': mean_gc12,
                'mean_gc3': mean_gc3
            }
            self.results_data.append(result)
            
            # Add to treeview
            self.tree.insert('', tk.END, values=(
                species_name,
                data['count'],
                f"{reg_results['slope']:.4f}",
                f"{reg_results['intercept']:.4f}",
                f"{reg_results['r_squared']:.4f}",
                f"{reg_results['p_value']:.2e}",
                f"{reg_results['std_err']:.4f}",
                f"{reg_results['pearson_r']:.4f}",
                f"{mean_gc12:.2f}",
                f"{mean_gc3:.2f}"
            ))
        
        self.status_bar.config(text=f"Analysis complete - {len(self.results_data)} species analyzed")
    
    def save_results(self):
        """Export results to CSV file"""
        if not self.results_data:
            messagebox.showwarning("Warning", "No results to save. Please analyze data first!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialdir=self.input_folder if self.input_folder else None
        )
        
        if not filename:
            return
        
        try:
            with open(filename, 'w', encoding='utf-8', newline='') as f:
                writer = csv.writer(f)
                
                # Header
                writer.writerow([
                    'Species', 'Gene Count', 'Slope', 'Intercept', 'R²', 
                    'P-value', 'Std Error', 'Pearson r', 'Mean GC12', 'Mean GC3',
                    'Regression Equation'
                ])
                
                # Data
                for r in self.results_data:
                    equation = f"GC12 = {r['slope']:.4f} * GC3 + {r['intercept']:.4f}"
                    writer.writerow([
                        r['species'],
                        r['count'],
                        f"{r['slope']:.4f}",
                        f"{r['intercept']:.4f}",
                        f"{r['r_squared']:.4f}",
                        f"{r['p_value']:.2e}",
                        f"{r['std_err']:.4f}",
                        f"{r['pearson_r']:.4f}",
                        f"{r['mean_gc12']:.2f}",
                        f"{r['mean_gc3']:.2f}",
                        equation
                    ])
            
            messagebox.showinfo("Success", f"Results exported to:\n{filename}")
            self.status_bar.config(text=f"Exported: {os.path.basename(filename)}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {str(e)}")


def main():
    root = tk.Tk()
    app = GC12GC3Analyzer(root)
    root.mainloop()


if __name__ == "__main__":
    main()
