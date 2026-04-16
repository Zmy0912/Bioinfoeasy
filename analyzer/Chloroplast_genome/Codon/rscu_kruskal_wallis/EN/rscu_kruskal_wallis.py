# -*- coding: utf-8 -*-
"""
RSCU Kruskal-Wallis Analysis Software
用于对不同物种间密码子RSCU值进行Kruskal-Wallis非参数检验
Kruskal-Wallis non-parametric test for RSCU values across species
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import numpy as np
from scipy.stats import kruskal
import os
import warnings
warnings.filterwarnings('ignore')


class RSCUAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("RSCU Kruskal-Wallis Analysis Tool")
        self.root.geometry("1000x800")
        self.root.configure(bg="#f0f0f0")
        
        self.data = None
        self.results = None
        self.amino_acids = {}
        
        self.setup_styles()
        self.create_widgets()
        
    def setup_styles(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure("Title.TLabel", font=("Arial", 16, "bold"), background="#f0f0f0")
        style.configure("Normal.TLabel", font=("Arial", 10), background="#f0f0f0")
        style.configure("Header.TLabel", font=("Arial", 11, "bold"), background="#f0f0f0")
        
    def create_widgets(self):
        # Title
        title_frame = tk.Frame(self.root, bg="#2c3e50")
        title_frame.pack(fill=tk.X)
        ttk.Label(title_frame, text="RSCU Kruskal-Wallis Analysis Tool", 
                  style="Title.TLabel", foreground="white", background="#2c3e50").pack(pady=10)
        
        # Info section
        info_frame = tk.LabelFrame(self.root, text="Info", bg="#f0f0f0", font=("Arial", 10))
        info_frame.pack(fill=tk.X, padx=20, pady=10)
        
        info_text = """This software performs Kruskal-Wallis tests on RSCU (Relative Synonymous Codon Usage) 
data to analyze differences in codon usage preferences between species groups.

- Kruskal-Wallis is a non-parametric test for comparing multiple groups
- For each codon, the test compares RSCU values across species groups
- Significant results indicate different codon usage preferences between groups"""
        ttk.Label(info_frame, text=info_text, background="#f0f0f0", justify=tk.LEFT, 
                  font=("Arial", 9)).pack(anchor=tk.W, padx=10, pady=5)
        
        # File selection
        file_frame = tk.LabelFrame(self.root, text="Data File", bg="#f0f0f0", font=("Arial", 10))
        file_frame.pack(fill=tk.X, padx=20, pady=10)
        
        self.file_path_var = tk.StringVar(value="No file selected")
        ttk.Label(file_frame, text="RSCU Data File:", background="#f0f0f0").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        ttk.Entry(file_frame, textvariable=self.file_path_var, width=60).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(file_frame, text="Browse...", command=self.browse_file).grid(row=0, column=2, padx=5, pady=5)
        ttk.Button(file_frame, text="Load Data", command=self.load_data).grid(row=0, column=3, padx=5, pady=5)
        
        # Group settings
        group_frame = tk.LabelFrame(self.root, text="Species Grouping", bg="#f0f0f0", font=("Arial", 10))
        group_frame.pack(fill=tk.X, padx=20, pady=10)
        
        self.group_method = tk.StringVar(value="auto")
        
        ttk.Label(group_frame, text="Grouping Method:", background="#f0f0f0").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(group_frame, text="Auto (by species prefix)", 
                        variable=self.group_method, value="auto").grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Radiobutton(group_frame, text="Manual (load group file)", 
                        variable=self.group_method, value="manual").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        self.group_file_var = tk.StringVar(value="No group file selected")
        ttk.Label(group_frame, text="Group File:", background="#f0f0f0").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        ttk.Entry(group_frame, textvariable=self.group_file_var, width=40).grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(group_frame, text="Browse...", command=self.browse_group_file).grid(row=1, column=2, padx=5, pady=5)
        
        # Group preview
        self.group_preview_text = tk.Text(group_frame, height=5, width=80, font=("Consolas", 9))
        self.group_preview_text.grid(row=2, column=0, columnspan=3, padx=10, pady=5)
        
        # Analysis options
        option_frame = tk.LabelFrame(self.root, text="Analysis Options", bg="#f0f0f0", font=("Arial", 10))
        option_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Label(option_frame, text="Analysis Scope:", background="#f0f0f0").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.analysis_scope = tk.StringVar(value="all")
        ttk.Radiobutton(option_frame, text="All Codons", variable=self.analysis_scope, value="all").grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Radiobutton(option_frame, text="By Amino Acid", variable=self.analysis_scope, value="by_amino").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        # Amino acid selection
        self.amino_frame = tk.Frame(option_frame, bg="#f0f0f0")
        self.amino_frame.grid(row=1, column=0, columnspan=4, sticky=tk.W, padx=10, pady=5)
        
        self.amino_vars = {}
        amino_list = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 
                      'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
        for i, amino in enumerate(amino_list):
            var = tk.BooleanVar(value=True)
            self.amino_vars[amino] = var
            cb = ttk.Checkbutton(self.amino_frame, text=amino, variable=var)
            cb.grid(row=i // 5, column=i % 5, sticky=tk.W, padx=2)
        
        # Significance level
        self.significance_level = tk.DoubleVar(value=0.05)
        ttk.Label(option_frame, text="Significance Level (alpha):", background="#f0f0f0").grid(row=0, column=3, sticky=tk.W, padx=(20, 5))
        ttk.Entry(option_frame, textvariable=self.significance_level, width=8).grid(row=0, column=4, padx=5)
        
        # Buttons
        button_frame = tk.Frame(self.root, bg="#f0f0f0")
        button_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Button(button_frame, text="Run Kruskal-Wallis Test", 
                   command=self.run_analysis).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Plot Boxplots", 
                   command=self.plot_boxplot).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Export Full Results", 
                   command=self.export_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Export Significant", 
                   command=self.export_significant).pack(side=tk.LEFT, padx=5)
        
        # Results display
        result_frame = tk.LabelFrame(self.root, text="Results", bg="#f0f0f0", font=("Arial", 10))
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        result_scroll = ttk.Scrollbar(result_frame)
        result_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.result_text = tk.Text(result_frame, yscrollcommand=result_scroll.set,
                                    font=("Consolas", 9), bg="#ffffff", wrap=tk.NONE)
        self.result_text.pack(fill=tk.BOTH, expand=True)
        result_scroll.config(command=self.result_text.yview)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = tk.Label(self.root, textvariable=self.status_var, bd=1, relief=tk.SUNKEN, 
                              anchor=tk.W, bg="#e0e0e0", font=("Arial", 9))
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Select RSCU Data File",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if filename:
            self.file_path_var.set(filename)
            
    def browse_group_file(self):
        filename = filedialog.askopenfilename(
            title="Select Group File",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("CSV files", "*.csv"), ("Text files", "*.txt")]
        )
        if filename:
            self.group_file_var.set(filename)
            
    def load_data(self):
        filepath = self.file_path_var.get()
        if not filepath or filepath == "No file selected":
            messagebox.showwarning("Warning", "Please select a data file first")
            return
            
        try:
            self.data = pd.read_excel(filepath)
            
            if 'Species' not in self.data.columns:
                messagebox.showerror("Error", "Data file must contain 'Species' column")
                return
                
            self.status_var.set(f"Loaded: {len(self.data)} species, {len(self.data.columns)-1} codons")
            
            # Parse amino acid info
            self.parse_amino_acids()
            
            # Process groups
            self.process_groups()
            
            messagebox.showinfo("Success", f"Data loaded successfully!\n{len(self.data)} species\n{len(self.data.columns)-1} codons\n\nSpecies list:\n" + "\n".join(self.data['Species'].tolist()[:10]) + ("\n..." if len(self.data) > 10 else ""))
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load data: {str(e)}")
            
    def parse_amino_acids(self):
        """Parse codon column names"""
        self.amino_acids = {}
        codon_columns = [col for col in self.data.columns if col != 'Species']
        
        for col in codon_columns:
            try:
                codon = col.split('(')[0]
                amino = col.split('(')[1].rstrip(')')
                if amino not in self.amino_acids:
                    self.amino_acids[amino] = []
                self.amino_acids[amino].append(col)
            except:
                pass
                
    def process_groups(self):
        """Process species groups"""
        self.group_preview_text.delete(1.0, tk.END)
        method = self.group_method.get()
        
        if method == "auto":
            species = self.data['Species'].tolist()
            
            # Try to find natural species groupings
            groups = {}
            for sp in species:
                parts = sp.replace('_', ' ').split()
                if len(parts) >= 2:
                    key = parts[0]
                else:
                    key = sp[:8] if len(sp) > 8 else sp
                    
                if key not in groups:
                    groups[key] = []
                groups[key].append(sp)
            
            # If too few groups, provide demo grouping
            if len(groups) < 3:
                self.group_preview_text.insert(tk.END, "Note: All species may belong to the same group\n")
                self.group_preview_text.insert(tk.END, "Use 'Manual' mode to upload a group file\n")
                self.group_preview_text.insert(tk.END, "Or provide geographic/phyylogenetic grouping info\n\n")
                self.group_preview_text.insert(tk.END, "[Using demo grouping for analysis]\n\n")
                n = len(species)
                num_groups = min(5, max(2, n // 5))
                groups = {}
                for i, sp in enumerate(species):
                    group_idx = i % num_groups + 1
                    key = f"Demo_Group{group_idx}"
                    if key not in groups:
                        groups[key] = []
                    groups[key].append(sp)
            
            self.species_groups = groups
            self.group_method.set("auto")
            
        else:
            # Load groups from file
            group_file = self.group_file_var.get()
            if not group_file or group_file == "No group file selected":
                messagebox.showwarning("Warning", "Please select a group file, or switch to auto mode")
                self.group_method.set("auto")
                self.process_groups()
                return
            self.load_groups_from_file(group_file)
            
        # Display group preview
        self.group_preview_text.insert(tk.END, f"Total {len(self.species_groups)} groups:\n")
        self.group_preview_text.insert(tk.END, "-" * 50 + "\n")
        for group_name, members in self.species_groups.items():
            self.group_preview_text.insert(tk.END, f"[{group_name}] ({len(members)} species): {', '.join(members)}\n")
            
    def load_groups_from_file(self, filepath):
        """Load groups from file"""
        try:
            if filepath.endswith('.xlsx') or filepath.endswith('.xls'):
                group_df = pd.read_excel(filepath)
            elif filepath.endswith('.csv'):
                group_df = pd.read_csv(filepath)
            else:
                group_df = pd.read_csv(filepath, delim_whitespace=True)
                    
            self.species_groups = {}
            for _, row in group_df.iterrows():
                species = str(row.iloc[0])
                group = str(row.iloc[1])
                if group not in self.species_groups:
                    self.species_groups[group] = []
                self.species_groups[group].append(species)
                
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load group file: {str(e)}")
            
    def run_analysis(self):
        """Run Kruskal-Wallis test"""
        if self.data is None:
            messagebox.showwarning("Warning", "Please load data first")
            return
            
        self.status_var.set("Running analysis...")
        self.result_text.delete(1.0, tk.END)
        
        alpha = self.significance_level.get()
        scope = self.analysis_scope.get()
        
        # Determine codons to analyze
        if scope == "all":
            codons_to_analyze = [col for col in self.data.columns if col != 'Species']
        else:
            codons_to_analyze = []
            for amino, var in self.amino_vars.items():
                if var.get() and amino in self.amino_acids:
                    codons_to_analyze.extend(self.amino_acids[amino])
                    
        if not codons_to_analyze:
            messagebox.showwarning("Warning", "Please select at least one amino acid to analyze")
            return
            
        results_list = []
        
        # Header
        self.result_text.insert(tk.END, "=" * 90 + "\n")
        self.result_text.insert(tk.END, "RSCU Kruskal-Wallis Test Results\n")
        self.result_text.insert(tk.END, "=" * 90 + "\n\n")
        
        # Basic info
        self.result_text.insert(tk.END, "[Basic Information]\n")
        self.result_text.insert(tk.END, f"  Number of species: {len(self.data)}\n")
        self.result_text.insert(tk.END, f"  Number of groups: {len(self.species_groups)}\n")
        self.result_text.insert(tk.END, f"  Codons analyzed: {len(codons_to_analyze)}\n")
        self.result_text.insert(tk.END, f"  Significance level: alpha = {alpha}\n\n")
        
        # Group details
        self.result_text.insert(tk.END, "[Group Details]\n")
        for group_name, members in self.species_groups.items():
            self.result_text.insert(tk.END, f"  {group_name}: {len(members)} species\n")
        self.result_text.insert(tk.END, "\n")
        
        # Results header
        self.result_text.insert(tk.END, "-" * 90 + "\n")
        self.result_text.insert(tk.END, f"{'Codon':<15} {'AA':<6} {'H-stat':<12} {'p-value':<15} {'Sig':<8} {'Interpretation'}\n")
        self.result_text.insert(tk.END, "-" * 90 + "\n")
        
        for codon in codons_to_analyze:
            try:
                # Extract amino acid name
                amino = codon.split('(')[1].rstrip(')') if '(' in codon else ''
                
                # Prepare group data
                groups_data = []
                group_names = list(self.species_groups.keys())
                for group_name in group_names:
                    species_in_group = self.species_groups[group_name]
                    mask = self.data['Species'].isin(species_in_group)
                    values = self.data.loc[mask, codon].values
                    groups_data.append(values)
                
                # Check data
                if any(len(g) < 2 for g in groups_data):
                    continue
                    
                # Kruskal-Wallis test
                stat, p_value = kruskal(*groups_data)
                
                # Significance
                if p_value < 0.001:
                    sig_level = "***"
                    sig_text = "Highly sig."
                elif p_value < 0.01:
                    sig_level = "**"
                    sig_text = "Very sig."
                elif p_value < alpha:
                    sig_level = "*"
                    sig_text = "Significant"
                else:
                    sig_level = "ns"
                    sig_text = "Not sig."
                
                # Generate interpretation
                interpretation = self.generate_interpretation(codon, groups_data, group_names, stat, p_value)
                
                results_list.append({
                    'Codon': codon,
                    'AminoAcid': amino,
                    'H_statistic': stat,
                    'p_value': p_value,
                    'Significance': sig_text,
                    'p_sign': sig_level,
                    'Interpretation': interpretation,
                    'GroupMeans': {gn: np.mean(g) for gn, g in zip(group_names, groups_data)},
                    'GroupMedians': {gn: np.median(g) for gn, g in zip(group_names, groups_data)},
                    'GroupSDs': {gn: np.std(g) for gn, g in zip(group_names, groups_data)}
                })
                
                self.result_text.insert(tk.END, f"{codon:<15} {amino:<6} {stat:<12.4f} {p_value:<15.2e} {sig_level:<8} {interpretation}\n")
                
            except Exception as e:
                pass
                
        self.results = results_list
        
        # Summary
        sig_count = sum(1 for r in results_list if r['p_value'] < alpha)
        self.result_text.insert(tk.END, "-" * 90 + "\n")
        self.result_text.insert(tk.END, "\n[Summary]\n")
        self.result_text.insert(tk.END, f"  Total codons analyzed: {len(results_list)}\n")
        self.result_text.insert(tk.END, f"  Significant codons (p < {alpha}): {sig_count}\n")
        self.result_text.insert(tk.END, f"  Non-significant: {len(results_list) - sig_count}\n")
        self.result_text.insert(tk.END, f"  Significant rate: {sig_count/len(results_list)*100:.1f}%\n\n")
        
        # Detailed interpretation
        self.result_text.insert(tk.END, "[Interpretation Guide]\n")
        if sig_count > 0:
            self.result_text.insert(tk.END, f"For {sig_count} codons with significant differences:\n\n")
            
            # Group by amino acid
            by_amino = {}
            for r in results_list:
                if r['p_value'] < alpha:
                    amino = r['AminoAcid']
                    if amino not in by_amino:
                        by_amino[amino] = []
                    by_amino[amino].append(r)
            
            for amino, codon_results in sorted(by_amino.items()):
                self.result_text.insert(tk.END, f"  {amino} ({len(codon_results)} codons significant):\n")
                for r in codon_results:
                    max_group = max(r['GroupMeans'].items(), key=lambda x: x[1])
                    min_group = min(r['GroupMeans'].items(), key=lambda x: x[1])
                    self.result_text.insert(tk.END, f"    - {r['Codon']}: p={r['p_value']:.2e}\n")
                    self.result_text.insert(tk.END, f"      Interpretation: {r['Significance']} difference between species groups\n")
                    self.result_text.insert(tk.END, f"      Highest preference: {max_group[0]} (mean={max_group[1]:.2f})\n")
                    self.result_text.insert(tk.END, f"      Lowest preference: {min_group[0]} (mean={min_group[1]:.2f})\n")
                self.result_text.insert(tk.END, "\n")
        
        self.status_var.set(f"Analysis complete, significant: {sig_count}/{len(results_list)}")
        
    def generate_interpretation(self, codon, groups_data, group_names, stat, p_value):
        """Generate result interpretation"""
        if p_value >= 0.05:
            return "No significant difference in codon usage"
        elif p_value >= 0.01:
            return "Some difference between species"
        elif p_value >= 0.001:
            return "Obvious difference between species"
        else:
            return "Very significant difference between species"
            
    def plot_boxplot(self):
        """Plot boxplots"""
        if self.results is None or len(self.results) == 0:
            messagebox.showwarning("Warning", "Please run analysis first")
            return
            
        try:
            import matplotlib.pyplot as plt
            plt.rcParams['font.sans-serif'] = ['Arial']
            plt.rcParams['axes.unicode_minus'] = False
            
            # Create figure
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle('RSCU Kruskal-Wallis Test - Boxplot Analysis', fontsize=14)
            
            # Select top 4 codons with smallest p-values
            top_results = sorted(self.results, key=lambda x: x['p_value'])[:4]
            
            for idx, (ax, result) in enumerate(zip(axes.flat, top_results)):
                codon = result['Codon']
                
                # Prepare boxplot data
                box_data = []
                labels = []
                for group_name in sorted(self.species_groups.keys()):
                    species_in_group = self.species_groups[group_name]
                    mask = self.data['Species'].isin(species_in_group)
                    values = self.data.loc[mask, codon].values
                    box_data.append(values)
                    labels.append(group_name)
                
                bp = ax.boxplot(box_data, labels=labels, patch_artist=True)
                
                # Set colors
                colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                ax.set_title(f"{codon}\np = {result['p_value']:.2e} ({result['Significance']})", fontsize=10)
                ax.set_ylabel('RSCU Value')
                ax.tick_params(axis='x', rotation=45)
                
            plt.tight_layout()
            
            # Save figure
            filepath = filedialog.asksaveasfilename(
                title="Save Boxplot",
                defaultextension=".png",
                filetypes=[("PNG image", "*.png"), ("PDF image", "*.pdf")],
                initialfile="RSCU_Boxplot"
            )
            
            if filepath:
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Boxplot saved to:\n{filepath}")
                
            plt.close()
            
        except ImportError:
            messagebox.showerror("Error", "Please install matplotlib: pip install matplotlib")
        except Exception as e:
            messagebox.showerror("Error", f"Plotting failed: {str(e)}")
            
    def export_results(self):
        """Export full results"""
        if self.results is None or len(self.results) == 0:
            messagebox.showwarning("Warning", "Please run analysis first")
            return
            
        filepath = filedialog.asksaveasfilename(
            title="Export Full Results",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("CSV files", "*.csv")],
            initialfile="RSCU_KruskalWallis_Full_Results"
        )
        
        if filepath:
            try:
                # Main results
                results_df = pd.DataFrame([{
                    'Codon': r['Codon'],
                    'AminoAcid': r['AminoAcid'],
                    'H_statistic': r['H_statistic'],
                    'p_value': r['p_value'],
                    'Significance': r['Significance'],
                    'Interpretation': r['Interpretation']
                } for r in self.results])
                
                # Group statistics
                group_names = sorted(self.species_groups.keys())
                stats_data = []
                for r in self.results:
                    row = {'Codon': r['Codon']}
                    for gn in group_names:
                        row[f'{gn}_Mean'] = r['GroupMeans'].get(gn, np.nan)
                        row[f'{gn}_Median'] = r['GroupMedians'].get(gn, np.nan)
                        row[f'{gn}_SD'] = r['GroupSDs'].get(gn, np.nan)
                    stats_data.append(row)
                stats_df = pd.DataFrame(stats_data)
                
                if filepath.endswith('.csv'):
                    results_df.to_csv(filepath, index=False, encoding='utf-8-sig')
                    stats_path = filepath.replace('.csv', '_GroupStats.csv')
                    stats_df.to_csv(stats_path, index=False, encoding='utf-8-sig')
                else:
                    with pd.ExcelWriter(filepath) as writer:
                        results_df.to_excel(writer, sheet_name='Test Results', index=False)
                        stats_df.to_excel(writer, sheet_name='Group Statistics', index=False)
                        
                messagebox.showinfo("Success", f"Results saved to:\n{filepath}\n\nGroup statistics saved")
                
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {str(e)}")
                
    def export_significant(self):
        """Export significant results"""
        if self.results is None or len(self.results) == 0:
            messagebox.showwarning("Warning", "Please run analysis first")
            return
            
        alpha = self.significance_level.get()
        sig_results = [r for r in self.results if r['p_value'] < alpha]
        
        if len(sig_results) == 0:
            messagebox.showinfo("Info", f"No significant results found (p < {alpha})")
            return
            
        filepath = filedialog.asksaveasfilename(
            title="Export Significant Results",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("CSV files", "*.csv")],
            initialfile="RSCU_KruskalWallis_Significant"
        )
        
        if filepath:
            try:
                sig_df = pd.DataFrame([{
                    'Codon': r['Codon'],
                    'AminoAcid': r['AminoAcid'],
                    'H_statistic': r['H_statistic'],
                    'p_value': r['p_value'],
                    'Significance': r['Significance'],
                    'Interpretation': r['Interpretation'],
                    'Highest_Group': max(r['GroupMeans'].items(), key=lambda x: x[1])[0],
                    'Highest_Mean': max(r['GroupMeans'].items(), key=lambda x: x[1])[1],
                    'Lowest_Group': min(r['GroupMeans'].items(), key=lambda x: x[1])[0],
                    'Lowest_Mean': min(r['GroupMeans'].items(), key=lambda x: x[1])[1]
                } for r in sig_results])
                
                if filepath.endswith('.csv'):
                    sig_df.to_csv(filepath, index=False, encoding='utf-8-sig')
                else:
                    sig_df.to_excel(filepath, index=False)
                    
                messagebox.showinfo("Success", f"Significant results saved to:\n{filepath}\n({len(sig_results)} records)")
                
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {str(e)}")


def main():
    root = tk.Tk()
    app = RSCUAnalyzerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
