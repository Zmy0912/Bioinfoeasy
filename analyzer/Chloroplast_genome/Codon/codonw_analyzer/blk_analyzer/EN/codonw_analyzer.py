"""
CodonW Codon Usage Bias Analyzer
A tool for processing CodonW output results and generating tables suitable for heatmap visualization
"""

import os
import re
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path


# Codon to amino acid mapping
CODON_AA_MAP = {
    'UUU': 'Phe', 'UUC': 'Phe',
    'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
    'AUG': 'Met',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr',
    'UAA': 'TER', 'UAG': 'TER', 'UGA': 'TER',
    'CAU': 'His', 'CAC': 'His',
    'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn',
    'AAA': 'Lys', 'AAG': 'Lys',
    'GAU': 'Asp', 'GAC': 'Asp',
    'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys',
    'UGG': 'Trp',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# All 61 codons (excluding stop codons)
ALL_CODONS = [c for c in CODON_AA_MAP.keys() if c not in ['UAA', 'UAG', 'UGA']]
ALL_CODONS.sort()

# Amino acid ordering (standard biochemical order)
AA_ORDER = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
            'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']


class CodonWAnalyzer:
    """CodonW result analyzer"""
    
    def __init__(self):
        self.all_data = {}
        self.notes = []
        
    def parse_file(self, file_path):
        """Parse a single CodonW result file"""
        species_name = Path(file_path).stem
        codon_data = {}
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        for line in lines:
            # Regex pattern to match codon data
            # Format: Phe UUU  791 1.27 or     UUC  450 0.73
            pattern = r'(?:\w{3}\s+)?([A-Z]{3})\s+(\d+)\s+([\d.]+)'
            matches = re.findall(pattern, line)
            
            for match in matches:
                codon, count, rscu = match
                if codon in CODON_AA_MAP:
                    aa = CODON_AA_MAP[codon]
                    codon_data[codon] = {
                        'AA': aa,
                        'Count': int(count),
                        'RSCU': float(rscu)
                    }
        
        return species_name, codon_data
    
    def analyze_folder(self, folder_path):
        """Analyze all CodonW results in the folder"""
        self.all_data = {}
        self.notes = []
        
        # Find all .blk files
        blk_files = list(Path(folder_path).glob('*.blk'))
        
        if not blk_files:
            raise ValueError(f"No .blk files found in folder {folder_path}")
        
        # Parse all files
        for file_path in blk_files:
            species_name, codon_data = self.parse_file(file_path)
            self.all_data[species_name] = codon_data
        
        # Check data completeness and record issues that need cleaning
        self._check_data_quality()
        
        return len(blk_files)
    
    def _check_data_quality(self):
        """Check data quality and record issues that need attention"""
        self.notes.append("=== Data Quality Report ===\n")
        self.notes.append(f"Processed CodonW results from {len(self.all_data)} species\n\n")
        
        # Check codon completeness for each species
        missing_codons = []
        for species, data in self.all_data.items():
            for codon in ALL_CODONS:
                if codon not in data:
                    missing_codons.append(f"{species} missing codon {codon}")
        
        if missing_codons:
            self.notes.append("⚠️ Missing codons:\n")
            for note in missing_codons:
                self.notes.append(f"  - {note}\n")
            self.notes.append("\n")
        
        # Check zero-count codons
        zero_codons = []
        for species, data in self.all_data.items():
            for codon in ALL_CODONS:
                if codon in data and data[codon]['Count'] == 0:
                    zero_codons.append(f"{species}: {codon} (RSCU={data[codon]['RSCU']})")
        
        if zero_codons:
            self.notes.append("⚠️ Zero-count codons (may affect statistical analysis):\n")
            for note in zero_codons[:20]:  # Show first 20 only
                self.notes.append(f"  - {note}\n")
            if len(zero_codons) > 20:
                self.notes.append(f"  ... and {len(zero_codons) - 20} more\n")
            self.notes.append("\n")
        
        # Check RSCU outliers
        high_rscu = []
        for species, data in self.all_data.items():
            for codon in ALL_CODONS:
                if codon in data and data[codon]['RSCU'] > 2.0:
                    high_rscu.append(f"{species}: {codon} (RSCU={data[codon]['RSCU']})")
        
        if high_rscu:
            self.notes.append("⚠️ RSCU values > 2.0 (potential outliers):\n")
            for note in high_rscu[:10]:
                self.notes.append(f"  - {note}\n")
            if len(high_rscu) > 10:
                self.notes.append(f"  ... and {len(high_rscu) - 10} more\n")
            self.notes.append("\n")
        
        self.notes.append("=== Cleaning Recommendations ===\n")
        self.notes.append("1. Data with missing codons may need manual completion or exclusion\n")
        self.notes.append("2. Zero-count codons may cause bias in heatmap analysis; consider whether to keep them\n")
        self.notes.append("3. RSCU outliers should be checked in the original data or analysis workflow\n")
        self.notes.append("4. For specific codon preference analysis, filter to the data of interest\n")
    
    def create_heatmap_table(self, output_file, metric='RSCU', sort_by_aa=False):
        """Create a table suitable for heatmap visualization"""
        if metric not in ['RSCU', 'Count']:
            raise ValueError("metric must be 'RSCU' or 'Count'")
        
        # Create DataFrame with species as rows and codons as columns
        species_list = sorted(self.all_data.keys())
        df_data = []
        
        for species in species_list:
            row = {'Species': species}
            codon_data = self.all_data[species]
            
            for codon in ALL_CODONS:
                aa = CODON_AA_MAP[codon]
                column_name = f"{codon}({aa})"
                if codon in codon_data:
                    row[column_name] = codon_data[codon][metric]
                else:
                    row[column_name] = None  # Missing value
            
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        
        # Sort columns by amino acid
        if sort_by_aa:
            # Extract amino acid from column names and sort by AA_ORDER
            columns = ['Species']
            for aa in AA_ORDER:
                aa_cols = [col for col in df.columns if col.endswith(f'({aa})')]
                aa_cols.sort()  # Sort codons within the same amino acid alphabetically
                columns.extend(aa_cols)
            df = df[columns]
        
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        
        return df
    
    def get_notes(self):
        """Get the data quality report"""
        return ''.join(self.notes)


class MainApp:
    """Main application interface"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("CodonW Codon Usage Bias Analyzer")
        self.root.geometry("800x700")
        
        self.analyzer = CodonWAnalyzer()
        
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup user interface"""
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        main_frame.rowconfigure(3, weight=1)
        
        # Title
        title_label = ttk.Label(main_frame, text="CodonW Codon Usage Bias Analyzer",
                               font=('Arial', 14, 'bold'))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # Input folder selection
        ttk.Label(main_frame, text="Input Folder:").grid(row=1, column=0, sticky=tk.W)
        self.input_entry = ttk.Entry(main_frame, width=50)
        self.input_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(main_frame, text="Browse...", command=self._browse_input).grid(row=1, column=2, padx=5)
        
        # Output folder selection
        ttk.Label(main_frame, text="Output Folder:").grid(row=2, column=0, sticky=tk.W)
        self.output_entry = ttk.Entry(main_frame, width=50)
        self.output_entry.grid(row=2, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(main_frame, text="Browse...", command=self._browse_output).grid(row=2, column=2, padx=5)
        
        # Metric selection
        metric_frame = ttk.Frame(main_frame)
        metric_frame.grid(row=3, column=0, columnspan=3, pady=10, sticky=tk.W)
        
        ttk.Label(metric_frame, text="Metric:").pack(side=tk.LEFT)
        self.metric_var = tk.StringVar(value='RSCU')
        ttk.Radiobutton(metric_frame, text="RSCU (Relative Synonymous Codon Usage)",
                        variable=self.metric_var, value='RSCU').pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(metric_frame, text="Count (Codon Count)",
                        variable=self.metric_var, value='Count').pack(side=tk.LEFT, padx=10)
        
        # Sort option
        sort_frame = ttk.Frame(main_frame)
        sort_frame.grid(row=4, column=0, columnspan=3, pady=5, sticky=tk.W)
        
        self.sort_by_aa_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(sort_frame, text="Sort codon columns by amino acid",
                        variable=self.sort_by_aa_var).pack(side=tk.LEFT)
        
        # Analysis buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=5, column=0, columnspan=3, pady=20)
        
        ttk.Button(button_frame, text="Start Analysis", command=self._analyze,
                  style='Accent.TButton').pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear", command=self._clear).pack(side=tk.LEFT, padx=5)
        
        # Status label
        self.status_label = ttk.Label(main_frame, text="Waiting for input...", foreground='blue')
        self.status_label.grid(row=6, column=0, columnspan=3, pady=5)
        
        # Results preview
        ttk.Label(main_frame, text="Data Quality Report:").grid(row=7, column=0, sticky=tk.W)
        self.text_area = scrolledtext.ScrolledText(main_frame, height=15, wrap=tk.WORD)
        self.text_area.grid(row=8, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        main_frame.rowconfigure(8, weight=1)
    
    def _browse_input(self):
        """Browse for input folder"""
        folder = filedialog.askdirectory(title="Select folder containing CodonW results")
        if folder:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, folder)
    
    def _browse_output(self):
        """Browse for output folder"""
        folder = filedialog.askdirectory(title="Select output folder")
        if folder:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, folder)
    
    def _analyze(self):
        """Execute analysis"""
        input_folder = self.input_entry.get()
        output_folder = self.output_entry.get()
        metric = self.metric_var.get()
        sort_by_aa = self.sort_by_aa_var.get()
        
        # Validate input
        if not input_folder:
            messagebox.showerror("Error", "Please select an input folder!")
            return
        
        if not os.path.exists(input_folder):
            messagebox.showerror("Error", "Input folder does not exist!")
            return
        
        if not output_folder:
            output_folder = input_folder
            self.output_entry.insert(0, output_folder)
        
        try:
            self.status_label.config(text="Analyzing...", foreground='orange')
            self.root.update()
            
            # Analyze data
            file_count = self.analyzer.analyze_folder(input_folder)
            
            # Generate heatmap table
            suffix = f"_sorted_{metric.lower()}.csv" if sort_by_aa else f"_{metric.lower()}.csv"
            output_file = os.path.join(output_folder, f"codonw_heatmap{suffix}")
            df = self.analyzer.create_heatmap_table(output_file, metric, sort_by_aa)
            
            # Generate data quality report
            notes_file = os.path.join(output_folder, "codonw_data_quality_report.txt")
            with open(notes_file, 'w', encoding='utf-8') as f:
                f.write(self.analyzer.get_notes())
            
            # Display results
            self.text_area.delete(1.0, tk.END)
            self.text_area.insert(tk.END, self.analyzer.get_notes())
            
            sort_info = " (sorted by amino acid)" if sort_by_aa else ""
            self.status_label.config(
                text=f"Analysis complete! Processed {file_count} files, results saved to {output_folder}",
                foreground='green'
            )
            
            messagebox.showinfo(
                "Success",
                f"Analysis complete{sort_info}!\n\n"
                f"Files processed: {file_count}\n"
                f"Heatmap table: {output_file}\n"
                f"Data quality report: {notes_file}"
            )
            
        except Exception as e:
            self.status_label.config(text="Analysis failed!", foreground='red')
            messagebox.showerror("Error", f"An error occurred during analysis:\n{str(e)}")
    
    def _clear(self):
        """Clear the interface"""
        self.input_entry.delete(0, tk.END)
        self.output_entry.delete(0, tk.END)
        self.text_area.delete(1.0, tk.END)
        self.status_label.config(text="Waiting for input...", foreground='blue')


def main():
    """Main function"""
    root = tk.Tk()
    
    # Configure style
    style = ttk.Style()
    style.configure('Accent.TButton', font=('Arial', 10, 'bold'))
    
    app = MainApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
