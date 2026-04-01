"""
GC Content Analyzer
A graphical tool for analyzing GC1, GC2, and GC3 content in gene sequences from FASTA files
"""

import os
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path
from Bio import SeqIO


class GCAnalyzer:
    """GC content analyzer"""
    
    def __init__(self):
        self.results = []
        self.notes = []
        
    def calculate_gc_content(self, sequence, position):
        """
        Calculate GC content at specific codon position
        
        Args:
            sequence: DNA sequence string
            position: 1, 2, or 3 (codon position)
        
        Returns:
            GC content as percentage
        """
        if len(sequence) == 0:
            return 0.0
        
        # Extract nucleotides at specific codon position
        if position == 1:
            # First codon position: indices 0, 3, 6, ...
            codon_bases = sequence[position-1::3]
        elif position == 2:
            # Second codon position: indices 1, 4, 7, ...
            codon_bases = sequence[position-1::3]
        elif position == 3:
            # Third codon position: indices 2, 5, 8, ...
            codon_bases = sequence[position-1::3]
        else:
            return 0.0
        
        if len(codon_bases) == 0:
            return 0.0
        
        # Count G and C
        gc_count = codon_bases.upper().count('G') + codon_bases.upper().count('C')
        
        # Calculate percentage
        gc_content = (gc_count / len(codon_bases)) * 100
        
        return gc_content
    
    def calculate_overall_gc(self, sequence):
        """Calculate overall GC content"""
        if len(sequence) == 0:
            return 0.0
        
        upper_seq = sequence.upper()
        gc_count = upper_seq.count('G') + upper_seq.count('C')
        
        return (gc_count / len(upper_seq)) * 100
    
    def analyze_fasta(self, fasta_file, source_file_name):
        """Analyze a FASTA file"""
        self.results = []
        self.notes = []
        
        try:
            # Read FASTA file
            for record in SeqIO.parse(fasta_file, "fasta"):
                gene_id = record.id
                description = record.description
                
                # Extract gene name (first part before special characters)
                gene_name = gene_id.split('[')[0] if '[' in gene_id else gene_id.split('|')[0] if '|' in gene_id else gene_id
                gene_name = gene_name.replace('>', '').strip()
                
                # Get sequence
                sequence = str(record.seq).upper()
                
                # Remove any non-DNA characters
                sequence = ''.join([c for c in sequence if c in 'ATGCN'])
                
                # Check if sequence length is valid (multiple of 3)
                seq_length = len(sequence)
                if seq_length % 3 != 0:
                    self.notes.append(f"⚠️ Warning: {gene_name} sequence length ({seq_length}) is not a multiple of 3\n")
                
                # Calculate GC content
                gc1 = self.calculate_gc_content(sequence, 1)
                gc2 = self.calculate_gc_content(sequence, 2)
                gc3 = self.calculate_gc_content(sequence, 3)
                gc_overall = self.calculate_overall_gc(sequence)
                gc12 = (gc1 + gc2) / 2  # Average of GC1 and GC2
                
                result = {
                    'Source File': source_file_name,
                    'Gene ID': gene_id,
                    'Gene Name': gene_name,
                    'Sequence Length': seq_length,
                    'GC1 (%)': round(gc1, 2),
                    'GC2 (%)': round(gc2, 2),
                    'GC3 (%)': round(gc3, 2),
                    'GC12 (%)': round(gc12, 2),
                    'Overall GC (%)': round(gc_overall, 2)
                }
                
                self.results.append(result)
            
            # Add summary notes
            self.notes.insert(0, f"=== Analysis Report ===\n")
            self.notes.insert(1, f"Total genes analyzed: {len(self.results)}\n")
            self.notes.insert(2, f"Source file: {source_file_name}\n\n")
            
            # Calculate statistics
            if len(self.results) > 0:
                gc1_avg = sum(r['GC1 (%)'] for r in self.results) / len(self.results)
                gc2_avg = sum(r['GC2 (%)'] for r in self.results) / len(self.results)
                gc3_avg = sum(r['GC3 (%)'] for r in self.results) / len(self.results)
                
                self.notes.insert(3, f"=== Average GC Content ===\n")
                self.notes.insert(4, f"GC1 average: {gc1_avg:.2f}%\n")
                self.notes.insert(5, f"GC2 average: {gc2_avg:.2f}%\n")
                self.notes.insert(6, f"GC3 average: {gc3_avg:.2f}%\n")
                self.notes.insert(7, f"GC3/GC1 ratio: {gc3_avg/gc1_avg:.2f}\n\n")
            
            return len(self.results)
            
        except Exception as e:
            raise Exception(f"Error parsing FASTA file: {str(e)}")
    
    def save_results(self, output_file):
        """Save results to CSV file"""
        if not self.results:
            raise Exception("No results to save")
        
        df = pd.DataFrame(self.results)
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        
        return df
    
    def get_notes(self):
        """Get analysis notes"""
        return ''.join(self.notes)


class MainApp:
    """Main application interface"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("GC Content Analyzer - GC1, GC2, GC3")
        self.root.geometry("800x700")
        
        self.analyzer = GCAnalyzer()
        
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
        
        # Title
        title_label = ttk.Label(main_frame, text="GC Content Analyzer (GC1, GC2, GC3)",
                               font=('Arial', 14, 'bold'))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # Input file selection
        ttk.Label(main_frame, text="Input FASTA File:").grid(row=1, column=0, sticky=tk.W)
        self.input_entry = ttk.Entry(main_frame, width=50)
        self.input_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(main_frame, text="Browse...", command=self._browse_input).grid(row=1, column=2, padx=5)
        
        # Output file selection
        ttk.Label(main_frame, text="Output CSV File:").grid(row=2, column=0, sticky=tk.W)
        self.output_entry = ttk.Entry(main_frame, width=50)
        self.output_entry.grid(row=2, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(main_frame, text="Browse...", command=self._browse_output).grid(row=2, column=2, padx=5)
        
        # Analysis buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=3, column=0, columnspan=3, pady=20)
        
        ttk.Button(button_frame, text="Analyze", command=self._analyze,
                  style='Accent.TButton').pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear", command=self._clear).pack(side=tk.LEFT, padx=5)
        
        # Status label
        self.status_label = ttk.Label(main_frame, text="Waiting for input...", foreground='blue')
        self.status_label.grid(row=4, column=0, columnspan=3, pady=5)
        
        # Results preview
        ttk.Label(main_frame, text="Analysis Report:").grid(row=5, column=0, sticky=tk.W)
        self.text_area = scrolledtext.ScrolledText(main_frame, height=20, wrap=tk.WORD)
        self.text_area.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        main_frame.rowconfigure(6, weight=1)
    
    def _browse_input(self):
        """Browse for input FASTA file"""
        file = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        )
        if file:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, file)
            
            # Auto-generate output file name
            input_path = Path(file)
            output_name = f"{input_path.stem}_GC_analysis.csv"
            output_path = input_path.parent / output_name
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, str(output_path))
    
    def _browse_output(self):
        """Browse for output CSV file"""
        file = filedialog.asksaveasfilename(
            title="Save analysis results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if file:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, file)
    
    def _analyze(self):
        """Execute analysis"""
        input_file = self.input_entry.get()
        output_file = self.output_entry.get()
        
        # Validate input
        if not input_file:
            messagebox.showerror("Error", "Please select an input FASTA file!")
            return
        
        if not os.path.exists(input_file):
            messagebox.showerror("Error", "Input file does not exist!")
            return
        
        if not output_file:
            input_path = Path(input_file)
            output_name = f"{input_path.stem}_GC_analysis.csv"
            output_file = str(input_path.parent / output_name)
            self.output_entry.insert(0, output_file)
        
        try:
            self.status_label.config(text="Analyzing...", foreground='orange')
            self.root.update()
            
            # Get source file name
            source_file_name = Path(input_file).name
            
            # Analyze FASTA file
            gene_count = self.analyzer.analyze_fasta(input_file, source_file_name)
            
            # Save results
            df = self.analyzer.save_results(output_file)
            
            # Display results
            self.text_area.delete(1.0, tk.END)
            self.text_area.insert(tk.END, self.analyzer.get_notes())
            
            # Add data preview
            self.text_area.insert(tk.END, f"\n=== Data Preview (First 10 rows) ===\n")
            self.text_area.insert(tk.END, df.head(10).to_string(index=False))
            
            if len(df) > 10:
                self.text_area.insert(tk.END, f"\n... and {len(df) - 10} more rows\n")
            
            self.status_label.config(
                text=f"Analysis complete! Analyzed {gene_count} genes, results saved to {output_file}",
                foreground='green'
            )
            
            messagebox.showinfo(
                "Success",
                f"Analysis complete!\n\n"
                f"Genes analyzed: {gene_count}\n"
                f"Results saved to: {output_file}"
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
