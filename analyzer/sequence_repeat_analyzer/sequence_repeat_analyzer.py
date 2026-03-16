#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence Repeat Analyzer
A program for identifying and analyzing scattered repeats in sequence files with a graphical interface
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText
import os
from collections import defaultdict
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# Configure matplotlib font settings
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False  # Fix minus sign display


class RepeatAnalyzer:
    """Repeat sequence analyzer"""
    
    def __init__(self):
        self.repeat_types = ['P', 'F', 'R']
    
    def parse_file(self, filepath):
        """Parse repeat sequence file"""
        data = []
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#') or not line:
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            repeat_len = int(parts[0])
                            repeat_type = parts[2]
                            pos1 = int(parts[3])
                            pos2 = int(parts[4])
                            errors = int(parts[5])
                            evalue = float(parts[6])
                            
                            data.append({
                                'length': repeat_len,
                                'type': repeat_type,
                                'pos1': pos1,
                                'pos2': pos2,
                                'errors': errors,
                                'evalue': evalue
                            })
                        except (ValueError, IndexError):
                            continue
            return data
        except Exception as e:
            print(f"Error parsing file: {e}")
            return []
    
    def analyze_directory(self, directory):
        """Analyze all text files in directory"""
        results = {}
        
        for filename in os.listdir(directory):
            filepath = os.path.join(directory, filename)
            
            if os.path.isfile(filepath) and filename.endswith('.txt'):
                data = self.parse_file(filepath)
                if data:
                    results[filename] = self._analyze_data(data)
        
        return results
    
    def _analyze_data(self, data):
        """Analyze repeat data"""
        # Group by length and count different types
        length_stats = defaultdict(lambda: defaultdict(int))
        
        for item in data:
            length = item['length']
            repeat_type = item['type']
            length_stats[length][repeat_type] += 1
        
        # Count by type
        type_stats = defaultdict(int)
        for item in data:
            type_stats[item['type']] += 1
        
        # Count by length range
        length_ranges = {
            '30-39': 0,
            '40-49': 0,
            '50-59': 0,
            '>=60': 0
        }
        
        for item in data:
            length = item['length']
            if 30 <= length <= 39:
                length_ranges['30-39'] += 1
            elif 40 <= length <= 49:
                length_ranges['40-49'] += 1
            elif 50 <= length <= 59:
                length_ranges['50-59'] += 1
            else:
                length_ranges['>=60'] += 1
        
        return {
            'total': len(data),
            'type_stats': dict(type_stats),
            'length_stats': {k: dict(v) for k, v in length_stats.items()},
            'length_ranges': length_ranges,
            'data': data
        }


class SequenceRepeatAnalyzerGUI:
    """Graphical interface main program"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Sequence Repeat Analyzer")
        self.root.geometry("1200x800")
        
        self.analyzer = RepeatAnalyzer()
        self.current_directory = ""
        self.analysis_results = {}
        
        self.setup_ui()
    
    def setup_ui(self):
        """Setup user interface"""
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=2)
        main_frame.rowconfigure(1, weight=1)
        
        # Left control panel
        control_frame = ttk.LabelFrame(main_frame, text="Control Panel", padding="10")
        control_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        # Folder selection
        ttk.Label(control_frame, text="Select Folder:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.dir_path_var = tk.StringVar()
        dir_entry = ttk.Entry(control_frame, textvariable=self.dir_path_var, width=40)
        dir_entry.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        browse_btn = ttk.Button(control_frame, text="Browse...", command=self.browse_directory)
        browse_btn.grid(row=1, column=1, padx=5)
        
        # Analyze button
        analyze_btn = ttk.Button(control_frame, text="Analyze Folder", command=self.analyze_directory)
        analyze_btn.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        
        # File list
        ttk.Label(control_frame, text="File List:").grid(row=3, column=0, sticky=tk.W, pady=5)
        
        self.file_listbox = tk.Listbox(control_frame, height=15, width=45)
        self.file_listbox.grid(row=4, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        self.file_listbox.bind('<<ListboxSelect>>', self.on_file_select)
        
        # Statistics info
        stats_frame = ttk.LabelFrame(control_frame, text="Statistics", padding="10")
        stats_frame.grid(row=5, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N), pady=10)
        stats_frame.columnconfigure(0, weight=1)
        stats_frame.rowconfigure(0, weight=1)
        
        self.stats_text = ScrolledText(stats_frame, height=10, state='disabled', wrap=tk.WORD)
        self.stats_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N))
        
        # Right chart area - with scroll support
        chart_container = ttk.LabelFrame(main_frame, text="Charts", padding="10")
        chart_container.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(5, 0))
        
        chart_container.columnconfigure(0, weight=1)
        chart_container.rowconfigure(0, weight=1)
        
        # Create canvas and scrollbars
        chart_canvas = tk.Canvas(chart_container, bg='white')
        chart_scrollbar_y = ttk.Scrollbar(chart_container, orient="vertical", command=chart_canvas.yview)
        chart_scrollbar_x = ttk.Scrollbar(chart_container, orient="horizontal", command=chart_canvas.xview)
        
        # Configure canvas scroll
        chart_canvas.configure(yscrollcommand=chart_scrollbar_y.set, xscrollcommand=chart_scrollbar_x.set)
        
        # Place canvas and scrollbars
        chart_canvas.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        chart_scrollbar_y.grid(row=0, column=1, sticky=(tk.N, tk.S))
        chart_scrollbar_x.grid(row=1, column=0, sticky=(tk.W, tk.E))
        
        # Create scrollable frame
        chart_frame = ttk.Frame(chart_canvas)
        chart_canvas.create_window((0, 0), window=chart_frame, anchor="nw")
        
        # Bind scroll event
        chart_frame.bind("<Configure>", lambda e: chart_canvas.configure(scrollregion=chart_canvas.bbox("all")))
        
        # Matplotlib charts
        self.figure = Figure(figsize=(10, 12), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, chart_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, chart_frame)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Bottom tabs
        tab_frame = ttk.Frame(main_frame)
        tab_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 0))
        
        tab_frame.columnconfigure(0, weight=1)
        tab_frame.rowconfigure(0, weight=1)
        
        self.notebook = ttk.Notebook(tab_frame)
        self.notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Detail data tab
        detail_frame = ttk.Frame(self.notebook)
        self.notebook.add(detail_frame, text="Details")
        
        detail_frame.columnconfigure(0, weight=1)
        detail_frame.rowconfigure(0, weight=1)
        
        self.detail_text = ScrolledText(detail_frame, height=15, state='disabled', wrap=tk.WORD)
        self.detail_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N))
        
        # Comparison analysis tab
        compare_frame = ttk.Frame(self.notebook)
        self.notebook.add(compare_frame, text="Comparison")
        
        compare_frame.columnconfigure(0, weight=1)
        compare_frame.rowconfigure(0, weight=1)
        
        self.compare_text = ScrolledText(compare_frame, height=15, state='disabled', wrap=tk.WORD)
        self.compare_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N))
        
        # File type statistics tab
        type_frame = ttk.Frame(self.notebook)
        self.notebook.add(type_frame, text="Type Statistics")
        
        type_frame.columnconfigure(0, weight=1)
        type_frame.rowconfigure(0, weight=1)
        
        self.type_text = ScrolledText(type_frame, height=15, state='disabled', wrap=tk.WORD)
        self.type_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N))
        
        # Initialize charts
        self.init_charts()
    
    def init_charts(self):
        """Initialize charts"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5, 'Please select a folder and analyze', 
                ha='center', va='center', fontsize=14, color='gray')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        self.canvas.draw()
    
    def browse_directory(self):
        """Browse folder"""
        directory = filedialog.askdirectory()
        if directory:
            self.dir_path_var.set(directory)
            self.update_file_list(directory)
    
    def update_file_list(self, directory):
        """Update file list"""
        self.file_listbox.delete(0, tk.END)
        
        if not os.path.exists(directory):
            return
        
        for filename in sorted(os.listdir(directory)):
            if filename.endswith('.txt'):
                self.file_listbox.insert(tk.END, filename)
    
    def analyze_directory(self):
        """Analyze folder"""
        directory = self.dir_path_var.get()
        
        if not directory or not os.path.exists(directory):
            messagebox.showwarning("Warning", "Please select a valid folder!")
            return
        
        self.current_directory = directory
        self.analysis_results = self.analyzer.analyze_directory(directory)
        
        if not self.analysis_results:
            messagebox.showwarning("Warning", "No analyzable files found!")
            return
        
        # Update file list
        self.update_file_list(directory)
        
        # Auto-select first file to display statistics
        if self.file_listbox.size() > 0:
            self.file_listbox.selection_set(0)
            self.on_file_select()
        
        # Display comparison analysis
        self.display_comparison()
        
        # Display file type statistics
        self.display_type_stats()
        
        # Update charts
        self.update_charts()
        
        messagebox.showinfo("Complete", f"Analysis complete! Analyzed {len(self.analysis_results)} files.")
    
    def on_file_select(self, event=None):
        """File selection event"""
        selection = self.file_listbox.curselection()
        
        if not selection:
            return
        
        filename = self.file_listbox.get(selection[0])
        
        if filename in self.analysis_results:
            self.display_stats(filename)
            self.display_details(filename)
    
    def display_stats(self, filename):
        """Display statistics"""
        result = self.analysis_results[filename]
        
        self.stats_text.config(state='normal')
        self.stats_text.delete(1.0, tk.END)
        
        text = f"File: {filename}\n"
        text += "=" * 40 + "\n\n"
        text += f"Total repeats: {result['total']}\n\n"
        
        text += "By type:\n"
        for rtype, count in sorted(result['type_stats'].items()):
            text += f"  {rtype}: {count}\n"
        
        text += "\nBy length range:\n"
        for range_name, count in result['length_ranges'].items():
            text += f"  {range_name}: {count}\n"
        
        self.stats_text.insert(tk.END, text)
        self.stats_text.see(tk.END)  # Scroll to end
        self.stats_text.config(state='disabled')
    
    def display_details(self, filename):
        """Display detailed data"""
        result = self.analysis_results[filename]
        
        self.detail_text.config(state='normal')
        self.detail_text.delete(1.0, tk.END)
        
        text = f"File: {filename}\n"
        text += "=" * 40 + "\n\n"
        
        text += "By length and type:\n"
        for length in sorted(result['length_stats'].keys()):
            text += f"\nLength {length}:\n"
            for rtype, count in sorted(result['length_stats'][length].items()):
                text += f"  {rtype}: {count}\n"
        
        self.detail_text.insert(tk.END, text)
        self.detail_text.see(tk.END)  # Scroll to end
        self.detail_text.config(state='disabled')
    
    def display_comparison(self):
        """Display comparison analysis"""
        self.compare_text.config(state='normal')
        self.compare_text.delete(1.0, tk.END)
        
        text = "Multi-file Comparison Analysis\n"
        text += "=" * 50 + "\n\n"
        
        # Create comparison table
        if self.analysis_results:
            # Summary comparison
            text += "Total repeats comparison:\n"
            text += "-" * 50 + "\n"
            for filename, result in self.analysis_results.items():
                text += f"{filename}: {result['total']}\n"
            
            text += "\n\nBy type comparison:\n"
            text += "-" * 50 + "\n"
            
            all_types = set()
            for result in self.analysis_results.values():
                all_types.update(result['type_stats'].keys())
            
            for rtype in sorted(all_types):
                text += f"\n{rtype} type:\n"
                for filename, result in self.analysis_results.items():
                    count = result['type_stats'].get(rtype, 0)
                    text += f"  {filename}: {count}\n"
            
            text += "\n\nBy length range comparison:\n"
            text += "-" * 50 + "\n"
            
            ranges = ['30-39', '40-49', '50-59', '>=60']
            for range_name in ranges:
                text += f"\n{range_name}:\n"
                for filename, result in self.analysis_results.items():
                    count = result['length_ranges'].get(range_name, 0)
                    text += f"  {filename}: {count}\n"
        
        self.compare_text.insert(tk.END, text)
        self.compare_text.see(tk.END)  # Scroll to end
        self.compare_text.config(state='disabled')
    
    def display_type_stats(self):
        """Display file type statistics"""
        self.type_text.config(state='normal')
        self.type_text.delete(1.0, tk.END)
        
        text = "File Type Statistics\n"
        text += "=" * 60 + "\n\n"
        
        if self.analysis_results:
            # Type description
            text += "Type Description:\n"
            text += "-" * 60 + "\n"
            text += "  P (Perfect): Perfect match (no errors)\n"
            text += "  F (Forward): Forward match\n"
            text += "  R (Reverse): Reverse match\n"
            text += "\n\n"
            
            # Detailed statistics for each file
            for filename, result in self.analysis_results.items():
                text += f"{'=' * 60}\n"
                text += f"File: {filename}\n"
                text += f"{'=' * 60}\n\n"
                
                text += f"Total repeats: {result['total']}\n\n"
                
                # By type
                text += "By type:\n"
                text += "-" * 60 + "\n"
                type_stats = result['type_stats']
                
                # Calculate percentage
                total = result['total']
                for rtype in sorted(type_stats.keys()):
                    count = type_stats[rtype]
                    percentage = (count / total * 100) if total > 0 else 0
                    text += f"  {rtype}: {count} ({percentage:.2f}%)\n"
                
                text += "\n\n"
                
                # By length range and type
                text += "By length range and type:\n"
                text += "-" * 60 + "\n"
                
                ranges = ['30-39', '40-49', '50-59', '>=60']
                for range_name in ranges:
                    range_count = result['length_ranges'].get(range_name, 0)
                    if range_count > 0:
                        text += f"\n  {range_name} bp (total {range_count}):\n"
                        
                        # Count types within this length range
                        range_type_stats = defaultdict(int)
                        for length, type_dict in result['length_stats'].items():
                            for rtype, count in type_dict.items():
                                if range_name == '30-39' and 30 <= length <= 39:
                                    range_type_stats[rtype] += count
                                elif range_name == '40-49' and 40 <= length <= 49:
                                    range_type_stats[rtype] += count
                                elif range_name == '50-59' and 50 <= length <= 59:
                                    range_type_stats[rtype] += count
                                elif range_name == '>=60' and length >= 60:
                                    range_type_stats[rtype] += count
                        
                        for rtype in sorted(range_type_stats.keys()):
                            count = range_type_stats[rtype]
                            percentage = (count / range_count * 100) if range_count > 0 else 0
                            text += f"    {rtype}: {count} ({percentage:.2f}%)\n"
                
                text += "\n\n"
                
                # By specific length
                text += "By specific length:\n"
                text += "-" * 60 + "\n"
                
                for length in sorted(result['length_stats'].keys()):
                    type_dict = result['length_stats'][length]
                    length_total = sum(type_dict.values())
                    
                    text += f"\n  Length {length} bp (total {length_total}):\n"
                    
                    for rtype in sorted(type_dict.keys()):
                        count = type_dict[rtype]
                        percentage = (count / length_total * 100) if length_total > 0 else 0
                        text += f"    {rtype}: {count} ({percentage:.2f}%)\n"
                
                text += "\n\n"
        
        self.type_text.insert(tk.END, text)
        self.type_text.see(tk.END)  # Scroll to end
        self.type_text.config(state='disabled')
    
    def update_charts(self):
        """Update charts"""
        self.figure.clear()
        
        if not self.analysis_results:
            self.init_charts()
            return
        
        # Create subplots
        gs = self.figure.add_gridspec(2, 2, hspace=0.4, wspace=0.3)
        
        # 1. Total repeats comparison bar chart
        ax1 = self.figure.add_subplot(gs[0, 0])
        filenames = list(self.analysis_results.keys())
        totals = [self.analysis_results[f]['total'] for f in filenames]
        
        bars = ax1.bar(filenames, totals, color='skyblue', alpha=0.8)
        ax1.set_title('Total Repeats Comparison', fontsize=10, fontweight='bold')
        ax1.set_ylabel('Repeats')
        ax1.tick_params(axis='x', rotation=45)
        
        # Display values on bars
        for bar, val in zip(bars, totals):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val}', ha='center', va='bottom', fontsize=8)
        
        # 2. By type comparison stacked bar chart
        ax2 = self.figure.add_subplot(gs[0, 1])
        
        all_types = sorted(set().union(*[r['type_stats'].keys() for r in self.analysis_results.values()]))
        colors = {'P': 'lightblue', 'F': 'lightcoral', 'R': 'lightgreen'}
        
        bottom = [0] * len(filenames)
        
        for rtype in all_types:
            values = [self.analysis_results[f]['type_stats'].get(rtype, 0) for f in filenames]
            ax2.bar(filenames, values, label=rtype, color=colors.get(rtype, 'gray'), bottom=bottom)
            bottom = [b + v for b, v in zip(bottom, values)]
        
        ax2.set_title('Comparison by Type', fontsize=10, fontweight='bold')
        ax2.set_ylabel('Repeats')
        ax2.legend(fontsize=8)
        ax2.tick_params(axis='x', rotation=45)
        
        # 3. By length range comparison stacked bar chart
        ax3 = self.figure.add_subplot(gs[1, 0])
        
        ranges = ['30-39', '40-49', '50-59', '>=60']
        range_colors = ['lightblue', 'lightcoral', 'lightgreen', 'orange']
        
        bottom = [0] * len(filenames)
        
        for i, range_name in enumerate(ranges):
            values = [self.analysis_results[f]['length_ranges'].get(range_name, 0) for f in filenames]
            ax3.bar(filenames, values, label=range_name, color=range_colors[i], bottom=bottom)
            bottom = [b + v for b, v in zip(bottom, values)]
        
        ax3.set_title('Comparison by Length Range', fontsize=10, fontweight='bold')
        ax3.set_ylabel('Repeats')
        ax3.legend(fontsize=8, ncol=2)
        ax3.tick_params(axis='x', rotation=45)
        
        # 4. File type distribution pie chart (first file)
        ax4 = self.figure.add_subplot(gs[1, 1])
        
        if filenames:
            first_file = filenames[0]
            type_stats = self.analysis_results[first_file]['type_stats']
            
            if type_stats:
                labels = list(type_stats.keys())
                sizes = list(type_stats.values())
                pie_colors = [colors.get(label, 'gray') for label in labels]
                
                ax4.pie(sizes, labels=labels, colors=pie_colors, autopct='%1.1f%%', startangle=90)
                ax4.set_title(f'{first_file}\nType Distribution', fontsize=10, fontweight='bold')
            else:
                ax4.text(0.5, 0.5, 'No data', ha='center', va='center')
                ax4.set_title('Type Distribution', fontsize=10)
        
        self.canvas.draw()


def main():
    """Main function"""
    root = tk.Tk()
    
    # Set icon and style
    try:
        # Try to set modern style
        style = ttk.Style()
        style.theme_use('clam')
    except:
        pass
    
    app = SequenceRepeatAnalyzerGUI(root)
    
    # Auto-detect test folder on Desktop
    test_dir = os.path.join(os.path.expanduser('~/Desktop'), 'test')
    if os.path.exists(test_dir):
        app.dir_path_var.set(test_dir)
        app.update_file_list(test_dir)
    
    root.mainloop()


if __name__ == "__main__":
    main()
