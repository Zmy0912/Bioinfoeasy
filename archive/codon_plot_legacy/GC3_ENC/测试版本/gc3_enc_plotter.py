import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import os


class GC3ENCPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("GC3-ENC Plotter")
        self.root.geometry("1100x900")
        
        # Style settings - default values
        self.style_settings = {
            'title': {
                'fontsize': 14,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xlabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'ylabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'yticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'curve': {
                'color': 'red',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'scatter': {
                'color': 'blue',
                'edgecolor': 'darkblue',
                'size': 50,
                'alpha': 0.6,
                'linewidth': 0.5,
                'visible': True
            },
            'legend': {
                'fontsize': 10,
                'color': 'black',
                'visible': True
            }
        }
        
        # Create interface
        self.create_widgets()
        
    def create_widgets(self):
        # Top control panel
        control_frame = ttk.LabelFrame(self.root, text="File Selection", padding=10)
        control_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # File path display
        ttk.Label(control_frame, text="Input file (.out):").grid(row=0, column=0, sticky=tk.W, pady=5)
        
        self.file_path = tk.StringVar()
        file_entry = ttk.Entry(control_frame, textvariable=self.file_path, width=50)
        file_entry.grid(row=0, column=1, padx=5, pady=5)
        
        # Browse button
        browse_btn = ttk.Button(control_frame, text="Browse...", command=self.browse_file)
        browse_btn.grid(row=0, column=2, padx=5, pady=5)
        
        # Species name input
        ttk.Label(control_frame, text="Species name:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.species_name = tk.StringVar(value="Species")
        name_entry = ttk.Entry(control_frame, textvariable=self.species_name, width=30)
        name_entry.grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)
        
        # Plot button
        plot_btn = ttk.Button(control_frame, text="Plot Graph", command=self.plot_graph)
        plot_btn.grid(row=1, column=2, padx=5, pady=5)
        
        # Style settings button
        style_btn = ttk.Button(control_frame, text="Style Settings", command=self.open_style_settings)
        style_btn.grid(row=1, column=3, padx=5, pady=5)
        
        # Save button
        save_btn = ttk.Button(control_frame, text="Save Image", command=self.save_image)
        save_btn.grid(row=1, column=4, padx=5, pady=5)
        
        # Plot area with fixed-size container
        plot_frame = ttk.LabelFrame(self.root, text="Graph Display", padding=10)
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create a fixed-size container frame
        self.canvas_container = tk.Frame(plot_frame, width=600, height=600, bg='white')
        self.canvas_container.pack_propagate(False)  # Prevent the container from shrinking
        self.canvas_container.pack()
        
        # Create matplotlib figure (square size, fixed)
        self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_container)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Status bar
        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(fill=tk.X, padx=5, pady=5)
        
        # Data storage
        self.gc3_data = []
        self.enc_data = []
        self.gene_names = []
        
        # Annotation settings
        self.annotation_settings = {
            'enabled': True,
            'n_genes': 5,
            'font_size': 8,
            'color': 'darkgreen'
        }
        
    def browse_file(self):
        """Browse and select file"""
        filename = filedialog.askopenfilename(
            title="Select .out file",
            filetypes=[("OUT files", "*.out"), ("All files", "*.*")]
        )
        if filename:
            self.file_path.set(filename)
            self.status_bar.config(text=f"Selected: {os.path.basename(filename)}")
            
    def load_data(self):
        """Load data from .out file"""
        filepath = self.file_path.get()
        
        if not filepath or not os.path.exists(filepath):
            messagebox.showerror("Error", "Please select a valid input file!")
            return False
        
        try:
            # Read file
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # Skip header line
            data_lines = lines[1:]
            
            self.gc3_data = []
            self.enc_data = []
            self.gene_names = []
            
            for line in data_lines:
                line = line.strip()
                if not line:
                    continue
                
                # Split data by whitespace
                parts = line.split()
                
                # Check if there are enough data columns (need at least 10 columns, indices 0-9)
                if len(parts) < 10:
                    continue
                
                try:
                    # Extract gene name (column 0), Nc (column 8, ENC) and GC3s (column 9, GC3)
                    gene_name = parts[0]
                    enc_str = parts[8]
                    gc3_str = parts[9]
                    
                    # Skip invalid data
                    if enc_str == '*****' or gc3_str == '*****':
                        continue
                    
                    enc = float(enc_str)
                    gc3 = float(gc3_str)
                    
                    self.gc3_data.append(gc3)
                    self.enc_data.append(enc)
                    self.gene_names.append(gene_name)
                    
                except (ValueError, IndexError) as e:
                    print(f"Error parsing line: {e}")
                    continue
            
            if len(self.gc3_data) == 0:
                messagebox.showerror("Error", "Could not extract valid data from file!")
                return False
            
            print(f"Loaded {len(self.gc3_data)} data points")
            print(f"GC3 range: {min(self.gc3_data):.3f} to {max(self.gc3_data):.3f}")
            print(f"ENC range: {min(self.enc_data):.2f} to {max(self.enc_data):.2f}")
            
            self.status_bar.config(text=f"Successfully loaded {len(self.gc3_data)} data points")
            return True
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file: {str(e)}")
            print(f"Error: {e}")
            return False
    
    def calculate_enc_curve(self):
        """Calculate theoretical ENC curve (Wright 1990)"""
        gc3_values = []
        enc_values = []
        
        # Generate GC3 values from 0.05 to 0.95
        for i in range(5, 96):
            gc3 = i / 100.0
            
            # Theoretical ENC calculation formula (Wright 1990)
            # ENC = 2 + GC3 + 29 / (GC3^2 + (1-GC3)^2)
            enc = 2 + gc3 + 29 / (gc3**2 + (1-gc3)**2)
            
            gc3_values.append(gc3)
            enc_values.append(enc)
        
        print(f"Theoretical curve: {len(gc3_values)} points, ENC range: {min(enc_values):.2f} to {max(enc_values):.2f}")
        
        return np.array(gc3_values), np.array(enc_values)
    
    def calculate_theoretical_enc(self, gc3):
        """Calculate theoretical ENC for a given GC3 value"""
        return 2 + gc3 + 29 / (gc3**2 + (1-gc3)**2)
    
    def find_genes_below_curve(self):
        """Find genes that deviate most below the theoretical curve"""
        if len(self.gc3_data) == 0:
            return []
        
        deviations = []
        for i, (gc3, enc) in enumerate(zip(self.gc3_data, self.enc_data)):
            theoretical_enc = self.calculate_theoretical_enc(gc3)
            # Deviation: negative means below curve, positive means above
            deviation = enc - theoretical_enc
            deviations.append((i, deviation, self.gene_names[i], gc3, enc))
        
        # Sort by deviation (ascending) to get genes below curve first
        deviations.sort(key=lambda x: x[1])
        
        # Get top n genes that are below curve (negative deviation)
        below_curve_genes = [(i, dev, name, gc3, enc) for i, dev, name, gc3, enc in deviations if dev < 0]
        
        n = self.annotation_settings['n_genes']
        return below_curve_genes[:n]
    
    def plot_graph(self):
        """Plot GC3-ENC relationship graph"""
        if not self.load_data():
            return
        
        # Clear old graph
        self.ax.clear()
        
        # Calculate and plot theoretical curve FIRST
        gc3_theory, enc_theory = self.calculate_enc_curve()
        print(f"Plotting theoretical curve with {len(gc3_theory)} points")
        if self.style_settings['curve']['visible']:
            self.ax.plot(gc3_theory, enc_theory, 
                        color=self.style_settings['curve']['color'],
                        linestyle=self.style_settings['curve']['linestyle'],
                        linewidth=self.style_settings['curve']['linewidth'],
                        zorder=1,  # Draw below data points
                        label='Theoretical ENC curve (no selection)')
        
        # Plot actual data points SECOND (on top of curve)
        print(f"Plotting {len(self.gc3_data)} data points")
        if self.style_settings['scatter']['visible']:
            self.ax.scatter(self.gc3_data, self.enc_data, 
                           c=self.style_settings['scatter']['color'],
                           alpha=self.style_settings['scatter']['alpha'],
                           s=self.style_settings['scatter']['size'],
                           edgecolors=self.style_settings['scatter']['edgecolor'],
                           linewidth=self.style_settings['scatter']['linewidth'],
                           zorder=2,  # Draw on top of curve
                           label=f'Observed data (n={len(self.gc3_data)})')
        
        # Set graph properties
        if self.style_settings['xlabel']['visible']:
            self.ax.set_xlabel('GC3 (GC content at third codon position)',
                              fontsize=self.style_settings['xlabel']['fontsize'],
                              fontweight=self.style_settings['xlabel']['fontweight'],
                              color=self.style_settings['xlabel']['color'])

        if self.style_settings['ylabel']['visible']:
            self.ax.set_ylabel('ENC (Effective Number of Codons)',
                              fontsize=self.style_settings['ylabel']['fontsize'],
                              fontweight=self.style_settings['ylabel']['fontweight'],
                              color=self.style_settings['ylabel']['color'])

        # Set tick labels size and color
        self.ax.tick_params(axis='x', labelsize=self.style_settings['xticklabels']['fontsize'],
                           labelcolor=self.style_settings['xticklabels']['color'])
        self.ax.tick_params(axis='y', labelsize=self.style_settings['yticklabels']['fontsize'],
                           labelcolor=self.style_settings['yticklabels']['color'])

        if self.style_settings['title']['visible']:
            self.ax.set_title(f'{self.species_name.get()}',
                             fontsize=self.style_settings['title']['fontsize'],
                             fontweight=self.style_settings['title']['fontweight'],
                             color=self.style_settings['title']['color'],
                             pad=20)
        
        # Set axis range - Y axis from 20 to 61 (full ENC range) + 5% top padding
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(20, 64.05)  # 61 + 5% of 61 = 61 + 3.05 = 64.05
        
        # Add grid
        self.ax.grid(True, alpha=0.3, linestyle=':')
        
        # Add legend
        if self.style_settings['legend']['visible']:
            self.ax.legend(loc='best',
                          fontsize=self.style_settings['legend']['fontsize'],
                          labelcolor=self.style_settings['legend']['color'])
        
        # Annotate genes below curve
        if self.annotation_settings['enabled'] and len(self.gc3_data) > 0:
            genes_below = self.find_genes_below_curve()
            if genes_below:
                for i, deviation, name, gc3, enc in genes_below:
                    self.ax.annotate(
                        name,
                        xy=(gc3, enc),
                        xytext=(5, -10),
                        textcoords='offset points',
                        fontsize=self.annotation_settings['font_size'],
                        color=self.annotation_settings['color'],
                        fontweight='bold',
                        arrowprops=dict(arrowstyle='->', color=self.annotation_settings['color'], lw=0.5)
                    )
        
        # Adjust layout
        self.fig.tight_layout()
        
        # Refresh canvas
        self.canvas.draw()
        
        self.status_bar.config(text=f"Graph plotted - {len(self.gc3_data)} data points")
    
    def save_image(self):
        """Save image"""
        if len(self.gc3_data) == 0:
            messagebox.showwarning("Warning", "Please plot the graph first!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Save Image",
            defaultextension=".pdf",
            filetypes=[
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("PNG files", "*.png"),
                ("EPS files", "*.eps"),
                ("All files", "*.*")
            ]
        )
        
        if filename:
            try:
                self.fig.savefig(filename, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Image saved to:\n{filename}")
                self.status_bar.config(text=f"Image saved: {os.path.basename(filename)}")
            except Exception as e:
                messagebox.showerror("Error", f"Save failed: {str(e)}")
    
    def open_style_settings(self):
        """Open style settings window"""
        style_window = tk.Toplevel(self.root)
        style_window.title("Style Settings")
        style_window.geometry("500x700")
        
        # Create notebook for tabs
        notebook = ttk.Notebook(style_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Tab 1: Title and Labels
        tab_labels = ttk.Frame(notebook)
        notebook.add(tab_labels, text="Title & Labels")

        # Title settings
        ttk.Label(tab_labels, text="Title", font=('Arial', 10, 'bold')).pack(pady=(10, 5))

        title_frame = ttk.Frame(tab_labels)
        title_frame.pack(fill=tk.X, padx=10)

        show_title_var = tk.BooleanVar(value=self.style_settings['title']['visible'])
        ttk.Checkbutton(title_frame, text="Show title",
                       variable=show_title_var,
                       command=lambda: self.toggle_visibility('title', show_title_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(title_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        title_size = ttk.Scale(title_frame, from_=8, to=24, orient=tk.HORIZONTAL)
        title_size.set(self.style_settings['title']['fontsize'])
        title_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(title_frame, text="Color",
                  command=lambda: self.choose_color('title', title_color_btn)).grid(row=2, column=0, pady=5)
        title_color_btn = ttk.Button(title_frame, text="", width=15,
                                    command=lambda: self.choose_color('title', title_color_btn))
        title_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(title_color_btn, self.style_settings['title']['color'])

        # X-axis label settings
        ttk.Label(tab_labels, text="X-axis Label", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        xlabel_frame = ttk.Frame(tab_labels)
        xlabel_frame.pack(fill=tk.X, padx=10)

        show_xlabel_var = tk.BooleanVar(value=self.style_settings['xlabel']['visible'])
        ttk.Checkbutton(xlabel_frame, text="Show X-axis label",
                       variable=show_xlabel_var,
                       command=lambda: self.toggle_visibility('xlabel', show_xlabel_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(xlabel_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        xlabel_size = ttk.Scale(xlabel_frame, from_=8, to=20, orient=tk.HORIZONTAL)
        xlabel_size.set(self.style_settings['xlabel']['fontsize'])
        xlabel_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(xlabel_frame, text="Color",
                  command=lambda: self.choose_color('xlabel', xlabel_color_btn)).grid(row=2, column=0, pady=5)
        xlabel_color_btn = ttk.Button(xlabel_frame, text="", width=15,
                                    command=lambda: self.choose_color('xlabel', xlabel_color_btn))
        xlabel_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(xlabel_color_btn, self.style_settings['xlabel']['color'])

        # X-axis tick labels settings
        ttk.Label(tab_labels, text="X-axis Tick Labels", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        xtick_frame = ttk.Frame(tab_labels)
        xtick_frame.pack(fill=tk.X, padx=10)

        ttk.Label(xtick_frame, text="Font size:").grid(row=0, column=0, sticky=tk.W, pady=5)
        xtick_size = ttk.Scale(xtick_frame, from_=6, to=18, orient=tk.HORIZONTAL)
        xtick_size.set(self.style_settings['xticklabels']['fontsize'])
        xtick_size.grid(row=0, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(xtick_frame, text="Color",
                  command=lambda: self.choose_color('xticklabels', xtick_color_btn)).grid(row=1, column=0, pady=5)
        xtick_color_btn = ttk.Button(xtick_frame, text="", width=15,
                                    command=lambda: self.choose_color('xticklabels', xtick_color_btn))
        xtick_color_btn.grid(row=1, column=1, padx=5, pady=5)
        self.update_color_button(xtick_color_btn, self.style_settings['xticklabels']['color'])

        # Y-axis label settings
        ttk.Label(tab_labels, text="Y-axis Label", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        ylabel_frame = ttk.Frame(tab_labels)
        ylabel_frame.pack(fill=tk.X, padx=10)

        show_ylabel_var = tk.BooleanVar(value=self.style_settings['ylabel']['visible'])
        ttk.Checkbutton(ylabel_frame, text="Show Y-axis label",
                       variable=show_ylabel_var,
                       command=lambda: self.toggle_visibility('ylabel', show_ylabel_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(ylabel_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ylabel_size = ttk.Scale(ylabel_frame, from_=8, to=20, orient=tk.HORIZONTAL)
        ylabel_size.set(self.style_settings['ylabel']['fontsize'])
        ylabel_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(ylabel_frame, text="Color",
                  command=lambda: self.choose_color('ylabel', ylabel_color_btn)).grid(row=2, column=0, pady=5)
        ylabel_color_btn = ttk.Button(ylabel_frame, text="", width=15,
                                    command=lambda: self.choose_color('ylabel', ylabel_color_btn))
        ylabel_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(ylabel_color_btn, self.style_settings['ylabel']['color'])

        # Y-axis tick labels settings
        ttk.Label(tab_labels, text="Y-axis Tick Labels", font=('Arial', 10, 'bold')).pack(pady=(20, 5))

        ytick_frame = ttk.Frame(tab_labels)
        ytick_frame.pack(fill=tk.X, padx=10)

        ttk.Label(ytick_frame, text="Font size:").grid(row=0, column=0, sticky=tk.W, pady=5)
        ytick_size = ttk.Scale(ytick_frame, from_=6, to=18, orient=tk.HORIZONTAL)
        ytick_size.set(self.style_settings['yticklabels']['fontsize'])
        ytick_size.grid(row=0, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(ytick_frame, text="Color",
                  command=lambda: self.choose_color('yticklabels', ytick_color_btn)).grid(row=1, column=0, pady=5)
        ytick_color_btn = ttk.Button(ytick_frame, text="", width=15,
                                    command=lambda: self.choose_color('yticklabels', ytick_color_btn))
        ytick_color_btn.grid(row=1, column=1, padx=5, pady=5)
        self.update_color_button(ytick_color_btn, self.style_settings['yticklabels']['color'])
        
        # Tab 2: Curve
        tab_curve = ttk.Frame(notebook)
        notebook.add(tab_curve, text="Theoretical Curve")
        
        ttk.Label(tab_curve, text="Theoretical ENC Curve", font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        curve_frame = ttk.Frame(tab_curve)
        curve_frame.pack(fill=tk.X, padx=10)
        
        show_curve_var = tk.BooleanVar(value=self.style_settings['curve']['visible'])
        ttk.Checkbutton(curve_frame, text="Show theoretical curve",
                       variable=show_curve_var,
                       command=lambda: self.toggle_visibility('curve', show_curve_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        ttk.Label(curve_frame, text="Line width:").grid(row=1, column=0, sticky=tk.W, pady=5)
        curve_width = ttk.Scale(curve_frame, from_=0.5, to=5, orient=tk.HORIZONTAL)
        curve_width.set(self.style_settings['curve']['linewidth'])
        curve_width.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Label(curve_frame, text="Line style:").grid(row=2, column=0, sticky=tk.W, pady=5)
        curve_style = ttk.Combobox(curve_frame, values=['-', '--', ':', '-.'], state='readonly', width=18)
        curve_style.set(self.style_settings['curve']['linestyle'])
        curve_style.grid(row=2, column=1, padx=5, pady=5)
        
        ttk.Button(curve_frame, text="Color",
                  command=lambda: self.choose_color('curve', curve_color_btn)).grid(row=3, column=0, pady=5)
        curve_color_btn = ttk.Button(curve_frame, text="", width=15,
                                   command=lambda: self.choose_color('curve', curve_color_btn))
        curve_color_btn.grid(row=3, column=1, padx=5, pady=5)
        self.update_color_button(curve_color_btn, self.style_settings['curve']['color'])
        
        # Tab 3: Scatter
        tab_scatter = ttk.Frame(notebook)
        notebook.add(tab_scatter, text="Scatter Points")
        
        ttk.Label(tab_scatter, text="Scatter Points", font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        scatter_frame = ttk.Frame(tab_scatter)
        scatter_frame.pack(fill=tk.X, padx=10)
        
        show_scatter_var = tk.BooleanVar(value=self.style_settings['scatter']['visible'])
        ttk.Checkbutton(scatter_frame, text="Show scatter points",
                       variable=show_scatter_var,
                       command=lambda: self.toggle_visibility('scatter', show_scatter_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        ttk.Label(scatter_frame, text="Point size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        scatter_size = ttk.Scale(scatter_frame, from_=10, to=200, orient=tk.HORIZONTAL)
        scatter_size.set(self.style_settings['scatter']['size'])
        scatter_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Label(scatter_frame, text="Transparency:").grid(row=2, column=0, sticky=tk.W, pady=5)
        scatter_alpha = ttk.Scale(scatter_frame, from_=0.1, to=1.0, orient=tk.HORIZONTAL)
        scatter_alpha.set(self.style_settings['scatter']['alpha'])
        scatter_alpha.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Label(scatter_frame, text="Edge width:").grid(row=3, column=0, sticky=tk.W, pady=5)
        scatter_edge = ttk.Scale(scatter_frame, from_=0, to=3, orient=tk.HORIZONTAL)
        scatter_edge.set(self.style_settings['scatter']['linewidth'])
        scatter_edge.grid(row=3, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Button(scatter_frame, text="Fill color",
                  command=lambda: self.choose_color('scatter', scatter_fill_btn)).grid(row=4, column=0, pady=5)
        scatter_fill_btn = ttk.Button(scatter_frame, text="", width=15,
                                    command=lambda: self.choose_color('scatter', scatter_fill_btn))
        scatter_fill_btn.grid(row=4, column=1, padx=5, pady=5)
        self.update_color_button(scatter_fill_btn, self.style_settings['scatter']['color'])
        
        ttk.Button(scatter_frame, text="Edge color",
                  command=lambda: self.choose_color('scatter_edge', scatter_edge_btn)).grid(row=5, column=0, pady=5)
        scatter_edge_btn = ttk.Button(scatter_frame, text="", width=15,
                                    command=lambda: self.choose_color('scatter_edge', scatter_edge_btn))
        scatter_edge_btn.grid(row=5, column=1, padx=5, pady=5)
        self.update_color_button(scatter_edge_btn, self.style_settings['scatter']['edgecolor'])
        
        # Tab 4: Legend
        tab_legend = ttk.Frame(notebook)
        notebook.add(tab_legend, text="Legend")
        
        ttk.Label(tab_legend, text="Legend", font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        legend_frame = ttk.Frame(tab_legend)
        legend_frame.pack(fill=tk.X, padx=10)
        
        show_legend_var = tk.BooleanVar(value=self.style_settings['legend']['visible'])
        ttk.Checkbutton(legend_frame, text="Show legend",
                       variable=show_legend_var,
                       command=lambda: self.toggle_visibility('legend', show_legend_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        ttk.Label(legend_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        legend_size = ttk.Scale(legend_frame, from_=6, to=16, orient=tk.HORIZONTAL)
        legend_size.set(self.style_settings['legend']['fontsize'])
        legend_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Button(legend_frame, text="Text color",
                  command=lambda: self.choose_color('legend', legend_color_btn)).grid(row=2, column=0, pady=5)
        legend_color_btn = ttk.Button(legend_frame, text="", width=15,
                                    command=lambda: self.choose_color('legend', legend_color_btn))
        legend_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(legend_color_btn, self.style_settings['legend']['color'])
        
        # Tab 5: Gene Annotation
        tab_annotation = ttk.Frame(notebook)
        notebook.add(tab_annotation, text="Gene Labels")
        
        ttk.Label(tab_annotation, text="Genes Below Curve Labels", font=('Arial', 10, 'bold')).pack(pady=(10, 5))
        
        annotation_frame = ttk.Frame(tab_annotation)
        annotation_frame.pack(fill=tk.X, padx=10)
        
        show_annotation_var = tk.BooleanVar(value=self.annotation_settings['enabled'])
        ttk.Checkbutton(annotation_frame, text="Show gene labels below curve",
                       variable=show_annotation_var,
                       command=lambda: self.toggle_annotation_visibility(show_annotation_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        ttk.Label(annotation_frame, text="Number of genes:").grid(row=1, column=0, sticky=tk.W, pady=5)
        annotation_count = ttk.Scale(annotation_frame, from_=3, to=20, orient=tk.HORIZONTAL)
        annotation_count.set(self.annotation_settings['n_genes'])
        annotation_count.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Label(annotation_frame, text="Font size:").grid(row=2, column=0, sticky=tk.W, pady=5)
        annotation_size = ttk.Scale(annotation_frame, from_=6, to=14, orient=tk.HORIZONTAL)
        annotation_size.set(self.annotation_settings['font_size'])
        annotation_size.grid(row=2, column=1, padx=5, pady=5, sticky=tk.EW)
        
        ttk.Button(annotation_frame, text="Label color",
                  command=lambda: self.choose_annotation_color(annotation_color_btn)).grid(row=3, column=0, pady=5)
        annotation_color_btn = ttk.Button(annotation_frame, text="", width=15,
                                        command=lambda: self.choose_annotation_color(annotation_color_btn))
        annotation_color_btn.grid(row=3, column=1, padx=5, pady=5)
        self.update_color_button(annotation_color_btn, self.annotation_settings['color'])
        
        # Apply and Reset buttons
        button_frame = ttk.Frame(style_window)
        button_frame.pack(fill=tk.X, padx=10, pady=10)

        apply_btn = ttk.Button(button_frame, text="Apply",
                              command=lambda: self.apply_settings({
                                  'title_size': title_size.get(),
                                  'xlabel_size': xlabel_size.get(),
                                  'ylabel_size': ylabel_size.get(),
                                  'xtick_size': xtick_size.get(),
                                  'ytick_size': ytick_size.get(),
                                  'curve_width': curve_width.get(),
                                  'curve_style': curve_style.get(),
                                  'scatter_size': scatter_size.get(),
                                  'scatter_alpha': scatter_alpha.get(),
                                  'scatter_edge': scatter_edge.get(),
                                  'legend_size': legend_size.get(),
                                  'annotation_count': annotation_count.get(),
                                  'annotation_size': annotation_size.get()
                              }, style_window))
        apply_btn.pack(side=tk.LEFT, padx=5)
        
        reset_btn = ttk.Button(button_frame, text="Reset to Default",
                              command=lambda: self.reset_settings(style_window))
        reset_btn.pack(side=tk.RIGHT, padx=5)
        
        # Store references
        self.style_refs = {
            'title_color_btn': title_color_btn,
            'xlabel_color_btn': xlabel_color_btn,
            'ylabel_color_btn': ylabel_color_btn,
            'curve_color_btn': curve_color_btn,
            'curve_style': curve_style,
            'scatter_fill_btn': scatter_fill_btn,
            'scatter_edge_btn': scatter_edge_btn,
            'legend_color_btn': legend_color_btn
        }
    
    def choose_color(self, element, button):
        """Choose color for an element"""
        color = colorchooser.askcolor(title="Choose Color")[1]
        if color:
            if element == 'curve':
                self.style_settings['curve']['color'] = color
            elif element == 'scatter':
                self.style_settings['scatter']['color'] = color
            elif element == 'scatter_edge':
                self.style_settings['scatter']['edgecolor'] = color
            elif element == 'title':
                self.style_settings['title']['color'] = color
            elif element == 'xlabel':
                self.style_settings['xlabel']['color'] = color
            elif element == 'ylabel':
                self.style_settings['ylabel']['color'] = color
            elif element == 'xticklabels':
                self.style_settings['xticklabels']['color'] = color
            elif element == 'yticklabels':
                self.style_settings['yticklabels']['color'] = color
            elif element == 'legend':
                self.style_settings['legend']['color'] = color

            self.update_color_button(button, color)
    
    def update_color_button(self, button, color):
        """Update color button appearance"""
        button.configure(text=color)
    
    def toggle_visibility(self, element, var):
        """Toggle visibility of element"""
        self.style_settings[element]['visible'] = var.get()
        print(f"Toggle {element} visibility: {self.style_settings[element]['visible']}")
    
    def apply_settings(self, settings, window):
        """Apply style settings"""
        self.style_settings['title']['fontsize'] = settings['title_size']
        self.style_settings['xlabel']['fontsize'] = settings['xlabel_size']
        self.style_settings['ylabel']['fontsize'] = settings['ylabel_size']
        self.style_settings['xticklabels']['fontsize'] = settings['xtick_size']
        self.style_settings['yticklabels']['fontsize'] = settings['ytick_size']
        self.style_settings['curve']['linewidth'] = settings['curve_width']
        self.style_settings['curve']['linestyle'] = settings['curve_style']
        self.style_settings['scatter']['size'] = settings['scatter_size']
        self.style_settings['scatter']['alpha'] = settings['scatter_alpha']
        self.style_settings['scatter']['linewidth'] = settings['scatter_edge']
        self.style_settings['legend']['fontsize'] = settings['legend_size']
        
        # Apply annotation settings
        self.annotation_settings['n_genes'] = int(settings['annotation_count'])
        self.annotation_settings['font_size'] = settings['annotation_size']

        # Replot if data exists
        if len(self.gc3_data) > 0:
            self.plot_graph()

        messagebox.showinfo("Success", "Style settings applied.")
        window.destroy()
    
    def toggle_annotation_visibility(self, var):
        """Toggle annotation visibility"""
        self.annotation_settings['enabled'] = var.get()
        if len(self.gc3_data) > 0:
            self.plot_graph()
    
    def choose_annotation_color(self, button):
        """Choose color for gene labels"""
        color = colorchooser.askcolor(title="Choose Label Color")[1]
        if color:
            self.annotation_settings['color'] = color
            self.update_color_button(button, color)
            if len(self.gc3_data) > 0:
                self.plot_graph()
    
    def reset_settings(self, window):
        """Reset to default settings"""
        self.style_settings = {
            'title': {
                'fontsize': 14,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xlabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'ylabel': {
                'fontsize': 12,
                'fontweight': 'bold',
                'color': 'black',
                'visible': True
            },
            'xticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'yticklabels': {
                'fontsize': 10,
                'color': 'black'
            },
            'curve': {
                'color': 'red',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'scatter': {
                'color': 'blue',
                'edgecolor': 'darkblue',
                'size': 50,
                'alpha': 0.6,
                'linewidth': 0.5,
                'visible': True
            },
            'legend': {
                'fontsize': 10,
                'color': 'black',
                'visible': True
            }
        }
        # Reset annotation settings
        self.annotation_settings = {
            'enabled': True,
            'n_genes': 5,
            'font_size': 8,
            'color': 'darkgreen'
        }
        messagebox.showinfo("Success", "Settings reset to default. Please click 'Plot Graph' to see changes.")
        window.destroy()


def main():
    root = tk.Tk()
    app = GC3ENCPlotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
