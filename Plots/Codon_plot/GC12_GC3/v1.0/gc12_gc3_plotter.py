import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import os


class GC12GC3Plotter:
    def __init__(self, root):
        self.root = root
        self.root.title("GC12-GC3 Plotter")
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
            'neutral_line': {
                'color': 'red',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'regression_line': {
                'color': 'blue',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'scatter': {
                'color': '#1f77b4',
                'edgecolor': 'black',
                'size': 50,
                'alpha': 0.6,
                'linewidth': 0.5,
                'visible': True
            },
            'annotation': {
                'fontsize': 12,
                'color': 'black',
                'visible': True
            },
            'legend': {
                'fontsize': 10,
                'color': 'black',
                'visible': True
            }
        }

        # Data storage
        self.gc3_data = []       # GC3 (%)  - X axis
        self.gc12_data = []      # GC12 (%) - Y axis
        self.gene_names = []
        self.regression_stats = {}

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
        """Load data from .out file and compute GC3 and GC12 values"""
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
            self.gc12_data = []
            self.gene_names = []

            for line in data_lines:
                line = line.strip()
                if not line:
                    continue

                # Split data by whitespace
                parts = line.split()

                # Need at least 11 columns (indices 0-10): name ... GC3s (9), GC (10)
                if len(parts) < 11:
                    continue

                try:
                    gene_name = parts[0]
                    gc3_str = parts[9]   # GC3s
                    gc_str = parts[10]   # overall GC

                    # Skip invalid data
                    if gc3_str == '*****' or gc_str == '*****':
                        continue

                    gc3 = float(gc3_str)
                    gc = float(gc_str)

                    # Accept both fractions (0-1) and percentages (0-100)
                    if gc3 > 1.0:
                        gc3 /= 100.0
                    if gc > 1.0:
                        gc /= 100.0

                    if not (0.0 <= gc3 <= 1.0 and 0.0 <= gc <= 1.0):
                        continue

                    # GC overall = (GC1 + GC2 + GC3) / 3
                    # => GC1 + GC2 = 3 * GC - GC3
                    # => GC12 = (GC1 + GC2) / 2 = (3 * GC - GC3) / 2
                    gc12 = (3.0 * gc - gc3) / 2.0

                    if not (0.0 <= gc12 <= 1.0):
                        continue

                    self.gc3_data.append(gc3 * 100.0)
                    self.gc12_data.append(gc12 * 100.0)
                    self.gene_names.append(gene_name)

                except (ValueError, IndexError) as e:
                    print(f"Error parsing line: {e}")
                    continue

            if len(self.gc3_data) == 0:
                messagebox.showerror("Error", "Could not extract valid data from file!")
                return False

            print(f"Loaded {len(self.gc3_data)} data points")
            print(f"GC3 (%) range: {min(self.gc3_data):.2f} to {max(self.gc3_data):.2f}")
            print(f"GC12 (%) range: {min(self.gc12_data):.2f} to {max(self.gc12_data):.2f}")

            self.status_bar.config(text=f"Successfully loaded {len(self.gc3_data)} data points")
            return True

        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file: {str(e)}")
            print(f"Error: {e}")
            return False

    def calculate_regression(self):
        """Linear regression of GC12 on GC3 (GC12 = slope * GC3 + intercept)"""
        x = np.array(self.gc3_data, dtype=float)
        y = np.array(self.gc12_data, dtype=float)

        if len(x) < 2 or np.std(x) == 0:
            slope, intercept, r = 0.0, float(y[0]) if len(y) else 0.0, 0.0
        else:
            slope, intercept = np.polyfit(x, y, 1)
            if np.std(y) == 0:
                r = 0.0
            else:
                r = float(np.corrcoef(x, y)[0, 1])

        self.regression_stats = {
            'slope': float(slope),
            'intercept': float(intercept),
            'r': r,
            'r_squared': r ** 2
        }
        print(f"Regression: GC12 = {slope:.4f} * GC3 + {intercept:.4f}, R^2 = {r ** 2:.4f}")
        return self.regression_stats

    def plot_graph(self):
        """Plot GC12-GC3 neutrality graph with regression"""
        if not self.load_data():
            return

        # Clear old graph
        self.ax.clear()

        stats = self.calculate_regression()
        slope = stats['slope']
        intercept = stats['intercept']
        r_squared = stats['r_squared']

        # Neutral expectation line y = x (FIRST, drawn below)
        if self.style_settings['neutral_line']['visible']:
            self.ax.plot([0, 100], [0, 100],
                         color=self.style_settings['neutral_line']['color'],
                         linestyle=self.style_settings['neutral_line']['linestyle'],
                         linewidth=self.style_settings['neutral_line']['linewidth'],
                         zorder=1,
                         label='Neutral expectation (y = x)')

        # Regression line
        if self.style_settings['regression_line']['visible']:
            reg_x = np.array([0.0, 100.0])
            reg_y = slope * reg_x + intercept
            self.ax.plot(reg_x, reg_y,
                         color=self.style_settings['regression_line']['color'],
                         linestyle=self.style_settings['regression_line']['linestyle'],
                         linewidth=self.style_settings['regression_line']['linewidth'],
                         zorder=2,
                         label='Linear regression')

        # Scatter points (on top)
        if self.style_settings['scatter']['visible']:
            self.ax.scatter(self.gc3_data, self.gc12_data,
                            c=self.style_settings['scatter']['color'],
                            alpha=self.style_settings['scatter']['alpha'],
                            s=self.style_settings['scatter']['size'],
                            edgecolors=self.style_settings['scatter']['edgecolor'],
                            linewidths=self.style_settings['scatter']['linewidth'],
                            zorder=3,
                            label=f'Observed data (n={len(self.gc3_data)})')

        # Set graph properties
        if self.style_settings['xlabel']['visible']:
            self.ax.set_xlabel('GC3 (%)',
                               fontsize=self.style_settings['xlabel']['fontsize'],
                               fontweight=self.style_settings['xlabel']['fontweight'],
                               color=self.style_settings['xlabel']['color'])

        if self.style_settings['ylabel']['visible']:
            self.ax.set_ylabel('GC12 (%)',
                               fontsize=self.style_settings['ylabel']['fontsize'],
                               fontweight=self.style_settings['ylabel']['fontweight'],
                               color=self.style_settings['ylabel']['color'])

        # Set tick labels size and color
        self.ax.tick_params(axis='x', labelsize=self.style_settings['xticklabels']['fontsize'],
                            labelcolor=self.style_settings['xticklabels']['color'])
        self.ax.tick_params(axis='y', labelsize=self.style_settings['yticklabels']['fontsize'],
                            labelcolor=self.style_settings['yticklabels']['color'])

        if self.style_settings['title']['visible']:
            self.ax.set_title(f'{self.species_name.get()} GC12-GC3 Plot',
                              fontsize=self.style_settings['title']['fontsize'],
                              fontweight=self.style_settings['title']['fontweight'],
                              color=self.style_settings['title']['color'],
                              pad=20)

        # Set axis range 0-100 (%) and keep the plot square
        self.ax.set_xlim(0, 100)
        self.ax.set_ylim(0, 100)
        self.ax.set_aspect('equal', adjustable='box')

        # Add grid
        self.ax.grid(True, alpha=0.3, linestyle=':')

        # Annotate regression equation and R^2
        if self.style_settings['annotation']['visible']:
            reg_text = f'y = {slope:.4f}x + {intercept:.4f}\nR² = {r_squared:.4f}'
            self.ax.text(0.05, 0.95, reg_text,
                         transform=self.ax.transAxes,
                         fontsize=self.style_settings['annotation']['fontsize'],
                         color=self.style_settings['annotation']['color'],
                         fontweight='bold',
                         verticalalignment='top')

        # Add legend
        if self.style_settings['legend']['visible']:
            self.ax.legend(loc='best',
                           fontsize=self.style_settings['legend']['fontsize'],
                           labelcolor=self.style_settings['legend']['color'])

        # Adjust layout
        self.fig.tight_layout()

        # Refresh canvas
        self.canvas.draw()

        self.status_bar.config(text=f"Graph plotted - {len(self.gc3_data)} data points, R² = {r_squared:.4f}")

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

        # Tab 2: Neutral Line
        tab_neutral = ttk.Frame(notebook)
        notebook.add(tab_neutral, text="Neutral Line")

        ttk.Label(tab_neutral, text="Neutral Expectation Line (y = x)", font=('Arial', 10, 'bold')).pack(pady=(10, 5))

        neutral_frame = ttk.Frame(tab_neutral)
        neutral_frame.pack(fill=tk.X, padx=10)

        show_neutral_var = tk.BooleanVar(value=self.style_settings['neutral_line']['visible'])
        ttk.Checkbutton(neutral_frame, text="Show neutral line",
                        variable=show_neutral_var,
                        command=lambda: self.toggle_visibility('neutral_line', show_neutral_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(neutral_frame, text="Line width:").grid(row=1, column=0, sticky=tk.W, pady=5)
        neutral_width = ttk.Scale(neutral_frame, from_=0.5, to=5, orient=tk.HORIZONTAL)
        neutral_width.set(self.style_settings['neutral_line']['linewidth'])
        neutral_width.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Label(neutral_frame, text="Line style:").grid(row=2, column=0, sticky=tk.W, pady=5)
        neutral_style = ttk.Combobox(neutral_frame, values=['-', '--', ':', '-.'], state='readonly', width=18)
        neutral_style.set(self.style_settings['neutral_line']['linestyle'])
        neutral_style.grid(row=2, column=1, padx=5, pady=5)

        ttk.Button(neutral_frame, text="Color",
                   command=lambda: self.choose_color('neutral_line', neutral_color_btn)).grid(row=3, column=0, pady=5)
        neutral_color_btn = ttk.Button(neutral_frame, text="", width=15,
                                       command=lambda: self.choose_color('neutral_line', neutral_color_btn))
        neutral_color_btn.grid(row=3, column=1, padx=5, pady=5)
        self.update_color_button(neutral_color_btn, self.style_settings['neutral_line']['color'])

        # Tab 3: Regression Line
        tab_reg = ttk.Frame(notebook)
        notebook.add(tab_reg, text="Regression Line")

        ttk.Label(tab_reg, text="Linear Regression Line", font=('Arial', 10, 'bold')).pack(pady=(10, 5))

        reg_frame = ttk.Frame(tab_reg)
        reg_frame.pack(fill=tk.X, padx=10)

        show_reg_var = tk.BooleanVar(value=self.style_settings['regression_line']['visible'])
        ttk.Checkbutton(reg_frame, text="Show regression line",
                        variable=show_reg_var,
                        command=lambda: self.toggle_visibility('regression_line', show_reg_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(reg_frame, text="Line width:").grid(row=1, column=0, sticky=tk.W, pady=5)
        reg_width = ttk.Scale(reg_frame, from_=0.5, to=5, orient=tk.HORIZONTAL)
        reg_width.set(self.style_settings['regression_line']['linewidth'])
        reg_width.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Label(reg_frame, text="Line style:").grid(row=2, column=0, sticky=tk.W, pady=5)
        reg_style = ttk.Combobox(reg_frame, values=['-', '--', ':', '-.'], state='readonly', width=18)
        reg_style.set(self.style_settings['regression_line']['linestyle'])
        reg_style.grid(row=2, column=1, padx=5, pady=5)

        ttk.Button(reg_frame, text="Color",
                   command=lambda: self.choose_color('regression_line', reg_color_btn)).grid(row=3, column=0, pady=5)
        reg_color_btn = ttk.Button(reg_frame, text="", width=15,
                                   command=lambda: self.choose_color('regression_line', reg_color_btn))
        reg_color_btn.grid(row=3, column=1, padx=5, pady=5)
        self.update_color_button(reg_color_btn, self.style_settings['regression_line']['color'])

        # Tab 4: Scatter
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

        # Tab 5: Annotation
        tab_annotation = ttk.Frame(notebook)
        notebook.add(tab_annotation, text="Annotation")

        ttk.Label(tab_annotation, text="Regression Statistics Text", font=('Arial', 10, 'bold')).pack(pady=(10, 5))

        annotation_frame = ttk.Frame(tab_annotation)
        annotation_frame.pack(fill=tk.X, padx=10)

        show_annotation_var = tk.BooleanVar(value=self.style_settings['annotation']['visible'])
        ttk.Checkbutton(annotation_frame, text="Show regression equation and R²",
                        variable=show_annotation_var,
                        command=lambda: self.toggle_visibility('annotation', show_annotation_var)).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)

        ttk.Label(annotation_frame, text="Font size:").grid(row=1, column=0, sticky=tk.W, pady=5)
        annotation_size = ttk.Scale(annotation_frame, from_=8, to=20, orient=tk.HORIZONTAL)
        annotation_size.set(self.style_settings['annotation']['fontsize'])
        annotation_size.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Button(annotation_frame, text="Text color",
                   command=lambda: self.choose_color('annotation', annotation_color_btn)).grid(row=2, column=0, pady=5)
        annotation_color_btn = ttk.Button(annotation_frame, text="", width=15,
                                          command=lambda: self.choose_color('annotation', annotation_color_btn))
        annotation_color_btn.grid(row=2, column=1, padx=5, pady=5)
        self.update_color_button(annotation_color_btn, self.style_settings['annotation']['color'])

        # Tab 6: Legend
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
                                   'neutral_width': neutral_width.get(),
                                   'neutral_style': neutral_style.get(),
                                   'reg_width': reg_width.get(),
                                   'reg_style': reg_style.get(),
                                   'scatter_size': scatter_size.get(),
                                   'scatter_alpha': scatter_alpha.get(),
                                   'scatter_edge': scatter_edge.get(),
                                   'annotation_size': annotation_size.get(),
                                   'legend_size': legend_size.get()
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
            'neutral_color_btn': neutral_color_btn,
            'neutral_style': neutral_style,
            'reg_color_btn': reg_color_btn,
            'reg_style': reg_style,
            'scatter_fill_btn': scatter_fill_btn,
            'scatter_edge_btn': scatter_edge_btn,
            'annotation_color_btn': annotation_color_btn,
            'legend_color_btn': legend_color_btn
        }

    def choose_color(self, element, button):
        """Choose color for an element"""
        color = colorchooser.askcolor(title="Choose Color")[1]
        if color:
            if element == 'neutral_line':
                self.style_settings['neutral_line']['color'] = color
            elif element == 'regression_line':
                self.style_settings['regression_line']['color'] = color
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
            elif element == 'annotation':
                self.style_settings['annotation']['color'] = color
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
        self.style_settings['neutral_line']['linewidth'] = settings['neutral_width']
        self.style_settings['neutral_line']['linestyle'] = settings['neutral_style']
        self.style_settings['regression_line']['linewidth'] = settings['reg_width']
        self.style_settings['regression_line']['linestyle'] = settings['reg_style']
        self.style_settings['scatter']['size'] = settings['scatter_size']
        self.style_settings['scatter']['alpha'] = settings['scatter_alpha']
        self.style_settings['scatter']['linewidth'] = settings['scatter_edge']
        self.style_settings['annotation']['fontsize'] = settings['annotation_size']
        self.style_settings['legend']['fontsize'] = settings['legend_size']

        # Replot if data exists
        if len(self.gc3_data) > 0:
            self.plot_graph()

        messagebox.showinfo("Success", "Style settings applied.")
        window.destroy()

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
            'neutral_line': {
                'color': 'red',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'regression_line': {
                'color': 'blue',
                'linestyle': '--',
                'linewidth': 2,
                'visible': True
            },
            'scatter': {
                'color': '#1f77b4',
                'edgecolor': 'black',
                'size': 50,
                'alpha': 0.6,
                'linewidth': 0.5,
                'visible': True
            },
            'annotation': {
                'fontsize': 12,
                'color': 'black',
                'visible': True
            },
            'legend': {
                'fontsize': 10,
                'color': 'black',
                'visible': True
            }
        }
        messagebox.showinfo("Success", "Settings reset to default. Please click 'Plot Graph' to see changes.")
        window.destroy()


def main():
    root = tk.Tk()
    app = GC12GC3Plotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
