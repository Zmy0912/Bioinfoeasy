import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import stats
import numpy as np

class GCAnalysisPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("GC12-GC3 Distribution Plotter")
        self.root.geometry("1200x800")
        
        self.csv_file = None
        self.df = None
        self.style_window = None
        
        # Style configuration variables
        self.scatter_color = tk.StringVar(value="#1f77b4")
        self.scatter_size = tk.IntVar(value=50)
        self.scatter_alpha = tk.DoubleVar(value=0.6)
        self.neutral_line_color = tk.StringVar(value="red")
        self.neutral_line_style = tk.StringVar(value="--")
        self.neutral_line_width = tk.IntVar(value=2)
        self.regression_line_color = tk.StringVar(value="blue")
        self.regression_line_style = tk.StringVar(value="--")
        self.regression_line_width = tk.IntVar(value=2)
        self.annotation_fontsize = tk.IntVar(value=15)
        self.annotation_color = tk.StringVar(value="black")
        self.use_custom_title = tk.BooleanVar(value=False)
        self.custom_title = tk.StringVar(value="")
        self.title_fontsize = tk.IntVar(value=16)
        self.title_color = tk.StringVar(value="black")
        self.axis_tick_fontsize = tk.IntVar(value=12)
        
        self.create_widgets()
        self.create_style_panel()
        
    def create_widgets(self):
        # Create main scrollable container
        self.canvas_scroll = tk.Canvas(self.root)
        self.scrollbar_y = tk.Scrollbar(self.root, orient="vertical", command=self.canvas_scroll.yview)
        self.scrollbar_x = tk.Scrollbar(self.root, orient="horizontal", command=self.canvas_scroll.xview)
        self.scrollable_frame = tk.Frame(self.canvas_scroll)
        
        # Configure scrolling
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas_scroll.configure(scrollregion=self.canvas_scroll.bbox("all"))
        )
        
        self.canvas_scroll.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas_scroll.configure(yscrollcommand=self.scrollbar_y.set)
        self.canvas_scroll.configure(xscrollcommand=self.scrollbar_x.set)
        
        # Pack scrollable container
        self.canvas_scroll.pack(side="left", fill="both", expand=True)
        self.scrollbar_y.pack(side="right", fill="y")
        self.scrollbar_x.pack(side="bottom", fill="x")
        
        # Control panel in scrollable frame
        control_frame = tk.Frame(self.scrollable_frame, padx=10, pady=10, relief=tk.RIDGE, borderwidth=2)
        control_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # File name display area
        file_frame = tk.Frame(control_frame, relief=tk.GROOVE, borderwidth=1)
        file_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(file_frame, text="Current File:", font=("Arial", 10, "bold")).pack(side=tk.LEFT, padx=5)
        self.file_label = tk.Label(file_frame, text="No file selected", 
                                  font=("Arial", 10), fg="gray", anchor="w")
        self.file_label.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        
        # Create row 1 for basic controls
        row1 = tk.Frame(control_frame)
        row1.pack(fill=tk.X, pady=5)
        
        # Select file button
        self.select_btn = tk.Button(row1, text="Select CSV File", 
                                   command=self.select_file, 
                                   bg="#4CAF50", fg="white", 
                                   font=("Arial", 12, "bold"),
                                   padx=20)
        self.select_btn.pack(side=tk.LEFT, padx=5)
        
        # Species name input
        tk.Label(row1, text="Species Name:", font=("Arial", 11)).pack(side=tk.LEFT, padx=5)
        self.species_name = tk.StringVar(value="NC01")
        self.species_entry = tk.Entry(row1, textvariable=self.species_name, 
                                      font=("Arial", 11), width=20)
        self.species_entry.pack(side=tk.LEFT, padx=5)
        
        # Custom title input
        tk.Checkbutton(row1, text="Custom Title:", variable=self.use_custom_title,
                    font=("Arial", 11)).pack(side=tk.LEFT, padx=5)
        self.custom_title_entry = tk.Entry(row1, textvariable=self.custom_title, 
                                        font=("Arial", 11), width=25)
        self.custom_title_entry.pack(side=tk.LEFT, padx=5)
        
        # Create row 2 for action buttons
        row2 = tk.Frame(control_frame)
        row2.pack(fill=tk.X, pady=5)
        
        # Plot button
        self.plot_btn = tk.Button(row2, text="Plot Chart", 
                                 command=self.plot_chart,
                                 bg="#2196F3", fg="white",
                                 font=("Arial", 12, "bold"),
                                 padx=20)
        self.plot_btn.pack(side=tk.LEFT, padx=5)
        
        # Save button
        self.save_btn = tk.Button(row2, text="Save Plot", 
                                 command=self.save_plot,
                                 bg="#FF9800", fg="white",
                                 font=("Arial", 12, "bold"),
                                 padx=20)
        self.save_btn.pack(side=tk.LEFT, padx=5)
        
        # Clear button
        self.clear_btn = tk.Button(row2, text="Clear", 
                                  command=self.clear_plot,
                                  bg="#f44336", fg="white",
                                  font=("Arial", 12, "bold"),
                                  padx=20)
        self.clear_btn.pack(side=tk.LEFT, padx=5)
        
        # Style settings button
        style_btn = tk.Button(row2, text="Style Settings", 
                             command=self.show_style_window,
                             bg="#9C27B0", fg="white",
                             font=("Arial", 12, "bold"),
                             padx=20)
        style_btn.pack(side=tk.LEFT, padx=5)
        
        # Status label
        self.status_label = tk.Label(row2, text="Please select a CSV file", 
                                    font=("Arial", 10), fg="gray")
        self.status_label.pack(side=tk.LEFT, padx=20)
        
        # Plot display area - use a fixed size container
        plot_frame_wrapper = tk.Frame(self.scrollable_frame, relief=tk.RIDGE, borderwidth=2)
        plot_frame_wrapper.pack(padx=10, pady=10)
        
        self.plot_container = tk.Frame(plot_frame_wrapper, width=850, height=850, bg='lightgray')
        self.plot_container.pack_propagate(False)  # Prevent container from resizing
        self.plot_container.pack(padx=10, pady=10)
        
        self.canvas = None
    
    def create_style_panel(self):
        # Create style settings window
        self.style_window = tk.Toplevel(self.root)
        self.style_window.title("Style Settings")
        self.style_window.geometry("400x680")
        self.style_window.withdraw()  # Initially hidden
        
        # Style settings panel
        style_frame = tk.Frame(self.style_window, padx=10, pady=10)
        style_frame.pack(fill=tk.BOTH, expand=True)
        
        # Scatter points settings
        tk.Label(style_frame, text="Scatter Points Settings", font=("Arial", 12, "bold")).pack(pady=5)
        
        scatter_frame = tk.Frame(style_frame)
        scatter_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(scatter_frame, text="Color:").grid(row=0, column=0, sticky=tk.W)
        tk.Entry(scatter_frame, textvariable=self.scatter_color, width=10).grid(row=0, column=1, padx=5)
        
        tk.Label(scatter_frame, text="Size:").grid(row=1, column=0, sticky=tk.W)
        tk.Scale(scatter_frame, from_=10, to=200, orient=tk.HORIZONTAL, 
                variable=self.scatter_size).grid(row=1, column=1, padx=5, sticky=tk.W+tk.E)
        
        tk.Label(scatter_frame, text="Alpha:").grid(row=2, column=0, sticky=tk.W)
        tk.Scale(scatter_frame, from_=0.1, to=1.0, resolution=0.1, orient=tk.HORIZONTAL, 
                variable=self.scatter_alpha).grid(row=2, column=1, padx=5, sticky=tk.W+tk.E)
        
        # Neutral expectation line settings
        tk.Label(style_frame, text="Neutral Expectation Line", font=("Arial", 12, "bold")).pack(pady=5)
        
        neutral_frame = tk.Frame(style_frame)
        neutral_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(neutral_frame, text="Color:").grid(row=0, column=0, sticky=tk.W)
        tk.Entry(neutral_frame, textvariable=self.neutral_line_color, width=10).grid(row=0, column=1, padx=5)
        
        tk.Label(neutral_frame, text="Style:").grid(row=1, column=0, sticky=tk.W)
        style_combo = ttk.Combobox(neutral_frame, textvariable=self.neutral_line_style, 
                                  values=["-", "--", ":", "-."], width=8)
        style_combo.grid(row=1, column=1, padx=5, sticky=tk.W)
        
        tk.Label(neutral_frame, text="Width:").grid(row=2, column=0, sticky=tk.W)
        tk.Scale(neutral_frame, from_=1, to=5, orient=tk.HORIZONTAL, 
                variable=self.neutral_line_width).grid(row=2, column=1, padx=5, sticky=tk.W+tk.E)
        
        # Regression line settings
        tk.Label(style_frame, text="Regression Line", font=("Arial", 12, "bold")).pack(pady=5)
        
        regression_frame = tk.Frame(style_frame)
        regression_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(regression_frame, text="Color:").grid(row=0, column=0, sticky=tk.W)
        tk.Entry(regression_frame, textvariable=self.regression_line_color, width=10).grid(row=0, column=1, padx=5)
        
        tk.Label(regression_frame, text="Style:").grid(row=1, column=0, sticky=tk.W)
        style_combo2 = ttk.Combobox(regression_frame, textvariable=self.regression_line_style, 
                                   values=["-", "--", ":", "-."], width=8)
        style_combo2.grid(row=1, column=1, padx=5, sticky=tk.W)
        
        tk.Label(regression_frame, text="Width:").grid(row=2, column=0, sticky=tk.W)
        tk.Scale(regression_frame, from_=1, to=5, orient=tk.HORIZONTAL, 
                variable=self.regression_line_width).grid(row=2, column=1, padx=5, sticky=tk.W+tk.E)
        
        # Annotation settings
        tk.Label(style_frame, text="Annotation Settings", font=("Arial", 12, "bold")).pack(pady=5)
        
        annotation_frame = tk.Frame(style_frame)
        annotation_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(annotation_frame, text="Font Size:").grid(row=0, column=0, sticky=tk.W)
        tk.Scale(annotation_frame, from_=10, to=25, orient=tk.HORIZONTAL, 
                variable=self.annotation_fontsize).grid(row=0, column=1, padx=5, sticky=tk.W+tk.E)
        
        tk.Label(annotation_frame, text="Color:").grid(row=1, column=0, sticky=tk.W)
        tk.Entry(annotation_frame, textvariable=self.annotation_color, width=10).grid(row=1, column=1, padx=5)
        
        # Title settings
        tk.Label(style_frame, text="Title Settings", font=("Arial", 12, "bold")).pack(pady=5)
        
        title_frame = tk.Frame(style_frame)
        title_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(title_frame, text="Font Size:").grid(row=0, column=0, sticky=tk.W)
        tk.Scale(title_frame, from_=10, to=30, orient=tk.HORIZONTAL, 
                variable=self.title_fontsize).grid(row=0, column=1, padx=5, sticky=tk.W+tk.E)
        
        tk.Label(title_frame, text="Color:").grid(row=1, column=0, sticky=tk.W)
        tk.Entry(title_frame, textvariable=self.title_color, width=10).grid(row=1, column=1, padx=5)
        
        # Axis settings
        tk.Label(style_frame, text="Axis Settings", font=("Arial", 12, "bold")).pack(pady=5)
        
        axis_frame = tk.Frame(style_frame)
        axis_frame.pack(fill=tk.X, pady=5)
        
        tk.Label(axis_frame, text="Tick Label Font Size:").grid(row=0, column=0, sticky=tk.W)
        tk.Scale(axis_frame, from_=8, to=20, orient=tk.HORIZONTAL, 
                variable=self.axis_tick_fontsize).grid(row=0, column=1, padx=5, sticky=tk.W+tk.E)
        
    def show_style_window(self):
        self.style_window.deiconify()
        
    def select_file(self):
        file_path = filedialog.askopenfilename(
            title="Select CSV File",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if file_path:
            self.csv_file = file_path
            try:
                self.df = pd.read_csv(file_path)
                file_name = file_path.split('\\')[-1] if '\\' in file_path else file_path.split('/')[-1]
                self.file_label.config(text=file_name, fg="blue")
                self.status_label.config(text=f"Loaded: {file_name}", fg="green")
                messagebox.showinfo("Success", f"Successfully loaded file: {file_name}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {str(e)}")
                self.file_label.config(text="Load failed", fg="red")
                self.status_label.config(text="Load failed", fg="red")
    
    def plot_chart(self):
        if self.df is None:
            messagebox.showwarning("Warning", "Please select a CSV file first!")
            return
        
        try:
            # Get GC12 and GC3 data
            gc3 = self.df['GC3 (%)'].values
            gc12 = self.df['GC12 (%)'].values
            
            # Calculate linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(gc3, gc12)
            r_squared = r_value ** 2
            
            # Create plot with fixed square size
            fig, ax = plt.subplots(figsize=(8, 8))
            
            # Draw scatter plot (GC3 on x-axis, GC12 on y-axis) - using style variables
            ax.scatter(gc3, gc12, 
                      alpha=self.scatter_alpha.get(), 
                      c=self.scatter_color.get(), 
                      s=self.scatter_size.get(), 
                      edgecolors='black', linewidth=0.5)
            
            # Draw neutral expectation line - using style variables
            ax.plot([0, 100], [0, 100], 
                   color=self.neutral_line_color.get(), 
                   linestyle=self.neutral_line_style.get(), 
                   linewidth=self.neutral_line_width.get(), 
                   label='Neutral Expectation', alpha=0.8)
            
            # Draw linear regression line - using style variables
            regression_x = np.array([0, 100])
            regression_y = slope * regression_x + intercept
            ax.plot(regression_x, regression_y, 
                   color=self.regression_line_color.get(), 
                   linestyle=self.regression_line_style.get(), 
                   linewidth=self.regression_line_width.get(), 
                   label='Linear Regression', alpha=0.8)
            
            # Set axis labels
            ax.set_xlabel('GC3 (%)', fontsize=14, fontweight='bold')
            ax.set_ylabel('GC12 (%)', fontsize=14, fontweight='bold')
            
            # Set axis tick label font size
            ax.tick_params(axis='both', which='major', labelsize=self.axis_tick_fontsize.get())
            
            # Set axis range to 0-100
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 100)
            
            # Add grid
            ax.grid(True, alpha=0.3, linestyle='--')
            
            # Add legend
            ax.legend(loc='best', fontsize=11, framealpha=0.9)
            
            # Annotate R-squared and regression equation - using style variables
            regression_text = f'R² = {r_squared:.4f}\ny = {slope:.4f}x + {intercept:.4f}'
            ax.text(0.05, 0.95, regression_text,
                   transform=ax.transAxes,
                   fontsize=self.annotation_fontsize.get(),
                   verticalalignment='top',
                   fontweight='bold',
                   color=self.annotation_color.get())
            
            # Add title (custom or default) with custom styling
            if self.use_custom_title.get():
                title_text = self.custom_title.get()
            else:
                species = self.species_name.get()
                title_text = f'{species} GC12-GC3 Distribution Plot'
            ax.set_title(title_text, fontsize=self.title_fontsize.get(), 
                       fontweight='bold', pad=20, color=self.title_color.get())
            
            plt.tight_layout()
            
            # Clear old plot
            if self.canvas:
                self.canvas.get_tk_widget().destroy()
            
            # Display plot in fixed size container
            self.canvas = FigureCanvasTkAgg(fig, master=self.plot_container)
            self.canvas.draw()
            # Center the canvas in the container
            self.canvas.get_tk_widget().place(relx=0.5, rely=0.5, anchor='center')
            
            self.status_label.config(text="Plot created successfully", fg="green")
            
        except KeyError as e:
            messagebox.showerror("Error", f"Missing required columns in CSV file: {str(e)}")
        except Exception as e:
            messagebox.showerror("Error", f"Error creating plot: {str(e)}")
            self.status_label.config(text="Plot failed", fg="red")
    
    def save_plot(self):
        if self.canvas is None:
            messagebox.showwarning("Warning", "Please create a plot first!")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save Plot",
            defaultextension=".pdf",
            filetypes=[
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("EPS files", "*.eps"),
                ("PNG files", "*.png"),
                ("All files", "*.*")
            ]
        )
        
        if file_path:
            try:
                # Determine save parameters based on file extension
                file_ext = file_path.split('.')[-1].lower()
                if file_ext in ['pdf', 'svg', 'eps']:
                    # Vector format, no dpi needed
                    self.canvas.figure.savefig(file_path, bbox_inches='tight', format=file_ext)
                else:
                    # Raster format, use high dpi
                    self.canvas.figure.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Plot saved to: {file_path}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save: {str(e)}")
    
    def clear_plot(self):
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        self.df = None
        self.csv_file = None
        self.file_label.config(text="No file selected", fg="gray")
        self.status_label.config(text="Please select a CSV file", fg="gray")


def main():
    root = tk.Tk()
    app = GCAnalysisPlotter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
