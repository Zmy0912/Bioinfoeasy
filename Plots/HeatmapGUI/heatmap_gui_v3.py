"""
Heatmap Visualization Tool - GUI Version V3
Fully rewritten version ensuring matplotlib figures display correctly
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from heatmap_visualizer import HeatmapVisualizer


class HeatmapGUI:
    """Heatmap Visualization Tool GUI Class"""

    COLORMAPS = [
        'viridis', 'plasma', 'inferno', 'magma', 'coolwarm',
        'RdBu', 'RdYlBu', 'PiYG', 'PRGn', 'Greens', 'Blues', 'Reds'
    ]

    COLORMAP_NAMES = {
        'viridis': 'viridis (Purple-Blue-Green-Yellow)',
        'plasma': 'plasma (Blue-Red-Yellow)',
        'inferno': 'inferno (Black-Red-Yellow)',
        'magma': 'magma (Black-Purple-Yellow)',
        'coolwarm': 'coolwarm (Blue-White-Red)',
        'RdBu': 'RdBu (Red-White-Blue)',
        'RdYlBu': 'RdYlBu (Red-Yellow-Blue)',
        'PiYG': 'PiYG (Purple-White-Green)',
        'PRGn': 'PRGn (Purple-White-Green)',
        'Greens': 'Greens (Green Gradient)',
        'Blues': 'Blues (Blue Gradient)',
        'Reds': 'Reds (Red Gradient)',
    }

    def __init__(self, root):
        self.root = root
        self.root.title("Heatmap Visualization Tool v3.1 - Enhanced")
        self.root.geometry("1300x900")

        # Set fonts
        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'sans-serif']
        plt.rcParams['axes.unicode_minus'] = False

        # Variables
        self.data_file = tk.StringVar()
        self.title_var = tk.StringVar(value="Heatmap Analysis")
        self.xlabel_var = tk.StringVar(value="Columns")
        self.ylabel_var = tk.StringVar(value="Rows")
        self.cmap_var = tk.StringVar(value="viridis")
        self.dpi_var = tk.IntVar(value=150)
        self.cluster_rows = tk.BooleanVar(value=True)
        self.cluster_cols = tk.BooleanVar(value=True)
        self.save_path_var = tk.StringVar(value="")

        # New: Color range
        self.vmin_var = tk.StringVar(value="")
        self.vmax_var = tk.StringVar(value="")

        # New: Dendrogram position
        self.dendro_row_width = tk.DoubleVar(value=0.12)
        self.dendro_col_height = tk.DoubleVar(value=0.15)

        # New: Color bar settings
        self.cbar_width_var = tk.DoubleVar(value=0.15)
        self.cbar_height_var = tk.StringVar(value="")

        # New: Font sizes
        self.title_size_var = tk.IntVar(value=16)
        self.label_size_var = tk.IntVar(value=12)
        self.tick_size_var = tk.IntVar(value=8)

        self.visualizer = None
        self.current_figure = None
        self.canvas_widget = None

        self._setup_ui()

    def _setup_ui(self):
        """Setup UI"""
        # Left panel - with scroll support
        left_panel = ttk.Frame(self.root, padding=10)
        left_panel.pack(side=tk.LEFT, fill=tk.Y)

        # Create scroll container
        canvas = tk.Canvas(left_panel, highlightthickness=0)
        scrollbar = ttk.Scrollbar(left_panel, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        # Layout scroll container
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Mouse wheel event
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        # Data selection
        ttk.Label(scrollable_frame, text="Data File:").pack(anchor=tk.W, pady=2)
        entry = ttk.Entry(scrollable_frame, textvariable=self.data_file, width=35)
        entry.pack(fill=tk.X, pady=2)
        ttk.Button(scrollable_frame, text="Browse...", command=self._browse_file).pack(fill=tk.X, pady=2)
        ttk.Button(scrollable_frame, text="Generate Sample Data", command=self._generate_sample).pack(fill=tk.X, pady=5)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Labels
        ttk.Label(scrollable_frame, text="Title:").pack(anchor=tk.W, pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.title_var, width=35).pack(fill=tk.X, pady=2)

        ttk.Label(scrollable_frame, text="X-axis Label:").pack(anchor=tk.W, pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.xlabel_var, width=35).pack(fill=tk.X, pady=2)

        ttk.Label(scrollable_frame, text="Y-axis Label:").pack(anchor=tk.W, pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.ylabel_var, width=35).pack(fill=tk.X, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Color scheme
        ttk.Label(scrollable_frame, text="Color Scheme:").pack(anchor=tk.W, pady=2)
        cmap_combo = ttk.Combobox(scrollable_frame, textvariable=self.cmap_var,
                                  values=[self.COLORMAP_NAMES.get(c, c) for c in self.COLORMAPS],
                                  state="readonly", width=32)
        cmap_combo.pack(fill=tk.X, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Clustering
        ttk.Checkbutton(scrollable_frame, text="Row Clustering", variable=self.cluster_rows).pack(anchor=tk.W, pady=2)
        ttk.Checkbutton(scrollable_frame, text="Column Clustering", variable=self.cluster_cols).pack(anchor=tk.W, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Color range
        ttk.Label(scrollable_frame, text="Color Range:").pack(anchor=tk.W, pady=2)
        ttk.Label(scrollable_frame, text="Minimum (leave blank for auto):").pack(anchor=tk.W, pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.vmin_var, width=10).pack(fill=tk.X, pady=2)
        ttk.Label(scrollable_frame, text="Maximum (leave blank for auto):").pack(anchor=tk.W, pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.vmax_var, width=10).pack(fill=tk.X, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Dendrogram position
        ttk.Label(scrollable_frame, text="Dendrogram Position:").pack(anchor=tk.W, pady=2)
        ttk.Label(scrollable_frame, text="Row Tree Width (0.05-0.25):").pack(anchor=tk.W, pady=2)
        ttk.Spinbox(scrollable_frame, from_=0.05, to=0.25, increment=0.01,
                    textvariable=self.dendro_row_width, width=10, format="%.2f").pack(fill=tk.X, pady=2)
        ttk.Label(scrollable_frame, text="Column Tree Height (0.05-0.25):").pack(anchor=tk.W, pady=2)
        ttk.Spinbox(scrollable_frame, from_=0.05, to=0.25, increment=0.01,
                    textvariable=self.dendro_col_height, width=10, format="%.2f").pack(fill=tk.X, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Color bar
        ttk.Label(scrollable_frame, text="Color Bar:").pack(anchor=tk.W, pady=2)
        ttk.Label(scrollable_frame, text="Width (0.05-0.3):").pack(anchor=tk.W, pady=2)
        ttk.Spinbox(scrollable_frame, from_=0.05, to=0.3, increment=0.01,
                    textvariable=self.cbar_width_var, width=10, format="%.2f").pack(fill=tk.X, pady=2)
        ttk.Label(scrollable_frame, text="Height (0.3-1.0, leave blank for auto):").pack(anchor=tk.W, pady=2)
        ttk.Entry(scrollable_frame, textvariable=self.cbar_height_var, width=10).pack(fill=tk.X, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Font sizes
        ttk.Label(scrollable_frame, text="Font Size:").pack(anchor=tk.W, pady=2)
        ttk.Label(scrollable_frame, text="Title:").pack(anchor=tk.W, pady=2)
        ttk.Spinbox(scrollable_frame, from_=10, to=32, textvariable=self.title_size_var, width=10).pack(fill=tk.X, pady=2)
        ttk.Label(scrollable_frame, text="Axis Labels:").pack(anchor=tk.W, pady=2)
        ttk.Spinbox(scrollable_frame, from_=8, to=20, textvariable=self.label_size_var, width=10).pack(fill=tk.X, pady=2)
        ttk.Label(scrollable_frame, text="Ticks:").pack(anchor=tk.W, pady=2)
        ttk.Spinbox(scrollable_frame, from_=4, to=16, textvariable=self.tick_size_var, width=10).pack(fill=tk.X, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # DPI
        ttk.Label(scrollable_frame, text="DPI (Resolution):").pack(anchor=tk.W, pady=2)
        ttk.Spinbox(scrollable_frame, from_=72, to=600, textvariable=self.dpi_var, width=10).pack(fill=tk.X, pady=2)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Buttons
        ttk.Button(scrollable_frame, text="Generate Heatmap", command=self._generate_heatmap).pack(fill=tk.X, pady=5)
        ttk.Button(scrollable_frame, text="Save Image", command=self._save_image).pack(fill=tk.X, pady=5)
        ttk.Button(scrollable_frame, text="Clear Preview", command=self._clear_preview).pack(fill=tk.X, pady=5)

        ttk.Separator(scrollable_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)

        # Status label
        self.status_label = ttk.Label(scrollable_frame, text="Ready", wraplength=300)
        self.status_label.pack(anchor=tk.W, pady=5)

        # Right canvas area
        right_panel = ttk.Frame(self.root)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.canvas_frame = ttk.Frame(right_panel)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Placeholder label
        self.placeholder_label = ttk.Label(self.canvas_frame, text="Heatmap preview will appear here",
                                           font=('Arial', 14))
        self.placeholder_label.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

    def _get_cmap_code(self):
        """Get colormap code"""
        selected_cmap = self.cmap_var.get()
        # If selected has parentheses, extract the English code
        if '(' in selected_cmap:
            cmap_code = selected_cmap.split('(')[1].rstrip(')')
        else:
            cmap_code = selected_cmap
        return cmap_code

    def _browse_file(self):
        """Browse file"""
        file_path = filedialog.askopenfilename(
            title="Select Data File",
            filetypes=[
                ("Excel Files", "*.xlsx *.xls"),
                ("CSV Files", "*.csv"),
                ("All Files", "*.*")
            ]
        )
        if file_path:
            self.data_file.set(file_path)
            self.status_label.config(text=f"Selected: {file_path}")

    def _generate_sample(self):
        """Generate sample data"""
        try:
            self.status_label.config(text="Generating sample data...")
            self.visualizer = HeatmapVisualizer(
                title=self.title_var.get(),
                xlabel=self.xlabel_var.get(),
                ylabel=self.ylabel_var.get(),
                cmap=self.cmap_var.get(),
                dpi=self.dpi_var.get(),
                cluster_rows=self.cluster_rows.get(),
                cluster_cols=self.cluster_cols.get()
            )
            self.visualizer.generate_sample_data(n_rows=15, n_cols=12)
            self.status_label.config(text="Sample data ready: 15 rows x 12 columns")
            self._generate_heatmap()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate sample data: {e}")
            self.status_label.config(text=f"Error: {e}")

    def _generate_heatmap(self):
        """Generate heatmap"""
        try:
            print("[DEBUG] Starting heatmap generation")

            # Clear old plot
            self._clear_preview()
            print("[DEBUG] Preview cleared")

            self.status_label.config(text="Generating heatmap...")

            # Parse color range
            try:
                vmin = float(self.vmin_var.get()) if self.vmin_var.get().strip() else None
                vmax = float(self.vmax_var.get()) if self.vmax_var.get().strip() else None
                print(f"[DEBUG] Color range: vmin={vmin}, vmax={vmax}")
            except ValueError:
                messagebox.showwarning("Warning", "Color range must be a number, or leave blank for auto")
                vmin = None
                vmax = None

            # Parse color bar height
            try:
                cbar_height = float(self.cbar_height_var.get()) if self.cbar_height_var.get().strip() else None
                print(f"[DEBUG] Color bar height: {cbar_height}")
            except ValueError:
                messagebox.showwarning("Warning", "Color bar height must be a number, or leave blank for auto")
                cbar_height = None

            # Get dendrogram position parameters
            dendro_row_width = self.dendro_row_width.get()
            dendro_col_height = self.dendro_col_height.get()
            print(f"[DEBUG] Dendrogram position: row width={dendro_row_width}, col height={dendro_col_height}")

            # Check data
            if self.visualizer is None or self.visualizer.data is None:
                print("[DEBUG] No existing data, loading...")
                if self.data_file.get():
                    # Load file
                    print(f"[DEBUG] Loading file: {self.data_file.get()}")
                    self.visualizer = HeatmapVisualizer(
                        data_file=self.data_file.get(),
                        title=self.title_var.get(),
                        xlabel=self.xlabel_var.get(),
                        ylabel=self.ylabel_var.get(),
                        cmap=self._get_cmap_code(),
                        dpi=self.dpi_var.get(),
                        cluster_rows=self.cluster_rows.get(),
                        cluster_cols=self.cluster_cols.get(),
                        vmin=vmin,
                        vmax=vmax,
                        title_size=self.title_size_var.get(),
                        label_size=self.label_size_var.get(),
                        tick_size=self.tick_size_var.get(),
                        cbar_width_ratio=self.cbar_width_var.get(),
                        cbar_height_ratio=cbar_height
                    )
                    # Set dendrogram position parameters
                    self.visualizer.dendro_row_width = dendro_row_width
                    self.visualizer.dendro_col_height = dendro_col_height
                    self.visualizer.load_data()
                    print(f"[DEBUG] Data loaded: {self.visualizer.data.shape}")
                else:
                    messagebox.showwarning("Warning", "Please select a data file or generate sample data")
                    return

            # Update parameters (if data exists)
            else:
                print("[DEBUG] Updating parameters with existing data...")
                visualizer = HeatmapVisualizer(
                    title=self.title_var.get(),
                    xlabel=self.xlabel_var.get(),
                    ylabel=self.ylabel_var.get(),
                    cmap=self._get_cmap_code(),
                    dpi=self.dpi_var.get(),
                    cluster_rows=self.cluster_rows.get(),
                    cluster_cols=self.cluster_cols.get(),
                    vmin=vmin,
                    vmax=vmax,
                    title_size=self.title_size_var.get(),
                    label_size=self.label_size_var.get(),
                    tick_size=self.tick_size_var.get(),
                    cbar_width_ratio=self.cbar_width_var.get(),
                    cbar_height_ratio=cbar_height
                )
                visualizer.data = self.visualizer.data
                # Set dendrogram position parameters
                visualizer.dendro_row_width = dendro_row_width
                visualizer.dendro_col_height = dendro_col_height
                self.visualizer = visualizer
                print(f"[DEBUG] Parameters updated")

            # Plot heatmap
            print("[DEBUG] Calling plot_heatmap()...")
            self.current_figure, ax = self.visualizer.plot_heatmap()
            print(f"[DEBUG] plot_heatmap() returned: {self.current_figure}")

            # Ensure figure is active
            if self.current_figure is None:
                raise Exception("Failed to create figure")

            # Draw immediately
            print("[DEBUG] Calling figure.canvas.draw()...")
            self.current_figure.canvas.draw()
            print("[DEBUG] figure.canvas.draw() completed")

            # Display in tkinter
            print("[DEBUG] Calling _show_figure()...")
            self._show_figure(self.current_figure)
            print("[DEBUG] _show_figure() completed")

            self.status_label.config(text="Heatmap generated successfully!")

        except Exception as e:
            import traceback
            error_msg = str(e)
            print(f"[ERROR] Heatmap generation failed: {error_msg}")
            traceback.print_exc()
            messagebox.showerror("Error", f"Heatmap generation failed:\n{error_msg}")
            self.status_label.config(text=f"Error: {error_msg}")

    def _show_figure(self, figure):
        """显示 matplotlib 图形到 tkinter"""
        try:
            print("[DEBUG] Starting _show_figure()")

            # 移除占位符
            self.placeholder_label.place_forget()
            print("[DEBUG] Placeholder removed")

            # 清理旧的 canvas
            if self.canvas_widget is not None:
                try:
                    self.canvas_widget.get_tk_widget().destroy()
                    print("[DEBUG] Old canvas destroyed")
                except Exception as e:
                    print(f"[DEBUG] Warning destroying old canvas: {e}")
                self.canvas_widget = None

            # 创建新的 canvas
            print("[DEBUG] Creating FigureCanvasTkAgg...")
            self.canvas_widget = FigureCanvasTkAgg(figure, master=self.canvas_frame)

            print("[DEBUG] Calling canvas.draw()...")
            self.canvas_widget.draw()

            # 获取 tkinter widget 并显示
            print("[DEBUG] Getting tkinter widget...")
            widget = self.canvas_widget.get_tk_widget()

            print("[DEBUG] Packing widget...")
            widget.pack(fill=tk.BOTH, expand=True)

            # 强制更新
            print("[DEBUG] Updating window...")
            self.root.update_idletasks()
            self.root.update()

            print("[DEBUG] Figure displayed successfully")

        except Exception as e:
            import traceback
            print(f"[ERROR] Failed to display figure: {e}")
            traceback.print_exc()
            raise Exception(f"Failed to display figure: {e}")

    def _clear_preview(self):
        """Clear preview"""
        # Destroy old canvas
        if self.canvas_widget is not None:
            try:
                self.canvas_widget.get_tk_widget().destroy()
            except:
                pass
            self.canvas_widget = None

        # Close matplotlib figure
        if self.current_figure is not None:
            try:
                plt.close(self.current_figure)
            except:
                pass
            self.current_figure = None

        # Show placeholder
        self.placeholder_label.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

        self.status_label.config(text="Preview cleared")

    def _save_image(self):
        """Save image"""
        if self.visualizer is None or self.visualizer.fig is None:
            messagebox.showwarning("Warning", "Please generate a heatmap first")
            return

        save_path = filedialog.asksaveasfilename(
            title="Save Heatmap",
            defaultextension=".png",
            filetypes=[
                ("PNG Image", "*.png"),
                ("All Files", "*.*")
            ]
        )

        if save_path:
            try:
                self.visualizer.save_heatmap(save_path)
                self.status_label.config(text=f"Saved to: {save_path}")
                messagebox.showinfo("Success", f"Image saved to:\n{save_path}")
            except Exception as e:
                messagebox.showerror("Error", f"Save failed:\n{e}")


def main():
    """Main function"""
    root = tk.Tk()
    app = HeatmapGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
