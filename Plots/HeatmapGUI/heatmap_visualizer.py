"""
Professional Heatmap Visualization Tool
Supports CSV/Excel data reading, clustering analysis and multiple color schemes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from pathlib import Path
import sys


class HeatmapVisualizer:
    """Heatmap Visualization Class"""

    def __init__(self, data_file=None, title="Heatmap Analysis",
                 xlabel="Columns", ylabel="Rows", cmap="viridis",
                 dpi=300, cluster_rows=True, cluster_cols=True,
                 vmin=None, vmax=None,
                 title_size=16, label_size=12, tick_size=8,
                 cbar_width_ratio=0.15, cbar_height_ratio=None,
                 dendro_row_width=0.12, dendro_col_height=0.15):
        """
        Initialize heatmap visualizer

        Parameters:
            data_file: Data file path (CSV or Excel)
            title: Chart title
            xlabel: X-axis label
            ylabel: Y-axis label
            cmap: Color scheme ('viridis', 'coolwarm', 'plasma', 'inferno',
                  'RdYlBu', 'RdBu', 'PiYG', 'PRGn', etc.)
            dpi: Output image resolution
            cluster_rows: Whether to cluster rows
            cluster_cols: Whether to cluster columns
            vmin: Color scale minimum (None for auto)
            vmax: Color scale maximum (None for auto)
            title_size: Title font size
            label_size: Axis label font size
            tick_size: Tick label font size
            cbar_width_ratio: Color bar width ratio to main plot (0.05-0.3)
            cbar_height_ratio: Color bar height ratio to main plot (None for auto)
            dendro_row_width: Row dendrogram width ratio to heatmap (0.05-0.25)
            dendro_col_height: Column dendrogram height ratio to heatmap (0.05-0.25)
        """
        self.data_file = data_file
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cmap = cmap
        self.dpi = dpi
        self.cluster_rows = cluster_rows
        self.cluster_cols = cluster_cols
        self.vmin = vmin
        self.vmax = vmax
        self.title_size = title_size
        self.label_size = label_size
        self.tick_size = tick_size
        self.cbar_width_ratio = cbar_width_ratio
        self.cbar_height_ratio = cbar_height_ratio
        self.dendro_row_width = dendro_row_width
        self.dendro_col_height = dendro_col_height
        self.data = None
        self.fig = None
        self.ax = None

        # Set font support (using English fonts)
        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'sans-serif']
        plt.rcParams['axes.unicode_minus'] = False

    def load_data(self, data_file=None):
        """
        Load data file

        Parameters:
            data_file: Data file path (CSV or Excel)
        """
        if data_file:
            self.data_file = data_file

        if not self.data_file:
            raise ValueError("Please provide a data file path")

        file_path = Path(self.data_file)

        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {self.data_file}")

        # Read based on file extension
        if file_path.suffix.lower() == '.csv':
            # Try multiple encodings
            encodings = ['utf-8', 'gbk', 'gb18030', 'gb2312', 'utf-16', 'latin1']

            for encoding in encodings:
                try:
                    # Try reading first
                    try:
                        self.data = pd.read_csv(self.data_file, index_col=0, encoding=encoding)
                        print(f"Successfully read file using encoding {encoding}")
                        break
                    except pd.errors.ParserError as e:
                        # Try other separators
                        print(f"Encoding {encoding} failed, trying other separators...")
                        for sep in [',', ';', '\t', '|']:
                            try:
                                self.data = pd.read_csv(
                                    self.data_file,
                                    index_col=0,
                                    encoding=encoding,
                                    sep=sep,
                                    on_bad_lines='warn'  # Warn but don't skip bad lines
                                )
                                print(f"Successfully read file using encoding {encoding} and separator '{sep}'")
                                break
                            except (pd.errors.ParserError, UnicodeDecodeError, UnicodeError):
                                continue
                        else:
                            raise ValueError(
                                f"CSV file format error:\n"
                                f"- File may contain inconsistent column counts\n"
                                f"- Check data format around line 4\n"
                                f"- Ensure all rows have the same number of columns\n"
                                f"- Suggest opening in Excel and saving as .xlsx format"
                            )
                except (UnicodeDecodeError, UnicodeError):
                    continue
            else:
                raise ValueError(f"Unable to recognize file encoding, please ensure file uses UTF-8, GBK or GB18030 encoding")
        elif file_path.suffix.lower() in ['.xlsx', '.xls']:
            self.data = pd.read_excel(self.data_file, index_col=0)
        else:
            raise ValueError("Unsupported file format, please use CSV or Excel file")

        # Ensure data is numeric
        self.data = self.data.astype(float)

        print(f"Data loaded successfully: {self.data.shape[0]} rows x {self.data.shape[1]} columns")
        return self.data

    def generate_sample_data(self, n_rows=15, n_cols=12, seed=42):
        """
        Generate sample data (for testing)

        Parameters:
            n_rows: Number of rows
            n_cols: Number of columns
            seed: Random seed
        """
        np.random.seed(seed)

        # Generate patterned sample data
        base_pattern = np.linspace(0, 10, n_cols)
        self.data = pd.DataFrame(
            np.random.randn(n_rows, n_cols) * 0.5 + np.tile(base_pattern, (n_rows, 1)),
            index=[f'Sample_{i+1}' for i in range(n_rows)],
            columns=[f'Gene_{j+1}' for j in range(n_cols)]
        )

        print(f"Sample data generated successfully: {self.data.shape[0]} rows x {self.data.shape[1]} columns")
        return self.data

    def _cluster_data(self):
        """
        Perform clustering analysis and return sorted data
        """
        clustered_data = self.data.copy()

        # Row clustering
        if self.cluster_rows:
            row_linkage = linkage(clustered_data.values, method='average')
            row_order = dendrogram(row_linkage, no_plot=True)['leaves']
            clustered_data = clustered_data.iloc[row_order, :]

        # Column clustering
        if self.cluster_cols:
            col_linkage = linkage(clustered_data.values.T, method='average')
            col_order = dendrogram(col_linkage, no_plot=True)['leaves']
            clustered_data = clustered_data.iloc[:, col_order]

        return clustered_data

    def plot_heatmap(self, save_path=None):
        """
        Plot heatmap

        Parameters:
            save_path: Save path (PNG format), if not specified then not saved
        """
        if self.data is None:
            raise ValueError("Please load data or generate sample data first")

        # Perform clustering
        clustered_data = self._cluster_data()

        # Calculate figure dimensions
        n_rows, n_cols = clustered_data.shape
        fig_width = max(8, n_cols * 0.6 + 2)
        fig_height = max(6, n_rows * 0.6 + 1.5)

        # 创建图形
        self.fig = plt.figure(figsize=(fig_width, fig_height))

        # Calculate main area ratio
        main_area_ratio = 1 - self.cbar_width_ratio
        main_area_ratio = max(0.7, min(0.95, main_area_ratio))

        # Set layout based on clustering
        if self.cluster_rows and self.cluster_cols:
            # Dual clustering: title at the top
            gs = self.fig.add_gridspec(
                nrows=3, ncols=3,
                width_ratios=[self.dendro_row_width, main_area_ratio, self.cbar_width_ratio],
                height_ratios=[0.08, self.dendro_col_height, 0.92 - self.dendro_col_height],
                wspace=0.01, hspace=0.01,
                left=0.02, right=0.98, top=0.98, bottom=0.02
            )
            # Row 1: Title space (no subplot, use suptitle)
            # Row 2: Column dendrogram
            ax_dendro_col = self.fig.add_subplot(gs[1, 1])
            # Row 3: Row dendrogram and heatmap
            ax_dendro_row = self.fig.add_subplot(gs[2, 0])
            ax_heatmap = self.fig.add_subplot(gs[2, 1])
            # Color bar on right side of rows 2-3
            cax = self.fig.add_subplot(gs[1:3, 2])

        elif self.cluster_rows:
            # Row clustering only
            gs = self.fig.add_gridspec(
                nrows=1, ncols=3,
                width_ratios=[self.dendro_row_width, main_area_ratio, self.cbar_width_ratio],
                wspace=0.01,
                left=0.02, right=0.98, top=0.90, bottom=0.05
            )
            ax_dendro_row = self.fig.add_subplot(gs[0, 0])
            ax_heatmap = self.fig.add_subplot(gs[0, 1])
            cax = self.fig.add_subplot(gs[0, 2])
            ax_dendro_col = None

        elif self.cluster_cols:
            # Column clustering only: title at the top
            gs = self.fig.add_gridspec(
                nrows=3, ncols=2,
                width_ratios=[main_area_ratio, self.cbar_width_ratio],
                height_ratios=[0.08, self.dendro_col_height, 0.92 - self.dendro_col_height],
                wspace=0.01, hspace=0.01,
                left=0.02, right=0.98, top=0.98, bottom=0.02
            )
            # Row 1: Title space
            # Row 2: Column dendrogram
            ax_dendro_col = self.fig.add_subplot(gs[1, 0])
            # Row 3: Heatmap
            ax_heatmap = self.fig.add_subplot(gs[2, 0])
            # Color bar on right side of rows 2-3
            cax = self.fig.add_subplot(gs[1:3, 1])
            ax_dendro_row = None

        else:
            # No clustering
            gs = self.fig.add_gridspec(
                nrows=1, ncols=2,
                width_ratios=[main_area_ratio, self.cbar_width_ratio],
                wspace=0.02,
                left=0.05, right=0.95, top=0.90, bottom=0.1
            )
            ax_heatmap = self.fig.add_subplot(gs[0, 0])
            cax = self.fig.add_subplot(gs[0, 1])
            ax_dendro_row = None
            ax_dendro_col = None

        # Draw column dendrogram
        if self.cluster_cols and ax_dendro_col is not None:
            col_linkage = linkage(clustered_data.values.T, method='average')
            dendrogram(col_linkage, ax=ax_dendro_col, orientation='top')
            ax_dendro_col.set_xticks([])
            ax_dendro_col.set_yticks([])
            for spine in ax_dendro_col.spines.values():
                spine.set_visible(False)

        # Draw row dendrogram
        if self.cluster_rows and ax_dendro_row is not None:
            row_linkage = linkage(clustered_data.values, method='average')
            dendrogram(row_linkage, ax=ax_dendro_row, orientation='left')
            ax_dendro_row.set_xticks([])
            ax_dendro_row.set_yticks([])
            for spine in ax_dendro_row.spines.values():
                spine.set_visible(False)

        # Draw heatmap
        sns.heatmap(
            clustered_data,
            ax=ax_heatmap,
            cmap=self.cmap,
            cbar_ax=cax,
            annot=False,
            fmt='.2f',
            linewidths=0.5,
            square=True,
            xticklabels=True,
            yticklabels=True,
            vmin=self.vmin,
            vmax=self.vmax
        )

        # Adjust color bar height (if specified)
        if self.cbar_height_ratio is not None and cax is not None:
            cbox = cax.get_position()
            new_bottom = cbox.y0 + (1 - self.cbar_height_ratio) / 2 * cbox.height
            new_height = cbox.height * self.cbar_height_ratio
            cax.set_position([cbox.x0, new_bottom, cbox.width, new_height])

        # Set title and labels
        if self.cluster_cols:
            # When column clustering exists, use fig.suptitle at the top
            self.fig.suptitle(self.title, fontsize=self.title_size, y=0.99)
        else:
            # When no column clustering, use traditional title
            ax_heatmap.set_title(self.title, fontsize=self.title_size, pad=20)

        ax_heatmap.set_xlabel(self.xlabel, fontsize=self.label_size)
        ax_heatmap.set_ylabel(self.ylabel, fontsize=self.label_size)

        # Adjust tick labels
        ax_heatmap.tick_params(axis='x', rotation=45, labelsize=self.tick_size)
        ax_heatmap.tick_params(axis='y', rotation=0, labelsize=self.tick_size)

        # Auto adjust layout (not needed for dual clustering)
        if not (self.cluster_rows and self.cluster_cols):
            plt.tight_layout()

        # Save image
        if save_path:
            self.save_heatmap(save_path)

        return self.fig, ax_heatmap

    def save_heatmap(self, save_path):
        """
        Save heatmap as PNG format

        Parameters:
            save_path: Save path
        """
        if self.fig is None:
            raise ValueError("Please plot heatmap first")

        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)

        if save_path.suffix.lower() != '.png':
            save_path = save_path.with_suffix('.png')

        self.fig.savefig(
            save_path,
            dpi=self.dpi,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none'
        )
        print(f"Heatmap saved to: {save_path}")

    def show(self):
        """Show heatmap"""
        if self.fig is not None:
            plt.show()
        else:
            print("Please plot heatmap first")


def main():
    """Main function example"""

    print("=" * 60)
    print("Professional Heatmap Visualization Tool")
    print("=" * 60)

    # Example: Generate heatmap using sample data
    print("\nExample: Generate sample data and plot heatmap (default settings)")
    print("-" * 60)

    visualizer = HeatmapVisualizer(
        title="Gene Expression Heatmap Analysis (viridis default)",
        xlabel="Genes",
        ylabel="Samples",
        cmap="viridis",
        dpi=300,
        cluster_rows=True,
        cluster_cols=True
    )

    # Generate sample data
    visualizer.generate_sample_data(n_rows=20, n_cols=15)

    # Plot and save heatmap
    fig, ax = visualizer.plot_heatmap(save_path="heatmap_example_viridis.png")

    print("\n" + "=" * 60)
    print("Program execution completed!")
    print("Generated heatmap saved as: heatmap_example.png")
    print("=" * 60)


if __name__ == "__main__":
    main()
