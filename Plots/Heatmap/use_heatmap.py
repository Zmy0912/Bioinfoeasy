from heatmap_visualizer import HeatmapVisualizer

# Create visualizer
visualizer = HeatmapVisualizer(
    data_file="HOt.xlsx",        # Your data file name
    title="Gene Expression Heatmap Analysis",  # Chart title
    xlabel="Genes",              # X-axis label
    ylabel="Samples",            # Y-axis label
    cmap="coolwarm",            # Use coolwarm colormap
    vmin=0,                     # Minimum value
    vmax=1.4,                   # Maximum value
    title_size=14,              # Small title
    label_size=10,              # Small axis labels
    tick_size=12,               # Tick labels
    cbar_width_ratio=0.02,      # Color bar width 20%
    cbar_height_ratio=0.48,     # Color bar height 80%
    dpi=300,                    # Output resolution
    cluster_rows=False,         # Disable row clustering
    cluster_cols=False          # Disable column clustering
)

# Load data
visualizer.load_data()

# Plot and save heatmap
visualizer.plot_heatmap(save_path="my_heatmap.png")

# Display image
visualizer.show()
