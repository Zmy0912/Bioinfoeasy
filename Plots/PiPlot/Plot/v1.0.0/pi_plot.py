#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pi Sliding Window Gene Labeling Program
Line plot with sequence length as x-axis and Pi value as y-axis
Labels genes at high Pi value positions

Author: Mingyuan Zhang
Email: myzhang0726@foxmail.com
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import re
import argparse
import os

# Set font support
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


def parse_pi_file(pi_file):
    """Parse Pi file, extract window positions and Pi values"""
    data = []
    with open(pi_file, 'r', encoding='utf-8') as f:
        # Skip first header line
        next(f)
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            # Extract window range
            window = parts[0]
            match = re.search(r'(\d+)-(\d+)', window)
            if match:
                window_start = int(match.group(1))
                window_end = int(match.group(2))
                midpoint = float(parts[1])
                pi = float(parts[2])
                theta = float(parts[3])

                data.append({
                    'window': window,
                    'window_start': window_start,
                    'window_end': window_end,
                    'midpoint': midpoint,
                    'pi': pi,
                    'theta': theta
                })

    return pd.DataFrame(data)


def parse_gff3_file(gff3_file):
    """Parse GFF3 file, extract gene position information"""
    genes = []

    with open(gff3_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            # Only extract gene rows
            if parts[2] != 'gene':
                continue

            seqid = parts[0]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Extract gene name
            gene_name = None
            match = re.search(r'Name=([^;]+)', attributes)
            if match:
                gene_name = match.group(1)
            else:
                match = re.search(r'gene=([^;]+)', attributes)
                if match:
                    gene_name = match.group(1)

            # Extract gene ID
            gene_id = None
            match = re.search(r'ID=gene-([^;]+)', attributes)
            if match:
                gene_id = match.group(1)

            if gene_name:
                genes.append({
                    'gene_name': gene_name,
                    'gene_id': gene_id,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'midpoint': (start + end) / 2
                })

    return pd.DataFrame(genes)


def find_high_pi_windows(pi_df, threshold_percentile=90):
    """Find high Pi value windows"""
    threshold = np.percentile(pi_df['pi'], threshold_percentile)
    high_pi_windows = pi_df[pi_df['pi'] >= threshold].copy()
    high_pi_windows = high_pi_windows.sort_values('pi', ascending=False)
    return high_pi_windows, threshold


def find_significant_genes(pi_df, genes_df, window_size=10000, significance_std=2.0,
                           absolute_threshold=None):
    """
    Find genes significant relative to their nearby regions

    Parameters:
        pi_df: Pi data frame
        genes_df: Gene data frame
        window_size: Window size for calculating local average (bp)
        significance_std: Significance criterion (standard deviation multiple)
        absolute_threshold: Absolute Pi value threshold, genes above this value are also considered significant

    Returns:
        List of significant genes, each element contains gene information and significance score
    """
    significant_genes = []

    for _, gene in genes_df.iterrows():
        gene_mid = gene['midpoint']
        gene_start = gene['start']
        gene_end = gene['end']

        # Calculate Pi values at gene position
        gene_pi_values = pi_df[
            (pi_df['window_start'] <= gene_end) &
            (pi_df['window_end'] >= gene_start)
        ]['pi']

        if len(gene_pi_values) == 0:
            continue

        # Use average Pi value of gene region
        gene_pi_mean = gene_pi_values.mean()

        # Calculate Pi values of surrounding region
        region_start = max(0, gene_mid - window_size)
        region_end = min(pi_df['midpoint'].max(), gene_mid + window_size)

        region_pi = pi_df[
            (pi_df['midpoint'] >= region_start) &
            (pi_df['midpoint'] <= region_end)
        ]['pi']

        if len(region_pi) == 0:
            continue

        region_mean = region_pi.mean()
        region_std = region_pi.std()

        # Calculate Z-score (standardized score relative to surrounding region)
        if region_std > 0:
            z_score = (gene_pi_mean - region_mean) / region_std
        else:
            z_score = 0

        # Determine if significant: meet any of the following conditions
        # 1. Z-score >= significance_std (relative significance)
        # 2. Pi value >= absolute_threshold (absolute significance)
        is_significant = False
        significance_type = "none"

        if z_score >= significance_std:
            is_significant = True
            significance_type = "relative"  # Relative significance

        if absolute_threshold is not None and gene_pi_mean >= absolute_threshold:
            is_significant = True
            if significance_type == "none":
                significance_type = "absolute"  # Absolute significance
            else:
                significance_type = "both"  # Both satisfied

        if is_significant:
            significant_genes.append({
                'gene_name': gene['gene_name'],
                'start': gene['start'],
                'end': gene['end'],
                'midpoint': gene['midpoint'],
                'pi_value': gene_pi_mean,
                'z_score': z_score,
                'region_mean': region_mean,
                'region_std': region_std,
                'significance_type': significance_type
            })

    # Sort by z-score in descending order, then by pi_value in descending order
    significant_genes.sort(key=lambda x: (x['z_score'], x['pi_value']), reverse=True)
    return significant_genes


def find_genes_near_windows(genes_df, window_start, window_end, extend=2000):
    """Find genes within a certain range of the window"""
    nearby_genes = genes_df[
        ((genes_df['start'] >= window_start - extend) & (genes_df['start'] <= window_end + extend)) |
        ((genes_df['end'] >= window_start - extend) & (genes_df['end'] <= window_end + extend)) |
        ((genes_df['start'] <= window_start) & (genes_df['end'] >= window_end))
    ].copy()

    # Calculate distance between gene center and window center
    window_center = (window_start + window_end) / 2
    nearby_genes['distance'] = abs(nearby_genes['midpoint'] - window_center)
    nearby_genes = nearby_genes.sort_values('distance')

    return nearby_genes


def plot_pi_with_genes(pi_file, gff3_file, output_file='pi_plot.png',
                       gene_fontsize=12, max_genes=15, window_size=10000,
                       significance_std=2.0, absolute_threshold=None):
    """
    Draw Pi sliding window plot and label genes significant relative to nearby regions

    Parameters:
        pi_file: Pi value file path
        gff3_file: GFF3 gene annotation file path
        output_file: Output image file name
        gene_fontsize: Gene label font size
        max_genes: Maximum number of genes to label
        window_size: Window size for calculating local average (bp)
        significance_std: Significance criterion (standard deviation multiple)
        absolute_threshold: Absolute Pi value threshold
    """

    # Read data
    print(f"Reading Pi file: {pi_file}")
    pi_df = parse_pi_file(pi_file)

    print(f"Reading GFF3 file: {gff3_file}")
    genes_df = parse_gff3_file(gff3_file)

    print(f"Found {len(pi_df)} sliding windows")
    print(f"Found {len(genes_df)} genes")

    # Find high Pi value windows (top 10% highest)
    high_pi_windows, pi_threshold = find_high_pi_windows(pi_df, threshold_percentile=90)

    print(f"\nHigh Pi threshold: {pi_threshold:.4f}")
    print(f"Found {len(high_pi_windows)} high Pi value windows")

    # Find genes significant relative to surrounding regions
    print(f"\nAnalyzing gene significance...")
    print(f"  Window size: {window_size} bp")
    print(f"  Relative significance criterion: {significance_std} standard deviations")
    if absolute_threshold is not None:
        print(f"  Absolute significance threshold: {absolute_threshold}")

    significant_genes = find_significant_genes(pi_df, genes_df,
                                            window_size=window_size,
                                            significance_std=significance_std,
                                            absolute_threshold=absolute_threshold)

    print(f"Found {len(significant_genes)} significant genes")

    # Select top N most significant genes to label
    genes_to_label = significant_genes[:max_genes]

    print(f"Labeling top {len(genes_to_label)} most significant genes\n")

    # Create figure
    fig, ax = plt.subplots(figsize=(20, 12))

    # Draw Pi line plot
    ax.plot(pi_df['midpoint'], pi_df['pi'], 'b-', linewidth=1.5, alpha=0.7, label='Pi value')

    # Draw high Pi value threshold line
    if absolute_threshold is not None:
        ax.axhline(y=absolute_threshold, color='purple', linestyle='--', linewidth=1.5,
                   label=f'Absolute Threshold ({absolute_threshold:.4f})')

    ax.axhline(y=pi_threshold, color='r', linestyle='--', linewidth=1,
               label=f'Relative Threshold (90th: {pi_threshold:.4f})')

    # Calculate labeling area above Y-axis
    max_pi = pi_df['pi'].max()
    arrow_length = max_pi * 0.01  # Arrow length (very short)
    label_spacing = max_pi * 0.15  # Vertical spacing between labels (reduced)
    min_horizontal_distance = max_pi * 0.05  # Minimum horizontal spacing

    # Sort genes by x-axis position first for easy overlap handling
    genes_sorted = sorted(genes_to_label, key=lambda x: x['midpoint'])

    # Define distance threshold for adjacent genes (15000bp)
    adjacent_gene_distance = 15000

    # Group adjacent genes (genes within 15000bp)
    gene_groups = []
    current_group = [genes_sorted[0]]

    for gene in genes_sorted[1:]:
        # Check distance between current gene and last gene in group
        if gene['midpoint'] - current_group[-1]['midpoint'] <= adjacent_gene_distance:
            # Within 15000bp range, add to same group
            current_group.append(gene)
        else:
            # Exceeds range, start new group
            gene_groups.append(current_group)
            current_group = [gene]

    # Add last group
    if current_group:
        gene_groups.append(current_group)

    # Pre-calculate all label positions to avoid overlaps
    label_positions = []

    for group in gene_groups:
        # Assign different height layers to genes in group
        group_size = len(group)
        base_height = max_pi * 0.03  # Base height (shorter)

        for i, gene_info in enumerate(group):
            gene_name = gene_info['gene_name']
            gene_mid = gene_info['midpoint']
            gene_pi = gene_info['pi_value']
            z_score = gene_info['z_score']
            sig_type = gene_info['significance_type']

            # Mark point at gene position
            marker_color = 'red' if sig_type == 'relative' else 'purple' if sig_type == 'absolute' else 'darkred'
            ax.plot(gene_mid, gene_pi, 'o', color=marker_color, markersize=5, zorder=5, alpha=0.8,
                   markeredgecolor='black', markeredgewidth=1)

            # Assign different height layers based on gene position in group
            layer_index = i % group_size
            layer_height = base_height + layer_index * label_spacing

            # Calculate label vertical position (directly above line plot, using short pointer)
            label_y = gene_pi + layer_height

            # Estimate text width (based on character count)
            label_width = len(gene_name) * max_pi * 0.008  # Rough estimate
            current_min_horizontal_distance = max(min_horizontal_distance, label_width + max_pi * 0.02)

            # Check and adjust position to avoid overlap
            final_y = label_y
            final_x = gene_mid

            # Multiple iterations to adjust position, ensure no overlap with any existing labels
            max_iterations = 50
            iteration = 0

            while iteration < max_iterations:
                has_overlap = False

                for (prev_x, prev_y, prev_width, prev_y_bottom) in label_positions:
                    # Check horizontal overlap
                    if abs(final_x - prev_x) < current_min_horizontal_distance:
                        # Check if vertical also overlaps
                        current_y_bottom = gene_pi  # Label bottom
                        current_y_top = final_y  # Label top

                        prev_y_bottom_prev = prev_y_bottom
                        prev_y_top = prev_y

                        # Check if vertical ranges overlap
                        if not (current_y_top < prev_y_bottom_prev or current_y_bottom > prev_y_top):
                            # Has overlap, move up
                            final_y = prev_y_top + label_spacing
                            has_overlap = True
                            break

                if not has_overlap:
                    break

                iteration += 1

            # Record label position, including bottom (for subsequent detection)
            label_positions.append((final_x, final_y, label_width, gene_pi))

            # Determine label color
            if sig_type == 'both':
                label_color = 'orange'
                edge_color = 'red'
                label_text = f'{gene_name}'
            elif sig_type == 'absolute':
                label_color = 'lightblue'
                edge_color = 'purple'
                label_text = f'{gene_name}'
            else:  # relative
                label_color = 'yellow'
                edge_color = 'red'
                label_text = f'{gene_name}'

            # Draw short pointer (vertical line)
            ax.annotate('', xy=(gene_mid, gene_pi), xytext=(gene_mid, final_y),
                       arrowprops=dict(arrowstyle='-',
                                     color=edge_color,
                                     linewidth=1.5,
                                     alpha=0.7),
                       zorder=5)

            # Add gene name label (above line plot, pointed by short pointer)
            ax.text(final_x, final_y, label_text,
                   fontsize=gene_fontsize,
                   fontweight='bold',
                   ha='center',
                   va='bottom',
                   bbox=dict(boxstyle='round,pad=0.4',
                           facecolor=label_color,
                           alpha=0.85,
                           edgecolor=edge_color,
                           linewidth=1.5),
                   zorder=6)

    # Set figure properties
    ax.set_xlabel('Sequence Position (bp)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Pi (Nucleotide Diversity)', fontsize=14, fontweight='bold')

    # Title away from figure body
    ax.set_title('Sliding Window Analysis: Nucleotide Diversity (Pi) Along Genome\n'
                 f'Markers indicate genes with elevated diversity',
                 fontsize=16, fontweight='bold', pad=30)

    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9, bbox_to_anchor=(0.98, 0.98))

    # Adjust X-axis and Y-axis range, leave space for labels
    # X-coordinate range is 40000bp larger than annotation data range
    ax.set_xlim(pi_df['midpoint'].min() - 20000, pi_df['midpoint'].max() + 20000)
    # Adjust Y-axis upper limit based on actual label height (add more space)
    if label_positions:
        max_label_y = max(pos[1] for pos in label_positions)
        ax.set_ylim(0, max_label_y * 1.25)
    else:
        ax.set_ylim(0, max_pi * 1.3)

    # Add significant gene information table outside right side of figure (doesn't affect main figure)
    if genes_to_label:
        info_text = "Significantly Elevated Genes:\n"
        info_text += "-" * 40 + "\n"

        for i, gene_info in enumerate(genes_to_label[:15], 1):
            gene_name = gene_info['gene_name']
            z_score = gene_info['z_score']
            gene_pi = gene_info['pi_value']
            sig_type = gene_info['significance_type']

            type_indicator = {'relative': '[R]', 'absolute': '[A]', 'both': '[R+A]'}
            info_text += f"{i}. {type_indicator.get(sig_type, '[?]')} {gene_name}\n"
            info_text += f"   Pi={gene_pi:.4f}, Z={z_score:.2f}\n"

        if len(genes_to_label) > 15:
            info_text += f"... and {len(genes_to_label) - 15} more\n"

        info_text += "\n[R]=Relative significant\n[A]=Absolute significant\n[R+A]=Both"

        # Use axes to position outside right side of figure
        fig.text(0.92, 0.5, info_text,
                fontsize=9,
                bbox=dict(boxstyle='round',
                         facecolor='lightyellow',
                         alpha=0.95,
                         edgecolor='orange',
                         linewidth=2),
                verticalalignment='center',
                family='monospace')

    # Adjust layout, leave space for right side information
    plt.tight_layout()
    plt.subplots_adjust(right=0.88)

    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nChart saved to: {output_file}")

    # Show figure
    plt.show()

    # Output significant gene details to file
    output_info = output_file.replace('.png', '_info.txt')
    with open(output_info, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("Significantly Elevated Genes Analysis\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Analysis Parameters:\n")
        f.write(f"  Window size for local comparison: {window_size} bp\n")
        f.write(f"  Relative significance threshold: {significance_std} standard deviations\n")
        if absolute_threshold is not None:
            f.write(f"  Absolute Pi threshold: {absolute_threshold}\n")
        f.write(f"  Number of genes analyzed: {len(genes_df)}\n")
        f.write(f"  Number of significantly elevated genes found: {len(significant_genes)}\n\n")

        f.write(f"Pi Threshold (90th percentile): {pi_threshold:.4f}\n\n")

        f.write("=" * 80 + "\n")
        f.write("Top Significantly Elevated Genes:\n")
        f.write("=" * 80 + "\n\n")

        for i, gene_info in enumerate(genes_to_label, 1):
            f.write(f"\n{i}. Gene: {gene_info['gene_name']}\n")
            f.write(f"   Position: {gene_info['start']}-{gene_info['end']} bp\n")
            f.write(f"   Gene midpoint: {gene_info['midpoint']:.0f} bp\n")
            f.write(f"   Gene Pi value: {gene_info['pi_value']:.4f}\n")
            f.write(f"   Region mean Pi: {gene_info['region_mean']:.4f}\n")
            f.write(f"   Region std Pi: {gene_info['region_std']:.4f}\n")
            f.write(f"   Z-score: {gene_info['z_score']:.4f}\n")
            f.write(f"   Significance type: {gene_info['significance_type']}\n")
            f.write(f"   Significance: {gene_info['z_score']:.2f}x above region average\n")

        if len(significant_genes) > max_genes:
            f.write(f"\n... and {len(significant_genes) - max_genes} more significantly elevated genes\n")

    print(f"Detailed information saved to: {output_info}")


def main():
    """Main function, handle command line arguments"""
    parser = argparse.ArgumentParser(
        description='Pi Sliding Window Gene Labeling Program - Label genes significant relative to surrounding regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default parameters (relative significance only)
  python pi_plot.py -p Pi.txt -g sequence.gff3

  # Use absolute threshold
  python pi_plot.py -p Pi.txt -g sequence.gff3 --absolute-threshold 0.05

  # Use both relative significance and absolute threshold
  python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 2.0 --absolute-threshold 0.05

  # Specify output file and font size
  python pi_plot.py -p Pi.txt -g sequence.gff3 -o my_plot.png --gene-fontsize 14

  # Adjust significance parameters
  python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 15000 --significance-std 2.5

  # Limit number of labeled genes
  python pi_plot.py -p Pi.txt -g sequence.gff3 --max-genes 10

Note:
  Relative significance: Gene Pi value is a certain multiple of standard deviations above surrounding region average
  Absolute significance: Gene Pi value is above a set absolute threshold
  Genes meeting either condition will be labeled
        """
    )

    parser.add_argument('-p', '--pi-file',
                       required=True,
                       help='Pi value file path (required)')

    parser.add_argument('-g', '--gff3-file',
                       required=True,
                       help='GFF3 gene annotation file path (required)')

    parser.add_argument('-o', '--output',
                       default='pi_sliding_window_plot.png',
                       help='Output image file name (default: pi_sliding_window_plot.png)')

    parser.add_argument('--gene-fontsize',
                       type=int,
                       default=12,
                       help='Gene label font size (default: 12)')

    parser.add_argument('--max-genes',
                       type=int,
                       default=15,
                       help='Maximum number of genes to label (default: 15)')

    parser.add_argument('--window-size',
                       type=int,
                       default=10000,
                       help='Window size for calculating local average, unit bp (default: 10000)')

    parser.add_argument('--significance-std',
                       type=float,
                       default=2.0,
                       help='Relative significance criterion, gene Pi value must be X standard deviations above surrounding region average (default: 2.0)')

    parser.add_argument('--absolute-threshold',
                       type=float,
                       default=None,
                       help='Absolute Pi value threshold, genes above this value are also considered significant (optional, recommended value such as 0.05)')

    args = parser.parse_args()

    # Check if files exist
    if not os.path.exists(args.pi_file):
        print(f"Error: Pi file does not exist: {args.pi_file}")
        return 1

    if not os.path.exists(args.gff3_file):
        print(f"Error: GFF3 file does not exist: {args.gff3_file}")
        return 1

    # Print parameter information
    print("=" * 70)
    print("Pi Sliding Window Gene Labeling Program")
    print("Author: Mingyuan Zhang")
    print("Email: myzhang0726@foxmail.com")
    print("=" * 70)
    print(f"Pi file: {args.pi_file}")
    print(f"GFF3 file: {args.gff3_file}")
    print(f"Output file: {args.output}")
    print(f"Gene label font size: {args.gene_fontsize}")
    print(f"Maximum number of labeled genes: {args.max_genes}")
    print(f"Comparison window size: {args.window_size} bp")
    print(f"Relative significance criterion: {args.significance_std} standard deviations")
    if args.absolute_threshold is not None:
        print(f"Absolute significance threshold: {args.absolute_threshold}")
    print("=" * 70)
    print()

    try:
        # Draw chart
        plot_pi_with_genes(
            pi_file=args.pi_file,
            gff3_file=args.gff3_file,
            output_file=args.output,
            gene_fontsize=args.gene_fontsize,
            max_genes=args.max_genes,
            window_size=args.window_size,
            significance_std=args.significance_std,
            absolute_threshold=args.absolute_threshold
        )

        print("\n" + "=" * 70)
        print("Program ran successfully!")
        print("=" * 70)
        return 0

    except Exception as e:
        print(f"\nError: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    exit(main())
