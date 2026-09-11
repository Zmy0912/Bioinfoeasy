#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pi Sliding Window Gene Annotation Program
Line plot with sequence position on x-axis and pi value on y-axis
Annotate genes at high pi value positions
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import re
import argparse
import os

# Set Chinese font support
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


def parse_pi_file(pi_file):
    """Parse Pi file and extract window positions and Pi values"""
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
    """Parse GFF3 file and extract gene position information"""
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
    """Find windows with high Pi values"""
    threshold = np.percentile(pi_df['pi'], threshold_percentile)
    high_pi_windows = pi_df[pi_df['pi'] >= threshold].copy()
    high_pi_windows = high_pi_windows.sort_values('pi', ascending=False)
    return high_pi_windows, threshold


def find_significant_genes(pi_df, genes_df, window_size=10000, significance_std=2.0,
                           absolute_threshold=None, lower_threshold=None,
                           compute_zscore=True):
    """
    Find genes that are significant relative to their surrounding regions

    Parameters:
        pi_df: Pi dataframe
        genes_df: Gene dataframe
        window_size: Window size for calculating local average (bp)
        significance_std: Significance standard (standard deviation multiplier)
        absolute_threshold: Absolute Pi value threshold, genes above this value are considered significant
        lower_threshold: Lower threshold, genes below this value are considered significant (lowest Pi genes)
        compute_zscore: Whether to calculate Z-score (False for absolute-only or lower-only mode)

    Returns:
        List of significant genes, each element contains gene information and significance score
    """
    significant_genes = []

    for _, gene in genes_df.iterrows():
        gene_mid = gene['midpoint']
        gene_start = gene['start']
        gene_end = gene['end']

        # Calculate Pi value at gene position
        gene_pi_values = pi_df[
            (pi_df['window_start'] <= gene_end) &
            (pi_df['window_end'] >= gene_start)
        ]['pi']

        if len(gene_pi_values) == 0:
            continue

        # Use average Pi value of gene region
        gene_pi_mean = gene_pi_values.mean()

        # Calculate Pi value of surrounding region (only when Z-score is needed)
        region_mean = 0
        region_std = 0
        z_score = 0

        if compute_zscore:
            region_start = max(0, gene_mid - window_size)
            region_end = min(pi_df['midpoint'].max(), gene_mid + window_size)

            region_pi = pi_df[
                (pi_df['midpoint'] >= region_start) &
                (pi_df['midpoint'] <= region_end)
            ]['pi']

            if len(region_pi) > 0:
                region_mean = region_pi.mean()
                region_std = region_pi.std()

                # Calculate Z-score (standardized score relative to surrounding region)
                if region_std > 0:
                    z_score = (gene_pi_mean - region_mean) / region_std

        # Determine significance: gene is significant if it meets any of the following conditions
        # 1. Z-score >= significance_std (relative significance)
        # 2. Pi value >= absolute_threshold (absolute significance)
        # 3. Pi value <= lower_threshold (lowest Pi value genes)
        is_significant = False
        significance_type = "none"

        # Record whether each condition is met
        meets_relative = compute_zscore and z_score >= significance_std
        meets_absolute = absolute_threshold is not None and gene_pi_mean >= absolute_threshold
        meets_lower = lower_threshold is not None and gene_pi_mean <= lower_threshold

        # Determine significance_type based on number and type of conditions met
        if meets_relative and meets_absolute and meets_lower:
            is_significant = True
            significance_type = "lower+both"  # Meets all three L+R+A
        elif meets_relative and meets_absolute:
            is_significant = True
            significance_type = "both"  # Meets relative and absolute R+A
        elif meets_relative and meets_lower:
            is_significant = True
            significance_type = "lower+relative"  # Meets lower and relative L+R
        elif meets_absolute and meets_lower:
            is_significant = True
            significance_type = "lower+absolute"  # Meets lower and absolute L+A
        elif meets_relative:
            is_significant = True
            significance_type = "relative"  # Relative only R
        elif meets_absolute:
            is_significant = True
            significance_type = "absolute"  # Absolute only A
        elif meets_lower:
            is_significant = True
            significance_type = "lower"  # Lower only L

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

    # Sort by Z-score descending, then by pi_value descending
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
                       significance_std=2.0, absolute_threshold=None, lower_threshold=None):
    """
    Draw Pi sliding window plot and annotate genes significant relative to surrounding regions

    Parameters:
        pi_file: Pi value file path
        gff3_file: GFF3 gene annotation file path
        output_file: Output image file name
        gene_fontsize: Gene annotation font size
        max_genes: Maximum number of genes to annotate
        window_size: Window size for calculating local average (bp)
        significance_std: Significance standard (standard deviation multiplier)
        absolute_threshold: Absolute Pi value threshold
        lower_threshold: Lower threshold, genes below this value are considered significant
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

    print(f"\nHigh Pi value threshold: {pi_threshold:.4f}")
    print(f"Found {len(high_pi_windows)} high Pi value windows")

    # Find genes significant relative to surrounding regions
    print(f"\nAnalyzing gene significance...")
    print(f"  Window size: {window_size} bp")
    if significance_std is not None:
        print(f"  Relative significance standard: {significance_std} standard deviations")
    else:
        print(f"  Relative significance: Disabled (no Z-score calculation)")
    if absolute_threshold is not None:
        print(f"  Absolute significance threshold: {absolute_threshold}")
    if lower_threshold is not None:
        print(f"  Lower threshold: {lower_threshold}")

    # Determine output logic based on user input parameters
    # Prioritize checking if significance_std is explicitly set (not None)
    # If significance_std is explicitly set, use Z-score mode (relative threshold priority)
    use_relative_explicit = (significance_std is not None)
    use_absolute_threshold = (absolute_threshold is not None)
    use_lower_threshold = (lower_threshold is not None)

    # Filter and sort genes based on requirements
    if use_relative_explicit:
        # significance_std explicitly specified: use Z-score analysis mode, sort by Z-score
        compute_zscore = True
        output_mode = "relative_priority"
        significant_genes = find_significant_genes(pi_df, genes_df,
                                            window_size=window_size,
                                            significance_std=significance_std,
                                            absolute_threshold=absolute_threshold,
                                            lower_threshold=lower_threshold,
                                            compute_zscore=compute_zscore)
        
        # In Z-score mode, prioritize genes meeting multiple conditions, then sort by Z-score
        # Priority: L+R+A (3 conditions) > R+A (2 conditions) > L+R (2 conditions) > L+A (2 conditions) > R (1 condition) > A (1 condition) > L (1 condition)
        def priority_score(gene):
            sig_type = gene['significance_type']
            if sig_type == 'lower+both':
                # L+R+A: Highest priority, sort by Z-score
                return (6, gene['z_score'], gene['pi_value'])
            elif sig_type == 'both':
                # R+A: Second priority, sort by Z-score
                return (5, gene['z_score'], gene['pi_value'])
            elif sig_type == 'lower+relative':
                # L+R: Third priority, sort by Z-score
                return (4, gene['z_score'], gene['pi_value'])
            elif sig_type == 'lower+absolute':
                # L+A: Fourth priority, sort by Z-score
                return (3, gene['z_score'], gene['pi_value'])
            elif sig_type == 'relative':
                # R: Fifth priority, sort by Z-score
                return (2, gene['z_score'], gene['pi_value'])
            elif sig_type == 'absolute':
                # A: Sixth priority, sort by Z-score
                return (1, gene['z_score'], gene['pi_value'])
            elif sig_type == 'lower':
                # L: Lowest priority
                return (0, gene['z_score'], gene['pi_value'])
            else:
                return (0, 0, 0)
        
        significant_genes.sort(key=priority_score, reverse=True)
        genes_to_label = significant_genes[:max_genes]

    elif use_lower_threshold and not use_absolute_threshold:
        # Only lower threshold: output only lower threshold results, sort by Pi value ascending, no Z-score calculation
        compute_zscore = False
        output_mode = "lower_only"
        significant_genes = find_significant_genes(pi_df, genes_df,
                                            window_size=window_size,
                                            significance_std=significance_std,
                                            absolute_threshold=absolute_threshold,
                                            lower_threshold=lower_threshold,
                                            compute_zscore=compute_zscore)
        lower_genes = [g for g in significant_genes if 'lower' in g['significance_type']]
        lower_genes.sort(key=lambda x: x['pi_value'])  # Sort ascending
        genes_to_label = lower_genes[:max_genes]

    elif use_absolute_threshold and not use_lower_threshold:
        # Only absolute threshold: directly filter all genes with Pi value greater than absolute threshold, sort by Pi value descending, no significance calculation
        compute_zscore = False
        output_mode = "absolute_only"

        # Directly filter genes with Pi value greater than absolute threshold
        significant_genes = []
        for _, gene in genes_df.iterrows():
            gene_start = gene['start']
            gene_end = gene['end']

            # Calculate Pi value at gene position
            gene_pi_values = pi_df[
                (pi_df['window_start'] <= gene_end) &
                (pi_df['window_end'] >= gene_start)
            ]['pi']

            if len(gene_pi_values) == 0:
                continue

            # Use average Pi value of gene region
            gene_pi_mean = gene_pi_values.mean()

            # Filter genes with Pi value greater than absolute threshold
            if gene_pi_mean > absolute_threshold:
                significant_genes.append({
                    'gene_name': gene['gene_name'],
                    'start': gene['start'],
                    'end': gene['end'],
                    'midpoint': gene['midpoint'],
                    'pi_value': gene_pi_mean,
                    'z_score': 0,
                    'region_mean': 0,
                    'region_std': 0,
                    'significance_type': "absolute"
                })

        # Sort by Pi value descending
        significant_genes.sort(key=lambda x: x['pi_value'], reverse=True)
        genes_to_label = significant_genes[:max_genes]

    elif use_lower_threshold and use_absolute_threshold:
        # Lower threshold and absolute threshold specified, and significance_std not explicitly set: distribute by gene count
        # For odd numbers, remainder assigned to absolute threshold
        # No Z-score calculation, sort by Pi value
        compute_zscore = False
        output_mode = "lower+absolute"
        significant_genes = find_significant_genes(pi_df, genes_df,
                                            window_size=window_size,
                                            significance_std=significance_std,
                                            absolute_threshold=absolute_threshold,
                                            lower_threshold=lower_threshold,
                                            compute_zscore=compute_zscore)
        lower_genes = [g for g in significant_genes if 'lower' in g['significance_type']]
        lower_genes.sort(key=lambda x: x['pi_value'])  # Ascending
        absolute_genes = [g for g in significant_genes if 'absolute' in g['significance_type']]
        absolute_genes.sort(key=lambda x: x['pi_value'], reverse=True)  # Descending

        # Limit allocation within max_genes range
        total_genes = len(significant_genes)
        genes_limit = min(total_genes, max_genes)
        absolute_count = (genes_limit + 1) // 2  # For odd numbers, remainder to absolute threshold
        lower_count = genes_limit // 2

        genes_to_label = absolute_genes[:absolute_count] + lower_genes[:lower_count]

    else:
        # No threshold parameters set: use default relative threshold mode
        compute_zscore = True
        output_mode = "relative_only"
        significant_genes = find_significant_genes(pi_df, genes_df,
                                            window_size=window_size,
                                            significance_std=significance_std if significance_std is not None else 2.0,
                                            absolute_threshold=absolute_threshold,
                                            lower_threshold=lower_threshold,
                                            compute_zscore=compute_zscore)
        relative_genes = [g for g in significant_genes if 'relative' in g['significance_type']]
        relative_genes.sort(key=lambda x: (x['z_score'], x['pi_value']), reverse=True)
        genes_to_label = relative_genes[:max_genes]

    # Output calculation mode information
    if compute_zscore:
        print(f"\nUsing Z-score analysis mode (sorted by Z-score)")
    else:
        print(f"\nUsing Pi value sorting mode (no Z-score calculation)")

    print(f"Found {len(significant_genes)} significant genes")

    if not genes_to_label:
        print("Warning: No genes met the criteria")
    else:
        print(f"Will annotate the top {len(genes_to_label)} most significant genes\n")

    # Create figure
    fig, ax = plt.subplots(figsize=(20, 12))

    # Draw Pi line plot
    ax.plot(pi_df['midpoint'], pi_df['pi'], 'b-', linewidth=1.5, alpha=0.7, label='Pi value')

    # Draw high Pi value threshold line
    if absolute_threshold is not None:
        ax.axhline(y=absolute_threshold, color='purple', linestyle='--', linewidth=1.5,
                   label=f'Absolute Threshold ({absolute_threshold:.4f})')

    if lower_threshold is not None:
        ax.axhline(y=lower_threshold, color='green', linestyle='--', linewidth=1.5,
                   label=f'Lower Threshold ({lower_threshold:.4f})')

    ax.axhline(y=pi_threshold, color='r', linestyle='--', linewidth=1,
               label=f'Relative Threshold (90th: {pi_threshold:.4f})')

    # Calculate annotation area above Y-axis
    max_pi = pi_df['pi'].max()
    arrow_length = max_pi * 0.05  # Arrow length (very short)
    label_spacing = max_pi * 0.15  # Vertical spacing between labels (reduced)
    min_horizontal_distance = max_pi * 0.05  # Minimum horizontal distance

    # Sort genes by x-axis position first for easier overlap handling
    genes_sorted = sorted(genes_to_label, key=lambda x: x['midpoint'])

    # Check if there are genes to annotate
    if not genes_sorted:
        # Set figure properties
        ax.set_xlabel('Sequence Position (bp)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Pi (Nucleotide Diversity)', fontsize=14, fontweight='bold')

        ax.set_title('Sliding Window Analysis: Nucleotide Diversity (Pi) Along Genome\n'
                     f'Markers indicate genes with elevated diversity',
                     fontsize=16, fontweight='bold', pad=30)

        ax.grid(True, alpha=0.3, linestyle=':')
        ax.legend(loc='upper right', fontsize=11, framealpha=0.9, bbox_to_anchor=(0.98, 0.98))

        # Adjust X-axis and Y-axis ranges
        ax.set_xlim(pi_df['midpoint'].min() - 20000, pi_df['midpoint'].max() + 20000)
        ax.set_ylim(0, max_pi * 1.3)

        # Save figure
        plt.tight_layout()
        plt.subplots_adjust(right=0.88)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nChart saved to: {output_file}")
        plt.show()

        # Output information to file
        output_info = output_file.replace('.png', '_info.txt')
        with open(output_info, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("No Significant Genes Found\n")
            f.write("=" * 80 + "\n\n")
            f.write("No genes met the specified criteria.\n\n")
            f.write(f"Analysis Parameters:\n")
            f.write(f"  Window size for local comparison: {window_size} bp\n")
            f.write(f"  Relative significance threshold: {significance_std} standard deviations\n")
            if absolute_threshold is not None:
                f.write(f"  Absolute Pi threshold: {absolute_threshold}\n")
            if lower_threshold is not None:
                f.write(f"  Lower Pi threshold: {lower_threshold}\n")
            f.write(f"  Number of genes analyzed: {len(genes_df)}\n")

        print(f"Detailed information saved to: {output_info}")
        return

    # Define adjacent gene distance threshold (15000bp)
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
            # Exceeded range, start new group
            gene_groups.append(current_group)
            current_group = [gene]

    # Add last group
    if current_group:
        gene_groups.append(current_group)

    # Pre-calculate all annotation positions to avoid overlaps
    label_positions = []

    for group in gene_groups:
        # Assign different height layers to genes within group
        group_size = len(group)
        base_height = max_pi * 0.03  # Base height (shorter)

        for i, gene_info in enumerate(group):
            gene_name = gene_info['gene_name']
            gene_mid = gene_info['midpoint']
            gene_pi = gene_info['pi_value']
            z_score = gene_info['z_score']
            sig_type = gene_info['significance_type']

            # Annotate point at gene position
            if 'lower' in sig_type:
                marker_color = 'green'
            elif sig_type == 'absolute':
                marker_color = 'purple'
            elif sig_type == 'both':
                marker_color = 'darkred'
            else:  # relative
                marker_color = 'red'
            ax.plot(gene_mid, gene_pi, 'o', color=marker_color, markersize=5, zorder=5, alpha=0.8,
                   markeredgecolor='black', markeredgewidth=1)

            # Assign different height layers based on gene's position in group
            layer_index = i % group_size
            layer_height = base_height + layer_index * label_spacing

            # Calculate vertical position of annotation (directly above line, with short pointer)
            label_y = gene_pi + layer_height

            # Estimate text width (based on character count)
            label_width = len(gene_name) * max_pi * 0.008  # Rough estimate
            current_min_horizontal_distance = max(min_horizontal_distance, label_width + max_pi * 0.02)

            # Check and adjust position to avoid overlap
            final_y = label_y
            final_x = gene_mid

            # Multiple iterations to adjust position, ensuring no overlap with any existing annotations
            max_iterations = 50
            iteration = 0

            while iteration < max_iterations:
                has_overlap = False

                for (prev_x, prev_y, prev_width, prev_y_bottom) in label_positions:
                    # Check horizontal direction overlap
                    if abs(final_x - prev_x) < current_min_horizontal_distance:
                        # Check if vertical direction also overlaps
                        current_y_bottom = gene_pi  # Annotation bottom
                        current_y_top = final_y  # Annotation top

                        prev_y_bottom_prev = prev_y_bottom
                        prev_y_top = prev_y

                        # Check if vertical ranges overlap
                        if not (current_y_top < prev_y_bottom_prev or current_y_bottom > prev_y_top):
                            # Overlap, move up
                            final_y = prev_y_top + label_spacing
                            has_overlap = True
                            break

                if not has_overlap:
                    break

                iteration += 1

            # Record annotation position, including bottom (for subsequent detection)
            label_positions.append((final_x, final_y, label_width, gene_pi))

            # Determine annotation color
            if 'lower' in sig_type:
                if 'both' in sig_type:
                    label_color = 'lightgreen'
                    edge_color = 'green'
                elif 'absolute' in sig_type:
                    label_color = 'lightcyan'
                    edge_color = 'green'
                elif 'relative' in sig_type:
                    label_color = 'lightyellow'
                    edge_color = 'green'
                else:
                    label_color = 'lightgreen'
                    edge_color = 'green'
            elif sig_type == 'both':
                label_color = 'orange'
                edge_color = 'red'
            elif sig_type == 'absolute':
                label_color = 'lightblue'
                edge_color = 'purple'
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

            # Add gene name annotation (above line, pointed to by short pointer)
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

    # Title away from main chart body
    ax.set_title('Sliding Window Analysis: Nucleotide Diversity (Pi) Along Genome\n'
                 f'Markers indicate genes with elevated diversity',
                 fontsize=16, fontweight='bold', pad=30)

    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9, bbox_to_anchor=(0.98, 0.98))

    # Adjust X-axis and Y-axis ranges, leave space for annotations
    # Horizontal coordinate range 40000bp larger than annotation data range
    ax.set_xlim(pi_df['midpoint'].min() - 20000, pi_df['midpoint'].max() + 20000)
    # Adjust Y-axis upper limit based on actual annotation height (add more space)
    if label_positions:
        max_label_y = max(pos[1] for pos in label_positions)
        ax.set_ylim(0, max_label_y * 1.25)
    else:
        ax.set_ylim(0, max_pi * 1.3)

    # Add significant genes information table outside right side of chart (doesn't affect main chart)
    if genes_to_label:
        if output_mode == "lower_only":
            info_text = "Lowest Pi Value Genes:\n"
        else:
            info_text = "Significantly Elevated Genes:\n"
        info_text += "-" * 40 + "\n"

        for i, gene_info in enumerate(genes_to_label[:15], 1):
            gene_name = gene_info['gene_name']
            z_score = gene_info['z_score']
            gene_pi = gene_info['pi_value']
            sig_type = gene_info['significance_type']

            type_indicator = {
                'relative': '[R]',
                'absolute': '[A]',
                'both': '[R+A]',
                'lower': '[L]',
                'lower+absolute': '[L+A]',
                'lower+relative': '[L+R]',
                'lower+both': '[L+R+A]'
            }
            info_text += f"{i}. {type_indicator.get(sig_type, '[?]')} {gene_name}\n"
            info_text += f"   Pi={gene_pi:.4f}, Z={z_score:.2f}\n"

        if len(genes_to_label) > 15:
            info_text += f"... and {len(genes_to_label) - 15} more\n"

        info_text += "\n\nSignificance Type Priority:\n[L+R+A]=Lower+Relative+Absolute (Highest)\n[R+A]=Relative+Absolute\n[L+R]=Lower+Relative\n[L+A]=Lower+Absolute\n[R]=Relative significant\n[A]=Absolute significant\n[L]=Lowest Pi value (Lowest)"

        # Use axes positioning outside right side of chart
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

    # Display figure
    plt.show()

    # Output significant genes detailed information to file
    output_info = output_file.replace('.png', '_info.txt')
    with open(output_info, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        if output_mode == "lower_only":
            f.write("Lowest Pi Value Genes Analysis\n")
        else:
            f.write("Significantly Elevated Genes Analysis\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Analysis Parameters:\n")
        f.write(f"  Window size for local comparison: {window_size} bp\n")
        f.write(f"  Relative significance threshold: {significance_std} standard deviations\n")
        if absolute_threshold is not None:
            f.write(f"  Absolute Pi threshold: {absolute_threshold}\n")
        if lower_threshold is not None:
            f.write(f"  Lower Pi threshold: {lower_threshold}\n")
        f.write(f"  Number of genes analyzed: {len(genes_df)}\n")
        f.write(f"  Number of significant genes found: {len(significant_genes)}\n")
        f.write(f"  Output mode: {output_mode}\n\n")

        f.write(f"Pi Threshold (90th percentile): {pi_threshold:.4f}\n\n")

        f.write("=" * 80 + "\n")
        if output_mode == "lower_only":
            f.write("Lowest Pi Value Genes:\n")
        else:
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
        description='Gene Diversity Sliding Window Plotting Program - Annotate genes significant relative to surrounding regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default parameters (relative significance only)
  python pi_plot.py -p Pi.txt -g sequence.gff3

  # Use lower threshold (detect lowest Pi value genes)
  python pi_plot.py -p Pi.txt -g sequence.gff3 --lower-threshold 0.01

  # Use absolute threshold
  python pi_plot.py -p Pi.txt -g sequence.gff3 --absolute-threshold 0.05

  # Use relative significance and absolute threshold simultaneously
  python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 2.0 --absolute-threshold 0.05

  # Use lower threshold and absolute threshold simultaneously (half to absolute threshold, half to lower threshold)
  python pi_plot.py -p Pi.txt -g sequence.gff3 --lower-threshold 0.01 --absolute-threshold 0.05

  # Use relative threshold and lower threshold (output relative threshold results only)
  python pi_plot.py -p Pi.txt -g sequence.gff3 --significance-std 2.0 --lower-threshold 0.01

  # Specify output file and font size
  python pi_plot.py -p Pi.txt -g sequence.gff3 -o my_plot.png --gene-fontsize 14

  # Adjust significance parameters
  python pi_plot.py -p Pi.txt -g sequence.gff3 --window-size 15000 --significance-std 2.5

  # Limit number of annotated genes
  python pi_plot.py -p Pi.txt -g sequence.gff3 --max-genes 10

Note:
  Relative significance: Gene Pi value is higher than surrounding region average by a certain multiple of standard deviation
  Absolute significance: Gene Pi value is higher than set absolute threshold
  Lower threshold: Gene Pi value is lower than set lower threshold (lowest Pi value genes)

  Parameter logic:
  - Lower threshold only: Output lower threshold results only, sorted by Pi value ascending (no Z-score calculation)
  - Absolute threshold only: Output absolute threshold results, sorted by Pi value descending (no Z-score calculation)
  - Relative threshold only (default): Output Z-score results only
  - Lower threshold + Absolute threshold: For odd numbers, remainder to absolute threshold, half to absolute threshold, half to lower threshold (no Z-score calculation)
  - Other combinations containing relative threshold: Output relative threshold results only (Z-score)

  Genes meeting any condition will be annotated
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
                       help='Gene annotation font size (default: 12)')

    parser.add_argument('--max-genes',
                       type=int,
                       default=15,
                       help='Maximum number of genes to annotate (default: 15)')

    parser.add_argument('--window-size',
                       type=int,
                       default=10000,
                       help='Window size for calculating local average, unit bp (default: 10000)')

    parser.add_argument('--significance-std',
                       type=float,
                       default=2.0,
                       help='Relative significance standard, gene Pi value needs to be higher than surrounding region average by standard deviation multiplier (default: 2.0)')

    parser.add_argument('--absolute-threshold',
                       type=float,
                       default=None,
                       help='Absolute Pi value threshold, genes above this value are considered significant (optional, suggested value like 0.05)')

    parser.add_argument('--lower-threshold',
                       type=float,
                       default=None,
                       help='Lower Pi value threshold, genes below this value are considered significant (optional, used to detect lowest Pi value genes)')

    args = parser.parse_args()

    # Check if significance_std parameter is explicitly set
    # Get original command line arguments to detect if user explicitly specified --significance-std
    import sys
    significance_std_explicit = '--significance-std' in sys.argv

    # If user explicitly set --significance-std, enable relative threshold logic
    # Otherwise ignore default value of significance_std
    if not significance_std_explicit:
        # User did not explicitly set, set significance_std to None, indicating no relative threshold
        args.significance_std = None

    # Check if files exist
    if not os.path.exists(args.pi_file):
        print(f"Error: Pi file does not exist: {args.pi_file}")
        return 1

    if not os.path.exists(args.gff3_file):
        print(f"Error: GFF3 file does not exist: {args.gff3_file}")
        return 1

    # Print parameter information
    print("=" * 70)
    print("Gene Diversity Sliding Window Plotting Program")
    print("=" * 70)
    print(f"Pi file: {args.pi_file}")
    print(f"GFF3 file: {args.gff3_file}")
    print(f"Output file: {args.output}")
    print(f"Gene annotation font size: {args.gene_fontsize}")
    print(f"Maximum genes to annotate: {args.max_genes}")
    print(f"Comparison window size: {args.window_size} bp")
    print(f"Relative significance standard: {args.significance_std} standard deviations")
    if args.absolute_threshold is not None:
        print(f"Absolute significance threshold: {args.absolute_threshold}")
    if args.lower_threshold is not None:
        print(f"Lower threshold: {args.lower_threshold}")
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
            absolute_threshold=args.absolute_threshold,
            lower_threshold=args.lower_threshold
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
