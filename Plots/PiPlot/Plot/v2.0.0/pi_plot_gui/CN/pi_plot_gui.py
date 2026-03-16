#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因动态性滑动窗口绘图程序 - 图形化界面版本
支持矢量和PDF格式输出
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import re
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False


class PiPlotProcessor:
    """数据处理和绘图核心类"""
    
    @staticmethod
    def parse_pi_file(pi_file):
        """解析Pi文件，提取窗口位置和Pi值"""
        data = []
        with open(pi_file, 'r', encoding='utf-8') as f:
            next(f)  # 跳过第一行标题
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < 4:
                    continue

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

    @staticmethod
    def parse_gff3_file(gff3_file):
        """解析GFF3文件，提取基因位置信息"""
        genes = []

        with open(gff3_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                if parts[2] != 'gene':
                    continue

                seqid = parts[0]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                gene_name = None
                match = re.search(r'Name=([^;]+)', attributes)
                if match:
                    gene_name = match.group(1)
                else:
                    match = re.search(r'gene=([^;]+)', attributes)
                    if match:
                        gene_name = match.group(1)

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

    @staticmethod
    def find_high_pi_windows(pi_df, threshold_percentile=90):
        """找出高Pi值的窗口"""
        threshold = np.percentile(pi_df['pi'], threshold_percentile)
        high_pi_windows = pi_df[pi_df['pi'] >= threshold].copy()
        high_pi_windows = high_pi_windows.sort_values('pi', ascending=False)
        return high_pi_windows, threshold

    @staticmethod
    def find_significant_genes(pi_df, genes_df, window_size=10000, significance_std=2.0,
                               absolute_threshold=None, lower_threshold=None,
                               compute_zscore=True):
        """找到相对于其附近区域显著的基因"""
        significant_genes = []

        for _, gene in genes_df.iterrows():
            gene_mid = gene['midpoint']
            gene_start = gene['start']
            gene_end = gene['end']

            gene_pi_values = pi_df[
                (pi_df['window_start'] <= gene_end) &
                (pi_df['window_end'] >= gene_start)
            ]['pi']

            if len(gene_pi_values) == 0:
                continue

            gene_pi_mean = gene_pi_values.mean()

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

                    if region_std > 0:
                        z_score = (gene_pi_mean - region_mean) / region_std

            is_significant = False
            significance_type = "none"

            meets_relative = compute_zscore and z_score >= significance_std
            meets_absolute = absolute_threshold is not None and gene_pi_mean >= absolute_threshold
            meets_lower = lower_threshold is not None and gene_pi_mean <= lower_threshold

            if meets_relative and meets_absolute and meets_lower:
                is_significant = True
                significance_type = "lower+both"
            elif meets_relative and meets_absolute:
                is_significant = True
                significance_type = "both"
            elif meets_relative and meets_lower:
                is_significant = True
                significance_type = "lower+relative"
            elif meets_absolute and meets_lower:
                is_significant = True
                significance_type = "lower+absolute"
            elif meets_relative:
                is_significant = True
                significance_type = "relative"
            elif meets_absolute:
                is_significant = True
                significance_type = "absolute"
            elif meets_lower:
                is_significant = True
                significance_type = "lower"

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

        significant_genes.sort(key=lambda x: (x['z_score'], x['pi_value']), reverse=True)
        return significant_genes

    @staticmethod
    def plot_pi_with_genes(pi_file, gff3_file, output_file='pi_plot',
                          gene_fontsize=12, max_genes=15, window_size=10000,
                          significance_std=2.0, absolute_threshold=None, 
                          lower_threshold=None, output_format='png'):
        """绘制Pi滑动窗口图并标注显著基因"""
        
        pi_df = PiPlotProcessor.parse_pi_file(pi_file)
        genes_df = PiPlotProcessor.parse_gff3_file(gff3_file)

        high_pi_windows, pi_threshold = PiPlotProcessor.find_high_pi_windows(pi_df, threshold_percentile=90)

        use_relative_explicit = (significance_std is not None)
        use_absolute_threshold = (absolute_threshold is not None)
        use_lower_threshold = (lower_threshold is not None)

        if use_relative_explicit:
            compute_zscore = True
            output_mode = "relative_priority"
            significant_genes = PiPlotProcessor.find_significant_genes(pi_df, genes_df,
                                                window_size=window_size,
                                                significance_std=significance_std,
                                                absolute_threshold=absolute_threshold,
                                                lower_threshold=lower_threshold,
                                                compute_zscore=compute_zscore)
            
            def priority_score(gene):
                sig_type = gene['significance_type']
                if sig_type == 'lower+both':
                    return (6, gene['z_score'], gene['pi_value'])
                elif sig_type == 'both':
                    return (5, gene['z_score'], gene['pi_value'])
                elif sig_type == 'lower+relative':
                    return (4, gene['z_score'], gene['pi_value'])
                elif sig_type == 'lower+absolute':
                    return (3, gene['z_score'], gene['pi_value'])
                elif sig_type == 'relative':
                    return (2, gene['z_score'], gene['pi_value'])
                elif sig_type == 'absolute':
                    return (1, gene['z_score'], gene['pi_value'])
                elif sig_type == 'lower':
                    return (0, gene['z_score'], gene['pi_value'])
                else:
                    return (0, 0, 0)
            
            significant_genes.sort(key=priority_score, reverse=True)
            genes_to_label = significant_genes[:max_genes]

        elif use_lower_threshold and not use_absolute_threshold:
            compute_zscore = False
            output_mode = "lower_only"
            significant_genes = PiPlotProcessor.find_significant_genes(pi_df, genes_df,
                                                window_size=window_size,
                                                significance_std=significance_std,
                                                absolute_threshold=absolute_threshold,
                                                lower_threshold=lower_threshold,
                                                compute_zscore=compute_zscore)
            lower_genes = [g for g in significant_genes if 'lower' in g['significance_type']]
            lower_genes.sort(key=lambda x: x['pi_value'])
            genes_to_label = lower_genes[:max_genes]

        elif use_absolute_threshold and not use_lower_threshold:
            compute_zscore = False
            output_mode = "absolute_only"

            significant_genes = []
            for _, gene in genes_df.iterrows():
                gene_start = gene['start']
                gene_end = gene['end']

                gene_pi_values = pi_df[
                    (pi_df['window_start'] <= gene_end) &
                    (pi_df['window_end'] >= gene_start)
                ]['pi']

                if len(gene_pi_values) == 0:
                    continue

                gene_pi_mean = gene_pi_values.mean()

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

            significant_genes.sort(key=lambda x: x['pi_value'], reverse=True)
            genes_to_label = significant_genes[:max_genes]

        elif use_lower_threshold and use_absolute_threshold:
            compute_zscore = False
            output_mode = "lower+absolute"
            significant_genes = PiPlotProcessor.find_significant_genes(pi_df, genes_df,
                                                window_size=window_size,
                                                significance_std=significance_std,
                                                absolute_threshold=absolute_threshold,
                                                lower_threshold=lower_threshold,
                                                compute_zscore=compute_zscore)
            lower_genes = [g for g in significant_genes if 'lower' in g['significance_type']]
            lower_genes.sort(key=lambda x: x['pi_value'])
            absolute_genes = [g for g in significant_genes if 'absolute' in g['significance_type']]
            absolute_genes.sort(key=lambda x: x['pi_value'], reverse=True)

            total_genes = len(significant_genes)
            genes_limit = min(total_genes, max_genes)
            absolute_count = (genes_limit + 1) // 2
            lower_count = genes_limit // 2

            genes_to_label = absolute_genes[:absolute_count] + lower_genes[:lower_count]

        else:
            compute_zscore = True
            output_mode = "relative_only"
            significant_genes = PiPlotProcessor.find_significant_genes(pi_df, genes_df,
                                                window_size=window_size,
                                                significance_std=significance_std if significance_std is not None else 2.0,
                                                absolute_threshold=absolute_threshold,
                                                lower_threshold=lower_threshold,
                                                compute_zscore=compute_zscore)
            relative_genes = [g for g in significant_genes if 'relative' in g['significance_type']]
            relative_genes.sort(key=lambda x: (x['z_score'], x['pi_value']), reverse=True)
            genes_to_label = relative_genes[:max_genes]

        fig, ax = plt.subplots(figsize=(20, 12))

        ax.plot(pi_df['midpoint'], pi_df['pi'], 'b-', linewidth=1.5, alpha=0.7, label='Pi value')

        if absolute_threshold is not None:
            ax.axhline(y=absolute_threshold, color='purple', linestyle='--', linewidth=1.5,
                       label=f'Absolute Threshold ({absolute_threshold:.4f})')

        if lower_threshold is not None:
            ax.axhline(y=lower_threshold, color='green', linestyle='--', linewidth=1.5,
                       label=f'Lower Threshold ({lower_threshold:.4f})')

        ax.axhline(y=pi_threshold, color='r', linestyle='--', linewidth=1,
                   label=f'Relative Threshold (90th: {pi_threshold:.4f})')

        max_pi = pi_df['pi'].max()
        arrow_length = max_pi * 0.05
        label_spacing = max_pi * 0.15
        min_horizontal_distance = max_pi * 0.05

        genes_sorted = sorted(genes_to_label, key=lambda x: x['midpoint'])

        if genes_sorted:
            adjacent_gene_distance = 15000

            gene_groups = []
            current_group = [genes_sorted[0]]

            for gene in genes_sorted[1:]:
                if gene['midpoint'] - current_group[-1]['midpoint'] <= adjacent_gene_distance:
                    current_group.append(gene)
                else:
                    gene_groups.append(current_group)
                    current_group = [gene]

            if current_group:
                gene_groups.append(current_group)

            label_positions = []

            for group in gene_groups:
                group_size = len(group)
                base_height = max_pi * 0.03

                for i, gene_info in enumerate(group):
                    gene_name = gene_info['gene_name']
                    gene_mid = gene_info['midpoint']
                    gene_pi = gene_info['pi_value']
                    z_score = gene_info['z_score']
                    sig_type = gene_info['significance_type']

                    if 'lower' in sig_type:
                        marker_color = 'green'
                    elif sig_type == 'absolute':
                        marker_color = 'purple'
                    elif sig_type == 'both':
                        marker_color = 'darkred'
                    else:
                        marker_color = 'red'
                    
                    ax.plot(gene_mid, gene_pi, 'o', color=marker_color, markersize=5, zorder=5, 
                           alpha=0.8, markeredgecolor='black', markeredgewidth=1)

                    layer_index = i % group_size
                    layer_height = base_height + layer_index * label_spacing
                    label_y = gene_pi + layer_height

                    label_width = len(gene_name) * max_pi * 0.008
                    current_min_horizontal_distance = max(min_horizontal_distance, label_width + max_pi * 0.02)

                    final_y = label_y
                    final_x = gene_mid

                    max_iterations = 50
                    iteration = 0

                    while iteration < max_iterations:
                        has_overlap = False

                        for (prev_x, prev_y, prev_width, prev_y_bottom) in label_positions:
                            if abs(final_x - prev_x) < current_min_horizontal_distance:
                                current_y_bottom = gene_pi
                                current_y_top = final_y
                                prev_y_bottom_prev = prev_y_bottom
                                prev_y_top = prev_y

                                if not (current_y_top < prev_y_bottom_prev or current_y_bottom > prev_y_top):
                                    final_y = prev_y_top + label_spacing
                                    has_overlap = True
                                    break

                        if not has_overlap:
                            break

                        iteration += 1

                    label_positions.append((final_x, final_y, label_width, gene_pi))

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
                    else:
                        label_color = 'yellow'
                        edge_color = 'red'
                    
                    label_text = f'{gene_name}'

                    ax.annotate('', xy=(gene_mid, gene_pi), xytext=(gene_mid, final_y),
                               arrowprops=dict(arrowstyle='-',
                                             color=edge_color,
                                             linewidth=1.5,
                                             alpha=0.7),
                               zorder=5)

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

        ax.set_xlabel('Sequence Position (bp)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Pi (Nucleotide Diversity)', fontsize=14, fontweight='bold')

        ax.set_title('Sliding Window Analysis: Nucleotide Diversity (Pi) Along Genome\n'
                     f'Markers indicate genes with elevated diversity',
                     fontsize=16, fontweight='bold', pad=30)

        ax.grid(True, alpha=0.3, linestyle=':')
        ax.legend(loc='upper right', fontsize=11, framealpha=0.9, bbox_to_anchor=(0.98, 0.98))

        ax.set_xlim(pi_df['midpoint'].min() - 20000, pi_df['midpoint'].max() + 20000)
        
        if label_positions:
            max_label_y = max(pos[1] for pos in label_positions)
            ax.set_ylim(0, max_label_y * 1.25)
        else:
            ax.set_ylim(0, max_pi * 1.3)

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

            fig.text(0.92, 0.5, info_text,
                    fontsize=9,
                    bbox=dict(boxstyle='round',
                             facecolor='lightyellow',
                             alpha=0.95,
                             edgecolor='orange',
                             linewidth=2),
                    verticalalignment='center',
                    family='monospace')

        plt.tight_layout()
        plt.subplots_adjust(right=0.88)

        output_files = []
        
        if output_format == 'png':
            png_file = f"{output_file}.png"
            plt.savefig(png_file, dpi=300, bbox_inches='tight')
            output_files.append(('PNG', png_file))
        elif output_format == 'svg':
            svg_file = f"{output_file}.svg"
            plt.savefig(svg_file, format='svg', bbox_inches='tight', dpi=300)
            output_files.append(('SVG', svg_file))
        elif output_format == 'pdf':
            pdf_file = f"{output_file}.pdf"
            plt.savefig(pdf_file, format='pdf', bbox_inches='tight')
            output_files.append(('PDF', pdf_file))
        elif output_format == 'all':
            png_file = f"{output_file}.png"
            plt.savefig(png_file, dpi=300, bbox_inches='tight')
            output_files.append(('PNG', png_file))
            
            svg_file = f"{output_file}.svg"
            plt.savefig(svg_file, format='svg', bbox_inches='tight', dpi=300)
            output_files.append(('SVG', svg_file))
            
            pdf_file = f"{output_file}.pdf"
            plt.savefig(pdf_file, format='pdf', bbox_inches='tight')
            output_files.append(('PDF', pdf_file))

        output_info = f"{output_file}_info.txt"
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

        plt.close()

        return output_files, len(significant_genes), genes_to_label


class PiPlotGUI:
    """图形化界面类"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("基因动态性滑动窗口绘图程序 v2.0")
        self.root.geometry("1100x750")
        
        self.pi_file_path = tk.StringVar()
        self.gff3_file_path = tk.StringVar()
        self.output_dir = tk.StringVar()
        
        self.setup_ui()
        
    def setup_ui(self):
        """设置界面布局"""
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        # 主选项卡
        main_frame = ttk.Frame(notebook)
        notebook.add(main_frame, text='主界面')
        self.setup_main_tab(main_frame)
        
        # 帮助选项卡
        help_frame = ttk.Frame(notebook)
        notebook.add(help_frame, text='使用帮助')
        self.setup_help_tab(help_frame)
        
    def setup_main_tab(self, parent):
        """设置主选项卡"""
        # 左侧控制面板
        control_frame = ttk.LabelFrame(parent, text="参数设置", padding="10")
        control_frame.pack(side='left', fill='both', expand=False, padx=5, pady=5)
        
        # 文件选择区域
        file_frame = ttk.LabelFrame(control_frame, text="文件选择", padding="10")
        file_frame.pack(fill='x', pady=5)
        
        ttk.Label(file_frame, text="Pi文件:").grid(row=0, column=0, sticky='w', pady=2)
        ttk.Entry(file_frame, textvariable=self.pi_file_path, width=35).grid(row=0, column=1, pady=2, padx=5)
        ttk.Button(file_frame, text="浏览...", command=self.browse_pi_file).grid(row=0, column=2, pady=2)
        
        ttk.Label(file_frame, text="GFF3文件:").grid(row=1, column=0, sticky='w', pady=2)
        ttk.Entry(file_frame, textvariable=self.gff3_file_path, width=35).grid(row=1, column=1, pady=2, padx=5)
        ttk.Button(file_frame, text="浏览...", command=self.browse_gff3_file).grid(row=1, column=2, pady=2)
        
        ttk.Label(file_frame, text="输出目录:").grid(row=2, column=0, sticky='w', pady=2)
        ttk.Entry(file_frame, textvariable=self.output_dir, width=35).grid(row=2, column=1, pady=2, padx=5)
        ttk.Button(file_frame, text="浏览...", command=self.browse_output_dir).grid(row=2, column=2, pady=2)
        
        # 输出格式选择
        format_frame = ttk.LabelFrame(control_frame, text="输出格式", padding="10")
        format_frame.pack(fill='x', pady=5)
        
        self.output_format = tk.StringVar(value='png')
        ttk.Radiobutton(format_frame, text="PNG格式", variable=self.output_format, value='png').pack(anchor='w', pady=2)
        ttk.Radiobutton(format_frame, text="SVG格式 (矢量图)", variable=self.output_format, value='svg').pack(anchor='w', pady=2)
        ttk.Radiobutton(format_frame, text="PDF格式", variable=self.output_format, value='pdf').pack(anchor='w', pady=2)
        ttk.Radiobutton(format_frame, text="全部格式 (PNG+SVG+PDF)", variable=self.output_format, value='all').pack(anchor='w', pady=2)
        
        # 参数设置区域
        param_frame = ttk.LabelFrame(control_frame, text="绘图参数", padding="10")
        param_frame.pack(fill='x', pady=5)
        
        # 输出文件名
        ttk.Label(param_frame, text="输出文件名:").grid(row=0, column=0, sticky='w', pady=2)
        self.output_name = ttk.Entry(param_frame, width=25)
        self.output_name.insert(0, 'pi_sliding_window_plot')
        self.output_name.grid(row=0, column=1, pady=2)
        
        # 字体大小
        ttk.Label(param_frame, text="字体大小:").grid(row=1, column=0, sticky='w', pady=2)
        self.gene_fontsize = ttk.Spinbox(param_frame, from_=8, to=18, width=10)
        self.gene_fontsize.set(12)
        self.gene_fontsize.grid(row=1, column=1, sticky='w', pady=2)
        
        # 最大基因数
        ttk.Label(param_frame, text="最大基因数:").grid(row=2, column=0, sticky='w', pady=2)
        self.max_genes = ttk.Spinbox(param_frame, from_=5, to=50, width=10)
        self.max_genes.set(15)
        self.max_genes.grid(row=2, column=1, sticky='w', pady=2)
        
        # 窗口大小
        ttk.Label(param_frame, text="窗口大小(bp):").grid(row=3, column=0, sticky='w', pady=2)
        self.window_size = ttk.Spinbox(param_frame, from_=1000, to=100000, increment=1000, width=10)
        self.window_size.set(10000)
        self.window_size.grid(row=3, column=1, sticky='w', pady=2)
        
        # 显著性标准
        ttk.Label(param_frame, text="相对显著性:").grid(row=4, column=0, sticky='w', pady=2)
        self.significance_std = ttk.Spinbox(param_frame, from_=0.5, to=5.0, increment=0.1, width=10)
        self.significance_std.set(2.0)
        self.significance_std.grid(row=4, column=1, sticky='w', pady=2)
        
        # 绝对阈值
        ttk.Label(param_frame, text="绝对阈值:").grid(row=5, column=0, sticky='w', pady=2)
        self.absolute_threshold = ttk.Entry(param_frame, width=10)
        self.absolute_threshold.insert(0, '')
        self.absolute_threshold.grid(row=5, column=1, sticky='w', pady=2)
        
        # 下限阈值
        ttk.Label(param_frame, text="下限阈值:").grid(row=6, column=0, sticky='w', pady=2)
        self.lower_threshold = ttk.Entry(param_frame, width=10)
        self.lower_threshold.insert(0, '')
        self.lower_threshold.grid(row=6, column=1, sticky='w', pady=2)
        
        # 操作按钮
        button_frame = ttk.Frame(control_frame)
        button_frame.pack(fill='x', pady=10)
        
        ttk.Button(button_frame, text="开始分析", command=self.start_analysis).pack(side='left', padx=5)
        ttk.Button(button_frame, text="重置参数", command=self.reset_params).pack(side='left', padx=5)
        
        # 进度条
        self.progress = ttk.Progressbar(control_frame, mode='indeterminate')
        self.progress.pack(fill='x', pady=10)
        
        # 状态标签
        self.status_label = ttk.Label(control_frame, text="准备就绪", wraplength=300)
        self.status_label.pack(fill='x', pady=5)
        
        # 右侧日志区域
        log_frame = ttk.LabelFrame(parent, text="运行日志", padding="10")
        log_frame.pack(side='right', fill='both', expand=True, padx=5, pady=5)
        
        self.log_text = scrolledtext.ScrolledText(log_frame, wrap='word', height=30)
        self.log_text.pack(fill='both', expand=True)
        
    def setup_help_tab(self, parent):
        """设置帮助选项卡"""
        help_text = """
基因动态性滑动窗口绘图程序 - 使用说明
========================================

一、程序功能
-----------
本程序用于分析基因组核苷酸多样性（Pi）的滑动窗口数据，并自动标注相对于其周围区域显著升高的基因。

主要功能：
• 解析 Pi 滑动窗口数据文件
• 解析 GFF3 基因注释文件
• 识别相对于周围区域显著升高的基因
• 支持三种显著性判断标准
• 智能标签布局，自动避免重叠
• 生成 PNG、SVG（矢量图）和 PDF 格式的图表

二、输入文件格式
---------------

1. Pi 文件格式：
   Window    Midpoint    Pi    Theta
   1-10000   5000        0.0234 0.0256
   10001-20000  15000    0.0287 0.0301

2. GFF3 文件格式：
   ##gff-version 3
   chr1    Ensembl    gene    1000    5000    .    +    .    ID=gene-001;Name=GENE1

三、参数说明
-----------

必填参数：
• Pi文件：Pi滑动窗口数据文件路径
• GFF3文件：基因注释文件路径
• 输出目录：结果文件保存位置

可选参数：
• 输出格式：可选择PNG、SVG（矢量图）、PDF或全部格式
• 输出文件名：默认为 pi_sliding_window_plot
• 字体大小：基因标注字体大小（8-18）
• 最大基因数：最多标注的基因数量（5-50）
• 窗口大小：计算局部平均值的窗口大小，单位bp（1000-100000）
• 相对显著性：相对显著性标准，标准差倍数（0.5-5.0）
• 绝对阈值：绝对Pi值阈值（可选）
• 下限阈值：下限Pi值阈值，用于检测最低Pi值基因（可选）

四、显著性标准说明
-----------------

1. 相对显著性：基因Pi值高于周围区域平均值一定倍数标准差
   • 值越小：标注的基因越多
   • 值越大：仅标注最显著的基因
   • 推荐值：1.5-3.0

2. 绝对显著性：基因Pi值高于设定的绝对阈值
   • 适合关注特定Pi值范围的基因
   • 建议取Pi值分布的75-90分位数

3. 下限阈值：基因Pi值低于设定的下限阈值
   • 用于检测基因组中Pi值最低的基因
   • 适合寻找保守区域

五、输出文件
-----------

程序将生成以下文件：
1. 图表文件（PNG/SVG/PDF）：显示Pi滑动窗口曲线和显著基因标注
2. 信息文本文件（_info.txt）：包含显著基因的详细分析结果

六、输出格式对比
---------------

• PNG格式：位图格式，适合屏幕查看
• SVG格式：矢量图格式，无损缩放，适合打印和出版
• PDF格式：矢量格式，适合文档嵌入和打印

建议：
• 屏幕查看：使用PNG格式
• 学术出版：使用SVG或PDF格式
• 多用途：选择"全部格式"

七、使用示例
-----------

场景1：快速浏览数据
• 使用默认参数运行
• 输出PNG格式即可

场景2：学术出版
• 调整参数获得最佳效果
• 选择SVG或PDF格式
• 适当调整字体大小和基因数量

场景3：检测局部高变异区域
• 使用较小的窗口（5000bp）
• 较低的显著性标准（1.5）
• 选择SVG格式获得高质量图表

八、常见问题
-----------

Q: 为什么有些高Pi值的基因没有被标注？
A: 程序基于相对于周围区域的显著性进行标注，而非绝对Pi值。可以：
   • 降低相对显著性参数
   • 增加绝对阈值参数
   • 减小窗口大小

Q: 基因标签重叠怎么办？
A: 可以：
   • 减少最大基因数
   • 减小字体大小
   • 程序已内置智能布局算法

Q: 如何获得更清晰的图表？
A: 建议：
   • 减少标注基因数量（10-15个）
   • 调整字体大小（10-12）
   • 选择SVG或PDF矢量格式

九、参数调优建议
---------------

窗口大小：
• 小窗口（5000-10000bp）：检测局部变化
• 中窗口（10000-20000bp）：平衡局部和全局特征（推荐）
• 大窗口（20000-50000bp）：检测大尺度趋势

显著性标准：
• 1.5：宽松标准，标注更多基因
• 2.0：标准平衡（默认）
• 2.5-3.0：严格标准，仅标注最显著基因

十、联系支持
-----------

如有问题或建议，请联系：myzhang0726@foxmail.com
        """
        
        text_widget = scrolledtext.ScrolledText(parent, wrap='word', height=40)
        text_widget.pack(fill='both', expand=True, padx=10, pady=10)
        text_widget.insert('1.0', help_text)
        text_widget.config(state='disabled')
        
    def browse_pi_file(self):
        """浏览选择Pi文件"""
        filename = filedialog.askopenfilename(
            title="选择Pi文件",
            filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")]
        )
        if filename:
            self.pi_file_path.set(filename)
            self.log(f"已选择Pi文件: {filename}")
            
    def browse_gff3_file(self):
        """浏览选择GFF3文件"""
        filename = filedialog.askopenfilename(
            title="选择GFF3文件",
            filetypes=[("GFF3文件", "*.gff3"), ("GFF文件", "*.gff"), ("所有文件", "*.*")]
        )
        if filename:
            self.gff3_file_path.set(filename)
            self.log(f"已选择GFF3文件: {filename}")
            
    def browse_output_dir(self):
        """浏览选择输出目录"""
        dirname = filedialog.askdirectory(title="选择输出目录")
        if dirname:
            self.output_dir.set(dirname)
            self.log(f"已选择输出目录: {dirname}")
            
    def log(self, message):
        """添加日志信息"""
        self.log_text.insert('end', message + '\n')
        self.log_text.see('end')
        self.root.update_idletasks()
        
    def reset_params(self):
        """重置参数"""
        self.output_name.delete(0, 'end')
        self.output_name.insert(0, 'pi_sliding_window_plot')
        self.gene_fontsize.set(12)
        self.max_genes.set(15)
        self.window_size.set(10000)
        self.significance_std.set(2.0)
        self.absolute_threshold.delete(0, 'end')
        self.lower_threshold.delete(0, 'end')
        self.output_format.set('png')
        self.log("参数已重置")
        
    def validate_params(self):
        """验证参数"""
        if not self.pi_file_path.get():
            messagebox.showerror("错误", "请选择Pi文件")
            return False
        
        if not self.gff3_file_path.get():
            messagebox.showerror("错误", "请选择GFF3文件")
            return False
            
        if not os.path.exists(self.pi_file_path.get()):
            messagebox.showerror("错误", "Pi文件不存在")
            return False
            
        if not os.path.exists(self.gff3_file_path.get()):
            messagebox.showerror("错误", "GFF3文件不存在")
            return False
            
        if not self.output_dir.get():
            messagebox.showerror("错误", "请选择输出目录")
            return False
            
        return True
        
    def start_analysis(self):
        """开始分析"""
        if not self.validate_params():
            return
            
        self.log_text.delete('1.0', 'end')
        self.log("开始分析...")
        self.progress.start()
        self.status_label.config(text="正在分析...")
        
        # 在新线程中运行分析，避免阻塞界面
        thread = threading.Thread(target=self.run_analysis)
        thread.daemon = True
        thread.start()
        
    def run_analysis(self):
        """运行分析（在新线程中）"""
        try:
            pi_file = self.pi_file_path.get()
            gff3_file = self.gff3_file_path.get()
            output_dir = self.output_dir.get()
            output_name = self.output_name.get()
            output_format = self.output_format.get()
            
            gene_fontsize = int(self.gene_fontsize.get())
            max_genes = int(self.max_genes.get())
            window_size = int(self.window_size.get())
            significance_std = float(self.significance_std.get())
            
            absolute_threshold = self.absolute_threshold.get()
            absolute_threshold = float(absolute_threshold) if absolute_threshold else None
            
            lower_threshold = self.lower_threshold.get()
            lower_threshold = float(lower_threshold) if lower_threshold else None
            
            self.log(f"分析参数:")
            self.log(f"  Pi文件: {pi_file}")
            self.log(f"  GFF3文件: {gff3_file}")
            self.log(f"  输出目录: {output_dir}")
            self.log(f"  输出格式: {output_format.upper()}")
            self.log(f"  字体大小: {gene_fontsize}")
            self.log(f"  最大基因数: {max_genes}")
            self.log(f"  窗口大小: {window_size} bp")
            self.log(f"  相对显著性: {significance_std}")
            if absolute_threshold:
                self.log(f"  绝对阈值: {absolute_threshold}")
            if lower_threshold:
                self.log(f"  下限阈值: {lower_threshold}")
            self.log("")
            
            output_file = os.path.join(output_dir, output_name)
            
            self.log("正在处理数据...")
            output_files, num_genes, genes_to_label = PiPlotProcessor.plot_pi_with_genes(
                pi_file, gff3_file, output_file,
                gene_fontsize, max_genes, window_size,
                significance_std, absolute_threshold, lower_threshold,
                output_format
            )
            
            self.log("")
            self.log("分析完成!")
            self.log(f"找到 {num_genes} 个显著基因")
            self.log(f"标注了 {len(genes_to_label)} 个基因")
            self.log("")
            
            self.log("生成的文件:")
            for fmt, filepath in output_files:
                self.log(f"  [{fmt}] {filepath}")
            
            info_file = f"{output_file}_info.txt"
            if os.path.exists(info_file):
                self.log(f"  [INFO] {info_file}")
            
            self.log("")
            self.log("显著基因列表:")
            for i, gene in enumerate(genes_to_label[:15], 1):
                self.log(f"  {i}. {gene['gene_name']} (Pi={gene['pi_value']:.4f}, Z={gene['z_score']:.2f})")
            
            self.root.after(0, lambda: self.status_label.config(text="分析完成!"))
            self.root.after(0, lambda: messagebox.showinfo("完成", f"分析完成!\n找到 {num_genes} 个显著基因\n生成 {len(output_files)} 个文件"))
            
        except Exception as e:
            error_msg = f"分析出错: {str(e)}"
            self.log(error_msg)
            self.root.after(0, lambda: self.status_label.config(text="分析失败"))
            self.root.after(0, lambda: messagebox.showerror("错误", error_msg))
            
        finally:
            self.progress.stop()


def main():
    """主函数"""
    root = tk.Tk()
    app = PiPlotGUI(root)
    root.mainloop()


if __name__ == '__main__':
    main()
