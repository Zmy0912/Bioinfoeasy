import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def plot_codon_usage_table_and_chart(excel_file_path, output_path='codon_usage_chart.png', colormap_name='tab20c', show_codon_labels=True, figsize=(18, 13), legend_block_height=0.12, legend_block_width=0.6, show_legend_blocks=True, codon_label_size=10, xlabel_size=14, ylabel_size=14, tick_label_size=12, ytick_label_size=12, legend_text_size=8, border_width=1.0, legend_bottom_offset=0.22):
    """
    绘制密码子使用频率的组合图：上方为表格，下方为堆积柱状图

    参数:
        excel_file_path: Excel文件路径
        output_path: 输出图片路径
        colormap_name: 配色方案名称，可选：'tab20c', 'tab20b', 'tab20', 'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'viridis', 'plasma', 'inferno'
        show_codon_labels: 是否在柱状图内显示密码子名称，True为显示，False为不显示
        figsize: 画布大小，格式为 (width, height)，默认为 (18, 13)
        legend_block_height: 密码子堆叠块的高度，默认为 0.12
        legend_block_width: 密码子堆叠块的宽度，默认为 0.6
        show_legend_blocks: 是否显示密码子堆叠块图例，True为显示，False为不显示
        codon_label_size: 柱状图内密码子名称字体大小，默认为 10
        xlabel_size: 横坐标标签字体大小，默认为 14
        ylabel_size: 纵坐标标签字体大小，默认为 14
        tick_label_size: 横坐标轴刻度标签字体大小，默认为 12
        ytick_label_size: 纵坐标轴刻度标签字体大小，默认为 12
        legend_text_size: 堆叠块图例文字大小，默认为 8
        border_width: 图表边框和刻度线粗细，默认为 1.0
        legend_bottom_offset: 堆叠块距离X坐标轴的距离，默认为 0.22
    """

    # 读取Excel文件，尝试不同的引擎
    try:
        df = pd.read_excel(excel_file_path, engine='openpyxl')
    except:
        try:
            df = pd.read_excel(excel_file_path, engine='xlrd')
        except:
            df = pd.read_excel(excel_file_path)
    
    # 获取列名（适应不同的列名格式）
    print("Excel列名:", df.columns.tolist())

    # 检查是否有数据
    if len(df) == 0:
        raise ValueError("Excel文件为空，没有数据")

    # 尝试识别密码子、频数、RSCU列
    codon_col = None
    freq_col = None
    rscu_col = None
    aa_col = None

    for col in df.columns:
        col_lower = str(col).lower()
        col_str = str(col)
        if 'codon' in col_lower or '密码子' in col_str:
            codon_col = col
        elif 'freq' in col_lower or '频数' in col_str or 'count' in col_lower or '数量' in col_str:
            freq_col = col
        elif 'rscu' in col_lower or 'rscu' in col_str:
            rscu_col = col
        elif 'amino' in col_lower or 'aa' in col_lower or '氨基酸' in col_str:
            aa_col = col
    
    if codon_col is None:
        codon_col = df.columns[0]
    if freq_col is None:
        freq_col = df.columns[1]
    if rscu_col is None:
        rscu_col = df.columns[2]
    
    # 如果没有氨基酸列，从密码子推导
    if aa_col is None:
        # 确保codon列的数据是字符串类型
        df[codon_col] = df[codon_col].astype(str)
        df['AminoAcid'] = df[codon_col].apply(codon_to_aa)
        aa_col = 'AminoAcid'
    
    # 重命名列以便使用
    rename_dict = {
        codon_col: 'Codon',
        freq_col: 'Frequency',
        rscu_col: 'RSCU'
    }
    if aa_col is not None:
        rename_dict[aa_col] = 'AminoAcid'

    df = df.rename(columns=rename_dict)

    # 确保必要列存在
    if 'Codon' not in df.columns:
        raise ValueError(f"未找到密码子列，可用列名: {df.columns.tolist()}")
    if 'RSCU' not in df.columns and 'Frequency' not in df.columns:
        raise ValueError(f"未找到RSCU或频数列，可用列名: {df.columns.tolist()}")
    
    # 按氨基酸分组
    amino_order = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Tyr',
                   'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly']
    
    # 确保数据按氨基酸正确排序
    df['AminoAcid'] = pd.Categorical(df['AminoAcid'], categories=amino_order, ordered=True)
    df = df.sort_values('AminoAcid')

    # 创建图形
    fig = plt.figure(figsize=figsize)

    # ===== 堆积柱状图 =====
    ax_chart = fig.add_axes([0.08, 0.22, 0.85, 0.68])
    
    # 为每个氨基酸准备堆积数据
    bottom_values = {aa: 0 for aa in amino_order}
    codon_colors = {}
    position_colors = {}  # 存储每个位置的固定颜色

    # 为每个位置分配固定颜色（位置1、2、3...使用不同的颜色）
    # 使用用户指定的色板
    try:
        colormap = plt.get_cmap(colormap_name)
    except:
        print(f"警告：未找到色板 '{colormap_name}'，使用默认色板 'tab20c'")
        colormap = plt.get_cmap('tab20c')

    # 按氨基酸分组，为每个氨基酸的密码子按位置分配颜色
    for aa in amino_order:
        aa_data = df[df['AminoAcid'] == aa].reset_index(drop=True)
        for idx, row in aa_data.iterrows():
            codon = row['Codon']
            position = idx  # 密码子在氨基酸中的位置
            if position not in position_colors:
                position_colors[position] = colormap(position % 20)
            codon_colors[codon] = position_colors[position]

    # 绘制堆积柱状图
    x_positions = np.arange(len(amino_order))
    bar_width = 0.8

    # 存储每个氨基酸的密码子信息，用于标注
    aa_codons_info = {}

    for aa in amino_order:
        aa_data = df[df['AminoAcid'] == aa]
        if len(aa_data) == 0:
            continue

        aa_codons = aa_data['Codon'].tolist()
        aa_rscu = aa_data['RSCU'].tolist()

        aa_codons_info[aa] = {
            'codons': aa_codons,
            'rscu': aa_rscu
        }

        for codon, rscu in zip(aa_codons, aa_rscu):
            ax_chart.bar(aa, rscu, bar_width,
                        bottom=bottom_values[aa],
                        color=codon_colors[codon],
                        edgecolor='white',
                        linewidth=0.5)
            bottom_values[aa] += rscu
    
    # 设置坐标轴
    ax_chart.set_xlabel('Amino Acid', fontsize=xlabel_size)
    ax_chart.set_ylabel('RSCU', fontsize=ylabel_size)
    ax_chart.set_xticks(x_positions)
    ax_chart.set_xticklabels(amino_order, fontsize=tick_label_size)
    # 设置Y轴刻度标签大小
    ax_chart.tick_params(axis='y', labelsize=ytick_label_size)

    # 设置边框和刻度线粗细
    for spine in ax_chart.spines.values():
        spine.set_linewidth(border_width)
    ax_chart.tick_params(width=border_width)

    # 计算每个氨基酸的RSCU总和（堆积后的总高度）
    aa_rscu_sums = df.groupby('AminoAcid')['RSCU'].sum()
    max_total_rscu = aa_rscu_sums.max()

    # 根据堆积后的最大总高度设置纵坐标上限
    if max_total_rscu > 3.0:
        # 如果堆积后的最大值大于3.0，设置稍大一点的值
        y_upper_limit = max_total_rscu * 1.05  # 增加5%余量
    else:
        # 如果最大值小于等于3.0，默认为3.0
        y_upper_limit = 3.0
    ax_chart.set_ylim(0, y_upper_limit)

    # 在每个氨基酸柱子下方标注对应的密码子
    if show_codon_labels:
        for idx, aa in enumerate(amino_order):
            if aa in aa_codons_info:
                codons = aa_codons_info[aa]['codons']
                rscu_values = aa_codons_info[aa]['rscu']

                # 计算每个色块的中心位置
                bottom = 0
                for codon, rscu in zip(codons, rscu_values):
                    # 色块中心位置
                    center_y = bottom + rscu / 2

                    # 标注密码子名称
                    ax_chart.text(idx, center_y, codon,
                                ha='center', va='center',
                                fontsize=codon_label_size,
                                color='black' if codon_colors[codon] > (0.5, 0.5, 0.5) else 'white')

                    bottom += rscu
    
    # 添加网格线
    ax_chart.grid(True, axis='y', alpha=0.3, linestyle='--')
    ax_chart.set_axisbelow(True)

    # 添加标题
    fig.text(0.5, 0.02, '密码子使用频率分析 / Codon Usage Analysis',
            ha='center', fontsize=16, fontweight='bold')

    # ===== 在X轴下方添加密码子堆叠图例 =====
    if show_legend_blocks:
        ax_legend = fig.add_axes([0.08, 0.06, 0.85, 0.14], sharex=ax_chart)
        ax_legend.axis('off')
        ax_legend.set_xlim(ax_chart.get_xlim())

        # 设置固定的图例高度范围
        ax_legend.set_ylim(0, 1.0)

        # 为每个氨基酸创建纵向堆叠方块图例
        block_width = legend_block_width  # 使用可调整的宽度
        block_height = legend_block_height  # 使用可调整的高度
        y_start = 0.02  # 从底部开始

        for idx, aa in enumerate(amino_order):
            if aa in aa_codons_info:
                codons = aa_codons_info[aa]['codons']

                x_center = idx
                x_left = x_center - block_width / 2

                for j, codon in enumerate(codons):
                    y_bottom = y_start + j * block_height
                    color = codon_colors[codon]

                    # 绘制方块（纵向堆叠，高度相同）
                    rect = plt.Rectangle((x_left, y_bottom), block_width, block_height,
                                       facecolor=color, edgecolor='white', linewidth=0.8)
                    ax_legend.add_patch(rect)

                    # 添加密码子名称
                    ax_legend.text(x_center, y_bottom + block_height / 2, codon,
                                 ha='center', va='center',
                                 fontsize=legend_text_size,
                                 color='black' if color > (0.5, 0.5, 0.5) else 'white')

            # 调整柱状图的X轴标签位置，上移以避免与图例重叠
            ax_chart.xaxis.set_label_coords(0.5, -legend_bottom_offset)
    else:
        # 不显示堆叠块时，调整柱状图位置
        ax_chart.xaxis.set_label_coords(0.5, -0.05)

    # 调整布局
    fig.subplots_adjust(left=0.08, right=0.95, top=0.98, bottom=0.08)
    return fig, df, output_path

def codon_to_aa(codon):
    """密码子到氨基酸的映射（使用RNA碱基U）"""
    genetic_code = {
        # Phe
        'UUU': 'Phe', 'UUC': 'Phe',
        # Leu
        'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        # Ile
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        # Met
        'AUG': 'Met',
        # Val
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        # Ser
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
        # Pro
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        # Thr
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        # Ala
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        # Tyr
        'UAU': 'Tyr', 'UAC': 'Tyr',
        # Stop
        'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
        # His
        'CAU': 'His', 'CAC': 'His',
        # Gln
        'CAA': 'Gln', 'CAG': 'Gln',
        # Asn
        'AAU': 'Asn', 'AAC': 'Asn',
        # Lys
        'AAA': 'Lys', 'AAG': 'Lys',
        # Asp
        'GAU': 'Asp', 'GAC': 'Asp',
        # Glu
        'GAA': 'Glu', 'GAG': 'Glu',
        # Cys
        'UGU': 'Cys', 'UGC': 'Cys',
        # Trp
        'UGG': 'Trp',
        # Arg
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
        # Gly
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }
    return genetic_code.get(codon.upper(), 'Unknown')


class CodonUsageApp:
    """交互式密码子使用频率分析应用程序"""
    
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("密码子使用频率分析工具 / Codon Usage Analysis Tool")
        self.root.geometry("1200x800")
        self.root.minsize(1000, 700)

        self.excel_file_path = None
        self.fig = None
        self.df = None
        self.canvas = None
        self.lbl_file = None
        self.lbl_hint = None
        self.colormap_var = tk.StringVar(value='tab20c')
        self.show_codon_labels_var = tk.BooleanVar(value=True)
        self.show_legend_blocks_var = tk.BooleanVar(value=True)
        self.fig_width_var = tk.StringVar(value='18')
        self.fig_height_var = tk.StringVar(value='13')
        self.legend_block_height_var = tk.StringVar(value='0.12')
        self.legend_block_width_var = tk.StringVar(value='0.6')

        # 文本大小设置（取消加粗功能）
        self.codon_label_size_var = tk.StringVar(value='10')
        self.xlabel_size_var = tk.StringVar(value='14')
        self.ylabel_size_var = tk.StringVar(value='14')
        self.tick_label_size_var = tk.StringVar(value='12')
        self.ytick_label_size_var = tk.StringVar(value='12')
        self.legend_text_size_var = tk.StringVar(value='8')
        self.border_width_var = tk.StringVar(value='1.0')
        self.legend_bottom_offset_var = tk.StringVar(value='0.22')  # 堆叠块距离X坐标轴的距离
        
    def create_widgets(self):
        """创建界面组件"""
        # 顶部工具栏 - 第一行
        toolbar_row1 = tk.Frame(self.root, bd=1, relief=tk.RAISED)
        toolbar_row1.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        # 选择Excel文件按钮
        btn_load = tk.Button(toolbar_row1, text="选择Excel文件", command=self.load_excel_file,
                            bg='#4CAF50', fg='white', font=('Arial', 10, 'bold'),
                            padx=15, pady=5)
        btn_load.pack(side=tk.LEFT, padx=5)

        # 文件路径显示
        if self.lbl_file is None:
            self.lbl_file = tk.Label(toolbar_row1, text="未选择文件", fg='gray')
        self.lbl_file.pack(side=tk.LEFT, padx=10)

        # 分隔符
        tk.Frame(toolbar_row1, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 保存按钮
        btn_pdf = tk.Button(toolbar_row1, text="保存为PDF", command=self.save_as_pdf,
                           bg='#2196F3', fg='white', font=('Arial', 10, 'bold'),
                           padx=15, pady=5)
        btn_pdf.pack(side=tk.LEFT, padx=5)

        btn_svg = tk.Button(toolbar_row1, text="保存为SVG矢量图", command=self.save_as_svg,
                           bg='#FF9800', fg='white', font=('Arial', 10, 'bold'),
                           padx=15, pady=5)
        btn_svg.pack(side=tk.LEFT, padx=5)

        btn_png = tk.Button(toolbar_row1, text="保存为PNG", command=self.save_as_png,
                           bg='#9C27B0', fg='white', font=('Arial', 10, 'bold'),
                           padx=15, pady=5)
        btn_png.pack(side=tk.LEFT, padx=5)

        # 顶部工具栏 - 第二行
        toolbar_row2 = tk.Frame(self.root, bd=1, relief=tk.RAISED)
        toolbar_row2.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        # 配色方案选择
        tk.Label(toolbar_row2, text="配色方案:", font=('Arial', 10)).pack(side=tk.LEFT, padx=5)
        colormap_options = ['tab20c', 'tab20b', 'tab20', 'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'viridis', 'plasma', 'inferno']
        self.colormap_combo = ttk.Combobox(toolbar_row2, textvariable=self.colormap_var,
                                           values=colormap_options, state='readonly', width=10)
        self.colormap_combo.pack(side=tk.LEFT, padx=5)

        # 密码子标签显示选项
        self.chk_codon_labels = tk.Checkbutton(toolbar_row2, text="显示密码子名称",
                                               variable=self.show_codon_labels_var,
                                               font=('Arial', 10),
                                               command=self.on_codon_label_toggle)
        self.chk_codon_labels.pack(side=tk.LEFT, padx=10)

        # 堆叠块显示选项
        self.chk_legend_blocks = tk.Checkbutton(toolbar_row2, text="显示堆叠块",
                                               variable=self.show_legend_blocks_var,
                                               font=('Arial', 10),
                                               command=self.on_legend_blocks_toggle)
        self.chk_legend_blocks.pack(side=tk.LEFT, padx=10)

        # 分隔符
        tk.Frame(toolbar_row2, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 画布大小调整
        tk.Label(toolbar_row2, text="画布大小:", font=('Arial', 10)).pack(side=tk.LEFT, padx=2)
        tk.Label(toolbar_row2, text="宽:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_width = tk.Entry(toolbar_row2, textvariable=self.fig_width_var, width=5)
        self.entry_width.pack(side=tk.LEFT, padx=2)
        tk.Label(toolbar_row2, text="高:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_height = tk.Entry(toolbar_row2, textvariable=self.fig_height_var, width=5)
        self.entry_height.pack(side=tk.LEFT, padx=2)

        # 密码子堆叠块高度调整
        tk.Frame(toolbar_row2, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        tk.Label(toolbar_row2, text="堆叠块高度:", font=('Arial', 10)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_height = tk.Entry(toolbar_row2, textvariable=self.legend_block_height_var, width=6)
        self.entry_legend_height.pack(side=tk.LEFT, padx=2)

        # 密码子堆叠块宽度调整
        tk.Label(toolbar_row2, text="宽度:", font=('Arial', 10)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_width = tk.Entry(toolbar_row2, textvariable=self.legend_block_width_var, width=6)
        self.entry_legend_width.pack(side=tk.LEFT, padx=2)

        # 分隔符
        tk.Frame(toolbar_row2, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 重新绘制按钮
        btn_redraw = tk.Button(toolbar_row2, text="重新绘制", command=self.redraw_chart,
                               bg='#607D8B', fg='white', font=('Arial', 10, 'bold'),
                               padx=15, pady=5)
        btn_redraw.pack(side=tk.LEFT, padx=5)

        # 顶部工具栏 - 第三行（文本设置）
        toolbar_row3 = tk.Frame(self.root, bd=1, relief=tk.RAISED)
        toolbar_row3.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)

        # 柱状图内密码子名称设置
        tk.Label(toolbar_row3, text="柱状图文字:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="字号:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_codon_label_size = tk.Entry(toolbar_row3, textvariable=self.codon_label_size_var, width=4)
        self.entry_codon_label_size.pack(side=tk.LEFT, padx=2)

        # 分隔符
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 坐标轴标签设置
        tk.Label(toolbar_row3, text="坐标轴:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="X轴字号:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_xlabel_size = tk.Entry(toolbar_row3, textvariable=self.xlabel_size_var, width=4)
        self.entry_xlabel_size.pack(side=tk.LEFT, padx=2)

        tk.Label(toolbar_row3, text="Y轴字号:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_ylabel_size = tk.Entry(toolbar_row3, textvariable=self.ylabel_size_var, width=4)
        self.entry_ylabel_size.pack(side=tk.LEFT, padx=2)

        # 分隔符
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 刻度标签设置
        tk.Label(toolbar_row3, text="刻度:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="X轴字号:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_tick_label_size = tk.Entry(toolbar_row3, textvariable=self.tick_label_size_var, width=4)
        self.entry_tick_label_size.pack(side=tk.LEFT, padx=2)

        tk.Label(toolbar_row3, text="Y轴字号:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_ytick_label_size = tk.Entry(toolbar_row3, textvariable=self.ytick_label_size_var, width=4)
        self.entry_ytick_label_size.pack(side=tk.LEFT, padx=2)

        # 分隔符
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 边框粗细设置
        tk.Label(toolbar_row3, text="边框:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="粗细:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_border_width = tk.Entry(toolbar_row3, textvariable=self.border_width_var, width=5)
        self.entry_border_width.pack(side=tk.LEFT, padx=2)

        # 分隔符
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 堆叠块位置设置
        tk.Label(toolbar_row3, text="堆叠块位置:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="距离X轴:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_bottom_offset = tk.Entry(toolbar_row3, textvariable=self.legend_bottom_offset_var, width=5)
        self.entry_legend_bottom_offset.pack(side=tk.LEFT, padx=2)

        # 分隔符
        tk.Frame(toolbar_row3, width=2, bd=1, relief=tk.SUNKEN).pack(side=tk.LEFT, fill=tk.Y, padx=10)

        # 堆叠块文字设置
        tk.Label(toolbar_row3, text="堆叠块文字:", font=('Arial', 10, 'bold')).pack(side=tk.LEFT, padx=5)
        tk.Label(toolbar_row3, text="字号:", font=('Arial', 9)).pack(side=tk.LEFT, padx=2)
        self.entry_legend_text_size = tk.Entry(toolbar_row3, textvariable=self.legend_text_size_var, width=4)
        self.entry_legend_text_size.pack(side=tk.LEFT, padx=2)

        # 主绘图区域容器（包含滚动条）
        self.plot_container = tk.Frame(self.root)
        self.plot_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # 创建垂直滚动条
        v_scrollbar = tk.Scrollbar(self.plot_container, orient=tk.VERTICAL, command=self._on_vscroll)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # 创建水平滚动条
        h_scrollbar = tk.Scrollbar(self.plot_container, orient=tk.HORIZONTAL, command=self._on_hscroll)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)

        # 创建可滚动的画布容器
        self.scroll_canvas = tk.Canvas(self.plot_container, bg='white')
        self.scroll_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # 配置滚动画布的滚动区域
        self.scroll_canvas.config(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)

        # 创建画布容器（固定大小，用于放置matplotlib画布）
        self.canvas_frame = tk.Frame(self.scroll_canvas, bg='white')
        self.canvas_window = self.scroll_canvas.create_window((0, 0), window=self.canvas_frame, anchor='nw')

        # 绑定画布大小变化事件
        self.canvas_frame.bind('<Configure>', self._on_frame_configure)

        # 提示信息
        if self.lbl_hint is None:
            self.lbl_hint = tk.Label(self.root,
                                   text="请先选择包含密码子数据的Excel文件，程序将自动生成表格和堆积柱状图",
                                   fg='blue', font=('Arial', 10))
        self.lbl_hint.pack(side=tk.BOTTOM, pady=5)
        
    def load_excel_file(self):
        """加载Excel文件"""
        file_path = filedialog.askopenfilename(
            title="选择Excel文件",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )

        if file_path:
            self.excel_file_path = file_path
            self.lbl_file.config(text=os.path.basename(file_path), fg='black')
            self.plot_chart()
            
    def plot_chart(self):
        """绘制图表"""
        try:
            # 清除旧图表
            if self.canvas:
                self.canvas.get_tk_widget().destroy()
            if hasattr(self, 'toolbar'):
                self.toolbar.destroy()

            # 生成图表
            colormap_name = self.colormap_var.get()
            show_codon_labels = self.show_codon_labels_var.get()

            # 获取画布大小
            try:
                fig_width = float(self.fig_width_var.get())
                fig_height = float(self.fig_height_var.get())
                if fig_width < 10 or fig_height < 8:
                    raise ValueError("画布大小过小")
                if fig_width > 30 or fig_height > 30:
                    raise ValueError("画布大小过大")
                figsize = (fig_width, fig_height)
            except ValueError as e:
                messagebox.showwarning("警告", f"画布大小格式错误：\n{str(e)}\n\n将使用默认大小 (18, 13)")
                figsize = (18, 13)
                self.fig_width_var.set('18')
                self.fig_height_var.set('13')

            # 获取密码子堆叠块高度
            try:
                legend_block_height = float(self.legend_block_height_var.get())
                if legend_block_height < 0.05 or legend_block_height > 0.3:
                    raise ValueError("堆叠块高度应在 0.05 到 0.3 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"堆叠块高度格式错误：\n{str(e)}\n\n将使用默认值 0.12")
                legend_block_height = 0.12
                self.legend_block_height_var.set('0.12')

            # 获取密码子堆叠块宽度
            try:
                legend_block_width = float(self.legend_block_width_var.get())
                if legend_block_width < 0.1 or legend_block_width > 1.0:
                    raise ValueError("堆叠块宽度应在 0.1 到 1.0 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"堆叠块宽度格式错误：\n{str(e)}\n\n将使用默认值 0.6")
                legend_block_width = 0.6
                self.legend_block_width_var.set('0.6')

            # 获取文本设置
            try:
                codon_label_size = float(self.codon_label_size_var.get())
                if codon_label_size < 6 or codon_label_size > 24:
                    raise ValueError("柱状图文字字号应在 6 到 24 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"柱状图文字字号格式错误：\n{str(e)}\n\n将使用默认值 10")
                codon_label_size = 10
                self.codon_label_size_var.set('10')

            try:
                xlabel_size = float(self.xlabel_size_var.get())
                ylabel_size = float(self.ylabel_size_var.get())
                if xlabel_size < 8 or xlabel_size > 30 or ylabel_size < 8 or ylabel_size > 30:
                    raise ValueError("坐标轴字号应在 8 到 30 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"坐标轴字号格式错误：\n{str(e)}\n\n将使用默认值 14")
                xlabel_size = 14
                ylabel_size = 14
                self.xlabel_size_var.set('14')
                self.ylabel_size_var.set('14')

            try:
                tick_label_size = float(self.tick_label_size_var.get())
                ytick_label_size = float(self.ytick_label_size_var.get())
                if tick_label_size < 6 or tick_label_size > 20 or ytick_label_size < 6 or ytick_label_size > 20:
                    raise ValueError("刻度标签字号应在 6 到 20 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"刻度标签字号格式错误：\n{str(e)}\n\n将使用默认值 12")
                tick_label_size = 12
                ytick_label_size = 12
                self.tick_label_size_var.set('12')
                self.ytick_label_size_var.set('12')

            try:
                border_width = float(self.border_width_var.get())
                if border_width < 0.5 or border_width > 5.0:
                    raise ValueError("边框粗细应在 0.5 到 5.0 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"边框粗细格式错误：\n{str(e)}\n\n将使用默认值 1.0")
                border_width = 1.0
                self.border_width_var.set('1.0')

            try:
                legend_text_size = float(self.legend_text_size_var.get())
                if legend_text_size < 4 or legend_text_size > 16:
                    raise ValueError("堆叠块文字字号应在 4 到 16 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"堆叠块文字字号格式错误：\n{str(e)}\n\n将使用默认值 8")
                legend_text_size = 8
                self.legend_text_size_var.set('8')

            # 获取堆叠块距离X坐标轴的距离
            try:
                legend_bottom_offset = float(self.legend_bottom_offset_var.get())
                if legend_bottom_offset < 0.05 or legend_bottom_offset > 0.5:
                    raise ValueError("堆叠块距离X轴应在 0.05 到 0.5 之间")
            except ValueError as e:
                messagebox.showwarning("警告", f"堆叠块距离X轴格式错误：\n{str(e)}\n\n将使用默认值 0.22")
                legend_bottom_offset = 0.22
                self.legend_bottom_offset_var.set('0.22')

            self.fig, self.df, _ = plot_codon_usage_table_and_chart(self.excel_file_path, colormap_name=colormap_name, show_codon_labels=show_codon_labels, figsize=figsize, legend_block_height=legend_block_height, legend_block_width=legend_block_width, show_legend_blocks=self.show_legend_blocks_var.get(), codon_label_size=codon_label_size, xlabel_size=xlabel_size, ylabel_size=ylabel_size, tick_label_size=tick_label_size, ytick_label_size=ytick_label_size, legend_text_size=legend_text_size, border_width=border_width, legend_bottom_offset=legend_bottom_offset)

            # 在Tkinter中显示（使用固定大小，不随窗口缩放）
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
            self.canvas.draw()

            # 将画布添加到frame，使用固定宽高，不随窗口缩放
            canvas_widget = self.canvas.get_tk_widget()
            canvas_widget.pack_propagate(False)  # 禁止frame自动调整大小
            canvas_widget.pack(fill=tk.BOTH, expand=False)

            # 设置画布的固定大小（基于figsize和dpi）
            dpi = 100
            canvas_width = int(figsize[0] * dpi)
            canvas_height = int(figsize[1] * dpi)
            canvas_widget.config(width=canvas_width, height=canvas_height)

            # 更新canvas_frame的大小以匹配画布
            self.canvas_frame.config(width=canvas_width, height=canvas_height)
            self.scroll_canvas.config(scrollregion=(0, 0, canvas_width, canvas_height))

            # 添加导航工具栏（支持缩放、平移等）
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.canvas_frame)
            self.toolbar.update()

            # 强制更新滚动区域
            self.root.update_idletasks()
            self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox("all"))

            self.lbl_hint.config(text="图表生成成功！使用下方工具栏可放大、缩小、平移图表", fg='green')

        except FileNotFoundError:
            messagebox.showerror("错误", f"文件不存在：\n{self.excel_file_path}")
            self.lbl_hint.config(text="错误：文件不存在", fg='red')
        except PermissionError:
            messagebox.showerror("错误", f"文件被占用或无权限访问：\n{self.excel_file_path}")
            self.lbl_hint.config(text="错误：文件被占用或无权限", fg='red')
        except Exception as e:
            messagebox.showerror("错误", f"加载文件失败：\n{str(e)}")
            self.lbl_hint.config(text=f"错误：{str(e)}", fg='red')

    def _on_frame_configure(self, event):
        """画布大小变化时更新滚动区域"""
        self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox("all"))

    def _on_vscroll(self, *args):
        """垂直滚动回调"""
        self.scroll_canvas.yview(*args)

    def _on_hscroll(self, *args):
        """水平滚动回调"""
        self.scroll_canvas.xview(*args)

    def redraw_chart(self):
        """使用新的配色方案重新绘制图表"""
        if self.excel_file_path is None:
            messagebox.showwarning("警告", "请先选择Excel文件")
            return
        self.plot_chart()

    def on_codon_label_toggle(self):
        """密码子标签显示/隐藏切换"""
        if self.excel_file_path is None:
            # 如果还没加载文件，只更新复选框状态
            return
        # 重新绘制图表
        self.redraw_chart()

    def on_legend_blocks_toggle(self):
        """堆叠块显示/隐藏切换"""
        if self.excel_file_path is None:
            # 如果还没加载文件，只更新复选框状态
            return
        # 重新绘制图表
        self.redraw_chart()
            
    def save_as_pdf(self):
        """保存为PDF"""
        if self.fig is None:
            messagebox.showwarning("警告", "请先加载Excel文件生成图表")
            return
            
        file_path = filedialog.asksaveasfilename(
            title="保存为PDF",
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")]
        )
        
        if file_path:
            self.fig.savefig(file_path, format='pdf', dpi=300, bbox_inches='tight')
            messagebox.showinfo("成功", f"PDF文件已保存至：\n{file_path}")
            
    def save_as_svg(self):
        """保存为SVG矢量图"""
        if self.fig is None:
            messagebox.showwarning("警告", "请先加载Excel文件生成图表")
            return
            
        file_path = filedialog.asksaveasfilename(
            title="保存为SVG矢量图",
            defaultextension=".svg",
            filetypes=[("SVG files", "*.svg"), ("All files", "*.*")]
        )
        
        if file_path:
            self.fig.savefig(file_path, format='svg', dpi=300, bbox_inches='tight')
            messagebox.showinfo("成功", f"SVG文件已保存至：\n{file_path}")
            
    def save_as_png(self):
        """保存为PNG"""
        if self.fig is None:
            messagebox.showwarning("警告", "请先加载Excel文件生成图表")
            return
            
        file_path = filedialog.asksaveasfilename(
            title="保存为PNG",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
        )
        
        if file_path:
            self.fig.savefig(file_path, format='png', dpi=300, bbox_inches='tight')
            messagebox.showinfo("成功", f"PNG文件已保存至：\n{file_path}")
            
    def run(self):
        """运行应用程序"""
        self.root.mainloop()


# 命令行模式（保留向后兼容）
def run_cli():
    """命令行模式运行"""
    import sys
    if len(sys.argv) > 1:
        excel_path = sys.argv[1]
    else:
        excel_path = '新建文件夹/A1.xlsx'
    
    if len(sys.argv) > 2:
        output_path = sys.argv[2]
    else:
        output_path = 'codon_usage_chart.pdf'
    
    fig, _, saved_path = plot_codon_usage_table_and_chart(excel_path, output_path)
    print(f"图表已保存至: {saved_path}")
    plt.show()


# 调试函数：查看Excel文件内容
def debug_excel_file(file_path):
    """调试函数：查看Excel文件内容"""
    try:
        print(f"\n正在读取文件: {file_path}")
        df = pd.read_excel(file_path)
        print(f"\n文件列名: {df.columns.tolist()}")
        print(f"数据类型: {df.dtypes}")
        print(f"\n前10行数据:")
        print(df.head(10))
        print(f"\n数据形状: {df.shape}")
        return df
    except Exception as e:
        print(f"读取Excel文件失败: {e}")
        import traceback
        traceback.print_exc()
        return None


# 使用示例
if __name__ == '__main__':
    import sys

    # 如果有命令行参数且是 --debug，则进入调试模式
    if len(sys.argv) > 1 and sys.argv[1] == '--debug':
        print("=== 调试模式 ===")
        debug_excel_file('新建文件夹/A1.xlsx')
        input("\n按Enter键退出...")
    else:
        # 启动交互式应用程序
        try:
            app = CodonUsageApp()
            app.create_widgets()
            app.run()
        except Exception as e:
            print(f"程序运行出错: {e}")
            import traceback
            traceback.print_exc()
            input("按Enter键退出...")

