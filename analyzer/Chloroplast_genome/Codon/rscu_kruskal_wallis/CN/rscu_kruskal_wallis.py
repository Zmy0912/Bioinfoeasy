# -*- coding: utf-8 -*-
"""
RSCU Kruskal-Wallis 检验分析软件
用于对不同物种间密码子RSCU值进行Kruskal-Wallis非参数检验
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import numpy as np
from scipy.stats import kruskal
import os
import warnings
warnings.filterwarnings('ignore')


class RSCUAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("RSCU Kruskal-Wallis 检验分析软件")
        self.root.geometry("1000x800")
        self.root.configure(bg="#f0f0f0")
        
        self.data = None
        self.results = None
        self.amino_acids = {}
        
        self.setup_styles()
        self.create_widgets()
        
    def setup_styles(self):
        style = ttk.Style()
        style.theme_use('clam')
        style.configure("Title.TLabel", font=("Microsoft YaHei", 16, "bold"), background="#f0f0f0")
        style.configure("Normal.TLabel", font=("Microsoft YaHei", 10), background="#f0f0f0")
        style.configure("Header.TLabel", font=("Microsoft YaHei", 11, "bold"), background="#f0f0f0")
        
    def create_widgets(self):
        # 标题
        title_frame = tk.Frame(self.root, bg="#2c3e50")
        title_frame.pack(fill=tk.X)
        ttk.Label(title_frame, text="🔬 密码子RSCU值 Kruskal-Wallis 检验分析软件", 
                  style="Title.TLabel", foreground="white", background="#2c3e50").pack(pady=10)
        
        # 说明区域
        info_frame = tk.LabelFrame(self.root, text="📊 分析说明", bg="#f0f0f0", font=("Microsoft YaHei", 10))
        info_frame.pack(fill=tk.X, padx=20, pady=10)
        
        info_text = """• 本软件用于比较不同物种间密码子RSCU使用偏好的差异
• Kruskal-Wallis检验是非参数检验，适用于多组样本比较
• 分析方法：对每个密码子，比较不同物种组在该密码子RSCU值上的差异
• 显著性结果表示该密码子的使用偏好在不同物种组间存在显著差异"""
        ttk.Label(info_frame, text=info_text, background="#f0f0f0", justify=tk.LEFT, 
                  font=("Microsoft YaHei", 9)).pack(anchor=tk.W, padx=10, pady=5)
        
        # 文件选择区域
        file_frame = tk.LabelFrame(self.root, text="📁 数据文件", bg="#f0f0f0", font=("Microsoft YaHei", 10))
        file_frame.pack(fill=tk.X, padx=20, pady=10)
        
        self.file_path_var = tk.StringVar(value="未选择文件")
        ttk.Label(file_frame, text="RSCU数据文件:", background="#f0f0f0").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        ttk.Entry(file_frame, textvariable=self.file_path_var, width=60).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(file_frame, text="浏览...", command=self.browse_file).grid(row=0, column=2, padx=5, pady=5)
        ttk.Button(file_frame, text="加载数据", command=self.load_data, style="Action.TButton").grid(row=0, column=3, padx=5, pady=5)
        
        # 分组设置区域
        group_frame = tk.LabelFrame(self.root, text="📋 物种分组设置", bg="#f0f0f0", font=("Microsoft YaHei", 10))
        group_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # 分组方式选择
        self.group_method = tk.StringVar(value="auto")
        
        ttk.Label(group_frame, text="分组方式:", background="#f0f0f0").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        ttk.Radiobutton(group_frame, text="按物种名自动分组（前缀相同为一组）", 
                        variable=self.group_method, value="auto").grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Radiobutton(group_frame, text="手动定义分组（加载分组文件）", 
                        variable=self.group_method, value="manual").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        self.group_file_var = tk.StringVar(value="未选择分组文件")
        ttk.Label(group_frame, text="分组文件:", background="#f0f0f0").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        ttk.Entry(group_frame, textvariable=self.group_file_var, width=40).grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(group_frame, text="浏览...", command=self.browse_group_file).grid(row=1, column=2, padx=5, pady=5)
        
        # 分组预览
        self.group_preview_text = tk.Text(group_frame, height=5, width=80, font=("Consolas", 9))
        self.group_preview_text.grid(row=2, column=0, columnspan=3, padx=10, pady=5)
        
        # 分析选项区域
        option_frame = tk.LabelFrame(self.root, text="⚙️ 分析选项", bg="#f0f0f0", font=("Microsoft YaHei", 10))
        option_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Label(option_frame, text="分析范围:", background="#f0f0f0").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.analysis_scope = tk.StringVar(value="all")
        ttk.Radiobutton(option_frame, text="全部密码子", variable=self.analysis_scope, value="all").grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Radiobutton(option_frame, text="按氨基酸选择", variable=self.analysis_scope, value="by_amino").grid(row=0, column=2, sticky=tk.W, padx=5)
        
        # 氨基酸选择框架
        self.amino_frame = tk.Frame(option_frame, bg="#f0f0f0")
        self.amino_frame.grid(row=1, column=0, columnspan=4, sticky=tk.W, padx=10, pady=5)
        
        self.amino_vars = {}
        amino_list = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 
                      'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
        for i, amino in enumerate(amino_list):
            var = tk.BooleanVar(value=True)
            self.amino_vars[amino] = var
            cb = ttk.Checkbutton(self.amino_frame, text=amino, variable=var)
            cb.grid(row=i // 5, column=i % 5, sticky=tk.W, padx=2)
        
        # 显著性水平
        self.significance_level = tk.DoubleVar(value=0.05)
        ttk.Label(option_frame, text="显著性水平 (α):", background="#f0f0f0").grid(row=0, column=3, sticky=tk.W, padx=(20, 5))
        ttk.Entry(option_frame, textvariable=self.significance_level, width=8).grid(row=0, column=4, padx=5)
        
        # 按钮区域
        button_frame = tk.Frame(self.root, bg="#f0f0f0")
        button_frame.pack(fill=tk.X, padx=20, pady=10)
        
        ttk.Button(button_frame, text="▶ 执行 Kruskal-Wallis 检验", style="Action.TButton", 
                   command=self.run_analysis).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="📊 绘制箱线图", style="Action.TButton", 
                   command=self.plot_boxplot).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="📥 导出完整结果", style="Action.TButton", 
                   command=self.export_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="📥 导出显著结果", style="Action.TButton", 
                   command=self.export_significant).pack(side=tk.LEFT, padx=5)
        
        # 结果显示区域
        result_frame = tk.LabelFrame(self.root, text="📈 检验结果", bg="#f0f0f0", font=("Microsoft YaHei", 10))
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # 创建Text和Scrollbar
        result_scroll = ttk.Scrollbar(result_frame)
        result_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.result_text = tk.Text(result_frame, yscrollcommand=result_scroll.set,
                                    font=("Consolas", 9), bg="#ffffff", wrap=tk.NONE)
        self.result_text.pack(fill=tk.BOTH, expand=True)
        result_scroll.config(command=self.result_text.yview)
        
        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = tk.Label(self.root, textvariable=self.status_var, bd=1, relief=tk.SUNKEN, 
                              anchor=tk.W, bg="#e0e0e0", font=("Microsoft YaHei", 9))
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="选择RSCU数据文件",
            filetypes=[("Excel文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        )
        if filename:
            self.file_path_var.set(filename)
            
    def browse_group_file(self):
        filename = filedialog.askopenfilename(
            title="选择分组文件",
            filetypes=[("Excel文件", "*.xlsx *.xls"), ("CSV文件", "*.csv"), ("文本文件", "*.txt")]
        )
        if filename:
            self.group_file_var.set(filename)
            
    def load_data(self):
        filepath = self.file_path_var.get()
        if not filepath or filepath == "未选择文件":
            messagebox.showwarning("警告", "请先选择数据文件")
            return
            
        try:
            self.data = pd.read_excel(filepath)
            
            if 'Species' not in self.data.columns:
                messagebox.showerror("错误", "数据文件必须包含'Species'列")
                return
                
            self.status_var.set(f"已加载: {len(self.data)} 个物种, {len(self.data.columns)-1} 个密码子")
            
            # 解析氨基酸信息
            self.parse_amino_acids()
            
            # 处理分组
            self.process_groups()
            
            messagebox.showinfo("成功", f"数据加载成功！\n共 {len(self.data)} 个物种\n{len(self.data.columns)-1} 个密码子\n\n物种列表:\n" + "\n".join(self.data['Species'].tolist()[:10]) + ("\n..." if len(self.data) > 10 else ""))
            
        except Exception as e:
            messagebox.showerror("错误", f"加载数据失败: {str(e)}")
            
    def parse_amino_acids(self):
        """解析密码子列名"""
        self.amino_acids = {}
        codon_columns = [col for col in self.data.columns if col != 'Species']
        
        for col in codon_columns:
            try:
                codon = col.split('(')[0]
                amino = col.split('(')[1].rstrip(')')
                if amino not in self.amino_acids:
                    self.amino_acids[amino] = []
                self.amino_acids[amino].append(col)
            except:
                pass
                
    def process_groups(self):
        """处理物种分组"""
        self.group_preview_text.delete(1.0, tk.END)
        method = self.group_method.get()
        
        if method == "auto":
            species = self.data['Species'].tolist()
            
            # 尝试找自然的物种分组 - 按物种名前几个字符相似性
            groups = {}
            for sp in species:
                # 尝试用空格或下划线分割后的第二部分作为组标识
                parts = sp.replace('_', ' ').split()
                if len(parts) >= 2:
                    # 使用种名的关键部分作为分组
                    key = parts[0]  # 属名作为第一组
                else:
                    key = sp[:8] if len(sp) > 8 else sp
                    
                if key not in groups:
                    groups[key] = []
                groups[key].append(sp)
            
            # 如果分组数太少（单一属或所有物种都不同），提供建议
            if len(groups) < 3:
                self.group_preview_text.insert(tk.END, "⚠️ 提示: 检测到所有物种可能属于同一类群\n")
                self.group_preview_text.insert(tk.END, "   请使用【手动定义分组】上传分组文件\n")
                self.group_preview_text.insert(tk.END, f"   或提供物种的地理分布/系统发育分组信息\n\n")
                self.group_preview_text.insert(tk.END, "   示例分组方案（仅演示）:\n")
                # 提供一个默认的均匀分组方案用于演示
                n = len(species)
                num_groups = min(5, max(2, n // 5))  # 根据样本量确定组数
                groups = {}
                for i, sp in enumerate(species):
                    group_idx = i % num_groups + 1
                    key = f"Demo_Group{group_idx}"
                    if key not in groups:
                        groups[key] = []
                    groups[key].append(sp)
                self.group_preview_text.insert(tk.END, "   [使用演示分组继续分析]\n\n")
            
            self.species_groups = groups
            self.group_method.set("auto")
            
        else:
            # 从文件加载分组
            group_file = self.group_file_var.get()
            if not group_file or group_file == "未选择分组文件":
                messagebox.showwarning("提示", "请先选择分组文件，或切换到自动分组模式")
                self.group_method.set("auto")
                self.process_groups()
                return
            self.load_groups_from_file(group_file)
            
        # 显示分组预览
        self.group_preview_text.insert(tk.END, f"共 {len(self.species_groups)} 个分组:\n")
        self.group_preview_text.insert(tk.END, "-" * 50 + "\n")
        for group_name, members in self.species_groups.items():
            self.group_preview_text.insert(tk.END, f"【{group_name}】({len(members)}个物种): {', '.join(members)}\n")
            
    def load_groups_from_file(self, filepath):
        """从文件加载分组"""
        try:
            if filepath.endswith('.xlsx') or filepath.endswith('.xls'):
                group_df = pd.read_excel(filepath)
            elif filepath.endswith('.csv'):
                group_df = pd.read_csv(filepath)
            else:
                with open(filepath, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                    group_df = pd.read_csv(filepath, delim_whitespace=True)
                    
            self.species_groups = {}
            for _, row in group_df.iterrows():
                species = str(row.iloc[0])
                group = str(row.iloc[1])
                if group not in self.species_groups:
                    self.species_groups[group] = []
                self.species_groups[group].append(species)
                
        except Exception as e:
            messagebox.showerror("错误", f"加载分组文件失败: {str(e)}")
            
    def run_analysis(self):
        """执行Kruskal-Wallis检验"""
        if self.data is None:
            messagebox.showwarning("警告", "请先加载数据")
            return
            
        self.status_var.set("正在执行分析...")
        self.result_text.delete(1.0, tk.END)
        
        alpha = self.significance_level.get()
        scope = self.analysis_scope.get()
        
        # 确定要分析的密码子
        if scope == "all":
            codons_to_analyze = [col for col in self.data.columns if col != 'Species']
        else:
            codons_to_analyze = []
            for amino, var in self.amino_vars.items():
                if var.get() and amino in self.amino_acids:
                    codons_to_analyze.extend(self.amino_acids[amino])
                    
        if not codons_to_analyze:
            messagebox.showwarning("警告", "请至少选择一个氨基酸进行分析")
            return
            
        results_list = []
        
        # 标题
        self.result_text.insert(tk.END, "=" * 90 + "\n")
        self.result_text.insert(tk.END, "🔬 RSCU Kruskal-Wallis 检验分析报告\n")
        self.result_text.insert(tk.END, "=" * 90 + "\n\n")
        
        # 基本信息
        self.result_text.insert(tk.END, "【基本信息】\n")
        self.result_text.insert(tk.END, f"  分析物种数: {len(self.data)} 个\n")
        self.result_text.insert(tk.END, f"  物种分组数: {len(self.species_groups)} 组\n")
        self.result_text.insert(tk.END, f"  分析密码子数: {len(codons_to_analyze)} 个\n")
        self.result_text.insert(tk.END, f"  显著性水平: α = {alpha}\n\n")
        
        # 分组详情
        self.result_text.insert(tk.END, "【分组详情】\n")
        for group_name, members in self.species_groups.items():
            self.result_text.insert(tk.END, f"  {group_name}: {len(members)}个物种\n")
        self.result_text.insert(tk.END, "\n")
        
        # 检验结果表头
        self.result_text.insert(tk.END, "-" * 90 + "\n")
        self.result_text.insert(tk.END, f"{'密码子':<15} {'氨基酸':<8} {'H统计量':<12} {'p值':<15} {'显著性':<8} {'解读'}\n")
        self.result_text.insert(tk.END, "-" * 90 + "\n")
        
        for codon in codons_to_analyze:
            try:
                # 提取氨基酸名称
                amino = codon.split('(')[1].rstrip(')') if '(' in codon else ''
                
                # 准备各组数据
                groups_data = []
                group_names = list(self.species_groups.keys())
                for group_name in group_names:
                    species_in_group = self.species_groups[group_name]
                    mask = self.data['Species'].isin(species_in_group)
                    values = self.data.loc[mask, codon].values
                    groups_data.append(values)
                
                # 检查是否有足够数据
                if any(len(g) < 2 for g in groups_data):
                    continue
                    
                # Kruskal-Wallis检验
                stat, p_value = kruskal(*groups_data)
                
                # 判断显著性
                if p_value < 0.001:
                    sig_level = "***"
                    sig_text = "极显著"
                elif p_value < 0.01:
                    sig_level = "**"
                    sig_text = "高度显著"
                elif p_value < alpha:
                    sig_level = "*"
                    sig_text = "显著"
                else:
                    sig_level = "ns"
                    sig_text = "不显著"
                
                # 生成解读
                interpretation = self.generate_interpretation(codon, groups_data, group_names, stat, p_value)
                
                results_list.append({
                    'Codon': codon,
                    'AminoAcid': amino,
                    'H_statistic': stat,
                    'p_value': p_value,
                    'Significance': sig_text,
                    'p_sign': sig_level,
                    'Interpretation': interpretation,
                    'GroupMeans': {gn: np.mean(g) for gn, g in zip(group_names, groups_data)},
                    'GroupMedians': {gn: np.median(g) for gn, g in zip(group_names, groups_data)},
                    'GroupSDs': {gn: np.std(g) for gn, g in zip(group_names, groups_data)}
                })
                
                self.result_text.insert(tk.END, f"{codon:<15} {amino:<8} {stat:<12.4f} {p_value:<15.2e} {sig_level:<8} {interpretation}\n")
                
            except Exception as e:
                pass
                
        self.results = results_list
        
        # 统计摘要
        sig_count = sum(1 for r in results_list if r['p_value'] < alpha)
        self.result_text.insert(tk.END, "-" * 90 + "\n")
        self.result_text.insert(tk.END, "\n【统计摘要】\n")
        self.result_text.insert(tk.END, f"  总分析密码子数: {len(results_list)}\n")
        self.result_text.insert(tk.END, f"  显著差异密码子数 (p < {alpha}): {sig_count}\n")
        self.result_text.insert(tk.END, f"  无显著差异密码子数: {len(results_list) - sig_count}\n")
        self.result_text.insert(tk.END, f"  显著率: {sig_count/len(results_list)*100:.1f}%\n\n")
        
        # 详细解读
        self.result_text.insert(tk.END, "【结果解读】\n")
        if sig_count > 0:
            self.result_text.insert(tk.END, f"  在{p_sig_count}个显著差异的密码子中:\n\n")
            
            # 按氨基酸分组展示
            by_amino = {}
            for r in results_list:
                if r['p_value'] < alpha:
                    amino = r['AminoAcid']
                    if amino not in by_amino:
                        by_amino[amino] = []
                    by_amino[amino].append(r)
            
            for amino, codon_results in sorted(by_amino.items()):
                self.result_text.insert(tk.END, f"  ▶ {amino}氨基酸 ({len(codon_results)}个密码子显著差异):\n")
                for r in codon_results:
                    max_group = max(r['GroupMeans'].items(), key=lambda x: x[1])
                    min_group = min(r['GroupMeans'].items(), key=lambda x: x[1])
                    self.result_text.insert(tk.END, f"    • {r['Codon']}: p={r['p_value']:.2e}\n")
                    self.result_text.insert(tk.END, f"      解读: 该密码子使用偏好在不同物种间{r['Significance']}\n")
                    self.result_text.insert(tk.END, f"      {max_group[0]}组偏好性最高(均值={max_group[1]:.2f}), {min_group[0]}组偏好性最低(均值={min_group[1]:.2f})\n")
                self.result_text.insert(tk.END, "\n")
        
        self.status_var.set(f"分析完成，显著差异: {sig_count}/{len(results_list)}")
        
    def generate_interpretation(self, codon, groups_data, group_names, stat, p_value):
        """生成结果解读"""
        if p_value >= 0.05:
            return "密码子使用无显著差异"
        elif p_value >= 0.01:
            return "物种间存在一定差异"
        elif p_value >= 0.001:
            return "物种间存在明显差异"
        else:
            return "物种间存在极显著差异"
            
    def plot_boxplot(self):
        """绘制箱线图"""
        if self.results is None or len(self.results) == 0:
            messagebox.showwarning("警告", "请先执行分析")
            return
            
        try:
            import matplotlib.pyplot as plt
            plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
            plt.rcParams['axes.unicode_minus'] = False
            
            # 创建图形
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle('RSCU Kruskal-Wallis检验 - 箱线图分析', fontsize=14)
            
            # 选择p值最小的4个密码子绘图
            top_results = sorted(self.results, key=lambda x: x['p_value'])[:4]
            
            for idx, (ax, result) in enumerate(zip(axes.flat, top_results)):
                codon = result['Codon']
                
                # 准备箱线图数据
                box_data = []
                labels = []
                for group_name in sorted(self.species_groups.keys()):
                    species_in_group = self.species_groups[group_name]
                    mask = self.data['Species'].isin(species_in_group)
                    values = self.data.loc[mask, codon].values
                    box_data.append(values)
                    labels.append(group_name)
                
                bp = ax.boxplot(box_data, labels=labels, patch_artist=True)
                
                # 设置颜色
                colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                ax.set_title(f"{codon}\np = {result['p_value']:.2e} ({result['Significance']})", fontsize=10)
                ax.set_ylabel('RSCU值')
                ax.tick_params(axis='x', rotation=45)
                
            plt.tight_layout()
            
            # 保存图形
            filepath = filedialog.asksaveasfilename(
                title="保存箱线图",
                defaultextension=".png",
                filetypes=[("PNG图像", "*.png"), ("PDF图像", "*.pdf")],
                initialfile="RSCU_Boxplot"
            )
            
            if filepath:
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                messagebox.showinfo("成功", f"箱线图已保存至:\n{filepath}")
                
            plt.close()
            
        except ImportError:
            messagebox.showerror("错误", "请先安装matplotlib: pip install matplotlib")
        except Exception as e:
            messagebox.showerror("错误", f"绘图失败: {str(e)}")
            
    def export_results(self):
        """导出完整结果"""
        if self.results is None or len(self.results) == 0:
            messagebox.showwarning("警告", "请先执行分析")
            return
            
        filepath = filedialog.asksaveasfilename(
            title="导出完整结果",
            defaultextension=".xlsx",
            filetypes=[("Excel文件", "*.xlsx"), ("CSV文件", "*.csv")],
            initialfile="RSCU_KruskalWallis_Full_Results"
        )
        
        if filepath:
            try:
                # 主结果表
                results_df = pd.DataFrame([{
                    '密码子': r['Codon'],
                    '氨基酸': r['AminoAcid'],
                    'H统计量': r['H_statistic'],
                    'p值': r['p_value'],
                    '显著性': r['Significance'],
                    '结果解读': r['Interpretation']
                } for r in self.results])
                
                # 分组统计表
                group_names = sorted(self.species_groups.keys())
                stats_data = []
                for r in self.results:
                    row = {'密码子': r['Codon']}
                    for gn in group_names:
                        row[f'{gn}_均值'] = r['GroupMeans'].get(gn, np.nan)
                        row[f'{gn}_中位数'] = r['GroupMedians'].get(gn, np.nan)
                        row[f'{gn}_标准差'] = r['GroupSDs'].get(gn, np.nan)
                    stats_data.append(row)
                stats_df = pd.DataFrame(stats_data)
                
                if filepath.endswith('.csv'):
                    results_df.to_csv(filepath, index=False, encoding='utf-8-sig')
                    stats_path = filepath.replace('.csv', '_分组统计.csv')
                    stats_df.to_csv(stats_path, index=False, encoding='utf-8-sig')
                else:
                    with pd.ExcelWriter(filepath) as writer:
                        results_df.to_excel(writer, sheet_name='检验结果', index=False)
                        stats_df.to_excel(writer, sheet_name='分组统计', index=False)
                        
                messagebox.showinfo("成功", f"结果已保存至:\n{filepath}\n\n分组统计已保存")
                
            except Exception as e:
                messagebox.showerror("错误", f"导出失败: {str(e)}")
                
    def export_significant(self):
        """导出显著结果"""
        if self.results is None or len(self.results) == 0:
            messagebox.showwarning("警告", "请先执行分析")
            return
            
        alpha = self.significance_level.get()
        sig_results = [r for r in self.results if r['p_value'] < alpha]
        
        if len(sig_results) == 0:
            messagebox.showinfo("提示", f"没有发现显著差异 (p < {alpha}) 的结果")
            return
            
        filepath = filedialog.asksaveasfilename(
            title="导出显著结果",
            defaultextension=".xlsx",
            filetypes=[("Excel文件", "*.xlsx"), ("CSV文件", "*.csv")],
            initialfile="RSCU_KruskalWallis_Significant"
        )
        
        if filepath:
            try:
                sig_df = pd.DataFrame([{
                    '密码子': r['Codon'],
                    '氨基酸': r['AminoAcid'],
                    'H统计量': r['H_statistic'],
                    'p值': r['p_value'],
                    '显著性水平': r['Significance'],
                    '结果解读': r['Interpretation'],
                    '偏好最高组': max(r['GroupMeans'].items(), key=lambda x: x[1])[0],
                    '偏好最高组均值': max(r['GroupMeans'].items(), key=lambda x: x[1])[1],
                    '偏好最低组': min(r['GroupMeans'].items(), key=lambda x: x[1])[0],
                    '偏好最低组均值': min(r['GroupMeans'].items(), key=lambda x: x[1])[1]
                } for r in sig_results])
                
                if filepath.endswith('.csv'):
                    sig_df.to_csv(filepath, index=False, encoding='utf-8-sig')
                else:
                    sig_df.to_excel(filepath, index=False)
                    
                messagebox.showinfo("成功", f"显著结果已保存至:\n{filepath}\n(共 {len(sig_results)} 条记录)")
                
            except Exception as e:
                messagebox.showerror("错误", f"导出失败: {str(e)}")


def main():
    root = tk.Tk()
    app = RSCUAnalyzerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
