"""
CodonW密码子偏好性统计分析工具
用于处理CodonW输出结果，生成热图可直接使用的表格
"""

import os
import re
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path


# 密码子与氨基酸的对应关系
CODON_AA_MAP = {
    'UUU': 'Phe', 'UUC': 'Phe',
    'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
    'AUG': 'Met',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr',
    'UAA': 'TER', 'UAG': 'TER', 'UGA': 'TER',
    'CAU': 'His', 'CAC': 'His',
    'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn',
    'AAA': 'Lys', 'AAG': 'Lys',
    'GAU': 'Asp', 'GAC': 'Asp',
    'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys',
    'UGG': 'Trp',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# 所有61个密码子（排除终止密码子）
ALL_CODONS = [c for c in CODON_AA_MAP.keys() if c not in ['UAA', 'UAG', 'UGA']]
ALL_CODONS.sort()

# 氨基酸排序顺序（按标准生物化学顺序）
AA_ORDER = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
            'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']


class CodonWAnalyzer:
    """CodonW结果分析器"""
    
    def __init__(self):
        self.all_data = {}
        self.notes = []
        
    def parse_file(self, file_path):
        """解析单个CodonW结果文件"""
        species_name = Path(file_path).stem
        codon_data = {}
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        for line in lines:
            # 匹配密码子数据的正则表达式
            # 格式：Phe UUU  791 1.27 或     UUC  450 0.73
            pattern = r'(?:\w{3}\s+)?([A-Z]{3})\s+(\d+)\s+([\d.]+)'
            matches = re.findall(pattern, line)
            
            for match in matches:
                codon, count, rscu = match
                if codon in CODON_AA_MAP:
                    aa = CODON_AA_MAP[codon]
                    codon_data[codon] = {
                        'AA': aa,
                        'Count': int(count),
                        'RSCU': float(rscu)
                    }
        
        return species_name, codon_data
    
    def analyze_folder(self, folder_path):
        """分析文件夹中的所有CodonW结果"""
        self.all_data = {}
        self.notes = []
        
        # 查找所有.blk文件
        blk_files = list(Path(folder_path).glob('*.blk'))
        
        if not blk_files:
            raise ValueError(f"在文件夹 {folder_path} 中未找到 .blk 文件")
        
        # 解析所有文件
        for file_path in blk_files:
            species_name, codon_data = self.parse_file(file_path)
            self.all_data[species_name] = codon_data
        
        # 检查数据完整性并记录需要清理的内容
        self._check_data_quality()
        
        return len(blk_files)
    
    def _check_data_quality(self):
        """检查数据质量，记录需要注意的问题"""
        self.notes.append("=== 数据质量报告 ===\n")
        self.notes.append(f"共处理 {len(self.all_data)} 个物种的CodonW结果\n\n")
        
        # 检查每个物种的密码子完整性
        missing_codons = []
        for species, data in self.all_data.items():
            for codon in ALL_CODONS:
                if codon not in data:
                    missing_codons.append(f"{species} 缺少密码子 {codon}")
        
        if missing_codons:
            self.notes.append("⚠️ 缺失的密码子：\n")
            for note in missing_codons:
                self.notes.append(f"  - {note}\n")
            self.notes.append("\n")
        
        # 检查零值密码子
        zero_codons = []
        for species, data in self.all_data.items():
            for codon in ALL_CODONS:
                if codon in data and data[codon]['Count'] == 0:
                    zero_codons.append(f"{species}: {codon} (RSCU={data[codon]['RSCU']})")
        
        if zero_codons:
            self.notes.append("⚠️ 计数为0的密码子（可能影响统计分析）：\n")
            for note in zero_codons[:20]:  # 只显示前20个
                self.notes.append(f"  - {note}\n")
            if len(zero_codons) > 20:
                self.notes.append(f"  ... 还有 {len(zero_codons) - 20} 个\n")
            self.notes.append("\n")
        
        # 检查RSCU异常值
        high_rscu = []
        for species, data in self.all_data.items():
            for codon in ALL_CODONS:
                if codon in data and data[codon]['RSCU'] > 2.0:
                    high_rscu.append(f"{species}: {codon} (RSCU={data[codon]['RSCU']})")
        
        if high_rscu:
            self.notes.append("⚠️ RSCU值 > 2.0 的密码子（可能存在异常）：\n")
            for note in high_rscu[:10]:
                self.notes.append(f"  - {note}\n")
            if len(high_rscu) > 10:
                self.notes.append(f"  ... 还有 {len(high_rscu) - 10} 个\n")
            self.notes.append("\n")
        
        self.notes.append("=== 清理建议 ===\n")
        self.notes.append("1. 缺失密码子的数据可能需要手动补充或排除\n")
        self.notes.append("2. 计数为0的密码子在热图分析中可能导致偏差，考虑是否保留\n")
        self.notes.append("3. RSCU异常值建议检查原始数据或分析流程\n")
        self.notes.append("4. 对于特定密码子的偏好性分析，可以筛选感兴趣的数据\n")
    
    def create_heatmap_table(self, output_file, metric='RSCU', sort_by_aa=False):
        """创建用于热图绘制的表格"""
        if metric not in ['RSCU', 'Count']:
            raise ValueError("metric 必须是 'RSCU' 或 'Count'")
        
        # 创建DataFrame，行为物种，列为密码子
        species_list = sorted(self.all_data.keys())
        df_data = []
        
        for species in species_list:
            row = {'Species': species}
            codon_data = self.all_data[species]
            
            for codon in ALL_CODONS:
                aa = CODON_AA_MAP[codon]
                column_name = f"{codon}({aa})"
                if codon in codon_data:
                    row[column_name] = codon_data[codon][metric]
                else:
                    row[column_name] = None  # 缺失值
            
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        
        # 按氨基酸排序列
        if sort_by_aa:
            # 提取列名中的氨基酸部分，并按照AA_ORDER排序
            columns = ['Species']
            for aa in AA_ORDER:
                aa_cols = [col for col in df.columns if col.endswith(f'({aa})')]
                aa_cols.sort()  # 同一氨基酸内的密码子按字母序排列
                columns.extend(aa_cols)
            df = df[columns]
        
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        
        return df
    
    def get_notes(self):
        """获取数据质量报告"""
        return ''.join(self.notes)


class MainApp:
    """主应用程序界面"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("CodonW密码子偏好性统计分析工具")
        self.root.geometry("800x700")
        
        self.analyzer = CodonWAnalyzer()
        
        self._setup_ui()
    
    def _setup_ui(self):
        """设置用户界面"""
        # 主框架
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 配置网格权重
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        main_frame.rowconfigure(3, weight=1)
        
        # 标题
        title_label = ttk.Label(main_frame, text="CodonW密码子偏好性统计分析工具",
                               font=('Arial', 14, 'bold'))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # 输入文件夹选择
        ttk.Label(main_frame, text="输入文件夹：").grid(row=1, column=0, sticky=tk.W)
        self.input_entry = ttk.Entry(main_frame, width=50)
        self.input_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(main_frame, text="浏览...", command=self._browse_input).grid(row=1, column=2, padx=5)
        
        # 输出文件夹选择
        ttk.Label(main_frame, text="输出文件夹：").grid(row=2, column=0, sticky=tk.W)
        self.output_entry = ttk.Entry(main_frame, width=50)
        self.output_entry.grid(row=2, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Button(main_frame, text="浏览...", command=self._browse_output).grid(row=2, column=2, padx=5)
        
        # 度量指标选择
        metric_frame = ttk.Frame(main_frame)
        metric_frame.grid(row=3, column=0, columnspan=3, pady=10, sticky=tk.W)
        
        ttk.Label(metric_frame, text="统计指标：").pack(side=tk.LEFT)
        self.metric_var = tk.StringVar(value='RSCU')
        ttk.Radiobutton(metric_frame, text="RSCU (相对同义密码子使用度)",
                        variable=self.metric_var, value='RSCU').pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(metric_frame, text="Count (密码子计数)",
                        variable=self.metric_var, value='Count').pack(side=tk.LEFT, padx=10)
        
        # 排序选项
        sort_frame = ttk.Frame(main_frame)
        sort_frame.grid(row=4, column=0, columnspan=3, pady=5, sticky=tk.W)
        
        self.sort_by_aa_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(sort_frame, text="按氨基酸排序密码子列",
                        variable=self.sort_by_aa_var).pack(side=tk.LEFT)
        
        # 分析按钮
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=5, column=0, columnspan=3, pady=20)
        
        ttk.Button(button_frame, text="开始分析", command=self._analyze,
                  style='Accent.TButton').pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="清空", command=self._clear).pack(side=tk.LEFT, padx=5)
        
        # 进度标签
        self.status_label = ttk.Label(main_frame, text="等待输入...", foreground='blue')
        self.status_label.grid(row=6, column=0, columnspan=3, pady=5)
        
        # 结果预览
        ttk.Label(main_frame, text="数据质量报告：").grid(row=7, column=0, sticky=tk.W)
        self.text_area = scrolledtext.ScrolledText(main_frame, height=15, wrap=tk.WORD)
        self.text_area.grid(row=8, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        main_frame.rowconfigure(8, weight=1)
    
    def _browse_input(self):
        """浏览输入文件夹"""
        folder = filedialog.askdirectory(title="选择包含CodonW结果的文件夹")
        if folder:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, folder)
    
    def _browse_output(self):
        """浏览输出文件夹"""
        folder = filedialog.askdirectory(title="选择输出文件夹")
        if folder:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, folder)
    
    def _analyze(self):
        """执行分析"""
        input_folder = self.input_entry.get()
        output_folder = self.output_entry.get()
        metric = self.metric_var.get()
        sort_by_aa = self.sort_by_aa_var.get()
        
        # 验证输入
        if not input_folder:
            messagebox.showerror("错误", "请选择输入文件夹！")
            return
        
        if not os.path.exists(input_folder):
            messagebox.showerror("错误", "输入文件夹不存在！")
            return
        
        if not output_folder:
            output_folder = input_folder
            self.output_entry.insert(0, output_folder)
        
        try:
            self.status_label.config(text="正在分析...", foreground='orange')
            self.root.update()
            
            # 分析数据
            file_count = self.analyzer.analyze_folder(input_folder)
            
            # 生成热图表格
            suffix = f"_sorted_{metric.lower()}.csv" if sort_by_aa else f"_{metric.lower()}.csv"
            output_file = os.path.join(output_folder, f"codonw_heatmap{suffix}")
            df = self.analyzer.create_heatmap_table(output_file, metric, sort_by_aa)
            
            # 生成数据质量报告
            notes_file = os.path.join(output_folder, "codonw_data_quality_report.txt")
            with open(notes_file, 'w', encoding='utf-8') as f:
                f.write(self.analyzer.get_notes())
            
            # 显示结果
            self.text_area.delete(1.0, tk.END)
            self.text_area.insert(tk.END, self.analyzer.get_notes())
            
            sort_info = "（按氨基酸排序）" if sort_by_aa else ""
            self.status_label.config(
                text=f"分析完成！处理了 {file_count} 个文件，结果已保存到 {output_folder}",
                foreground='green'
            )
            
            messagebox.showinfo(
                "成功",
                f"分析完成{sort_info}！\n\n"
                f"处理文件数：{file_count}\n"
                f"热图表格：{output_file}\n"
                f"数据质量报告：{notes_file}"
            )
            
        except Exception as e:
            self.status_label.config(text="分析失败！", foreground='red')
            messagebox.showerror("错误", f"分析过程中发生错误：\n{str(e)}")
    
    def _clear(self):
        """清空界面"""
        self.input_entry.delete(0, tk.END)
        self.output_entry.delete(0, tk.END)
        self.text_area.delete(1.0, tk.END)
        self.status_label.config(text="等待输入...", foreground='blue')


def main():
    """主函数"""
    root = tk.Tk()
    
    # 设置样式
    style = ttk.Style()
    style.configure('Accent.TButton', font=('Arial', 10, 'bold'))
    
    app = MainApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
