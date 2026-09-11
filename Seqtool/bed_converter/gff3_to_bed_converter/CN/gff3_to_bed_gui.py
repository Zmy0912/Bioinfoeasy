#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GFF3 to BED Converter
一个具有图形化界面的GFF3转BED格式转换工具
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os


class GFF3ToBEDConverter:
    def __init__(self, root):
        self.root = root
        self.root.title("GFF3 to BED 转换器")
        self.root.geometry("800x750")
        self.gff3_file = ""
        self.bed_file = ""
        self.parsed_features = []
        self.gene_mapping = {}
        self.bed_preview_content = []
        
        self.setup_ui()
    
    def setup_ui(self):
        # 标题
        title_label = tk.Label(
            self.root,
            text="GFF3 to BED 格式转换器",
            font=("Microsoft YaHei", 16, "bold"),
            fg="#2c3e50",
            pady=10
        )
        title_label.pack()
        
        # 主框架
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # 文件选择区域
        file_frame = ttk.LabelFrame(main_frame, text="文件选择", padding="10")
        file_frame.pack(fill=tk.X, pady=5)
        
        # GFF3文件选择
        gff3_label = ttk.Label(file_frame, text="GFF3 文件:", font=("Microsoft YaHei", 10))
        gff3_label.grid(row=0, column=0, sticky=tk.W, pady=5)
        
        gff3_entry = ttk.Entry(file_frame, width=45)
        gff3_entry.grid(row=0, column=1, padx=5, pady=5)
        self.gff3_entry = gff3_entry
        
        gff3_btn = ttk.Button(
            file_frame,
            text="浏览",
            command=self.select_gff3_file
        )
        gff3_btn.grid(row=0, column=2, padx=5, pady=5)
        
        # BED文件保存
        bed_label = ttk.Label(file_frame, text="BED 文件:", font=("Microsoft YaHei", 10))
        bed_label.grid(row=1, column=0, sticky=tk.W, pady=5)
        
        bed_entry = ttk.Entry(file_frame, width=45)
        bed_entry.grid(row=1, column=1, padx=5, pady=5)
        self.bed_entry = bed_entry
        
        bed_btn = ttk.Button(
            file_frame,
            text="浏览",
            command=self.select_bed_file
        )
        bed_btn.grid(row=1, column=2, padx=5, pady=5)
        
        # 选项区域
        option_frame = ttk.LabelFrame(main_frame, text="转换选项", padding="10")
        option_frame.pack(fill=tk.X, pady=5)
        
        # 特征类型选择
        self.feature_types = {
            "gene": tk.BooleanVar(value=True),
            "mRNA": tk.BooleanVar(value=True),
            "exon": tk.BooleanVar(value=True),
            "CDS": tk.BooleanVar(value=True),
        }
        
        type_frame = ttk.Frame(option_frame)
        type_frame.pack(fill=tk.X, pady=5)
        
        row, col = 0, 0
        for feature, var in self.feature_types.items():
            checkbox = ttk.Checkbutton(type_frame, text=f"转换 {feature} 特征", variable=var)
            checkbox.grid(row=row, column=col, padx=10, pady=2, sticky=tk.W)
            col += 1
            if col > 1:
                col = 0
                row += 1
        
        # 基因名称选项
        name_frame = ttk.Frame(option_frame)
        name_frame.pack(fill=tk.X, pady=5)
        
        self.use_custom_name = tk.BooleanVar(value=False)
        custom_name_check = ttk.Checkbutton(
            name_frame,
            text="使用自定义基因名称(不使用GFF3中的名称)",
            variable=self.use_custom_name
        )
        custom_name_check.grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        # 转换按钮
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=10)
        
        preview_btn = ttk.Button(
            button_frame,
            text="预览输出",
            command=self.preview_output,
            width=15
        )
        preview_btn.pack(side=tk.LEFT, padx=5)
        
        convert_btn = ttk.Button(
            button_frame,
            text="开始转换",
            command=self.convert_gff3_to_bed,
            width=15
        )
        convert_btn.pack(side=tk.LEFT, padx=5)
        
        clear_btn = ttk.Button(
            button_frame,
            text="清空",
            command=self.clear_all,
            width=10
        )
        clear_btn.pack(side=tk.LEFT, padx=5)
        
        # 日志区域
        log_frame = ttk.LabelFrame(main_frame, text="转换日志", padding="10")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        self.log_text = scrolledtext.ScrolledText(
            log_frame,
            height=8,
            width=70,
            wrap=tk.WORD,
            font=("Consolas", 9)
        )
        self.log_text.pack(fill=tk.BOTH, expand=True)
        
        # 配置日志颜色
        self.log_text.tag_config("info", foreground="blue")
        self.log_text.tag_config("success", foreground="green")
        self.log_text.tag_config("error", foreground="red")
        self.log_text.tag_config("warning", foreground="orange")
    
    def select_gff3_file(self):
        file_path = filedialog.askopenfilename(
            title="选择GFF3文件",
            filetypes=[("GFF3文件", "*.gff3 *.gff"), ("所有文件", "*.*")]
        )
        if file_path:
            self.gff3_file = file_path
            self.gff3_entry.delete(0, tk.END)
            self.gff3_entry.insert(0, file_path)
            # 自动设置BED文件路径
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            bed_path = os.path.join(os.path.dirname(file_path), f"{base_name}.bed")
            self.bed_entry.delete(0, tk.END)
            self.bed_entry.insert(0, bed_path)
            self.bed_file = bed_path
            self.log(f"已选择GFF3文件: {file_path}", "info")
            # 自动解析文件以获取基因信息
            self.parse_gff3_file()
    
    def select_bed_file(self):
        file_path = filedialog.asksaveasfilename(
            title="保存BED文件",
            defaultextension=".bed",
            filetypes=[("BED文件", "*.bed"), ("所有文件", "*.*")]
        )
        if file_path:
            self.bed_file = file_path
            self.bed_entry.delete(0, tk.END)
            self.bed_entry.insert(0, file_path)
            self.log(f"已设置BED输出文件: {file_path}", "info")
    
    def parse_gff3_file(self):
        """解析GFF3文件"""
        self.parsed_features = []
        self.gene_mapping = {}
        
        try:
            with open(self.gff3_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    # 跳过注释行和空行
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        seqid = parts[0]
                        source = parts[1]
                        feature_type = parts[2]
                        start = int(parts[3])
                        end = int(parts[4])
                        score = parts[5]
                        strand = parts[6]
                        phase = parts[7]
                        attributes = parts[8]
                        
                        feature = {
                            'line_num': line_num,
                            'seqid': seqid,
                            'source': source,
                            'type': feature_type,
                            'start': start,
                            'end': end,
                            'score': score,
                            'strand': strand,
                            'phase': phase,
                            'attributes': attributes
                        }
                        
                        self.parsed_features.append(feature)
                        
                        # 提取基因ID
                        gene_id = self.extract_gene_id(feature)
                        if gene_id and feature_type == 'gene':
                            self.gene_mapping[gene_id] = gene_id
            
            self.log(f"成功解析 {len(self.parsed_features)} 个特征", "success")
            
            # 显示基因统计
            if self.gene_mapping:
                self.log(f"检测到 {len(self.gene_mapping)} 个基因:", "info")
                gene_list = sorted(self.gene_mapping.keys())
                for i, gene in enumerate(gene_list):
                    if i < 5:  # 只显示前5个
                        self.log(f"  - {gene}", "info")
                    elif i == 5:
                        self.log(f"  ... 还有 {len(gene_list) - 5} 个基因", "info")
                        break
            
        except Exception as e:
            self.log(f"解析GFF3文件出错: {str(e)}", "error")
            self.parsed_features = []
            self.gene_mapping = {}
    
    def extract_gene_id(self, feature):
        """从特征中提取基因ID"""
        attributes = feature['attributes']
        gene_id = ""
        
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('ID='):
                gene_id = attr.split('=', 1)[1]
                break
            elif attr.startswith('Name='):
                gene_id = attr.split('=', 1)[1]
                break
            elif attr.startswith('gene='):
                gene_id = attr.split('=', 1)[1]
                break
        
        return gene_id
    
    def find_gene_for_feature(self, feature, all_features):
        """为非基因特征找到对应的基因"""
        gene_id = ""
        gene_name = ""
        
        # 首先从Parent属性中找
        attributes = feature['attributes']
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('Parent=gene-'):
                gene_name = attr.split('=', 1)[1].replace('gene-', '')
                break
            elif attr.startswith('Parent='):
                parent_id = attr.split('=', 1)[1]
                # 查找这个父ID是否是基因
                for f in all_features:
                    if f['type'] == 'gene':
                        f_gene_id = self.extract_gene_id(f)
                        if f_gene_id == parent_id or parent_id.startswith('gene-' + f_gene_id):
                            gene_name = f_gene_id
                            break
                if gene_name:
                    break
        
        # 如果没找到,尝试从gene属性找
        if not gene_name:
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene='):
                    gene_name = attr.split('=', 1)[1]
                    break
        
        # 如果还是没找到,在所有基因中查找包含该特征的基因
        if not gene_name:
            feature_start = feature['start']
            feature_end = feature['end']
            feature_seqid = feature['seqid']
            
            for f in all_features:
                if f['type'] == 'gene' and f['seqid'] == feature_seqid:
                    # 检查特征是否在基因范围内
                    if (f['start'] <= feature_start <= f['end'] or
                        f['start'] <= feature_end <= f['end'] or
                        (feature_start <= f['start'] and feature_end >= f['end'])):
                        gene_name = self.extract_gene_id(f)
                        break
        
        return gene_name if gene_name else feature['seqid']
    
    def generate_bed_entries(self):
        """生成BED条目"""
        bed_entries = []
        
        selected_types = [ft for ft, var in self.feature_types.items() if var.get()]
        
        if not selected_types:
            self.log("错误: 请至少选择一种特征类型!", "error")
            return None
        
        # 为每个特征确定基因名称
        for feature in self.parsed_features:
            if feature['type'] not in selected_types:
                continue
            
            # 确定基因名称
            if feature['type'] == 'gene':
                gene_name = self.extract_gene_id(feature)
            else:
                gene_name = self.find_gene_for_feature(feature, self.parsed_features)
            
            # 起始和结束位置
            start = feature['start']
            end = feature['end']
            
            # 特征类型
            feature_type = feature['type']
            
            # phase (可选列)
            phase = feature['phase'] if feature['phase'] != '.' else ''
            
            # 构建BED条目
            if phase:
                bed_entry = f"{gene_name}\t{start}\t{end}\t{feature_type}\t{phase}"
            else:
                bed_entry = f"{gene_name}\t{start}\t{end}\t{feature_type}"
            
            bed_entries.append(bed_entry)
        
        return bed_entries
    
    def show_preview_window(self):
        """显示预览窗口"""
        preview_window = tk.Toplevel(self.root)
        preview_window.title("BED输出预览")
        preview_window.geometry("900x600")
        
        # 工具栏
        toolbar = ttk.Frame(preview_window, padding="5")
        toolbar.pack(fill=tk.X)
        
        save_btn = ttk.Button(
            toolbar,
            text="保存到文件",
            command=lambda: self.save_from_preview(preview_window)
        )
        save_btn.pack(side=tk.LEFT, padx=5)
        
        close_btn = ttk.Button(
            toolbar,
            text="关闭",
            command=preview_window.destroy
        )
        close_btn.pack(side=tk.LEFT, padx=5)
        
        # 信息标签
        info_label = ttk.Label(
            toolbar,
            text=f"共 {len(self.bed_preview_content)} 行数据 (显示前500行)",
            font=("Microsoft YaHei", 9)
        )
        info_label.pack(side=tk.LEFT, padx=20)
        
        # 预览文本框
        text_frame = ttk.Frame(preview_window)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        preview_text = scrolledtext.ScrolledText(
            text_frame,
            wrap=tk.NONE,
            font=("Consolas", 10)
        )
        preview_text.pack(fill=tk.BOTH, expand=True)
        
        # 显示内容
        preview_text.insert(tk.END, "# BED format file converted from GFF3\n")
        preview_text.insert(tk.END, "# Columns: geneID/transcriptID  start  end  featureType  phase(optional)\n")
        preview_text.insert(tk.END, "\n")
        
        display_lines = self.bed_preview_content[:500]  # 最多显示500行
        for line in display_lines:
            preview_text.insert(tk.END, line + "\n")
        
        if len(self.bed_preview_content) > 500:
            preview_text.insert(tk.END, f"\n... 还有 {len(self.bed_preview_content) - 500} 行未显示\n")
        
        preview_text.config(state=tk.DISABLED)
    
    def save_from_preview(self, preview_window):
        """从预览窗口保存"""
        if not self.bed_file:
            messagebox.showwarning("警告", "请先设置BED输出文件!")
            return
        
        if self.save_bed_file(self.bed_preview_content, self.bed_file):
            messagebox.showinfo("成功", f"转换成功!\nBED文件已保存到:\n{self.bed_file}")
            preview_window.destroy()
    
    def preview_output(self):
        """预览输出内容"""
        if not self.gff3_file:
            messagebox.showwarning("警告", "请先选择GFF3文件!")
            return
        
        if not self.parsed_features:
            self.log("正在解析GFF3文件...", "info")
            self.parse_gff3_file()
        
        if not self.parsed_features:
            messagebox.showerror("错误", "无法解析GFF3文件!")
            return
        
        self.log("正在生成BED预览...", "info")
        self.bed_preview_content = self.generate_bed_entries()
        
        if self.bed_preview_content is None:
            return
        
        self.log(f"生成了 {len(self.bed_preview_content)} 个BED条目", "success")
        
        # 显示预览窗口
        self.show_preview_window()
    
    def save_bed_file(self, bed_entries, file_path):
        """保存BED文件"""
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                # 写入BED文件头
                f.write("# BED format file converted from GFF3\n")
                f.write("# Columns: geneID/transcriptID  start  end  featureType  phase(optional)\n")
                f.write("\n")
                
                for entry in bed_entries:
                    f.write(entry + '\n')
            
            return True
        except Exception as e:
            self.log(f"保存BED文件出错: {str(e)}", "error")
            return False
    
    def convert_gff3_to_bed(self):
        """执行转换"""
        if not self.gff3_file:
            messagebox.showwarning("警告", "请先选择GFF3文件!")
            return
        
        if not self.bed_file:
            messagebox.showwarning("警告", "请设置BED输出文件!")
            return
        
        if not os.path.exists(self.gff3_file):
            messagebox.showerror("错误", f"GFF3文件不存在: {self.gff3_file}")
            return
        
        if not self.parsed_features:
            self.log("正在解析GFF3文件...", "info")
            self.parse_gff3_file()
        
        if not self.parsed_features:
            messagebox.showerror("错误", "无法解析GFF3文件!")
            return
        
        self.log("=" * 60, "info")
        self.log("开始转换...", "info")
        self.log(f"GFF3文件: {self.gff3_file}", "info")
        self.log(f"BED文件: {self.bed_file}", "info")
        
        # 统计特征类型
        type_count = {}
        for f in self.parsed_features:
            ftype = f['type']
            type_count[ftype] = type_count.get(ftype, 0) + 1
        
        self.log("特征统计:", "info")
        for ftype, count in sorted(type_count.items()):
            self.log(f"  {ftype}: {count}", "info")
        
        # 生成BED条目
        self.log("正在生成BED条目...", "info")
        bed_entries = self.generate_bed_entries()
        
        if bed_entries is None:
            return
        
        self.log(f"生成了 {len(bed_entries)} 个BED条目", "success")
        
        # 显示基因名称统计
        gene_names = set()
        for entry in bed_entries:
            gene_name = entry.split('\t')[0]
            gene_names.add(gene_name)
        
        self.log(f"共涉及 {len(gene_names)} 个不同的基因名称", "success")
        
        # 保存BED文件
        self.log(f"正在保存BED文件...", "info")
        if self.save_bed_file(bed_entries, self.bed_file):
            self.log(f"转换成功! BED文件已保存到:", "success")
            self.log(f"  {self.bed_file}", "success")
            messagebox.showinfo("成功", f"转换成功!\nBED文件已保存到:\n{self.bed_file}")
        else:
            messagebox.showerror("错误", "保存BED文件失败!")
        
        self.log("=" * 60, "info")
    
    def clear_all(self):
        """清空所有输入"""
        self.gff3_file = ""
        self.bed_file = ""
        self.parsed_features = []
        self.gene_mapping = {}
        self.bed_preview_content = []
        self.gff3_entry.delete(0, tk.END)
        self.bed_entry.delete(0, tk.END)
        self.log_text.delete(1.0, tk.END)
        self.log("已清空所有内容", "info")
    
    def log(self, message, level="info"):
        """添加日志"""
        self.log_text.insert(tk.END, message + "\n", level)
        self.log_text.see(tk.END)
        self.root.update()


def main():
    root = tk.Tk()
    app = GFF3ToBEDConverter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
