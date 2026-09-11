import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
import re


class BEDConverter:
    def __init__(self, root):
        self.root = root
        self.root.title("BED/GFF3文件坐标转换工具")
        self.root.geometry("900x700")
        
        self.input_file_path = ""
        self.output_data = []
        self.file_type = "bed"  # bed or gff3
        
        self.create_widgets()
    
    def create_widgets(self):
        # 标题
        title_label = tk.Label(
            self.root,
            text="BED/GFF3文件坐标转换工具 - 将每个基因内部坐标改为从0开始",
            font=("Arial", 14, "bold")
        )
        title_label.pack(pady=10)
        
        # 文件类型选择
        type_frame = ttk.LabelFrame(self.root, text="文件类型")
        type_frame.pack(fill="x", padx=10, pady=5)
        
        self.file_type_var = tk.StringVar(value="bed")
        ttk.Radiobutton(type_frame, text="BED格式", variable=self.file_type_var, 
                        value="bed", command=self.on_file_type_change).pack(side="left", padx=20, pady=5)
        ttk.Radiobutton(type_frame, text="GFF3格式", variable=self.file_type_var, 
                        value="gff3", command=self.on_file_type_change).pack(side="left", padx=20, pady=5)
        
        # 文件选择框架
        file_frame = ttk.LabelFrame(self.root, text="文件选择")
        file_frame.pack(fill="x", padx=10, pady=5)
        
        # 输入文件
        input_frame = ttk.Frame(file_frame)
        input_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(input_frame, text="输入文件:").pack(side="left", padx=5)
        self.input_entry = ttk.Entry(input_frame)
        self.input_entry.pack(side="left", fill="x", expand=True, padx=5)
        
        ttk.Button(input_frame, text="浏览", command=self.browse_input).pack(side="left", padx=5)
        
        # 操作按钮
        button_frame = ttk.Frame(self.root)
        button_frame.pack(pady=10)
        
        ttk.Button(button_frame, text="加载文件", command=self.load_file).pack(side="left", padx=5)
        ttk.Button(button_frame, text="转换", command=self.convert).pack(side="left", padx=5)
        ttk.Button(button_frame, text="保存文件", command=self.save_file).pack(side="left", padx=5)
        ttk.Button(button_frame, text="清空", command=self.clear_all).pack(side="left", padx=5)
        
        # 统计信息
        self.stats_label = ttk.Label(self.root, text="")
        self.stats_label.pack(pady=5)
        
        # 预览区域
        preview_frame = ttk.LabelFrame(self.root, text="转换结果预览")
        preview_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.preview_text = scrolledtext.ScrolledText(preview_frame, wrap=tk.NONE)
        self.preview_text.pack(fill="both", expand=True, padx=5, pady=5)
        
        # 添加水平和垂直滚动条
        h_scroll = ttk.Scrollbar(self.preview_text, orient="horizontal", command=self.preview_text.xview)
        h_scroll.pack(side="bottom", fill="x")
        self.preview_text.configure(xscrollcommand=h_scroll.set)
    
    def on_file_type_change(self):
        self.file_type = self.file_type_var.get()
        # 清空当前数据
        self.input_entry.delete(0, tk.END)
        self.output_data = []
        self.preview_text.delete(1.0, tk.END)
        self.stats_label.config(text="")
    
    def browse_input(self):
        if self.file_type == "bed":
            filetypes = [("BED文件", "*.bed *.BED"), ("文本文件", "*.txt"), ("所有文件", "*.*")]
            title = "选择BED文件"
        else:
            filetypes = [("GFF3文件", "*.gff *.gff3"), ("文本文件", "*.txt"), ("所有文件", "*.*")]
            title = "选择GFF3文件"
        
        file_path = filedialog.askopenfilename(
            title=title,
            filetypes=filetypes
        )
        if file_path:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, file_path)
            self.input_file_path = file_path
    
    def load_file(self):
        self.input_file_path = self.input_entry.get().strip()
        
        if not self.input_file_path:
            messagebox.showwarning("警告", "请先选择输入文件！")
            return
        
        if not os.path.exists(self.input_file_path):
            messagebox.showerror("错误", "文件不存在！")
            return
        
        try:
            with open(self.input_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(tk.END, content)
            file_type_name = "BED" if self.file_type == "bed" else "GFF3"
            self.stats_label.config(text=f"原始{file_type_name}文件已加载")
            
        except Exception as e:
            messagebox.showerror("错误", f"加载文件失败: {str(e)}")
    
    def convert(self):
        if self.file_type == "bed":
            self.convert_bed()
        else:
            self.convert_gff3()
    
    def convert_bed(self):
        self.input_file_path = self.input_entry.get().strip()
        
        if not self.input_file_path:
            messagebox.showwarning("警告", "请先选择输入文件！")
            return
        
        if not os.path.exists(self.input_file_path):
            messagebox.showerror("错误", "文件不存在！")
            return
        
        try:
            # 读取文件
            with open(self.input_file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # 解析数据
            gene_data = {}
            comments = []
            
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    comments.append(line)
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                gene_id = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                feature_type = parts[3]
                phase = parts[4] if len(parts) > 4 else ''
                
                if gene_id not in gene_data:
                    gene_data[gene_id] = {
                        'gene_start': None,
                        'gene_end': None,
                        'features': []
                    }
                
                # 更新基因的起始和结束位置
                if gene_data[gene_id]['gene_start'] is None:
                    gene_data[gene_id]['gene_start'] = start
                    gene_data[gene_id]['gene_end'] = end
                else:
                    gene_data[gene_id]['gene_start'] = min(gene_data[gene_id]['gene_start'], start)
                    gene_data[gene_id]['gene_end'] = max(gene_data[gene_id]['gene_end'], end)
                
                # 保存特征信息
                gene_data[gene_id]['features'].append({
                    'start': start,
                    'end': end,
                    'type': feature_type,
                    'phase': phase,
                    'original_line': parts
                })
            
            # 转换坐标
            self.output_data = []
            gene_count = 0
            feature_count = 0
            
            for gene_id, data in gene_data.items():
                gene_start = data['gene_start']
                gene_count += 1
                
                # 按原始顺序处理特征
                for feature in data['features']:
                    new_start = feature['start'] - gene_start
                    new_end = feature['end'] - gene_start
                    
                    new_line = f"{gene_id}\t{new_start}\t{new_end}\t{feature['type']}"
                    if feature['phase']:
                        new_line += f"\t{feature['phase']}"
                    
                    self.output_data.append(new_line)
                    feature_count += 1
            
            # 显示预览
            self.preview_text.delete(1.0, tk.END)
            
            # 添加注释行
            for comment in comments:
                self.preview_text.insert(tk.END, comment + "\n")
            
            # 添加转换后的数据
            for line in self.output_data:
                self.preview_text.insert(tk.END, line + "\n")
            
            self.stats_label.config(text=f"转换完成: {gene_count} 个基因, {feature_count} 个特征")
            messagebox.showinfo("成功", f"转换完成！\n基因数量: {gene_count}\n特征数量: {feature_count}")
            
        except Exception as e:
            messagebox.showerror("错误", f"转换失败: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def convert_gff3(self):
        self.input_file_path = self.input_entry.get().strip()
        
        if not self.input_file_path:
            messagebox.showwarning("警告", "请先选择输入文件！")
            return
        
        if not os.path.exists(self.input_file_path):
            messagebox.showerror("错误", "文件不存在！")
            return
        
        try:
            # 读取文件
            with open(self.input_file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # 解析数据
            gene_data = {}
            gene_id_map = {}  # 存储gene行，方便后续查找
            transcript_map = {}  # 存储转录本与基因的映射关系
            
            for line_idx, line in enumerate(lines):
                line = line.rstrip('\n')
                
                # 跳过空行和注释行（但保存注释行）
                if not line.strip():
                    continue
                
                if line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 8:
                    continue
                
                seqid = parts[0]
                source = parts[1]
                feature_type = parts[2]
                start = int(parts[3]) - 1  # GFF3是1-based，转为0-based
                end = int(parts[4])
                score = parts[5]
                strand = parts[6]
                phase = parts[7]
                attributes = parts[8] if len(parts) > 8 else ''
                
                # 解析属性
                attrs_dict = self.parse_gff3_attributes(attributes)
                
                gene_id = attrs_dict.get('gene_id') or attrs_dict.get('gene', '')
                transcript_id = attrs_dict.get('transcript_id') or attrs_dict.get('Parent', '')
                
                # 存储gene行
                if feature_type == 'gene' and gene_id:
                    if gene_id not in gene_id_map:
                        gene_id_map[gene_id] = []
                    gene_id_map[gene_id].append({
                        'line_idx': line_idx,
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'original_line': line
                    })
                
                # 建立转录本与基因的映射
                if feature_type == 'mRNA' or feature_type == 'transcript':
                    if transcript_id and gene_id:
                        transcript_map[transcript_id] = gene_id
                
                # 处理特征数据
                if gene_id:
                    if gene_id not in gene_data:
                        gene_data[gene_id] = {
                            'gene_start': None,
                            'gene_end': None,
                            'seqid': seqid,
                            'strand': strand,
                            'features': []
                        }
                    
                    # 更新基因的起始和结束位置
                    if gene_data[gene_id]['gene_start'] is None:
                        gene_data[gene_id]['gene_start'] = start
                        gene_data[gene_id]['gene_end'] = end
                    else:
                        gene_data[gene_id]['gene_start'] = min(gene_data[gene_id]['gene_start'], start)
                        gene_data[gene_id]['gene_end'] = max(gene_data[gene_id]['gene_end'], end)
                    
                    # 保存特征信息
                    gene_data[gene_id]['features'].append({
                        'line_idx': line_idx,
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'type': feature_type,
                        'score': score,
                        'strand': strand,
                        'phase': phase,
                        'attributes': attributes,
                        'original_line': line
                    })
            
            # 转换坐标
            self.output_data = []
            gene_count = 0
            feature_count = 0
            
            # 按原始行号排序
            for gene_id, data in gene_data.items():
                gene_start = data['gene_start']
                gene_count += 1
                
                # 按原始顺序处理特征
                data['features'].sort(key=lambda x: x['line_idx'])
                
                for feature in data['features']:
                    new_start = feature['start'] - gene_start
                    new_end = feature['end'] - gene_start
                    
                    # 转换为1-based (GFF3标准)
                    new_start_gff = new_start + 1
                    new_end_gff = new_end
                    
                    # 构建新的行
                    new_line = "\t".join([
                        feature['seqid'],
                        feature['score'],  # source列，使用原来的score值或'.'
                        feature['type'],
                        str(new_start_gff),
                        str(new_end_gff),
                        feature['score'],
                        feature['strand'],
                        feature['phase'],
                        feature['attributes']
                    ])
                    
                    self.output_data.append({
                        'line_idx': feature['line_idx'],
                        'content': new_line
                    })
                    feature_count += 1
            
            # 按原始行号排序输出
            self.output_data.sort(key=lambda x: x['line_idx'])
            
            # 重建完整文件内容（包括注释行和空行）
            converted_lines = []
            output_dict = {item['line_idx']: item['content'] for item in self.output_data}
            
            for line_idx, line in enumerate(lines):
                line_stripped = line.strip()
                
                # 保留注释行
                if line_stripped.startswith('#'):
                    converted_lines.append(line.rstrip('\n'))
                # 保留空行
                elif not line_stripped:
                    converted_lines.append('')
                # 转换数据行
                elif line_idx in output_dict:
                    converted_lines.append(output_dict[line_idx])
            
            self.output_content = converted_lines
            
            # 显示预览
            self.preview_text.delete(1.0, tk.END)
            for line in converted_lines:
                self.preview_text.insert(tk.END, line + "\n")
            
            self.stats_label.config(text=f"GFF3转换完成: {gene_count} 个基因, {feature_count} 个特征")
            messagebox.showinfo("成功", f"GFF3转换完成！\n基因数量: {gene_count}\n特征数量: {feature_count}")
            
        except Exception as e:
            messagebox.showerror("错误", f"转换失败: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def parse_gff3_attributes(self, attributes_str):
        """解析GFF3属性字段"""
        attrs_dict = {}
        if not attributes_str or attributes_str == '.':
            return attrs_dict
        
        parts = attributes_str.split(';')
        for part in parts:
            part = part.strip()
            if '=' in part:
                key, value = part.split('=', 1)
                attrs_dict[key] = value
        
        return attrs_dict
    
    def save_file(self):
        if self.file_type == "bed":
            if not self.output_data:
                messagebox.showwarning("警告", "没有可保存的数据，请先进行转换！")
                return
        else:
            if not hasattr(self, 'output_content') or not self.output_content:
                messagebox.showwarning("警告", "没有可保存的数据，请先进行转换！")
                return
        
        # 根据文件类型设置默认扩展名
        if self.file_type == "bed":
            default_extension = ".bed"
            initialfile = "基因成员_转换后.bed"
            filetypes = [("BED文件", "*.bed"), ("文本文件", "*.txt"), ("所有文件", "*.*")]
            title = "保存转换后的BED文件"
        else:
            default_extension = ".gff3"
            initialfile = "转换后.gff3"
            filetypes = [("GFF3文件", "*.gff3 *.gff"), ("文本文件", "*.txt"), ("所有文件", "*.*")]
            title = "保存转换后的GFF3文件"
        
        file_path = filedialog.asksaveasfilename(
            title=title,
            defaultextension=default_extension,
            initialfile=initialfile,
            filetypes=filetypes
        )
        
        if not file_path:
            return
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                if self.file_type == "bed":
                    # BED格式保存
                    with open(self.input_file_path, 'r', encoding='utf-8') as orig_f:
                        original_lines = orig_f.readlines()
                    
                    # 提取注释行
                    comments = []
                    for line in original_lines:
                        if line.strip().startswith('#'):
                            comments.append(line)
                    
                    # 写入注释行
                    for comment in comments:
                        f.write(comment)
                    
                    # 写入转换后的数据
                    for line in self.output_data:
                        f.write(line + "\n")
                else:
                    # GFF3格式保存
                    for line in self.output_content:
                        if line == '':
                            f.write("\n")
                        else:
                            f.write(line + "\n")
            
            messagebox.showinfo("成功", f"文件已保存到:\n{file_path}")
            
        except Exception as e:
            messagebox.showerror("错误", f"保存文件失败: {str(e)}")
    
    def clear_all(self):
        self.input_entry.delete(0, tk.END)
        self.input_file_path = ""
        self.output_data = []
        if hasattr(self, 'output_content'):
            delattr(self, 'output_content')
        self.preview_text.delete(1.0, tk.END)
        self.stats_label.config(text="")


def main():
    root = tk.Tk()
    app = BEDConverter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
