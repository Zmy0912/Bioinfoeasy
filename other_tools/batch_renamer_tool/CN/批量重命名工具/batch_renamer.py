import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
from typing import Dict, Tuple


class BatchRenamerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("批量重命名工具")
        self.root.geometry("1000x800")
        
        # 文件路径
        self.name_mapping_file = tk.StringVar()
        self.target_file = tk.StringVar()
        self.output_file = tk.StringVar()
        
        # 名称映射字典
        self.name_mapping = {}
        
        # 创建界面
        self.create_widgets()
        
    def create_widgets(self):
        # 标题
        title_frame = tk.Frame(self.root, pady=10)
        title_frame.pack(fill=tk.X)
        tk.Label(title_frame, text="批量重命名工具", font=("Arial", 16, "bold")).pack()
        
        # 文件选择区域
        file_frame = tk.LabelFrame(self.root, text="文件选择", padx=10, pady=10)
        file_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # 名称映射文件选择
        tk.Label(file_frame, text="名称映射文件(TXT):").grid(row=0, column=0, sticky="w", pady=5)
        tk.Entry(file_frame, textvariable=self.name_mapping_file, width=50).grid(row=0, column=1, sticky="ew", padx=5)
        tk.Button(file_frame, text="浏览...", command=self.select_name_mapping_file).grid(row=0, column=2, padx=5)
        
        # 目标文件选择
        tk.Label(file_frame, text="目标文件:").grid(row=1, column=0, sticky="w", pady=5)
        tk.Entry(file_frame, textvariable=self.target_file, width=50).grid(row=1, column=1, sticky="ew", padx=5)
        tk.Button(file_frame, text="浏览...", command=self.select_target_file).grid(row=1, column=2, padx=5)
        
        # 输出文件选择
        tk.Label(file_frame, text="输出文件:").grid(row=2, column=0, sticky="w", pady=5)
        tk.Entry(file_frame, textvariable=self.output_file, width=50).grid(row=2, column=1, sticky="ew", padx=5)
        tk.Button(file_frame, text="浏览...", command=self.select_output_file).grid(row=2, column=2, padx=5)
        
        file_frame.columnconfigure(1, weight=1)
        
        # 操作按钮
        button_frame = tk.Frame(self.root, pady=10)
        button_frame.pack(fill=tk.X)
        tk.Button(button_frame, text="加载名称映射", command=self.load_name_mapping, 
                 bg="#4CAF50", fg="white", font=("Arial", 10, "bold"), padx=20).pack(side=tk.LEFT, padx=10)
        tk.Button(button_frame, text="预览替换", command=self.preview_replacement,
                 bg="#2196F3", fg="white", font=("Arial", 10, "bold"), padx=20).pack(side=tk.LEFT, padx=10)
        tk.Button(button_frame, text="执行替换", command=self.execute_replacement,
                 bg="#FF9800", fg="white", font=("Arial", 10, "bold"), padx=20).pack(side=tk.LEFT, padx=10)
        
        # 预览区域
        preview_frame = tk.LabelFrame(self.root, text="预览区域", padx=10, pady=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # 创建分屏显示
        paned_window = tk.PanedWindow(preview_frame, orient=tk.HORIZONTAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        
        # 左侧:原始文件预览
        left_frame = tk.Frame(paned_window)
        tk.Label(left_frame, text="原始文件预览", font=("Arial", 10, "bold")).pack(anchor="w")
        self.original_text = scrolledtext.ScrolledText(left_frame, wrap=tk.WORD, height=20)
        self.original_text.pack(fill=tk.BOTH, expand=True)
        paned_window.add(left_frame, minsize=400)
        
        # 右侧:替换后预览
        right_frame = tk.Frame(paned_window)
        tk.Label(right_frame, text="替换后预览", font=("Arial", 10, "bold")).pack(anchor="w")
        self.preview_text = scrolledtext.ScrolledText(right_frame, wrap=tk.WORD, height=20)
        self.preview_text.pack(fill=tk.BOTH, expand=True)
        paned_window.add(right_frame, minsize=400)
        
        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = tk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor="w")
        status_bar.pack(fill=tk.X, side=tk.BOTTOM)
        
        # 统计信息
        self.stats_var = tk.StringVar(value="")
        stats_frame = tk.Frame(self.root, pady=5)
        stats_frame.pack(fill=tk.X)
        tk.Label(stats_frame, textvariable=self.stats_var, fg="blue").pack()
    
    def select_name_mapping_file(self):
        filepath = filedialog.askopenfilename(
            title="选择名称映射文件",
            filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")]
        )
        if filepath:
            self.name_mapping_file.set(filepath)
            # 自动加载
            self.load_name_mapping()
    
    def select_target_file(self):
        filepath = filedialog.askopenfilename(
            title="选择目标文件",
            filetypes=[("所有文件", "*.*")]
        )
        if filepath:
            self.target_file.set(filepath)
            # 自动预览
            self.preview_target_file()
            # 如果已加载映射,自动预览替换
            if self.name_mapping:
                self.preview_replacement()
    
    def select_output_file(self):
        filepath = filedialog.asksaveasfilename(
            title="选择输出文件",
            filetypes=[("文本文件", "*.txt"), ("所有文件", "*.*")],
            defaultextension=".txt"
        )
        if filepath:
            self.output_file.set(filepath)
    
    def load_name_mapping(self):
        mapping_file = self.name_mapping_file.get()
        if not mapping_file:
            messagebox.showwarning("警告", "请先选择名称映射文件")
            return
        
        if not os.path.exists(mapping_file):
            messagebox.showerror("错误", "文件不存在")
            return
        
        try:
            self.name_mapping = {}
            with open(mapping_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    old_name = parts[0].strip()
                    new_name = parts[1].strip()
                    if old_name and new_name:
                        self.name_mapping[old_name] = new_name
            
            if self.name_mapping:
                self.status_var.set(f"成功加载 {len(self.name_mapping)} 条名称映射")
                messagebox.showinfo("成功", f"成功加载 {len(self.name_mapping)} 条名称映射")
                
                # 如果已选择目标文件,自动预览
                if self.target_file.get():
                    self.preview_replacement()
            else:
                messagebox.showwarning("警告", "未找到有效的名称映射")
                
        except Exception as e:
            messagebox.showerror("错误", f"加载文件失败: {str(e)}")
    
    def preview_target_file(self):
        target_file = self.target_file.get()
        if not target_file or not os.path.exists(target_file):
            return
        
        try:
            with open(target_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            self.original_text.delete(1.0, tk.END)
            self.original_text.insert(1.0, content)
            
        except Exception as e:
            messagebox.showerror("错误", f"读取目标文件失败: {str(e)}")
    
    def preview_replacement(self):
        if not self.name_mapping:
            messagebox.showwarning("警告", "请先加载名称映射")
            return
        
        target_file = self.target_file.get()
        if not target_file:
            messagebox.showwarning("警告", "请先选择目标文件")
            return
        
        if not os.path.exists(target_file):
            messagebox.showerror("错误", "目标文件不存在")
            return
        
        try:
            with open(target_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 显示原始内容
            self.original_text.delete(1.0, tk.END)
            self.original_text.insert(1.0, content)
            
            # 执行替换预览
            preview_content, replaced_count = self.replace_names(content)
            
            # 显示预览结果
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, preview_content)
            
            # 更新统计信息
            self.stats_var.set(f"预览: 共替换 {replaced_count} 处"
                             f" | 名称映射: {len(self.name_mapping)} 条")
            self.status_var.set("预览完成")
            
        except Exception as e:
            messagebox.showerror("错误", f"预览替换失败: {str(e)}")
    
    def replace_names(self, content: str) -> Tuple[str, int]:
        """
        替换名称,返回替换后的内容和替换次数
        """
        replaced_count = 0
        result = content
        
        # 按照旧名称长度排序,从长到短,避免部分匹配问题
        sorted_names = sorted(self.name_mapping.keys(), key=len, reverse=True)
        
        for old_name in sorted_names:
            new_name = self.name_mapping[old_name]
            if old_name in result:
                count = result.count(old_name)
                result = result.replace(old_name, new_name)
                replaced_count += count
        
        return result, replaced_count
    
    def execute_replacement(self):
        if not self.name_mapping:
            messagebox.showwarning("警告", "请先加载名称映射")
            return
        
        target_file = self.target_file.get()
        if not target_file:
            messagebox.showwarning("警告", "请先选择目标文件")
            return
        
        output_file = self.output_file.get()
        if not output_file:
            messagebox.showwarning("警告", "请先选择输出文件")
            return
        
        if not os.path.exists(target_file):
            messagebox.showerror("错误", "目标文件不存在")
            return
        
        try:
            # 确认对话框
            result = messagebox.askyesno(
                "确认", 
                f"即将执行替换操作:\n"
                f"输入文件: {target_file}\n"
                f"输出文件: {output_file}\n"
                f"名称映射: {len(self.name_mapping)} 条\n\n"
                f"是否继续?"
            )
            
            if not result:
                return
            
            # 读取目标文件
            with open(target_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 执行替换
            new_content, replaced_count = self.replace_names(content)
            
            # 写入输出文件
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(new_content)
            
            # 更新预览
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, new_content)
            
            # 更新统计信息
            self.stats_var.set(f"完成: 共替换 {replaced_count} 处"
                             f" | 名称映射: {len(self.name_mapping)} 条")
            self.status_var.set(f"替换完成! 文件已保存到: {output_file}")
            
            messagebox.showinfo("成功", f"替换完成!\n\n共替换 {replaced_count} 处\n文件已保存到:\n{output_file}")
            
        except Exception as e:
            messagebox.showerror("错误", f"执行替换失败: {str(e)}")


def main():
    root = tk.Tk()
    app = BatchRenamerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
