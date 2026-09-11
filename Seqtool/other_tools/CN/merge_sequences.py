#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from pathlib import Path


def get_sequence_files(folder_path):
    sequence_extensions = ['.txt', '.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']
    sequence_files = []
    
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_ext = Path(file).suffix.lower()
            if file_ext in sequence_extensions:
                sequence_files.append(os.path.join(root, file))
    
    return sequence_files


def merge_sequences(input_folder, output_file, output_format='fasta'):
    sequence_files = get_sequence_files(input_folder)
    
    if not sequence_files:
        print(f"警告: 在 {input_folder} 中未找到任何序列文件")
        return False
    
    print(f"找到 {len(sequence_files)} 个序列文件:")
    for f in sequence_files:
        print(f"  - {os.path.basename(f)}")
    
    with open(output_file, 'w', encoding='utf-8') as out_f:
        for seq_file in sequence_files:
            print(f"正在处理: {os.path.basename(seq_file)}")
            
            try:
                with open(seq_file, 'r', encoding='utf-8') as in_f:
                    content = in_f.read()
                    
                    if output_format.lower() in ['fasta', 'fa']:
                        file_name = Path(seq_file).stem
                        
                        lines = content.strip().split('\n')
                        if lines and lines[0].startswith('>'):
                            out_f.write(content + '\n\n')
                        else:
                            sequence = ''.join(line.strip() for line in lines if line.strip())
                            if sequence:
                                out_f.write(f">{file_name}\n{sequence}\n\n")
                    
                    elif output_format.lower() == 'txt':
                        out_f.write(f"=== 文件: {os.path.basename(seq_file)} ===\n")
                        out_f.write(content + '\n\n')
                    
                    else:
                        out_f.write(content + '\n\n')
                        
            except Exception as e:
                print(f"  错误: 无法读取文件 - {e}")
                continue
    
    print(f"\n成功合并序列到: {output_file}")
    return True


def main():
    print("=" * 50)
    print("基因序列合并工具")
    print("=" * 50)
    
    input_folder = input("\n请输入包含序列文件的文件夹路径: ").strip()
    input_folder = input_folder.strip('"').strip("'")
    
    if not os.path.exists(input_folder):
        print(f"错误: 文件夹不存在 - {input_folder}")
        return
    
    default_output = os.path.join(input_folder, "merged_sequences.fasta")
    output_file = input(f"请输入输出文件路径 (默认: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    
    print("\n请选择输出格式:")
    print("  1. FASTA 格式 (.fasta/.fa) - 推荐")
    print("  2. 纯文本格式 (.txt)")
    print("  3. 原始格式 (保持原样)")
    
    format_choice = input("\n请输入选择 (1-3, 默认1): ").strip()
    
    format_map = {
        '1': 'fasta',
        '2': 'txt',
        '3': 'raw'
    }
    
    output_format = format_map.get(format_choice, 'fasta')
    
    if output_format == 'txt' and not Path(output_file).suffix:
        output_file += '.txt'
    elif output_format == 'fasta' and not Path(output_file).suffix:
        output_file += '.fasta'
    
    print("\n" + "=" * 50)
    print("开始合并序列...")
    print("=" * 50)
    
    success = merge_sequences(input_folder, output_file, output_format)
    
    if success:
        print("\n" + "=" * 50)
        print("合并完成!")
        print("=" * 50)
    else:
        print("\n合并失败，请检查输入。")


if __name__ == "__main__":
    main()
