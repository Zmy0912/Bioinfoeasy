#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
from pathlib import Path


def read_fasta_file(file_path):
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    if line:
                        current_seq.append(line)
            
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
                
    except Exception as e:
        print(f"读取文件错误: {e}")
        return {}
    
    return sequences


def read_txt_file(file_path):
    sequences = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
        if '>' in content:
            return read_fasta_file(file_path)
        
        sequence = re.sub(r'[^ATCGN\s]', '', content.upper())
        sequence = ''.join(sequence.split())
        
        if sequence:
            file_name = Path(file_path).stem
            sequences[file_name] = sequence
            
    except Exception as e:
        print(f"读取文件错误: {e}")
    
    return sequences


def get_sequences_from_file(file_path):
    file_ext = Path(file_path).suffix.lower()
    
    if file_ext in ['.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']:
        return read_fasta_file(file_path)
    else:
        return read_txt_file(file_path)


def list_all_sequences(sequences):
    print("\n" + "=" * 60)
    print(f"当前文件中的序列列表 (共 {len(sequences)} 个):")
    print("=" * 60)
    for i, (seq_id, seq) in enumerate(sequences.items(), 1):
        seq_preview = seq[:50] + '...' if len(seq) > 50 else seq
        print(f"{i}. {seq_id}")
        print(f"   长度: {len(seq)} bp")
        print(f"   预览: {seq_preview}")
        print()


def remove_sequences(input_file, seq_ids_to_remove, output_file, output_format='fasta'):
    all_sequences = get_sequences_from_file(input_file)
    
    if not all_sequences:
        print("错误: 文件中没有找到有效的序列")
        return False
    
    list_all_sequences(all_sequences)
    
    if seq_ids_to_remove:
        ids_to_remove = seq_ids_to_remove
    else:
        print("\n请输入要删除的序列编号（用空格或逗号分隔，如: 1 3 5）:")
        user_input = input("> ").strip()
        
        nums = re.split(r'[,，\s]+', user_input)
        ids_to_remove = []
        for num in nums:
            if num.strip():
                try:
                    idx = int(num.strip()) - 1
                    if 0 <= idx < len(all_sequences):
                        ids_to_remove.append(list(all_sequences.keys())[idx])
                except ValueError:
                    if num.strip() in all_sequences:
                        ids_to_remove.append(num.strip())
    
    if not ids_to_remove:
        print("错误: 没有选择要删除的序列")
        return False
    
    print("\n" + "=" * 60)
    print(f"准备删除以下 {len(ids_to_remove)} 个序列:")
    print("=" * 60)
    for seq_id in ids_to_remove:
        if seq_id in all_sequences:
            print(f"  - {seq_id} (长度: {len(all_sequences[seq_id])} bp)")
    
    confirm = input("\n确认删除这些序列吗？(y/n): ").strip().lower()
    if confirm != 'y':
        print("操作已取消")
        return False
    
    remaining_sequences = {k: v for k, v in all_sequences.items() if k not in ids_to_remove}
    
    if not remaining_sequences:
        print("警告: 删除后没有剩余序列，将创建空文件")
    else:
        print(f"\n剩余序列数量: {len(remaining_sequences)}")
    
    try:
        with open(output_file, 'w', encoding='utf-8') as out_f:
            for seq_id, sequence in remaining_sequences.items():
                if output_format.lower() in ['fasta', 'fa']:
                    out_f.write(f">{seq_id}\n{sequence}\n\n")
                elif output_format.lower() == 'txt':
                    out_f.write(f"=== 序列: {seq_id} ===\n")
                    out_f.write(f"长度: {len(sequence)} bp\n")
                    for i in range(0, len(sequence), 80):
                        out_f.write(sequence[i:i+80] + '\n')
                    out_f.write('\n')
                else:
                    out_f.write(f">{seq_id}\n{sequence}\n\n")
        
        print(f"\n成功保存剩余序列到: {output_file}")
        return True
        
    except Exception as e:
        print(f"写入文件错误: {e}")
        return False


def main():
    print("=" * 60)
    print("基因序列删除工具")
    print("=" * 60)
    
    input_file = input("\n请输入序列文件路径: ").strip()
    input_file = input_file.strip('"').strip("'")
    
    if not os.path.exists(input_file):
        print(f"错误: 文件不存在 - {input_file}")
        return
    
    file_ext = Path(input_file).suffix.lower()
    supported_ext = ['.txt', '.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']
    
    if file_ext not in supported_ext:
        print(f"警告: 文件类型 '{file_ext}' 可能不是序列文件，尝试继续...")
    
    default_output = os.path.join(os.path.dirname(input_file), "remaining_sequences.fasta")
    output_file = input(f"请输入输出文件路径 (默认: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    
    print("\n请选择输出格式:")
    print("  1. FASTA 格式 (.fasta/.fa) - 推荐")
    print("  2. 纯文本格式 (.txt) - 便于查看")
    print("  3. FASTA 格式（每行80字符）- 适合进一步分析")
    
    format_choice = input("\n请输入选择 (1-3, 默认1): ").strip()
    
    format_map = {
        '1': 'fasta',
        '2': 'txt',
        '3': 'fasta'
    }
    
    output_format = format_map.get(format_choice, 'fasta')
    
    if output_format == 'txt' and not Path(output_file).suffix:
        output_file += '.txt'
    elif output_format == 'fasta' and not Path(output_file).suffix:
        output_file += '.fasta'
    
    print("\n" + "=" * 60)
    print("开始读取文件...")
    print("=" * 60)
    
    success = remove_sequences(input_file, [], output_file, output_format)
    
    if success:
        print("\n" + "=" * 60)
        print("操作完成!")
        print("=" * 60)
    else:
        print("\n操作失败或已取消。")


if __name__ == "__main__":
    main()
