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


def complement_sequence(sequence):
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'N': 'N', 'n': 'n'
    }
    
    complement = []
    for base in sequence:
        complement.append(complement_map.get(base, base))
    
    return ''.join(complement)


def reverse_complement(sequence):
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'N': 'N', 'n': 'n'
    }
    
    reverse_seq = sequence[::-1]
    complement = []
    for base in reverse_seq:
        complement.append(complement_map.get(base, base))
    
    return ''.join(complement)


def list_all_sequences(sequences):
    print("\n" + "=" * 60)
    print(f"当前文件中的序列列表 (共 {len(sequences)} 个):")
    print("=" * 60)
    for i, (seq_id, seq) in enumerate(sequences.items(), 1):
        seq_preview = seq[:40] + '...' if len(seq) > 40 else seq
        print(f"{i}. {seq_id}")
        print(f"   长度: {len(seq)} bp, 预览: {seq_preview}")
        print()


def complement_sequences(input_file, seq_ids, output_file, output_format='fasta', complement_type='rc'):
    all_sequences = get_sequences_from_file(input_file)
    
    if not all_sequences:
        print("错误: 文件中没有找到有效的序列")
        return False
    
    list_all_sequences(all_sequences)
    
    if seq_ids:
        selected_ids = seq_ids
    else:
        print("\n请输入要互补的序列编号（用空格或逗号分隔，如: 1 3 5 或 all）:")
        user_input = input("> ").strip()
        
        if user_input.lower() == 'all':
            selected_ids = list(all_sequences.keys())
        else:
            nums = re.split(r'[,，\s]+', user_input)
            selected_ids = []
            for num in nums:
                if num.strip():
                    try:
                        idx = int(num.strip()) - 1
                        if 0 <= idx < len(all_sequences):
                            selected_ids.append(list(all_sequences.keys())[idx])
                    except ValueError:
                        if num.strip() in all_sequences:
                            selected_ids.append(num.strip())
    
    if not selected_ids:
        print("错误: 没有选择任何序列")
        return False
    
    print("\n请选择互补类型:")
    print("  1. 仅互补 (A<->T, G<->C)")
    print("  2. 反向互补 (推荐，生成反向互补链)")
    
    comp_choice = input("\n请输入选择 (1-2, 默认2): ").strip()
    
    if comp_choice == '1':
        complement_type = 'complement'
        comp_name = "互补"
    else:
        complement_type = 'rc'
        comp_name = "反向互补"
    
    keep_original = input("\n是否在输出文件中保留原始序列？(y/n, 默认n): ").strip().lower()
    keep_original = (keep_original == 'y')
    
    try:
        with open(output_file, 'w', encoding='utf-8') as out_f:
            if keep_original:
                for seq_id in selected_ids:
                    if seq_id not in all_sequences:
                        continue
                    sequence = all_sequences[seq_id]
                    out_f.write(f">{seq_id}\n{sequence}\n\n")
            
            for seq_id in selected_ids:
                if seq_id not in all_sequences:
                    print(f"警告: 序列 '{seq_id}' 不存在，已跳过")
                    continue
                
                original_seq = all_sequences[seq_id]
                
                if complement_type == 'complement':
                    complement_seq = complement_sequence(original_seq)
                else:
                    complement_seq = reverse_complement(original_seq)
                
                new_id = seq_id + "CDNA"
                
                if output_format.lower() in ['fasta', 'fa']:
                    out_f.write(f">{new_id}\n{complement_seq}\n\n")
                elif output_format.lower() == 'txt':
                    out_f.write(f"=== 序列: {new_id} ===\n")
                    out_f.write(f"原序列ID: {seq_id}\n")
                    out_f.write(f"处理类型: {comp_name}\n")
                    out_f.write(f"原序列长度: {len(original_seq)} bp\n")
                    out_f.write(f"互补序列长度: {len(complement_seq)} bp\n")
                    for i in range(0, len(complement_seq), 80):
                        out_f.write(complement_seq[i:i+80] + '\n')
                    out_f.write('\n')
                else:
                    out_f.write(f">{new_id}\n{complement_seq}\n\n")
        
        print(f"\n成功处理 {len(selected_ids)} 个序列")
        print(f"处理类型: {comp_name}")
        if keep_original:
            print(f"原始序列: 已保留")
        else:
            print(f"原始序列: 未保留")
        print(f"输出文件: {output_file}")
        return True
        
    except Exception as e:
        print(f"写入文件错误: {e}")
        return False


def main():
    print("=" * 60)
    print("基因序列互补工具")
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
    
    default_output = os.path.join(os.path.dirname(input_file), "complemented_sequences.fasta")
    output_file = input(f"请输入输出文件路径 (默认: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    
    print("\n请选择输出格式:")
    print("  1. FASTA 格式 (.fasta/.fa) - 推荐")
    print("  2. 纯文本格式 (.txt) - 便于查看（包含原序列信息）")
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
    
    success = complement_sequences(input_file, [], output_file, output_format)
    
    if success:
        print("\n" + "=" * 60)
        print("互补处理完成!")
        print("=" * 60)
    else:
        print("\n处理失败，请检查输入。")


if __name__ == "__main__":
    main()
