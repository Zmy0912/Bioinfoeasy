#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
from pathlib import Path


def parse_fasta_headers(file_path):
    sequences_info = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            current_header = None
            current_seq = []
            
            for line in f:
                line_stripped = line.strip()
                if line_stripped.startswith('>'):
                    if current_header is not None:
                        sequences_info.append({
                            'header': current_header,
                            'sequence': ''.join(current_seq),
                            'seq_length': len(''.join(current_seq))
                        })
                    current_header = line_stripped[1:]
                    current_seq = []
                else:
                    if line_stripped:
                        current_seq.append(line_stripped)
            
            if current_header is not None:
                sequences_info.append({
                    'header': current_header,
                    'sequence': ''.join(current_seq),
                    'seq_length': len(''.join(current_seq))
                })
                
    except Exception as e:
        print(f"读取文件错误: {e}")
        return []
    
    return sequences_info


def parse_txt_file(file_path):
    sequences_info = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
        if '>' in content:
            return parse_fasta_headers(file_path)
        
        file_name = Path(file_path).stem
        sequence = re.sub(r'[^ATCGN\s]', '', content.upper())
        sequence = ''.join(sequence.split())
        
        if sequence:
            sequences_info.append({
                'header': file_name,
                'sequence': sequence,
                'seq_length': len(sequence)
            })
            
    except Exception as e:
        print(f"读取文件错误: {e}")
    
    return sequences_info


def get_sequence_info(file_path):
    file_ext = Path(file_path).suffix.lower()
    
    if file_ext in ['.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']:
        return parse_fasta_headers(file_path)
    else:
        return parse_txt_file(file_path)


def extract_id_and_description(header):
    parts = re.split(r'\s+', header, maxsplit=1)
    
    if len(parts) >= 2:
        seq_id = parts[0]
        description = parts[1]
    else:
        seq_id = header
        description = ""
    
    return seq_id, description


def analyze_sequences(sequences_info):
    print("\n" + "=" * 60)
    print("序列信息统计:")
    print("=" * 60)
    print(f"  总序列数: {len(sequences_info)}")
    
    if sequences_info:
        lengths = [seq['seq_length'] for seq in sequences_info]
        print(f"  最短序列: {min(lengths)} bp")
        print(f"  最长序列: {max(lengths)} bp")
        print(f"  平均长度: {sum(lengths) // len(lengths)} bp")
        
        has_desc = sum(1 for seq in sequences_info
                      if extract_id_and_description(seq['header'])[1])
        print(f"  有描述的序列: {has_desc} / {len(sequences_info)}")


def display_sequences_preview(sequences_info, limit=10):
    print("\n" + "=" * 60)
    print(f"序列预览 (显示前 {min(limit, len(sequences_info))} 个):")
    print("=" * 60)
    
    for i, seq_info in enumerate(sequences_info[:limit], 1):
        header = seq_info['header']
        seq_id, description = extract_id_and_description(header)
        
        print(f"\n序列 {i}:")
        print(f"  ID: {seq_id}")
        if description:
            print(f"  描述: {description[:60]}{'...' if len(description) > 60 else ''}")
        print(f"  序列长度: {seq_info['seq_length']} bp")
        print(f"  完整头部: {header[:80]}{'...' if len(header) > 80 else ''}")
    
    if len(sequences_info) > limit:
        print(f"\n... 还有 {len(sequences_info) - limit} 个序列")


def extract_to_text_file(sequences_info, output_file, options):
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("序列标识符及补充信息提取结果\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("统计信息:\n")
            f.write("-" * 40 + "\n")
            f.write(f"总序列数: {len(sequences_info)}\n")
            
            if sequences_info:
                lengths = [seq['seq_length'] for seq in sequences_info]
                f.write(f"最短序列: {min(lengths)} bp\n")
                f.write(f"最长序列: {max(lengths)} bp\n")
                f.write(f"平均长度: {sum(lengths) // len(lengths)} bp\n")
            
            f.write("\n")
            
            f.write("=" * 80 + "\n")
            f.write("序列详细信息:\n")
            f.write("=" * 80 + "\n\n")
            
            for i, seq_info in enumerate(sequences_info, 1):
                header = seq_info['header']
                seq_id, description = extract_id_and_description(header)
                
                f.write(f"【序列 {i}】\n")
                f.write(f"序列ID: {seq_id}\n")
                
                if description:
                    f.write(f"描述信息: {description}\n")
                else:
                    f.write(f"描述信息: (无)\n")
                
                f.write(f"序列长度: {seq_info['seq_length']} bp\n")
                
                if options['show_full_header']:
                    f.write(f"完整头部: {header}\n")
                
                if options['show_sequence_preview']:
                    seq_preview = seq_info['sequence']
                    if len(seq_preview) > 100:
                        seq_preview = seq_preview[:100] + "..."
                    f.write(f"序列预览: {seq_preview}\n")
                
                f.write("\n")
            
            if options['include_table']:
                f.write("\n" + "=" * 80 + "\n")
                f.write("ID与描述对照表:\n")
                f.write("=" * 80 + "\n")
                f.write(f"{'序号':<8}{'ID':<30}{'描述':<30}\n")
                f.write("-" * 80 + "\n")
                for i, seq_info in enumerate(sequences_info, 1):
                    seq_id, description = extract_id_and_description(seq_info['header'])
                    desc_display = description[:28] + '..' if len(description) > 28 else description
                    f.write(f"{i:<8}{seq_id:<30}{desc_display:<30}\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("序列ID列表:\n")
            f.write("=" * 80 + "\n")
            for seq_info in sequences_info:
                seq_id = extract_id_and_description(seq_info['header'])[0]
                f.write(f"{seq_id}\n")
        
        return True
        
    except Exception as e:
        print(f"写入文件错误: {e}")
        return False


def main():
    print("=" * 60)
    print("序列标识符提取工具")
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
    
    print("\n正在读取文件...")
    sequences_info = get_sequence_info(input_file)
    
    if not sequences_info:
        print("错误: 文件中没有找到有效的序列")
        return
    
    analyze_sequences(sequences_info)
    
    display_sequences_preview(sequences_info)
    
    print("\n" + "=" * 60)
    print("请选择提取范围:")
    print("  1. 提取所有序列")
    print("  2. 指定序号提取")
    
    choice = input("\n请输入选择 (1-2, 默认1): ").strip()
    
    selected_sequences = sequences_info
    
    if choice == '2':
        user_input = input("请输入序号（用空格或逗号分隔，如: 1 3 5）: ").strip()
        
        nums = re.split(r'[,，\s]+', user_input)
        selected_indices = []
        for num in nums:
            if num.strip():
                try:
                    idx = int(num.strip()) - 1
                    if 0 <= idx < len(sequences_info):
                        selected_indices.append(idx)
                except ValueError:
                    pass
        
        if selected_indices:
            selected_sequences = [sequences_info[i] for i in selected_indices]
        else:
            print("未选择有效序号，将提取所有序列")
    
    print("\n" + "=" * 60)
    print("请选择输出内容:")
    print("  1. 仅ID和描述")
    print("  2. ID、描述和序列长度")
    print("  3. 完整信息（包含序列预览）")
    print("  4. 完整信息 + 对照表")
    print("  注: 序列ID列表将始终在文件末尾输出")
    
    output_choice = input("\n请输入选择 (1-4, 默认2): ").strip()
    
    options = {
        'show_full_header': False,
        'show_sequence_preview': False,
        'include_table': False
    }
    
    if output_choice == '1':
        options['show_full_header'] = False
        options['show_sequence_preview'] = False
    elif output_choice == '2':
        options['show_full_header'] = False
        options['show_sequence_preview'] = False
    elif output_choice == '3':
        options['show_full_header'] = True
        options['show_sequence_preview'] = True
    elif output_choice == '4':
        options['show_full_header'] = True
        options['show_sequence_preview'] = True
        options['include_table'] = True
    
    default_output = os.path.join(os.path.dirname(input_file), "sequence_info.txt")
    output_file = input(f"\n请输入输出文件路径 (默认: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    
    if not output_file.endswith('.txt'):
        output_file += '.txt'
    
    print("\n" + "=" * 60)
    print("正在提取信息...")
    print("=" * 60)
    
    success = extract_to_text_file(selected_sequences, output_file, options)
    
    if success:
        print("\n" + "=" * 60)
        print("提取完成!")
        print("=" * 60)
        print(f"\n输出文件: {output_file}")
        print(f"提取序列数: {len(selected_sequences)}")
    else:
        print("\n提取失败，请检查输入。")


if __name__ == "__main__":
    main()
