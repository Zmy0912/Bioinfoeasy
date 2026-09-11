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
    current_full_header = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line_stripped = line.strip()
                if line_stripped.startswith('>'):
                    if current_id is not None:
                        sequences[current_id] = {
                            'sequence': ''.join(current_seq),
                            'full_header': ''.join(current_full_header)
                        }
                    current_full_header = [line]
                    current_id = line_stripped[1:].split()[0]
                    current_seq = []
                else:
                    if line_stripped:
                        current_seq.append(line_stripped)
            
            if current_id is not None:
                sequences[current_id] = {
                    'sequence': ''.join(current_seq),
                    'full_header': ''.join(current_full_header)
                }
                
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
            sequences[file_name] = {
                'sequence': sequence,
                'full_header': f">{file_name}\n"
            }
            
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
    for i, seq_id in enumerate(sequences.keys(), 1):
        seq_data = sequences[seq_id]
        seq_preview = seq_data['sequence'][:40] + '...' if len(seq_data['sequence']) > 40 else seq_data['sequence']
        print(f"{i}. {seq_id}")
        print(f"   长度: {len(seq_data['sequence'])} bp, 预览: {seq_preview}")


def analyze_id_format(sequences):
    ids = list(sequences.keys())
    
    if not ids:
        return None
    
    first_id = ids[0]
    
    patterns = {
        'has_numbers': bool(re.search(r'\d', first_id)),
        'has_underscores': '_' in first_id,
        'has_pipe': '|' in first_id,
        'has_space': ' ' in first_id,
        'has_colon': ':' in first_id,
        'has_hyphen': '-' in first_id,
        'length': len(first_id),
        'starts_with_upper': first_id[0].isupper() if first_id else False,
        'all_upper': first_id.isupper(),
        'all_alpha': first_id.replace('_', '').replace('-', '').replace('|', '').isalpha(),
    }
    
    lengths = set(len(id) for id in ids)
    has_numbers = all(any(c.isdigit() for c in id) for id in ids)
    
    print("\n" + "=" * 60)
    print("序列ID格式分析:")
    print("=" * 60)
    print(f"  序列总数: {len(ids)}")
    print(f"  ID长度: {lengths}")
    print(f"  包含数字: {has_numbers}")
    print(f"  包含下划线: {patterns['has_underscores']}")
    print(f"  包含竖线: {patterns['has_pipe']}")
    print(f"  包含冒号: {patterns['has_colon']}")
    print(f"  全大写: {patterns['all_upper']}")
    print(f"  首字母大写: {patterns['starts_with_upper']}")
    
    return patterns


def generate_consistent_names(sequences, name_mapping):
    old_ids = list(sequences.keys())
    new_ids = {}
    
    for i, old_id in enumerate(old_ids):
        new_name = name_mapping.get(old_id, old_id)
        
        old_length = len(old_id)
        
        if len(new_name) != old_length:
            if len(new_name) > old_length:
                new_name = new_name[:old_length]
            elif len(new_name) < old_length:
                suffix = str(i + 1)
                new_name = new_name + suffix[:old_length - len(new_name)]
        
        new_ids[old_id] = new_name
    
    return new_ids


def rename_with_prefix_length(sequences, prefix):
    old_ids = list(sequences.keys())
    
    lengths = set(len(id) for id in old_ids)
    
    if len(lengths) == 1:
        target_length = list(lengths)[0]
    else:
        from collections import Counter
        length_counter = Counter(len(id) for id in old_ids)
        target_length = length_counter.most_common(1)[0][0]
        print(f"\n注意: 原ID长度不一致，将使用最常见长度: {target_length}")
    
    new_ids = {}
    for i, old_id in enumerate(old_ids):
        if len(prefix) >= target_length:
            new_name = prefix[:target_length]
        else:
            num_str = str(i + 1)
            new_name = prefix + num_str
            new_name = new_name[:target_length]
        
        new_ids[old_id] = new_name
    
    return new_ids


def rename_interactive_with_length_check(sequences):
    old_ids = list(sequences.keys())
    
    lengths = [len(id) for id in old_ids]
    unique_lengths = set(lengths)
    
    print(f"\n原ID长度信息:")
    print(f"  最小长度: {min(lengths)}")
    print(f"  最大长度: {max(lengths)}")
    print(f"  唯一长度值: {unique_lengths}")
    
    keep_length = input("\n是否保持原ID长度一致？(y/n, 默认y): ").strip().lower()
    if keep_length != 'n':
        keep_length = True
    else:
        keep_length = False
    
    new_ids = {}
    
    for i, old_id in enumerate(old_ids):
        seq_data = sequences[old_id]
        seq_preview = seq_data['sequence'][:30] + '...' if len(seq_data['sequence']) > 30 else seq_data['sequence']
        
        print(f"\n序列 {i+1}/{len(old_ids)}:")
        print(f"  原ID: {old_id} (长度: {len(old_id)})")
        print(f"  序列长度: {len(seq_data['sequence'])} bp")
        print(f"  预览: {seq_preview}")
        
        new_name = input(f"  新ID: ").strip()
        
        if not new_name:
            new_name = old_id
            print(f"  保持原ID: {new_name}")
        elif keep_length and len(new_name) != len(old_id):
            if len(new_name) > len(old_id):
                new_name = new_name[:len(old_id)]
                print(f"  截断至原长度: {new_name}")
            else:
                suffix = str(i + 1)
                new_name = new_name + suffix[:len(old_id) - len(new_name)]
                print(f"  填充至原长度: {new_name}")
        
        new_ids[old_id] = new_name
    
    return new_ids


def write_renamed_sequences(sequences, new_ids, output_file, output_format='fasta'):
    try:
        with open(output_file, 'w', encoding='utf-8') as out_f:
            for old_id, new_id in new_ids.items():
                seq_data = sequences[old_id]
                sequence = seq_data['sequence']
                
                if output_format.lower() in ['fasta', 'fa']:
                    out_f.write(f">{new_id}\n{sequence}\n\n")
                elif output_format.lower() == 'txt':
                    out_f.write(f"=== 序列: {new_id} ===\n")
                    out_f.write(f"原ID: {old_id}\n")
                    out_f.write(f"长度: {len(sequence)} bp\n")
                    for i in range(0, len(sequence), 80):
                        out_f.write(sequence[i:i+80] + '\n')
                    out_f.write('\n')
                else:
                    out_f.write(f">{new_id}\n{sequence}\n\n")
        
        return True
        
    except Exception as e:
        print(f"写入文件错误: {e}")
        return False


def main():
    print("=" * 60)
    print("基因序列重命名工具（位置格式保持）")
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
    sequences = get_sequences_from_file(input_file)
    
    if not sequences:
        print("错误: 文件中没有找到有效的序列")
        return
    
    list_all_sequences(sequences)
    
    patterns = analyze_id_format(sequences)
    print("\n" + "=" * 60)
    print("请选择重命名方式:")
    print("=" * 60)
    print("  1. 逐个手动输入新名称（保持原长度）")
    print("  2. 使用前缀 + 数字（自动匹配原长度）")
    print("  3. 使用前缀 + 填充零数字（自动匹配原长度）")
    print("  4. 使用指定名称列表（保持原长度）")
    print("  5. 保持原名不变")
    
    choice = input("\n请输入选择 (1-5): ").strip()
    
    new_ids = {}
    
    if choice == '1':
        new_ids = rename_interactive_with_length_check(sequences)
    
    elif choice == '2':
        prefix = input("请输入前缀 (如: ABC): ").strip().upper()
        if not prefix:
            prefix = 'SEQ'
        new_ids = rename_with_prefix_length(sequences, prefix)
    
    elif choice == '3':
        prefix = input("请输入前缀 (如: XYZ): ").strip().upper()
        if not prefix:
            prefix = 'SEQ'
        
        old_ids = list(sequences.keys())
        lengths = set(len(id) for id in old_ids)
        if len(lengths) == 1:
            target_length = list(lengths)[0]
        else:
            from collections import Counter
            target_length = Counter(len(id) for id in old_ids).most_common(1)[0][0]
        
        new_ids = {}
        for i, old_id in enumerate(old_ids):
            num_str = str(i + 1).zfill(target_length - len(prefix))
            new_ids[old_id] = prefix + num_str
    
    elif choice == '4':
        print(f"\n请输入 {len(sequences)} 个新名称（每行一个）:")
        print("直接回车结束输入:")
        
        new_names = []
        for i in range(len(sequences)):
            name = input(f"  名称 {i+1}/{len(sequences)}: ").strip().upper()
            if name:
                new_names.append(name)
            else:
                break
        
        old_ids = list(sequences.keys())
        new_ids = {}
        for i, old_id in enumerate(old_ids):
            if i < len(new_names):
                new_name = new_names[i]
                if len(new_name) > len(old_id):
                    new_name = new_name[:len(old_id)]
                elif len(new_name) < len(old_id):
                    suffix = str(i + 1)
                    new_name = new_name + suffix[:len(old_id) - len(new_name)]
                new_ids[old_id] = new_name
            else:
                new_ids[old_id] = old_id
    
    elif choice == '5':
        new_ids = {k: k for k in sequences.keys()}
    
    else:
        print("无效选择，保持原名不变")
        new_ids = {k: k for k in sequences.keys()}
    
    # 显示重命名预览
    print("\n" + "=" * 60)
    print("重命名预览（检查长度一致性）:")
    print("=" * 60)
    length_match_count = 0
    for i, (old_id, new_id) in enumerate(new_ids.items(), 1):
        if old_id != new_id:
            length_match = len(old_id) == len(new_id)
            if length_match:
                length_match_count += 1
            match_symbol = "✓" if length_match else "✗"
            print(f"{i}. {old_id} (长度:{len(old_id)}) -> {new_id} (长度:{len(new_id)}) {match_symbol}")
        else:
            print(f"{i}. {old_id} (未修改)")
    
    print(f"\n长度一致性: {length_match_count}/{len(new_ids)} 个ID保持原长度")
    
    confirm = input("\n确认应用这些重命名吗？(y/n): ").strip().lower()
    if confirm != 'y':
        print("操作已取消")
        return
    default_output = os.path.join(os.path.dirname(input_file), "renamed_sequences.fasta")
    output_file = input(f"\n请输入输出文件路径 (默认: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    print("\n请选择输出格式:")
    print("  1. FASTA 格式 (.fasta/.fa) - 推荐")
    print("  2. 纯文本格式 (.txt) - 便于查看（包含原ID信息）")
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
    print("正在保存重命名后的序列...")
    print("=" * 60)
    
    success = write_renamed_sequences(sequences, new_ids, output_file, output_format)
    
    if success:
        print("\n" + "=" * 60)
        print("重命名完成!")
        print("=" * 60)
        print(f"\n输出文件: {output_file}")
        print(f"处理序列数: {len(sequences)}")
        print(f"长度一致性: {length_match_count}/{len(new_ids)}")
    else:
        print("\n重命名失败，请检查输入。")


if __name__ == "__main__":
    main()
