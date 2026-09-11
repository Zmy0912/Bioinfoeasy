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
        print(f"Error reading file: {e}")
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
        print(f"Error reading file: {e}")
    
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
    print("Sequence Information Statistics:")
    print("=" * 60)
    print(f"  Total sequences: {len(sequences_info)}")
    
    if sequences_info:
        lengths = [seq['seq_length'] for seq in sequences_info]
        print(f"  Shortest sequence: {min(lengths)} bp")
        print(f"  Longest sequence: {max(lengths)} bp")
        print(f"  Average length: {sum(lengths) // len(lengths)} bp")
        
        has_desc = sum(1 for seq in sequences_info
                      if extract_id_and_description(seq['header'])[1])
        print(f"  Sequences with description: {has_desc} / {len(sequences_info)}")


def display_sequences_preview(sequences_info, limit=10):
    print("\n" + "=" * 60)
    print(f"Sequence Preview (showing first {min(limit, len(sequences_info))}):")
    print("=" * 60)
    
    for i, seq_info in enumerate(sequences_info[:limit], 1):
        header = seq_info['header']
        seq_id, description = extract_id_and_description(header)
        
        print(f"\nSequence {i}:")
        print(f"  ID: {seq_id}")
        if description:
            print(f"  Description: {description[:60]}{'...' if len(description) > 60 else ''}")
        print(f"  Sequence length: {seq_info['seq_length']} bp")
        print(f"  Full header: {header[:80]}{'...' if len(header) > 80 else ''}")
    
    if len(sequences_info) > limit:
        print(f"\n... and {len(sequences_info) - limit} more sequences")


def extract_to_text_file(sequences_info, output_file, options):
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("Sequence Identifier and Additional Information Extraction Results\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("Statistics:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total sequences: {len(sequences_info)}\n")
            
            if sequences_info:
                lengths = [seq['seq_length'] for seq in sequences_info]
                f.write(f"Shortest sequence: {min(lengths)} bp\n")
                f.write(f"Longest sequence: {max(lengths)} bp\n")
                f.write(f"Average length: {sum(lengths) // len(lengths)} bp\n")
            
            f.write("\n")
            
            f.write("=" * 80 + "\n")
            f.write("Detailed Sequence Information:\n")
            f.write("=" * 80 + "\n\n")
            
            for i, seq_info in enumerate(sequences_info, 1):
                header = seq_info['header']
                seq_id, description = extract_id_and_description(header)
                
                f.write(f"[Sequence {i}]\n")
                f.write(f"Sequence ID: {seq_id}\n")
                
                if description:
                    f.write(f"Description: {description}\n")
                else:
                    f.write(f"Description: (None)\n")
                
                f.write(f"Sequence length: {seq_info['seq_length']} bp\n")
                
                if options['show_full_header']:
                    f.write(f"Full header: {header}\n")
                
                if options['show_sequence_preview']:
                    seq_preview = seq_info['sequence']
                    if len(seq_preview) > 100:
                        seq_preview = seq_preview[:100] + "..."
                    f.write(f"Sequence preview: {seq_preview}\n")
                
                f.write("\n")
            
            if options['include_table']:
                f.write("\n" + "=" * 80 + "\n")
                f.write("ID and Description Reference Table:\n")
                f.write("=" * 80 + "\n")
                f.write(f"{'No.':<8}{'ID':<30}{'Description':<30}\n")
                f.write("-" * 80 + "\n")
                for i, seq_info in enumerate(sequences_info, 1):
                    seq_id, description = extract_id_and_description(seq_info['header'])
                    desc_display = description[:28] + '..' if len(description) > 28 else description
                    f.write(f"{i:<8}{seq_id:<30}{desc_display:<30}\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("Sequence ID List:\n")
            f.write("=" * 80 + "\n")
            for seq_info in sequences_info:
                seq_id = extract_id_and_description(seq_info['header'])[0]
                f.write(f"{seq_id}\n")
        
        return True
        
    except Exception as e:
        print(f"Error writing file: {e}")
        return False


def main():
    print("=" * 60)
    print("Sequence Identifier Extraction Tool")
    print("=" * 60)
    
    input_file = input("\nPlease enter sequence file path: ").strip()
    input_file = input_file.strip('"').strip("'")
    
    if not os.path.exists(input_file):
        print(f"Error: File does not exist - {input_file}")
        return
    
    file_ext = Path(input_file).suffix.lower()
    supported_ext = ['.txt', '.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']
    
    if file_ext not in supported_ext:
        print(f"Warning: File type '{file_ext}' may not be a sequence file, attempting to continue...")
    
    print("\nReading file...")
    sequences_info = get_sequence_info(input_file)
    
    if not sequences_info:
        print("Error: No valid sequences found in file")
        return
    
    analyze_sequences(sequences_info)
    
    display_sequences_preview(sequences_info)
    
    print("\n" + "=" * 60)
    print("Please select extraction range:")
    print("  1. Extract all sequences")
    print("  2. Specify sequence numbers")
    
    choice = input("\nPlease enter choice (1-2, default 1): ").strip()
    
    selected_sequences = sequences_info
    
    if choice == '2':
        user_input = input("Please enter sequence numbers (separated by space or comma, e.g., 1 3 5): ").strip()
        
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
            print("No valid sequence numbers selected, extracting all sequences")
    
    print("\n" + "=" * 60)
    print("Please select output content:")
    print("  1. ID and description only")
    print("  2. ID, description and sequence length")
    print("  3. Full information (including sequence preview)")
    print("  4. Full information + reference table")
    print("  Note: Sequence ID list will always be output at the end of the file")
    
    output_choice = input("\nPlease enter choice (1-4, default 2): ").strip()
    
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
    output_file = input(f"\nPlease enter output file path (default: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    
    if not output_file.endswith('.txt'):
        output_file += '.txt'
    
    print("\n" + "=" * 60)
    print("Extracting information...")
    print("=" * 60)
    
    success = extract_to_text_file(selected_sequences, output_file, options)
    
    if success:
        print("\n" + "=" * 60)
        print("Extraction complete!")
        print("=" * 60)
        print(f"\nOutput file: {output_file}")
        print(f"Extracted sequences: {len(selected_sequences)}")
    else:
        print("\nExtraction failed, please check input.")


if __name__ == "__main__":
    main()
