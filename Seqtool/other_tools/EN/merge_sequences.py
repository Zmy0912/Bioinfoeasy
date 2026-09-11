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
        print(f"Warning: No sequence files found in {input_folder}")
        return False
    
    print(f"Found {len(sequence_files)} sequence files:")
    for f in sequence_files:
        print(f"  - {os.path.basename(f)}")
    
    with open(output_file, 'w', encoding='utf-8') as out_f:
        for seq_file in sequence_files:
            print(f"Processing: {os.path.basename(seq_file)}")
            
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
                        out_f.write(f"=== File: {os.path.basename(seq_file)} ===\n")
                        out_f.write(content + '\n\n')
                    
                    else:
                        out_f.write(content + '\n\n')
                        
            except Exception as e:
                print(f"  Error: Unable to read file - {e}")
                continue
    
    print(f"\nSuccessfully merged sequences to: {output_file}")
    return True


def main():
    print("=" * 50)
    print("Gene Sequence Merger Tool")
    print("=" * 50)
    
    input_folder = input("\nEnter the folder path containing sequence files: ").strip()
    input_folder = input_folder.strip('"').strip("'")
    
    if not os.path.exists(input_folder):
        print(f"Error: Folder does not exist - {input_folder}")
        return
    
    default_output = os.path.join(input_folder, "merged_sequences.fasta")
    output_file = input(f"Enter output file path (default: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    
    print("\nPlease select output format:")
    print("  1. FASTA format (.fasta/.fa) - Recommended")
    print("  2. Plain text format (.txt)")
    print("  3. Raw format (keep original)")
    
    format_choice = input("\nEnter your choice (1-3, default 1): ").strip()
    
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
    print("Starting sequence merge...")
    print("=" * 50)
    
    success = merge_sequences(input_folder, output_file, output_format)
    
    if success:
        print("\n" + "=" * 50)
        print("Merge completed!")
        print("=" * 50)
    else:
        print("\nMerge failed. Please check your input.")


if __name__ == "__main__":
    main()
