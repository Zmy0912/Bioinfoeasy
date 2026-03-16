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
        print(f"Error reading file: {e}")
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
        print(f"Error reading file: {e}")
    
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
    print(f"Sequence list in current file (Total: {len(sequences)}):")
    print("=" * 60)
    for i, (seq_id, seq) in enumerate(sequences.items(), 1):
        seq_preview = seq[:40] + '...' if len(seq) > 40 else seq
        print(f"{i}. {seq_id}")
        print(f"   Length: {len(seq)} bp, Preview: {seq_preview}")
        print()


def complement_sequences(input_file, seq_ids, output_file, output_format='fasta', complement_type='rc'):
    all_sequences = get_sequences_from_file(input_file)
    
    if not all_sequences:
        print("Error: No valid sequences found in file")
        return False
    
    list_all_sequences(all_sequences)
    
    if seq_ids:
        selected_ids = seq_ids
    else:
        print("\nPlease enter the sequence numbers to complement (separated by space or comma, e.g., 1 3 5 or all):")
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
        print("Error: No sequences selected")
        return False
    
    print("\nPlease select complement type:")
    print("  1. Complement only (A<->T, G<->C)")
    print("  2. Reverse complement (Recommended, generates reverse complementary strand)")

    comp_choice = input("\nPlease enter your choice (1-2, default 2): ").strip()

    if comp_choice == '1':
        complement_type = 'complement'
        comp_name = "Complement"
    else:
        complement_type = 'rc'
        comp_name = "Reverse complement"
    
    keep_original = input("\nKeep original sequences in output file? (y/n, default n): ").strip().lower()
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
                    out_f.write(f"=== Sequence: {new_id} ===\n")
                    out_f.write(f"Original sequence ID: {seq_id}\n")
                    out_f.write(f"Processing type: {comp_name}\n")
                    out_f.write(f"Original sequence length: {len(original_seq)} bp\n")
                    out_f.write(f"Complement sequence length: {len(complement_seq)} bp\n")
                    for i in range(0, len(complement_seq), 80):
                        out_f.write(complement_seq[i:i+80] + '\n')
                    out_f.write('\n')
                else:
                    out_f.write(f">{new_id}\n{complement_seq}\n\n")
        
        print(f"\nSuccessfully processed {len(selected_ids)} sequences")
        print(f"Processing type: {comp_name}")
        if keep_original:
            print(f"Original sequences: Kept")
        else:
            print(f"Original sequences: Not kept")
        print(f"Output file: {output_file}")
        return True
        
    except Exception as e:
        print(f"Error writing file: {e}")
        return False


def main():
    print("=" * 60)
    print("DNA Sequence Complement Tool")
    print("=" * 60)
    
    input_file = input("\nPlease enter sequence file path: ").strip()
    input_file = input_file.strip('"').strip("'")
    
    if not os.path.exists(input_file):
        print(f"Error: File not found - {input_file}")
        return
    
    file_ext = Path(input_file).suffix.lower()
    supported_ext = ['.txt', '.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']
    
    if file_ext not in supported_ext:
        print(f"Warning: File type '{file_ext}' may not be a sequence file, attempting to continue...")
    
    default_output = os.path.join(os.path.dirname(input_file), "complemented_sequences.fasta")
    output_file = input(f"Please enter output file path (default: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    
    print("\nPlease select output format:")
    print("  1. FASTA format (.fasta/.fa) - Recommended")
    print("  2. Plain text format (.txt) - Easy to view (includes original sequence information)")
    print("  3. FASTA format (80 characters per line) - Suitable for further analysis")
    
    format_choice = input("\nPlease enter your choice (1-3, default 1): ").strip()
    
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
    print("Starting to read file...")
    print("=" * 60)
    
    success = complement_sequences(input_file, [], output_file, output_format)
    
    if success:
        print("\n" + "=" * 60)
        print("Complement processing completed!")
        print("=" * 60)
    else:
        print("\nProcessing failed, please check input.")


if __name__ == "__main__":
    main()
