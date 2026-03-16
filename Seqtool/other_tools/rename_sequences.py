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
            sequences[file_name] = {
                'sequence': sequence,
                'full_header': f">{file_name}\n"
            }
            
    except Exception as e:
        print(f"Error reading file: {e}")
    
    return sequences


def get_sequences_from_file(file_path):
    file_ext = Path(file_path).suffix.lower()
    
    if file_ext in ['.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']:
        return read_fasta_file(file_path)
    else:
        return read_txt_file(file_path)


def list_all_sequences(sequences):
    print("\n" + "=" * 60)
    print(f"List of sequences in current file (Total: {len(sequences)}):")
    print("=" * 60)
    for i, seq_id in enumerate(sequences.keys(), 1):
        seq_data = sequences[seq_id]
        seq_preview = seq_data['sequence'][:40] + '...' if len(seq_data['sequence']) > 40 else seq_data['sequence']
        print(f"{i}. {seq_id}")
        print(f"   Length: {len(seq_data['sequence'])} bp, Preview: {seq_preview}")


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
    print("Sequence ID Format Analysis:")
    print("=" * 60)
    print(f"  Total sequences: {len(ids)}")
    print(f"  ID lengths: {lengths}")
    print(f"  Contains numbers: {has_numbers}")
    print(f"  Contains underscores: {patterns['has_underscores']}")
    print(f"  Contains pipes: {patterns['has_pipe']}")
    print(f"  Contains colons: {patterns['has_colon']}")
    print(f"  All uppercase: {patterns['all_upper']}")
    print(f"  First letter uppercase: {patterns['starts_with_upper']}")
    
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
        print(f"\nNote: Original ID lengths are inconsistent, will use most common length: {target_length}")
    
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
    
    print(f"\nOriginal ID length information:")
    print(f"  Minimum length: {min(lengths)}")
    print(f"  Maximum length: {max(lengths)}")
    print(f"  Unique length values: {unique_lengths}")
    
    keep_length = input("\nKeep original ID length consistent? (y/n, default y): ").strip().lower()
    if keep_length != 'n':
        keep_length = True
    else:
        keep_length = False
    
    new_ids = {}
    
    for i, old_id in enumerate(old_ids):
        seq_data = sequences[old_id]
        seq_preview = seq_data['sequence'][:30] + '...' if len(seq_data['sequence']) > 30 else seq_data['sequence']
        
        print(f"\nSequence {i+1}/{len(old_ids)}:")
        print(f"  Original ID: {old_id} (Length: {len(old_id)})")
        print(f"  Sequence length: {len(seq_data['sequence'])} bp")
        print(f"  Preview: {seq_preview}")
        
        new_name = input(f"  New ID: ").strip()
        
        if not new_name:
            new_name = old_id
            print(f"  Keep original ID: {new_name}")
        elif keep_length and len(new_name) != len(old_id):
            if len(new_name) > len(old_id):
                new_name = new_name[:len(old_id)]
                print(f"  Truncated to original length: {new_name}")
            else:
                suffix = str(i + 1)
                new_name = new_name + suffix[:len(old_id) - len(new_name)]
                print(f"  Padded to original length: {new_name}")
        
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
                    out_f.write(f"=== Sequence: {new_id} ===\n")
                    out_f.write(f"Original ID: {old_id}\n")
                    out_f.write(f"Length: {len(sequence)} bp\n")
                    for i in range(0, len(sequence), 80):
                        out_f.write(sequence[i:i+80] + '\n')
                    out_f.write('\n')
                else:
                    out_f.write(f">{new_id}\n{sequence}\n\n")
        
        return True
        
    except Exception as e:
        print(f"Error writing file: {e}")
        return False


def main():
    print("=" * 60)
    print("Gene Sequence Renaming Tool (Position Format Preserved)")
    print("=" * 60)
    
    input_file = input("\nPlease enter the sequence file path: ").strip()
    input_file = input_file.strip('"').strip("'")
    
    if not os.path.exists(input_file):
        print(f"Error: File does not exist - {input_file}")
        return
    
    file_ext = Path(input_file).suffix.lower()
    supported_ext = ['.txt', '.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']
    
    if file_ext not in supported_ext:
        print(f"Warning: File type '{file_ext}' may not be a sequence file, attempting to continue...")
    
    print("\nReading file...")
    sequences = get_sequences_from_file(input_file)
    
    if not sequences:
        print("Error: No valid sequences found in the file")
        return
    
    list_all_sequences(sequences)
    
    patterns = analyze_id_format(sequences)
    print("\n" + "=" * 60)
    print("Please select a renaming method:")
    print("=" * 60)
    print("  1. Manually enter new names one by one (keep original length)")
    print("  2. Use prefix + numbers (auto-match original length)")
    print("  3. Use prefix + zero-padded numbers (auto-match original length)")
    print("  4. Use specified name list (keep original length)")
    print("  5. Keep original names unchanged")
    
    choice = input("\nPlease enter your choice (1-5): ").strip()
    
    new_ids = {}
    
    if choice == '1':
        new_ids = rename_interactive_with_length_check(sequences)
    
    elif choice == '2':
        prefix = input("Please enter prefix (e.g., ABC): ").strip().upper()
        if not prefix:
            prefix = 'SEQ'
        new_ids = rename_with_prefix_length(sequences, prefix)
    
    elif choice == '3':
        prefix = input("Please enter prefix (e.g., XYZ): ").strip().upper()
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
        print(f"\nPlease enter {len(sequences)} new names (one per line):")
        print("Press Enter to finish input:")
        
        new_names = []
        for i in range(len(sequences)):
            name = input(f"  Name {i+1}/{len(sequences)}: ").strip().upper()
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
        print("Invalid choice, keeping original names unchanged")
        new_ids = {k: k for k in sequences.keys()}
    
    # Show renaming preview
    print("\n" + "=" * 60)
    print("Renaming Preview (Check Length Consistency):")
    print("=" * 60)
    length_match_count = 0
    for i, (old_id, new_id) in enumerate(new_ids.items(), 1):
        if old_id != new_id:
            length_match = len(old_id) == len(new_id)
            if length_match:
                length_match_count += 1
            match_symbol = "✓" if length_match else "✗"
            print(f"{i}. {old_id} (Length:{len(old_id)}) -> {new_id} (Length:{len(new_id)}) {match_symbol}")
        else:
            print(f"{i}. {old_id} (Unchanged)")
    
    print(f"\nLength consistency: {length_match_count}/{len(new_ids)} IDs maintain original length")
    
    confirm = input("\nConfirm to apply these renames? (y/n): ").strip().lower()
    if confirm != 'y':
        print("Operation cancelled")
        return
    default_output = os.path.join(os.path.dirname(input_file), "renamed_sequences.fasta")
    output_file = input(f"\nPlease enter output file path (default: {default_output}): ").strip()
    output_file = output_file.strip('"').strip("'")
    
    if not output_file:
        output_file = default_output
    print("\nPlease select output format:")
    print("  1. FASTA format (.fasta/.fa) - Recommended")
    print("  2. Plain text format (.txt) - Easy to view (includes original ID info)")
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
    print("Saving renamed sequences...")
    print("=" * 60)
    
    success = write_renamed_sequences(sequences, new_ids, output_file, output_format)
    
    if success:
        print("\n" + "=" * 60)
        print("Renaming completed!")
        print("=" * 60)
        print(f"\nOutput file: {output_file}")
        print(f"Number of sequences processed: {len(sequences)}")
        print(f"Length consistency: {length_match_count}/{len(new_ids)}")
    else:
        print("\nRenaming failed, please check input.")


if __name__ == "__main__":
    main()
