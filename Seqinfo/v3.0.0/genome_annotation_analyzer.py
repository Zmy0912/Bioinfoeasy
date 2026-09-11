#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genome Annotation File Analysis Tool
Supported formats: GenBank (.gb, .gbk), GFF (.gff, .gff3), FASTA (.fasta, .fa)
Extracted information: Sequence ID, species name, description, sequence length, gene statistics
"""

import os
import sys
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

# Supported annotation file extensions
SUPPORTED_EXTENSIONS = {'.gb', '.gbk', '.gff', '.gff3', '.fasta', '.fa'}


def is_annotation_file(file_path):
    """Check whether the file is a supported annotation file."""
    return os.path.splitext(file_path)[1].lower() in SUPPORTED_EXTENSIONS


# File types shown in the open-file dialog
ANNOTATION_FILE_TYPES = [
    ('Annotation files', '*.gb *.gbk *.gff *.gff3 *.fasta *.fa'),
    ('GenBank files', '*.gb *.gbk'),
    ('GFF files', '*.gff *.gff3'),
    ('FASTA files', '*.fasta *.fa'),
    ('All files', '*.*'),
]


class GenomeAnnotationAnalyzer:
    """Genome annotation file analyzer"""
    
    def __init__(self, input_file):
        self.input_file = input_file
        self.file_ext = os.path.splitext(input_file)[1].lower()
        self.data = {}
        
    def analyze(self):
        """Analyze annotation file"""
        if self.file_ext in ['.gb', '.gbk']:
            return self._analyze_genbank()
        elif self.file_ext in ['.gff', '.gff3']:
            return self._analyze_gff()
        elif self.file_ext in ['.fasta', '.fa']:
            return self._analyze_fasta()
        else:
            raise ValueError(f"Unsupported file format: {self.file_ext}")
    
    def _analyze_genbank(self):
        """Analyze GenBank format file"""
        print(f"Analyzing GenBank file: {self.input_file}")
        
        for record in SeqIO.parse(self.input_file, "genbank"):
            data = {}
            
            # Extract basic information
            data['Sequence ID'] = record.id
            data['Sequence Name'] = record.name
            data['Species Name'] = self._extract_organism(record)
            data['Description'] = record.description
            data['Genome Length'] = len(record.seq)
            data['Topology'] = record.annotations.get('topology', 'linear')
            data['Molecule Type'] = record.annotations.get('molecule_type', 'DNA')
            data['Taxonomy'] = self._extract_taxonomy(record)
            
            # Gene statistics
            gene_counts = defaultdict(int)
            rrna_count = 0
            trna_count = 0
            other_rna_count = 0
            cds_count = 0
            gene_count = 0
            
            # Iterate through all features
            for feature in record.features:
                feature_type = feature.type.lower()
                
                # Count gene types
                if feature_type == 'gene':
                    gene_count += 1
                    gene_type = feature.qualifiers.get('gene_type', ['unknown'])[0]
                    gene_counts[gene_type] += 1
                
                elif feature_type == 'cds':
                    cds_count += 1
                
                elif feature_type == 'rrna':
                    rrna_count += 1
                
                elif feature_type == 'trna':
                    trna_count += 1
                
                elif feature_type in ['rrna_gene', 'ribosomal_rna']:
                    rrna_count += 1
                
                elif feature_type in ['trna_gene', 'transfer_rna']:
                    trna_count += 1
                
                elif feature_type in ['ncrna', 'ncrna_gene', 'non_coding_rna', 'misc_rna']:
                    other_rna_count += 1
            
            # Fill gene statistics
            data['Total Genes'] = gene_count
            data['CDS Count'] = cds_count
            data['rRNA Count'] = rrna_count
            data['tRNA Count'] = trna_count
            data['Other ncRNA Count'] = other_rna_count
            
            # Other gene types
            for gene_type, count in gene_counts.items():
                if gene_type != 'unknown':
                    data[f'{gene_type} Count'] = count
            
            # Other annotation information
            data['Annotation Date'] = record.annotations.get('date', 'N/A')
            data['Data Version'] = record.annotations.get('data_file_division', 'N/A')
            data['Annotation Source'] = record.annotations.get('source', 'N/A')
            
            # GC content
            gc_count = record.seq.upper().count('G') + record.seq.upper().count('C')
            data['GC Content (%)'] = f"{(gc_count / len(record.seq) * 100):.2f}" if len(record.seq) > 0 else '0'
            
            self.data[record.id] = data
        
        return self.data
    
    def _analyze_gff(self):
        """Analyze GFF/GFF3 format file"""
        print(f"Analyzing GFF file: {self.input_file}")
        
        try:
            from BCBio import GFF
            use_bcbio = True
        except ImportError:
            print("Warning: BCBio.GFF library is required, using simplified parsing...")
            use_bcbio = False
        
        if use_bcbio:
            for rec in GFF.parse(self.input_file):
                data = {}
                
                data['Sequence ID'] = rec.id
                data['Genome Length'] = len(rec.seq) if rec.seq else 'N/A'
                data['Species Name'] = self._extract_organism_from_gff(rec)
                data['Description'] = rec.description
                
                # Gene statistics
                gene_counts = defaultdict(int)
                rrna_count = 0
                trna_count = 0
                other_rna_count = 0
                cds_count = 0
                gene_count = 0
                
                for feature in rec.features:
                    feature_type = feature.type.lower()
                    
                    if feature_type == 'gene':
                        gene_count += 1
                    elif feature_type == 'cds':
                        cds_count += 1
                    elif feature_type in ['rrna', 'rrna_gene']:
                        rrna_count += 1
                    elif feature_type in ['trna', 'trna_gene']:
                        trna_count += 1
                    elif feature_type in ['ncrna', 'ncrna_gene', 'non_coding_rna', 'misc_rna']:
                        other_rna_count += 1
                
                data['Total Genes'] = gene_count
                data['CDS Count'] = cds_count
                data['rRNA Count'] = rrna_count
                data['tRNA Count'] = trna_count
                data['Other ncRNA Count'] = other_rna_count
                
                self.data[rec.id] = data
        else:
            return self._analyze_gff_simple()
        
        return self.data
    
    def _analyze_gff_simple(self):
        """Simplified GFF parsing (BCBio independent)"""
        data = {}
        gene_counts = defaultdict(int)
        rrna_count = 0
        trna_count = 0
        other_rna_count = 0
        cds_count = 0
        gene_count = 0
        seq_ids = set()
        seq_length = 0
        
        with open(self.input_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Parse sequence line
                if line.startswith('##sequence-region'):
                    parts = line.split()
                    if len(parts) >= 4:
                        seq_ids.add(parts[1])
                        try:
                            seq_length = int(parts[3])
                        except ValueError:
                            pass
                
                # Parse feature line
                if line.startswith('##'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                
                seq_id = parts[0]
                feature_type = parts[2].lower()
                
                seq_ids.add(seq_id)
                
                # Count genes
                if feature_type == 'gene':
                    gene_count += 1
                elif feature_type == 'cds':
                    cds_count += 1
                elif feature_type in ['rrna', 'rrna_gene']:
                    rrna_count += 1
                elif feature_type in ['trna', 'trna_gene']:
                    trna_count += 1
                elif feature_type in ['ncrna', 'ncrna_gene', 'non_coding_rna']:
                    other_rna_count += 1
        
        # Use filename as identifier
        seq_id = os.path.splitext(os.path.basename(self.input_file))[0]
        
        data['Sequence ID'] = seq_id
        data['Genome Length'] = seq_length if seq_length > 0 else 'N/A'
        data['Species Name'] = 'N/A'
        data['Description'] = f'GFF file: {os.path.basename(self.input_file)}'
        data['Total Genes'] = gene_count
        data['CDS Count'] = cds_count
        data['rRNA Count'] = rrna_count
        data['tRNA Count'] = trna_count
        data['Other ncRNA Count'] = other_rna_count
        
        self.data[seq_id] = data
        return self.data
    
    def _analyze_fasta(self):
        """Analyze FASTA format file"""
        print(f"Analyzing FASTA file: {self.input_file}")
        
        for record in SeqIO.parse(self.input_file, "fasta"):
            data = {}
            
            data['Sequence ID'] = record.id
            data['Sequence Name'] = record.name
            data['Description'] = record.description
            data['Genome Length'] = len(record.seq)
            data['Molecule Type'] = 'DNA/RNA'
            
            # Calculate GC content
            gc_count = record.seq.upper().count('G') + record.seq.upper().count('C')
            data['GC Content (%)'] = f"{(gc_count / len(record.seq) * 100):.2f}" if len(record.seq) > 0 else '0'
            
            # FASTA files typically have no gene annotations
            data['Total Genes'] = 'N/A'
            data['CDS Count'] = 'N/A'
            data['rRNA Count'] = 'N/A'
            data['tRNA Count'] = 'N/A'
            data['Other ncRNA Count'] = 'N/A'
            data['Species Name'] = 'N/A'
            
            self.data[record.id] = data
        
        return self.data
    
    def _extract_organism(self, record):
        """Extract species name from GenBank record"""
        if 'organism' in record.annotations:
            return record.annotations['organism']
        if 'source' in record.annotations:
            return record.annotations['source']
        return 'N/A'
    
    def _extract_taxonomy(self, record):
        """Extract taxonomy information"""
        if 'taxonomy' in record.annotations:
            return '; '.join(record.annotations['taxonomy'])
        return 'N/A'
    
    def _extract_organism_from_gff(self, rec):
        """Extract species name from GFF record"""
        return rec.annotations.get('organism', 'N/A')


def debug_genbank_features(input_file):
    """Debug function: Display all features in GenBank file"""
    print("\n" + "=" * 60)
    print("Debug Info - Display All Features:")
    print("=" * 60)
    
    for record in SeqIO.parse(input_file, "genbank"):
        print(f"\nRecord ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Number of Features: {len(record.features)}")
        
        for i, feature in enumerate(record.features):
            print(f"\nFeature {i + 1}:")
            print(f"  Type: {feature.type}")
            print(f"  Location: {feature.location}")
            if feature.location:
                print(f"  Length: {len(feature.location)}")
            
            print(f"  Qualifiers:")
            for key, value in feature.qualifiers.items():
                val_str = str(value)
                if len(val_str) > 100:
                    val_str = val_str[:100] + "..."
                print(f"    {key}: {val_str}")


def export_to_excel(data, output_file=None):
    """Export analysis results to Excel file"""
    import pandas as pd
    from datetime import datetime
    
    if not data:
        print("Error: No data to export")
        return None
    
    # Get all possible columns
    all_columns = set()
    for record_data in data.values():
        all_columns.update(record_data.keys())
    
    # Sort columns in reasonable order
    priority_columns = [
        'Source File', 'Sequence ID', 'Sequence Name', 'Species Name', 'Description',
        'Genome Length', 'Topology', 'Molecule Type', 'Taxonomy', 'GC Content (%)',
        'Total Genes', 'CDS Count', 'rRNA Count', 'tRNA Count', 'Other ncRNA Count',
        'Annotation Date', 'Data Version', 'Annotation Source'
    ]
    
    other_columns = sorted(all_columns - set(priority_columns))
    columns = [col for col in priority_columns if col in all_columns] + other_columns
    
    # Create DataFrame
    df = pd.DataFrame.from_dict(data, orient='index', columns=columns)
    df.reset_index(drop=True, inplace=True)
    
    # If output filename not specified, use first sequence ID
    if output_file is None:
        first_id = list(data.keys())[0]
        output_file = f"{first_id}.xlsx"
    elif os.path.isdir(output_file):
        output_file = os.path.join(output_file, 'annotation_summary.xlsx')
    
    # Ensure correct file extension
    if not output_file.lower().endswith('.xlsx'):
        output_file += '.xlsx'
    
    # Export to Excel (if the file is occupied by Excel, auto-save with timestamp suffix)
    try:
        df.to_excel(output_file, index=False, engine='openpyxl')
    except PermissionError:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        base, ext = os.path.splitext(output_file)
        output_file = f"{base}_{timestamp}{ext}"
        df.to_excel(output_file, index=False, engine='openpyxl')
    
    print(f"\nResults successfully exported to: {output_file}")
    print(f"\nNumber of records extracted: {len(data)}")
    print(f"Number of fields: {len(columns)}")
    
    return output_file


def select_paths_console():
    """Console fallback when tkinter is unavailable: choose input type and paths."""
    print('\nSelect input type:')
    print('  1 = a folder (analyze all annotation files inside it)')
    print('  2 = a single annotation file')
    print('  3 = multiple annotation files (separated by commas)')
    try:
        choice = input('Your choice (1/2/3, press Enter to cancel): ').strip()
    except Exception:
        return []
    if choice == '1':
        folder = input('Folder path (press Enter to cancel): ').strip().strip('"')
        return [folder] if folder else []
    if choice == '2':
        file_path = input('File path (press Enter to cancel): ').strip().strip('"')
        return [file_path] if file_path else []
    if choice == '3':
        raw = input('File paths, separated by commas (press Enter to cancel): ').strip()
        return [p.strip().strip('"') for p in raw.split(',') if p.strip()]
    return []


def select_inputs_gui():
    """Show a chooser dialog for a folder, a single file, or multiple files."""
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:
        return select_paths_console()

    selected = []

    def choose_folder():
        folder = filedialog.askdirectory(title='Select a folder containing annotation files')
        if folder:
            selected.append(folder)
            root.destroy()

    def choose_one_file():
        file_path = filedialog.askopenfilename(title='Select an annotation file',
                                               filetypes=ANNOTATION_FILE_TYPES)
        if file_path:
            selected.append(file_path)
            root.destroy()

    def choose_many_files():
        files = filedialog.askopenfilenames(title='Select annotation files',
                                            filetypes=ANNOTATION_FILE_TYPES)
        if files:
            selected.extend(files)
            root.destroy()

    root = tk.Tk()
    root.title('Genome Annotation Analyzer - Select Input')
    root.attributes('-topmost', True)
    tk.Label(root, text='Choose an input type:').pack(pady=(14, 6))
    tk.Button(root, text='Select a Folder (analyze all files in it)',
              width=48, command=choose_folder).pack(padx=20, pady=3)
    tk.Button(root, text='Select a Single File',
              width=48, command=choose_one_file).pack(padx=20, pady=3)
    tk.Button(root, text='Select Multiple Files (merge into one table)',
              width=48, command=choose_many_files).pack(padx=20, pady=3)
    tk.Button(root, text='Cancel',
              width=48, command=root.destroy).pack(padx=20, pady=(3, 14))
    root.mainloop()
    return selected


def collect_annotation_files(folder_path):
    """Recursively collect all supported annotation files under the folder."""
    files = []
    for root, dirs, names in os.walk(folder_path):
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        for name in sorted(names):
            if is_annotation_file(name):
                files.append(os.path.join(root, name))
    return files


def analyze_file_list(file_paths, source_labels=None, output_file=None):
    """Analyze multiple annotation files and merge them into one summary Excel file."""
    files = []
    seen = set()
    for file_path in file_paths:
        abs_path = os.path.abspath(file_path)
        if abs_path not in seen:
            seen.add(abs_path)
            files.append(abs_path)

    if not files:
        print('\nNo annotation files to analyze.')
        return None

    labels = source_labels or {}
    print(f'\nFound {len(files)} annotation file(s). Starting batch analysis...')
    combined = {}
    failed = []

    for i, file_path in enumerate(files, 1):
        label = labels.get(file_path, os.path.basename(file_path))
        print(f'\n[{i}/{len(files)}] Analyzing: {label}')
        try:
            analyzer = GenomeAnnotationAnalyzer(file_path)
            data = analyzer.analyze()
        except Exception as e:
            print(f'  Failed: {e}')
            failed.append((label, str(e)))
            continue

        for seq_id, record_data in data.items():
            if 'Source File' not in record_data:
                record_data['Source File'] = label
            combined[f'{i:04d}_{seq_id}'] = record_data

    if not combined:
        print('\nNo annotation file could be parsed successfully. No summary table generated.')
        return None

    # Default output: put the summary table next to the first input file
    if output_file is None:
        output_file = os.path.join(os.path.dirname(files[0]), 'annotation_summary.xlsx')

    out_path = export_to_excel(combined, output_file)

    print('\n' + '=' * 60)
    print('Batch Analysis Summary:')
    print('=' * 60)
    print(f'  Files processed   : {len(files)}')
    print(f'  Files succeeded   : {len(files) - len(failed)}')
    print(f'  Files failed      : {len(failed)}')
    print(f'  Records extracted : {len(combined)}')
    if failed:
        print('\nFailed file details:')
        for label, err in failed:
            print(f'  - {label} ({err})')
    return out_path


def analyze_folder(folder_path, output_file=None):
    """Batch analyze all annotation files in the folder and generate one summary table."""
    files = collect_annotation_files(folder_path)
    if not files:
        print('\nNo supported annotation files found in the folder:')
        print(f'  Folder: {folder_path}')
        print('  Supported formats: ' + ', '.join(sorted(SUPPORTED_EXTENSIONS)))
        return None

    labels = {os.path.abspath(f): os.path.relpath(f, folder_path) for f in files}

    # Default output: put the summary table into the analyzed folder
    if output_file is None:
        base_name = os.path.basename(os.path.normpath(folder_path)) or 'annotation'
        output_file = os.path.join(folder_path, f'{base_name}_annotation_summary.xlsx')

    return analyze_file_list(files, labels, output_file)


def analyze_single_file(input_file, output_file=None, debug_mode=False):
    """Analyze a single annotation file and export its result."""
    if not os.path.exists(input_file):
        print(f"Error: File not found: {input_file}")
        return None

    # If debug mode, display all features first
    if debug_mode and input_file.endswith(('.gb', '.gbk')):
        debug_genbank_features(input_file)

    # Analyze file
    analyzer = GenomeAnnotationAnalyzer(input_file)
    data = analyzer.analyze()

    for record_data in data.values():
        if 'Source File' not in record_data:
            record_data['Source File'] = os.path.basename(input_file)

    # Display analysis results summary
    print("\n" + "=" * 60)
    print("Analysis Results Summary:")
    print("=" * 60)
    for seq_id, record_data in data.items():
        print(f"\nSequence: {seq_id}")
        print("-" * 60)
        for key, value in sorted(record_data.items()):
            print(f"  {key}: {value}")

    # Export to Excel
    return export_to_excel(data, output_file)


def process_paths(input_paths, output_file=None):
    """Analyze mixed inputs (folders and/or files) selected by the user."""
    files = []
    labels = {}

    for path in input_paths:
        abs_path = os.path.abspath(path)
        if os.path.isdir(abs_path):
            for file_path in collect_annotation_files(abs_path):
                abs_file = os.path.abspath(file_path)
                if abs_file not in labels:
                    labels[abs_file] = os.path.relpath(abs_file, abs_path)
                files.append(abs_file)
        elif os.path.isfile(abs_path) and is_annotation_file(abs_path):
            if abs_path not in labels:
                labels[abs_path] = os.path.basename(abs_path)
            files.append(abs_path)

    # Remove duplicates while preserving order
    unique_files = []
    seen = set()
    for file_path in files:
        if file_path not in seen:
            seen.add(file_path)
            unique_files.append(file_path)

    if not unique_files:
        print('\nNo supported annotation files were found in the selected input.')
        print('  Supported formats: ' + ', '.join(sorted(SUPPORTED_EXTENSIONS)))
        return None

    if len(unique_files) == 1:
        return analyze_single_file(unique_files[0], output_file)

    return analyze_file_list(unique_files, labels, output_file)


def hold_console():
    """Keep the console window open (double-click / drag & drop launch)."""
    try:
        input('\nDone. Press Enter to exit...')
    except Exception:
        pass


def main():
    """Main function"""
    print("=" * 60)
    print("Genome Annotation File Analysis Tool")
    print("=" * 60)
    print()

    args = sys.argv[1:]
    debug_mode = '--debug' in args
    positional = [a for a in args if not a.startswith('--')]

    # Double-click launch (no argument): show the input chooser dialog
    if not positional:
        print('No file/folder argument detected.')
        print('A dialog will open for you to choose a folder or file(s).')
        print('(You can also drag & drop files or folders onto this script.)\n')
        selected = select_inputs_gui()
        if not selected:
            print('No input selected. Program exiting.')
            hold_console()
            return
        try:
            if len(selected) == 1 and os.path.isdir(selected[0]):
                analyze_folder(selected[0])
            else:
                process_paths(selected)
        except Exception as e:
            print(f"\nError: {str(e)}")
            import traceback
            traceback.print_exc()
        hold_console()
        return

    try:
        # Single command-line argument: folder batch or single file
        if len(positional) == 1:
            target = positional[0]
            if os.path.isdir(target):
                analyze_folder(target)
            else:
                analyze_single_file(target, None, debug_mode)

        # Multiple command-line arguments: distinguish inputs from the output file
        else:
            inputs = []
            output_file = None
            for arg in positional:
                if os.path.isdir(arg) or is_annotation_file(arg):
                    inputs.append(arg)
                elif output_file is None:
                    output_file = arg

            if not inputs:
                # Backward compatible: old usage <input_file> [output_file]
                analyze_single_file(positional[0],
                                    positional[1] if len(positional) > 1 else None,
                                    debug_mode)
            elif len(inputs) == 1 and os.path.isfile(inputs[0]):
                analyze_single_file(inputs[0], output_file, debug_mode)
            else:
                process_paths(inputs, output_file)

    except Exception as e:
        print(f"\nError: {str(e)}")
        import traceback
        traceback.print_exc()
        hold_console()
        return

    hold_console()


if __name__ == '__main__':
    main()
