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
    
    if not data:
        print("Error: No data to export")
        return None
    
    # Get all possible columns
    all_columns = set()
    for record_data in data.values():
        all_columns.update(record_data.keys())
    
    # Sort columns in reasonable order
    priority_columns = [
        'Sequence ID', 'Sequence Name', 'Species Name', 'Description', 'Genome Length',
        'Topology', 'Molecule Type', 'Taxonomy', 'GC Content (%)',
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
    
    # Ensure correct file extension
    if not output_file.endswith('.xlsx'):
        output_file += '.xlsx'
    
    # Export to Excel
    df.to_excel(output_file, index=False, engine='openpyxl')
    print(f"\nResults successfully exported to: {output_file}")
    print(f"\nNumber of records extracted: {len(data)}")
    print(f"Number of fields: {len(columns)}")
    
    return output_file


def main():
    """Main function"""
    print("=" * 60)
    print("Genome Annotation File Analysis Tool")
    print("=" * 60)
    print()
    
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python genome_annotation_analyzer.py <input_file> [output_file] [--debug]")
        print()
        print("Examples:")
        print("  python genome_annotation_analyzer.py annotation.gb")
        print("  python genome_annotation_analyzer.py annotation.gff result.xlsx")
        print("  python genome_annotation_analyzer.py annotation.gb --debug")
        print()
        print("Supported file formats:")
        print("  - GenBank (.gb, .gbk)")
        print("  - GFF/GFF3 (.gff, .gff3)")
        print("  - FASTA (.fasta, .fa)")
        print()
        print("Extracted information includes:")
        print("  - Sequence ID, species name, description, sequence length")
        print("  - Gene types and counts")
        print("  - rRNA count, tRNA count, total gene count")
        print("  - GC content, topology, taxonomy information, etc.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = None
    debug_mode = False
    
    # Parse arguments
    for arg in sys.argv[2:]:
        if arg == '--debug':
            debug_mode = True
        elif not arg.startswith('--'):
            output_file = arg
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: File not found: {input_file}")
        sys.exit(1)
    
    try:
        # If debug mode, display all features first
        if debug_mode and input_file.endswith(('.gb', '.gbk')):
            debug_genbank_features(input_file)
        
        # Analyze file
        analyzer = GenomeAnnotationAnalyzer(input_file)
        data = analyzer.analyze()
        
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
        export_to_excel(data, output_file)
        
    except Exception as e:
        print(f"\nError: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
