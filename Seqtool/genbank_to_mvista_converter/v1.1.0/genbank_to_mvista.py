#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GenBank/GFF3 to mVISTA Annotation File Converter
Batch convert GenBank or GFF3 format annotation files to mVISTA-compatible annotation files

mVISTA Format Description:
- Tab-delimited text file
- Contains fields: Feature, Start, End, Name, Direction, Annotation
- Feature: Feature type (gene, exon, CDS, rRNA, tRNA, mRNA, etc.)
- Start: Start position (1-based)
- End: End position (1-based)
- Name: Gene/feature name
- Direction: Direction (+/- or forward/reverse)
- Annotation: Annotation information

GFF3 Format Description:
- Column 1: seqid (sequence ID)
- Column 2: source (source)
- Column 3: type (feature type)
- Column 4: start (start position, 1-based)
- Column 5: end (end position, 1-based)
- Column 6: score (score)
- Column 7: strand (strand direction: +, -, .)
- Column 8: phase (phase: 0, 1, 2, .)
- Column 9: attributes (attributes, semicolon-separated key-value pairs)
"""

import os
import sys
from pathlib import Path
from typing import List, Dict
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QListWidget, QTextEdit,
    QFileDialog, QMessageBox, QGroupBox, QProgressBar, QSplitter,
    QRadioButton, QButtonGroup, QCheckBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


class ConversionThread(QThread):
    """Background conversion thread"""
    progress_updated = pyqtSignal(int)
    file_completed = pyqtSignal(str, str)
    error_occurred = pyqtSignal(str, str)
    all_completed = pyqtSignal(int, int)
    log_message = pyqtSignal(str)

    def __init__(self, input_files: List[str], output_dir: str,
                 feature_types: List[str], file_format: str = 'genbank',
                 rna_exon_as_utr: bool = False):
        super().__init__()
        self.input_files = input_files
        self.output_dir = output_dir
        self.feature_types = feature_types
        self.file_format = file_format
        self.rna_exon_as_utr = rna_exon_as_utr
        self.success_count = 0
        self.fail_count = 0

    def run(self):
        """Execute conversion"""
        total_files = len(self.input_files)
        
        for i, input_file in enumerate(self.input_files):
            try:
                self.log_message.emit(f"Processing: {Path(input_file).name}")
                
                if self.file_format == 'gff3':
                    # Read GFF3 file
                    self.convert_gff3_file(input_file, self.rna_exon_as_utr)
                else:
                    # Read GenBank file
                    records = list(SeqIO.parse(input_file, "genbank"))
                    
                    if not records:
                        raise ValueError("No valid GenBank records found")
                    
                    # Convert each record
                    for record in records:
                        output_file = self.convert_record(record, self.output_dir,
                                                         self.feature_types, self.rna_exon_as_utr)
                        self.file_completed.emit(input_file, output_file)
                        self.log_message.emit(f"✓ Success: {output_file}")
                
                self.success_count += 1
                
            except Exception as e:
                self.fail_count += 1
                error_msg = f"Conversion failed: {str(e)}"
                self.error_occurred.emit(input_file, error_msg)
                self.log_message.emit(f"✗ {error_msg}")
            
            # Update progress
            progress = int((i + 1) / total_files * 100)
            self.progress_updated.emit(progress)
        
        self.all_completed.emit(self.success_count, self.fail_count)

    def convert_record(self, record: SeqRecord, output_dir: str,
                      feature_types: List[str], rna_exon_as_utr: bool = False) -> str:
        """
        Convert GenBank record to mVISTA format

        Args:
            record: GenBank record
            output_dir: Output directory
            feature_types: List of feature types to extract
            rna_exon_as_utr: Whether to display RNA exons as UTR

        Returns:
            Output file path
        """
        # Create output filename
        record_id = record.id.replace(" ", "_")
        output_filename = f"{record_id}_mvista.txt"
        output_path = os.path.join(output_dir, output_filename)
        
        # Extract features and group by gene
        features_by_gene = {}  # {gene_name: {'strand': +1/-1, 'gene_range': (start, end), 'exons': [], 'utrs': []}}
        
        for feature in record.features:
            feature_type = feature.type.lower()

            # Process all relevant feature types
            # RNA-related features
            if feature_type in ['gene', 'exon', 'cds', '5utr', '3utr', 'utr', 'rna', 'rrna', 'trna', 'snrna', 'scrna', 'mrna', 'ncrna', 'mirna', 'lnc_rna']:
                # Get location information
                if feature.location is None:
                    continue

                # Get gene name as grouping identifier
                gene_name = self.get_gene_name(feature)

                # Check for compound location (e.g., multiple exons in CDS)
                from Bio.SeqFeature import CompoundLocation

                if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                    # Process compound location (e.g., multiple exon regions in CDS)
                    for part in feature.location.parts:
                        start = int(part.start) + 1
                        end = int(part.end)

                        # Get gene name as grouping identifier
                        gene_name_part = self.get_gene_name(feature)

                        # Process RNA features separately (as independent genes)
                        # Note: mRNA is not in this list because it should be processed as a regular gene
                        if feature_type in ['rna', 'rrna', 'trna', 'snrna', 'scrna', 'ncrna', 'mirna', 'lnc_rna',
                                           'rRNA', 'tRNA', 'snRNA', 'scRNA', 'ncRNA', 'miRNA', 'lnc_RNA']:
                            rna_type_lower = feature_type.lower()
                            # Create independent entry for RNA features
                            if gene_name_part not in features_by_gene:
                                features_by_gene[gene_name_part] = {
                                    'gene_start': start,
                                    'gene_end': end,
                                    'gene_strand': 1,  # Default to forward strand
                                    'exons': [],
                                    'utrs': [],
                                    'is_rna': True,
                                    'rna_type': rna_type_lower
                                }
                                # Get RNA strand direction
                                try:
                                    features_by_gene[gene_name_part]['gene_strand'] = feature.strand if feature.strand else 1
                                except (AttributeError, TypeError):
                                    features_by_gene[gene_name_part]['gene_strand'] = 1
                            else:
                                # Update RNA range
                                features_by_gene[gene_name_part]['gene_start'] = min(features_by_gene[gene_name_part]['gene_start'], start)
                                features_by_gene[gene_name_part]['gene_end'] = max(features_by_gene[gene_name_part]['gene_end'], end)
                                # Get RNA strand direction
                                try:
                                    features_by_gene[gene_name_part]['gene_strand'] = feature.strand if feature.strand else 1
                                except (AttributeError, TypeError):
                                    features_by_gene[gene_name_part]['gene_strand'] = 1
                            # Add RNA region as exon
                            features_by_gene[gene_name_part]['exons'].append({
                                'start': start,
                                'end': end
                            })
                        else:
                            # Processing of regular gene-related features
                            if gene_name_part not in features_by_gene:
                                features_by_gene[gene_name_part] = {
                                    'gene_start': None,
                                    'gene_end': None,
                                    'gene_strand': 1,
                                    'exons': [],
                                    'utrs': [],
                                    'cds_regions': [],
                                    'is_rna': False,
                                    'rna_type': None
                                }

                            # Record gene start and end positions
                            if feature_type == 'gene':
                                if features_by_gene[gene_name_part]['gene_start'] is None:
                                    features_by_gene[gene_name_part]['gene_start'] = start
                                if features_by_gene[gene_name_part]['gene_end'] is None:
                                    features_by_gene[gene_name_part]['gene_end'] = end

                            # Record mRNA
                            elif feature_type == 'mrna':
                                if features_by_gene[gene_name_part]['gene_start'] is None:
                                    features_by_gene[gene_name_part]['gene_start'] = start
                                if features_by_gene[gene_name_part]['gene_end'] is None:
                                    features_by_gene[gene_name_part]['gene_end'] = end

                            # Record each part of exon and CDS
                            elif feature_type in ['exon', 'cds']:
                                features_by_gene[gene_name_part]['exons'].append({
                                    'start': start,
                                    'end': end
                                })
                                if feature_type == 'cds':
                                    features_by_gene[gene_name_part]['cds_regions'].append({
                                        'start': start,
                                        'end': end
                                    })

                            # Record UTR
                            elif feature_type in ['5utr', '3utr', 'utr']:
                                features_by_gene[gene_name_part]['utrs'].append({
                                    'start': start,
                                    'end': end
                                })
                    continue  # Continue to next feature after processing compound location

                # Convert to 1-based coordinates (simple location)
                start = int(feature.location.start) + 1
                end = int(feature.location.end)

                # Process RNA features separately (as independent genes)
                # Note: mRNA is not in this list because it should be processed as a regular gene
                if feature_type in ['rna', 'rrna', 'trna', 'snrna', 'scrna', 'ncrna', 'mirna', 'lnc_rna',
                                   'rRNA', 'tRNA', 'snRNA', 'scRNA', 'ncRNA', 'miRNA', 'lnc_RNA']:
                    # Standardize RNA type to lowercase
                    rna_type_lower = feature_type.lower()
                    # Create independent entry for RNA features
                    if gene_name not in features_by_gene:
                        features_by_gene[gene_name] = {
                            'gene_start': start,
                            'gene_end': end,
                            'gene_strand': 1,  # Default to forward strand
                            'exons': [],
                            'utrs': [],
                            'is_rna': True,
                            'rna_type': rna_type_lower
                        }
                        # Add RNA itself as an exon
                        features_by_gene[gene_name]['exons'].append({
                            'start': start,
                            'end': end
                        })
                    else:
                        # Update RNA range
                        features_by_gene[gene_name]['gene_start'] = min(features_by_gene[gene_name]['gene_start'], start)
                        features_by_gene[gene_name]['gene_end'] = max(features_by_gene[gene_name]['gene_end'], end)
                        # Add RNA region as an exon
                        features_by_gene[gene_name]['exons'].append({
                            'start': start,
                            'end': end
                        })

                    # Get RNA strand direction
                    try:
                        features_by_gene[gene_name]['gene_strand'] = feature.strand if feature.strand else 1
                    except (AttributeError, TypeError):
                        features_by_gene[gene_name]['gene_strand'] = 1

                    continue

                # Processing of regular gene-related features (including gene, mRNA, exon, CDS, UTR, etc.)
                # If no entry for this gene exists, create one
                if gene_name not in features_by_gene:
                    features_by_gene[gene_name] = {
                        'gene_start': None,
                        'gene_end': None,
                        'gene_strand': 1,  # Default to forward strand
                        'exons': [],
                        'utrs': [],
                        'cds_regions': [],  # Add CDS region list
                        'is_rna': False,
                        'rna_type': None
                    }

                # Record gene start and end positions (from gene type features)
                if feature_type == 'gene':
                    features_by_gene[gene_name]['gene_start'] = start
                    features_by_gene[gene_name]['gene_end'] = end
                    # Get gene strand direction
                    try:
                        features_by_gene[gene_name]['gene_strand'] = feature.strand if feature.strand else 1
                    except (AttributeError, TypeError):
                        features_by_gene[gene_name]['gene_strand'] = 1

                # Record mRNA (in GenBank, mRNA also contains gene position information)
                elif feature_type == 'mrna':
                    # If gene range is not set, use mRNA range
                    if features_by_gene[gene_name]['gene_start'] is None:
                        features_by_gene[gene_name]['gene_start'] = start
                    if features_by_gene[gene_name]['gene_end'] is None:
                        features_by_gene[gene_name]['gene_end'] = end
                    # Get mRNA strand direction
                    try:
                        features_by_gene[gene_name]['gene_strand'] = feature.strand if feature.strand else 1
                    except (AttributeError, TypeError):
                        features_by_gene[gene_name]['gene_strand'] = 1

                # Record exon
                elif feature_type == 'exon':
                    features_by_gene[gene_name]['exons'].append({
                        'start': start,
                        'end': end
                    })

                # Record CDS (in GenBank, CDS usually represents exon regions)
                elif feature_type == 'cds':
                    features_by_gene[gene_name]['cds_regions'].append({
                        'start': start,
                        'end': end
                    })

                # Record UTR (including 5'UTR and 3'UTR)
                elif feature_type in ['5utr', '3utr', 'utr']:
                    features_by_gene[gene_name]['utrs'].append({
                        'start': start,
                        'end': end
                    })
        
        # If no exon but has cds, use cds as exon
        for gene_name in features_by_gene:
            if not features_by_gene[gene_name]['exons'] and features_by_gene[gene_name]['cds_regions']:
                features_by_gene[gene_name]['exons'] = features_by_gene[gene_name]['cds_regions']
        
        # Write mVISTA format file
        with open(output_path, 'w', encoding='utf-8') as f:
            # Sort by gene start position
            sorted_genes = sorted(
                [(name, data) for name, data in features_by_gene.items()],
                key=lambda x: (x[1]['gene_start'] or float('inf'))
            )
            
            # Write each gene and its features
            for gene_name, gene_data in sorted_genes:
                # Get gene start and end positions
                gene_start = gene_data['gene_start']
                gene_end = gene_data['gene_end']
                gene_strand = gene_data['gene_strand']
                is_rna = gene_data.get('is_rna', False)
                rna_type = gene_data.get('rna_type', None)

                # If no gene type feature found, use min and max of all features
                if gene_start is None or gene_end is None:
                    all_coords = []
                    all_coords.extend([(e['start'], e['end']) for e in gene_data['exons']])
                    all_coords.extend([(c['start'], c['end']) for c in gene_data.get('cds_regions', [])])
                    all_coords.extend([(u['start'], u['end']) for u in gene_data['utrs']])
                    if all_coords:
                        gene_start = min(all_coords)
                        gene_end = max([coord[1] for coord in all_coords])

                # Write gene line (use > for forward strand, < for reverse strand)
                if gene_start and gene_end:
                    strand_symbol = '>' if gene_strand == 1 else '<'

                    # Write RNA features directly
                    if is_rna and rna_type:
                        # Standardize RNA type name
                        rna_display = rna_type.replace('_', '').upper()
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name} ({rna_display})\n")

                        # Write RNA exons (sorted by position)
                        # Decide to display as exon or utr based on option
                        feature_label = "utr" if rna_exon_as_utr else "exon"
                        sorted_exons = sorted(gene_data['exons'], key=lambda x: x['start'])
                        for exon in sorted_exons:
                            f.write(f"{exon['start']} {exon['end']} {feature_label}\n")

                    # Regular gene features
                    else:
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name}\n")

                        # Write UTR (sorted by position)
                        sorted_utrs = sorted(gene_data['utrs'], key=lambda x: x['start'])
                        for utr in sorted_utrs:
                            f.write(f"{utr['start']} {utr['end']} utr\n")

                        # Write exons (sorted by position)
                        sorted_exons = sorted(gene_data['exons'], key=lambda x: x['start'])
                        for exon in sorted_exons:
                            f.write(f"{exon['start']} {exon['end']} exon\n")

                # Add blank line between genes
                f.write("\n")
        
        return output_path

    def get_gene_name(self, feature: SeqFeature) -> str:
        """Get gene name (for grouping)"""
        # Try to get name from different fields
        if 'gene' in feature.qualifiers:
            return feature.qualifiers['gene'][0]
        elif 'locus_tag' in feature.qualifiers:
            return feature.qualifiers['locus_tag'][0]
        elif 'protein_id' in feature.qualifiers:
            return feature.qualifiers['protein_id'][0]
        elif 'label' in feature.qualifiers:
            return feature.qualifiers['label'][0]
        else:
            # If no name found, use feature ID or generate one
            if hasattr(feature, 'id') and feature.id:
                return str(feature.id)
            else:
                # Use position information as identifier
                if feature.location:
                    pos = int(feature.location.start) + 1
                    return f"feature_{pos}"
                return "unknown"

    def get_feature_name(self, feature: SeqFeature) -> str:
        """Get feature name"""
        # Try to get name from different fields
        if 'gene' in feature.qualifiers:
            return feature.qualifiers['gene'][0]
        elif 'locus_tag' in feature.qualifiers:
            return feature.qualifiers['locus_tag'][0]
        elif 'protein_id' in feature.qualifiers:
            return feature.qualifiers['protein_id'][0]
        elif 'label' in feature.qualifiers:
            return feature.qualifiers['label'][0]
        else:
            return "."

    def get_feature_annotation(self, feature: SeqFeature) -> str:
        """Get feature annotation"""
        # Try to get annotation from product field
        if 'product' in feature.qualifiers:
            return feature.qualifiers['product'][0]
        elif 'note' in feature.qualifiers:
            return feature.qualifiers['note'][0]
        elif 'function' in feature.qualifiers:
            return feature.qualifiers['function'][0]
        else:
            return "."

    def convert_gff3_file(self, input_file: str, rna_exon_as_utr: bool = False):
        """
        Convert GFF3 file to mVISTA format

        Args:
            input_file: GFF3 file path
            rna_exon_as_utr: Whether to display RNA exons as UTR
        """
        # Create output filename
        input_path = Path(input_file)
        output_filename = f"{input_path.stem}_mvista.txt"
        output_path = os.path.join(self.output_dir, output_filename)

        # Read GFF3 file
        features_by_gene = {}  # {gene_name: {'strand': '+/-', 'gene_range': (start, end), 'exons': [], 'utrs': []}}
        mrna_to_gene = {}  # Store mRNA to gene mapping

        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()

                # Skip comment lines and empty lines
                if line.startswith('#') or not line:
                    continue

                # Parse GFF3 line
                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                seqid = parts[0]
                source = parts[1]
                feature_type = parts[2].lower()
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes_str = parts[8]

                # Process all relevant feature types (including RNA)
                if feature_type not in ['gene', 'mrna', 'rna', 'rrna', 'trna', 'snrna', 'scrna',
                                       'ncrna', 'mirna', 'lnc_rna', 'exon', 'cds',
                                       '5utr', '3utr', 'utr', 'five_prime_utr', 'three_prime_utr']:
                    continue

                # Parse attributes
                attributes = self.parse_gff3_attributes(attributes_str)

                # Process gene features - establish mRNA to gene mapping
                if feature_type == 'gene':
                    gene_id = attributes.get('ID', [f'gene_{start}'])[0]
                    # Temporarily use ID, will be overwritten by mRNA's Name later
                    if gene_id not in features_by_gene:
                        features_by_gene[gene_id] = {
                            'gene_start': start,
                            'gene_end': end,
                            'gene_strand': strand if strand in ['+', '-'] else '+',
                            'exons': [],
                            'utrs': [],
                            'cds_regions': [],
                            'gene_name': attributes.get('Name', [gene_id])[0],
                            'is_rna': False,
                            'rna_type': None
                        }

                # Process RNA-related features (as independent genes)
                elif feature_type in ['rna', 'rrna', 'trna', 'snrna', 'scrna', 'ncrna', 'mirna', 'lnc_rna']:
                    rna_id = attributes.get('ID', [f'{feature_type}_{start}'])[0]
                    rna_name = attributes.get('Name', [rna_id])[0]

                    if rna_id not in features_by_gene:
                        features_by_gene[rna_id] = {
                            'gene_start': start,
                            'gene_end': end,
                            'gene_strand': strand if strand in ['+', '-'] else '+',
                            'exons': [],
                            'utrs': [],
                            'cds_regions': [],
                            'gene_name': rna_name,
                            'is_rna': True,
                            'rna_type': feature_type
                        }
                        # Add RNA itself as an exon
                        features_by_gene[rna_id]['exons'].append({
                            'start': start,
                            'end': end
                        })
                    else:
                        # Update RNA range
                        features_by_gene[rna_id]['gene_start'] = min(features_by_gene[rna_id]['gene_start'], start)
                        features_by_gene[rna_id]['gene_end'] = max(features_by_gene[rna_id]['gene_end'], end)
                        # Add RNA region as an exon
                        features_by_gene[rna_id]['exons'].append({
                            'start': start,
                            'end': end
                        })

                    # Establish RNA to self mapping (for associating child features)
                    mrna_to_gene[rna_id] = rna_id

                # Process mRNA features - extract gene name and establish mapping
                elif feature_type == 'mrna':
                    mrna_id = attributes.get('ID', [f'mrna_{start}'])[0]
                    parent_id = attributes.get('Parent', [''])[0]
                    gene_name = attributes.get('Name', [parent_id])[0]

                    # If mRNA has parent, update that gene's name
                    if parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['gene_name'] = gene_name
                        features_by_gene[parent_id]['gene_start'] = start
                        features_by_gene[parent_id]['gene_end'] = end
                        features_by_gene[parent_id]['gene_strand'] = strand if strand in ['+', '-'] else '+'
                        # Establish mRNA to gene mapping
                        mrna_to_gene[mrna_id] = parent_id
                    else:
                        # If mRNA has no parent, process as independent gene
                        if mrna_id not in features_by_gene:
                            features_by_gene[mrna_id] = {
                                'gene_start': start,
                                'gene_end': end,
                                'gene_strand': strand if strand in ['+', '-'] else '+',
                                'exons': [],
                                'utrs': [],
                                'cds_regions': [],
                                'gene_name': gene_name,
                                'is_rna': True,
                                'rna_type': 'mrna'
                            }
                        else:
                            # Update mRNA range
                            features_by_gene[mrna_id]['gene_start'] = min(features_by_gene[mrna_id]['gene_start'], start)
                            features_by_gene[mrna_id]['gene_end'] = max(features_by_gene[mrna_id]['gene_end'], end)
                        # Establish mRNA to self mapping
                        mrna_to_gene[mrna_id] = mrna_id

                # Process exon
                elif feature_type == 'exon':
                    parent_id = attributes.get('Parent', [''])[0]
                    # Find corresponding gene through mRNA
                    gene_id = mrna_to_gene.get(parent_id, '')
                    if gene_id and gene_id in features_by_gene:
                        features_by_gene[gene_id]['exons'].append({
                            'start': start,
                            'end': end
                        })
                    # If no mRNA mapping, try to find gene directly through parent
                    elif parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['exons'].append({
                            'start': start,
                            'end': end
                        })

                # Process CDS (as supplement to exons)
                elif feature_type == 'cds':
                    parent_id = attributes.get('Parent', [''])[0]
                    # Find corresponding gene through mRNA
                    gene_id = mrna_to_gene.get(parent_id, '')
                    if gene_id and gene_id in features_by_gene:
                        features_by_gene[gene_id]['cds_regions'].append({
                            'start': start,
                            'end': end
                        })
                    # If no mRNA mapping, try to find gene directly through parent
                    elif parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['cds_regions'].append({
                            'start': start,
                            'end': end
                        })

                # Record UTR (including 5'UTR and 3'UTR)
                elif feature_type in ['5utr', '3utr', 'utr', 'five_prime_utr', 'three_prime_utr']:
                    parent_id = attributes.get('Parent', [''])[0]
                    gene_id = mrna_to_gene.get(parent_id, '')
                    if gene_id and gene_id in features_by_gene:
                        features_by_gene[gene_id]['utrs'].append({
                            'start': start,
                            'end': end
                        })
                    elif parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['utrs'].append({
                            'start': start,
                            'end': end
                        })

        # If no exon but has cds, use cds as exon
        for gene_id in features_by_gene:
            if not features_by_gene[gene_id]['exons'] and features_by_gene[gene_id]['cds_regions']:
                features_by_gene[gene_id]['exons'] = features_by_gene[gene_id]['cds_regions']

        # Write mVISTA format file
        with open(output_path, 'w', encoding='utf-8') as f:
            # Sort by gene start position
            sorted_genes = sorted(
                [(gene_id, data) for gene_id, data in features_by_gene.items()],
                key=lambda x: (x[1]['gene_start'] or float('inf'))
            )

            # Write each gene and its features
            for gene_id, gene_data in sorted_genes:
                # Get gene name (prefer mRNA's Name)
                gene_name = gene_data['gene_name']
                gene_start = gene_data['gene_start']
                gene_end = gene_data['gene_end']
                gene_strand = gene_data['gene_strand']
                is_rna = gene_data.get('is_rna', False)
                rna_type = gene_data.get('rna_type', None)

                # If no gene type feature found, use min and max of all features
                if gene_start is None or gene_end is None:
                    all_coords = []
                    all_coords.extend([(e['start'], e['end']) for e in gene_data['exons']])
                    all_coords.extend([(u['start'], u['end']) for u in gene_data['utrs']])
                    all_coords.extend([(c['start'], c['end']) for c in gene_data['cds_regions']])
                    if all_coords:
                        gene_start = min(all_coords)
                        gene_end = max([coord[1] for coord in all_coords])

                # Write gene line (use > for forward strand, < for reverse strand)
                if gene_start and gene_end:
                    strand_symbol = '>' if gene_strand == '+' else '<'

                    # Write RNA features directly
                    if is_rna and rna_type:
                        # Standardize RNA type name
                        rna_display = rna_type.replace('_', '').upper()
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name} ({rna_display})\n")

                        # Write RNA exons (sorted by position)
                        # Decide to display as exon or utr based on option
                        feature_label = "utr" if rna_exon_as_utr else "exon"
                        sorted_exons = sorted(gene_data['exons'], key=lambda x: x['start'])
                        for exon in sorted_exons:
                            f.write(f"{exon['start']} {exon['end']} {feature_label}\n")

                    # Regular gene features (must have exon or cds)
                    elif gene_data['exons'] or gene_data['cds_regions']:
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name}\n")

                        # Write UTR (sorted by position)
                        sorted_utrs = sorted(gene_data['utrs'], key=lambda x: x['start'])
                        for utr in sorted_utrs:
                            f.write(f"{utr['start']} {utr['end']} utr\n")

                        # Write exons (sorted by position)
                        sorted_exons = sorted(gene_data['exons'], key=lambda x: x['start'])
                        for exon in sorted_exons:
                            f.write(f"{exon['start']} {exon['end']} exon\n")

                    # Add blank line between genes
                    f.write("\n")

        self.file_completed.emit(input_file, output_path)

    def parse_gff3_attributes(self, attributes_str: str) -> Dict[str, List[str]]:
        """
        Parse GFF3 attributes field
        
        Args:
            attributes_str: Attributes string (column 9)
            
        Returns:
            Attributes dictionary {key: [value1, value2, ...]}
        """
        attributes = {}
        if not attributes_str or attributes_str == '.':
            return attributes
        
        # Split attribute key-value pairs
        pairs = attributes_str.split(';')
        for pair in pairs:
            if '=' in pair:
                key, value = pair.split('=', 1)
                if key not in attributes:
                    attributes[key] = []
                attributes[key].append(value)
        
        return attributes

    def get_gene_name_from_gff3(self, attributes: Dict[str, List[str]], 
                                 feature_type: str, position: int) -> str:
        """
        Get gene name from GFF3 attributes
        
        Args:
            attributes: Attributes dictionary
            feature_type: Feature type
            position: Feature position
            
        Returns:
            Gene name
        """
        # Try to get name from different fields
        if 'gene' in attributes:
            return attributes['gene'][0]
        elif 'gene_id' in attributes:
            return attributes['gene_id'][0]
        elif 'Name' in attributes:
            return attributes['Name'][0]
        elif 'name' in attributes:
            return attributes['name'][0]
        elif 'locus_tag' in attributes:
            return attributes['locus_tag'][0]
        elif 'ID' in attributes:
            # For genes, use ID directly
            if feature_type == 'gene':
                return attributes['ID'][0]
            # For other features, try to get gene name from Parent
            elif 'Parent' in attributes:
                return attributes['Parent'][0]
            return attributes['ID'][0]
        elif 'Parent' in attributes:
            return attributes['Parent'][0]
        else:
            # Use position information as identifier
            return f"feature_{position}"


class GenBankToMVISTAConverter(QMainWindow):
    """Main window"""

    def __init__(self):
        super().__init__()
        self.input_files = []
        self.file_format = 'genbank'  # Default to GenBank format
        self.init_ui()
        self.conversion_thread = None

    def init_ui(self):
        """Initialize UI"""
        self.setWindowTitle("GenBank/GFF3 to mVISTA Converter")
        self.setGeometry(100, 100, 1000, 700)

        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main layout
        main_layout = QVBoxLayout(central_widget)

        # Title
        title_label = QLabel("GenBank/GFF3 → mVISTA Batch Converter")
        title_label.setFont(QFont("Arial", 16, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)

        # Create splitter
        splitter = QSplitter(Qt.Horizontal)

        # Left panel - file selection and settings
        left_panel = self.create_left_panel()
        splitter.addWidget(left_panel)

        # Right panel - log display
        right_panel = self.create_right_panel()
        splitter.addWidget(right_panel)

        splitter.setSizes([500, 500])
        main_layout.addWidget(splitter)

        # Bottom progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(self.progress_bar)

    def create_left_panel(self) -> QWidget:
        """Create left control panel"""
        panel = QWidget()
        layout = QVBoxLayout(panel)

        # File format selection group
        format_group = QGroupBox("File Format")
        format_layout = QHBoxLayout()

        self.format_button_group = QButtonGroup()

        self.genbank_radio = QRadioButton("GenBank")
        self.genbank_radio.setChecked(True)
        self.genbank_radio.toggled.connect(self.on_format_changed)
        self.format_button_group.addButton(self.genbank_radio)
        format_layout.addWidget(self.genbank_radio)

        self.gff3_radio = QRadioButton("GFF3")
        self.gff3_radio.toggled.connect(self.on_format_changed)
        self.format_button_group.addButton(self.gff3_radio)
        format_layout.addWidget(self.gff3_radio)

        format_group.setLayout(format_layout)
        layout.addWidget(format_group)

        # Input file selection group
        input_group = QGroupBox("Input Files")
        input_layout = QVBoxLayout()

        # File list
        self.file_list = QListWidget()
        input_layout.addWidget(self.file_list)

        # Button layout
        btn_layout = QHBoxLayout()

        add_btn = QPushButton("Add Files")
        add_btn.clicked.connect(self.add_files)
        btn_layout.addWidget(add_btn)

        add_dir_btn = QPushButton("Add Folder")
        add_dir_btn.clicked.connect(self.add_directory)
        btn_layout.addWidget(add_dir_btn)

        clear_btn = QPushButton("Clear List")
        clear_btn.clicked.connect(self.clear_files)
        btn_layout.addWidget(clear_btn)

        input_layout.addLayout(btn_layout)

        # File count
        self.file_count_label = QLabel("Selected 0 files")
        input_layout.addWidget(self.file_count_label)

        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # Output directory selection group
        output_group = QGroupBox("Output Settings")
        output_layout = QVBoxLayout()
        
        dir_layout = QHBoxLayout()
        dir_layout.addWidget(QLabel("Output Directory:"))
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory...")
        dir_layout.addWidget(self.output_dir_edit)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.browse_output_dir)
        dir_layout.addWidget(browse_btn)
        
        output_layout.addLayout(dir_layout)
        output_group.setLayout(output_layout)
        layout.addWidget(output_group)
        
        # Feature type selection group (fixed to gene, exon, and UTR)
        feature_group = QGroupBox("Feature Type Information")
        feature_layout = QVBoxLayout()

        # Add info text
        info_label = QLabel("mVISTA format supports the following feature types:\n"
                            "• gene: Genes (with strand direction marker >/<)\n"
                            "• exon: Exons\n"
                            "• utr: Untranslated regions (5'UTR and 3'UTR)\n"
                            "• RNA: rRNA, tRNA, mRNA, snRNA, ncRNA, miRNA, lnc_RNA, etc.")
        info_label.setWordWrap(True)
        feature_layout.addWidget(info_label)

        feature_group.setLayout(feature_layout)
        layout.addWidget(feature_group)

        # RNA display options group
        rna_option_group = QGroupBox("RNA Display Options")
        rna_option_layout = QVBoxLayout()

        # Add checkbox
        self.rna_exon_as_utr_checkbox = QCheckBox("Display RNA exons as UTR")
        self.rna_exon_as_utr_checkbox.setToolTip("When checked, all RNA feature exons (rRNA, tRNA, etc.) will be displayed as UTR")
        rna_option_layout.addWidget(self.rna_exon_as_utr_checkbox)

        rna_option_group.setLayout(rna_option_layout)
        layout.addWidget(rna_option_group)
        
        # Convert button
        convert_btn = QPushButton("Start Conversion")
        convert_btn.setFont(QFont("Arial", 12, QFont.Bold))
        convert_btn.clicked.connect(self.start_conversion)
        layout.addWidget(convert_btn)
        
        return panel

    def create_right_panel(self) -> QWidget:
        """Create right log panel"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        
        # Log display
        log_label = QLabel("Conversion Log")
        log_label.setFont(QFont("Arial", 12, QFont.Bold))
        layout.addWidget(log_label)
        
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Consolas", 9))
        layout.addWidget(self.log_text)
        
        # Clear log button
        clear_log_btn = QPushButton("Clear Log")
        clear_log_btn.clicked.connect(self.clear_log)
        layout.addWidget(clear_log_btn)
        
        return panel

    def add_files(self):
        """Add files"""
        if self.file_format == 'genbank':
            files, _ = QFileDialog.getOpenFileNames(
                self, "Select GenBank files", "",
                "GenBank files (*.gb *.gbk *.gbff);;All files (*.*)"
            )
        else:
            files, _ = QFileDialog.getOpenFileNames(
                self, "Select GFF3 files", "",
                "GFF3 files (*.gff *.gff3);;All files (*.*)"
            )

        if files:
            for file in files:
                if file not in self.input_files:
                    self.input_files.append(file)
                    self.file_list.addItem(Path(file).name)
            self.update_file_count()

    def add_directory(self):
        """Add folder"""
        if self.file_format == 'genbank':
            directory = QFileDialog.getExistingDirectory(
                self, "Select folder containing GenBank files"
            )
            if directory:
                for file in Path(directory).glob("*.gb*"):
                    file_path = str(file)
                    if file_path not in self.input_files:
                        self.input_files.append(file_path)
                        self.file_list.addItem(file.name)
                self.update_file_count()
        else:
            directory = QFileDialog.getExistingDirectory(
                self, "Select folder containing GFF3 files"
            )
            if directory:
                for file in Path(directory).glob("*.gff*"):
                    file_path = str(file)
                    if file_path not in self.input_files:
                        self.input_files.append(file_path)
                        self.file_list.addItem(file.name)
                self.update_file_count()

    def clear_files(self):
        """Clear file list"""
        self.input_files.clear()
        self.file_list.clear()
        self.update_file_count()

    def update_file_count(self):
        """Update file count"""
        self.file_count_label.setText(f"Selected {len(self.input_files)} files")

    def on_format_changed(self):
        """File format changed"""
        if self.genbank_radio.isChecked():
            self.file_format = 'genbank'
        else:
            self.file_format = 'gff3'
        # Clear file list when format changes
        self.clear_files()

    def browse_output_dir(self):
        """Browse output directory"""
        directory = QFileDialog.getExistingDirectory(
            self, "Select output directory"
        )
        if directory:
            self.output_dir_edit.setText(directory)

    def start_conversion(self):
        """Start conversion"""
        # Validate input
        if not self.input_files:
            QMessageBox.warning(self, "Warning", "Please select input files first!")
            return

        output_dir = self.output_dir_edit.text()
        if not output_dir:
            QMessageBox.warning(self, "Warning", "Please select output directory!")
            return

        if not os.path.exists(output_dir):
            QMessageBox.warning(self, "Warning", "Output directory does not exist!")
            return

        # Clear log
        self.log_text.clear()
        self.log_message("=== Starting Conversion ===")
        self.log_message(f"Input files: {len(self.input_files)}")
        self.log_message(f"Output directory: {output_dir}")
        self.log_message(f"File format: {self.file_format.upper()}")
        self.log_message("Feature types: gene, exon, utr (fixed)")
        self.log_message("")

        # Reset progress bar
        self.progress_bar.setValue(0)

        # Create and start conversion thread
        rna_exon_as_utr = self.rna_exon_as_utr_checkbox.isChecked()
        self.conversion_thread = ConversionThread(
            self.input_files, output_dir, None, self.file_format, rna_exon_as_utr
        )
        # Record RNA display option in log
        if rna_exon_as_utr:
            self.log_message("RNA display option: Exons displayed as UTR")
        self.conversion_thread.progress_updated.connect(self.progress_bar.setValue)
        self.conversion_thread.file_completed.connect(self.on_file_completed)
        self.conversion_thread.error_occurred.connect(self.on_error_occurred)
        self.conversion_thread.all_completed.connect(self.on_all_completed)
        self.conversion_thread.log_message.connect(self.log_message)
        self.conversion_thread.start()

        # Disable buttons
        self.set_buttons_enabled(False)

    def on_file_completed(self, input_file: str, output_file: str):
        """File conversion completed callback"""
        self.log_message(f"Completed: {Path(input_file).name} → {Path(output_file).name}")

    def on_error_occurred(self, input_file: str, error_msg: str):
        """Error occurred callback"""
        self.log_message(f"Error: {Path(input_file).name}\n  {error_msg}")

    def on_all_completed(self, success_count: int, fail_count: int):
        """All files conversion completed callback"""
        self.log_message("")
        self.log_message("=== Conversion Complete ===")
        self.log_message(f"Success: {success_count} files")
        self.log_message(f"Failed: {fail_count} files")
        self.log_message(f"Total: {success_count + fail_count} files")
        
        # Show result dialog
        QMessageBox.information(
            self, "Conversion Complete",
            f"Conversion complete!\n\nSuccess: {success_count}\nFailed: {fail_count}"
        )
        
        # Enable buttons
        self.set_buttons_enabled(True)
        self.progress_bar.setValue(100)

    def log_message(self, message: str):
        """Add log message"""
        self.log_text.append(message)
        # Auto scroll to bottom
        scrollbar = self.log_text.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def clear_log(self):
        """Clear log"""
        self.log_text.clear()

    def set_buttons_enabled(self, enabled: bool):
        """Set button enabled state"""
        self.file_list.setEnabled(enabled)
        self.output_dir_edit.setEnabled(enabled)


def main():
    """Main function"""
    app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')
    
    # Create and show main window
    window = GenBankToMVISTAConverter()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
