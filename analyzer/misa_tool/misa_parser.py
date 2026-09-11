#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import pandas as pd
from collections import defaultdict
import argparse
import sys


class MISAParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.results = []
        self.seq_stats = defaultdict(lambda: {
            'Mononucleotide': 0,
            'Dinucleotide': 0,
            'Trinucleotide': 0,
            'Tetranucleotide': 0,
            'Pentanucleotide': 0,
            'Hexanucleotide': 0,
            'Total': 0
        })
        self.global_stats = defaultdict(int)
        
    def parse_file(self):
        print(f"Parsing file: {self.file_path}")
        
        with open(self.file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        data_start = 1
        
        for i, line in enumerate(lines[data_start:], start=data_start + 1):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            
            if len(parts) < 6:
                continue
            
            try:
                record = self._parse_record(parts)
                if record:
                    self.results.append(record)
                    self._update_stats(record)
            except Exception as e:
                print(f"Warning: Failed to parse line {i}: {e}")
        
        print(f"Successfully parsed {len(self.results)} SSR records")
        return self.results
    
    def _parse_record(self, parts):
        seq_id = parts[0].strip()
        ssr_nr = parts[1].strip()
        ssr_type = parts[2].strip()
        ssr = parts[3].strip()
        size = self._safe_int(parts[4])
        start = self._safe_int(parts[5])
        end = self._safe_int(parts[6]) if len(parts) > 6 else 0
        
        unit_size = self._extract_unit_size_from_ssr(ssr, ssr_type)
        
        record = {
            'Sequence_ID': seq_id,
            'SSR_Number': ssr_nr,
            'SSR_Type': ssr_type,
            'SSR': ssr,
            'SSR_Length': size,
            'Unit_Size': unit_size,
            'Start': start,
            'End': end
        }
        
        return record
    
    def _extract_unit_size_from_ssr(self, ssr, ssr_type):
        if ssr_type.startswith('p'):
            unit_size = int(ssr_type[1])
        else:
            matches = re.findall(r'\([A-Z]+\)\d+', ssr)
            unit_sizes = []
            for match in matches:
                inner = re.search(r'\(([A-Z]+)\)', match)
                if inner:
                    unit_sizes.append(len(inner.group(1)))
            
            if unit_sizes:
                unit_size = min(unit_sizes)
            else:
                unit_size = 0
        
        return unit_size
    
    def _safe_int(self, value):
        try:
            return int(value)
        except:
            return 0
    
    def _update_stats(self, record):
        seq_id = record['Sequence_ID']
        unit_size = record['Unit_Size']
        ssr_type = record['SSR_Type']
        
        self.seq_stats[seq_id]['Total'] += 1
        
        if unit_size == 1:
            self.seq_stats[seq_id]['Mononucleotide'] += 1
        elif unit_size == 2:
            self.seq_stats[seq_id]['Dinucleotide'] += 1
        elif unit_size == 3:
            self.seq_stats[seq_id]['Trinucleotide'] += 1
        elif unit_size == 4:
            self.seq_stats[seq_id]['Tetranucleotide'] += 1
        elif unit_size == 5:
            self.seq_stats[seq_id]['Pentanucleotide'] += 1
        elif unit_size == 6:
            self.seq_stats[seq_id]['Hexanucleotide'] += 1
        
        self.global_stats[ssr_type] += 1
        self.global_stats[f"{unit_size}bp"] = self.global_stats.get(f"{unit_size}bp", 0) + 1
    
    def generate_detailed_table(self):
        df = pd.DataFrame(self.results)
        return df
    
    def generate_sequence_table(self):
        seq_data = []
        
        for seq_id, stats in sorted(self.seq_stats.items()):
            seq_data.append({
                'Sequence_Name': seq_id,
                'Mononucleotide': stats['Mononucleotide'],
                'Dinucleotide': stats['Dinucleotide'],
                'Trinucleotide': stats['Trinucleotide'],
                'Tetranucleotide': stats['Tetranucleotide'],
                'Pentanucleotide': stats['Pentanucleotide'],
                'Hexanucleotide': stats['Hexanucleotide'],
                'Total_SSR': stats['Total']
            })
        
        df_seq = pd.DataFrame(seq_data)
        return df_seq
    
    def generate_summary_table(self):
        summary_data = []
        
        unit_sizes = ['1bp', '2bp', '3bp', '4bp', '5bp', '6bp']
        for unit_size in unit_sizes:
            count = self.global_stats.get(unit_size, 0)
            if count > 0:
                summary_data.append({
                    'Category': f"{unit_size} Unit",
                    'Count': count,
                    'Percentage(%)': f"{(count/len(self.results)*100):.2f}"
                })
        
        df_summary = pd.DataFrame(summary_data)
        return df_summary
    
    def export_results(self, output_prefix='misa_result'):
        print("\nExporting results...")
        
        df_detailed = self.generate_detailed_table()
        detailed_file = f"{output_prefix}_detailed.xlsx"
        df_detailed.to_excel(detailed_file, index=False, engine='openpyxl')
        print(f"Detailed table exported: {detailed_file}")
        
        df_seq = self.generate_sequence_table()
        seq_file = f"{output_prefix}_by_sequence.xlsx"
        df_seq.to_excel(seq_file, index=False, engine='openpyxl')
        print(f"Sequence statistics table exported: {seq_file}")
        
        df_summary = self.generate_summary_table()
        summary_file = f"{output_prefix}_summary.xlsx"
        df_summary.to_excel(summary_file, index=False, engine='openpyxl')
        print(f"Summary table exported: {summary_file}")
        
        csv_file = f"{output_prefix}_detailed.csv"
        df_detailed.to_csv(csv_file, index=False, encoding='utf-8-sig')
        print(f"Detailed CSV exported: {csv_file}")
        
        csv_seq_file = f"{output_prefix}_by_sequence.csv"
        df_seq.to_csv(csv_seq_file, index=False, encoding='utf-8-sig')
        print(f"Sequence statistics CSV exported: {csv_seq_file}")
        
        return detailed_file, seq_file, summary_file, csv_file, csv_seq_file
    
    def print_statistics(self):
        total_count = len(self.results)
        seq_count = len(self.seq_stats)
        
        print("\n" + "="*90)
        print("MISA SSR Statistics")
        print("="*90)
        print(f"Total SSR Count: {total_count}")
        print(f"Number of Sequences: {seq_count}")
        
        print("\nGlobal Statistics by Unit Size:")
        print("-" * 50)
        for unit_size in ['1bp', '2bp', '3bp', '4bp', '5bp', '6bp']:
            count = self.global_stats.get(unit_size, 0)
            if count > 0:
                percentage = (count / total_count * 100)
                print(f"  {unit_size}: {count} ({percentage:.2f}%)")
        
        print("\nStatistics by Sequence:")
        print("-" * 90)
        print(f"{'Sequence_Name':<16} {'Mono':<8} {'Di':<8} {'Tri':<8} {'Tetra':<8} {'Penta':<8} {'Hexa':<8} {'Total':<8}")
        print("-" * 90)
        for seq_id, stats_dict in sorted(self.seq_stats.items()):
            print(f"{seq_id:<16} {stats_dict['Mononucleotide']:<8} {stats_dict['Dinucleotide']:<8} {stats_dict['Trinucleotide']:<8} "
                  f"{stats_dict['Tetranucleotide']:<8} {stats_dict['Pentanucleotide']:<8} {stats_dict['Hexanucleotide']:<8} {stats_dict['Total']:<8}")
        print("="*90)


def main():
    parser = argparse.ArgumentParser(description='MISA SSR Parser Tool')
    parser.add_argument('-i', '--input', required=True, help='MISA output file path')
    parser.add_argument('-o', '--output', default='misa_result', help='Output file prefix')
    parser.add_argument('--no-export', action='store_true', help='Display statistics only, do not export files')
    
    args = parser.parse_args()
    
    misa_parser = MISAParser(args.input)
    
    misa_parser.parse_file()
    
    if not misa_parser.results:
        print("Error: No valid SSR records found, please check file format")
        sys.exit(1)
    
    misa_parser.print_statistics()
    
    if not args.no_export:
        misa_parser.export_results(args.output)
    else:
        print("\nDetailed table preview:")
        print(misa_parser.generate_detailed_table().head(20).to_string(index=False))
        
        print("\nSequence statistics table:")
        print(misa_parser.generate_sequence_table().to_string(index=False))
        
        print("\nSummary table:")
        print(misa_parser.generate_summary_table().to_string(index=False))


if __name__ == "__main__":
    main()
