# -*- coding: utf-8 -*-
"""
Chloroplast Genome Assembly Comparison Tool
Compares two assembled chloroplast genome sequences against a reference to determine
which assembly is more accurate.
Supports detection of quadripartite structure errors (LSC, IRa, SSC, IRb)
such as swap, inversion, and junction errors.

Dependencies: Python 3.x, Biopython, tkinter (built-in)
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import os
import time
from collections import Counter, defaultdict
from statistics import median

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


# ============================================================
# Sequence Analysis Engine
# ============================================================
class SequenceAnalyzer:
    """Core sequence analysis class, performs all alignment and scoring calculations"""

    def __init__(self, progress_callback=None):
        self.progress_callback = progress_callback

    def _report_progress(self, msg, pct):
        if self.progress_callback:
            self.progress_callback(msg, pct)

    def read_fasta(self, filepath):
        """Read a FASTA file and merge all records into a single sequence"""
        records = []
        # Try UTF-8 first, then fall back to other encodings
        for encoding in ['utf-8', 'utf-8-sig', 'gbk', 'gb18030', 'latin-1']:
            try:
                with open(filepath, 'r', encoding=encoding) as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                break
            except (UnicodeDecodeError, UnicodeError):
                continue
        if not records:
            raise ValueError(f"No valid sequences found in file {filepath}, or the sequence format is incorrect")
        if len(records) == 1:
            return str(records[0].seq).upper()
        combined = "".join(str(r.seq).upper() for r in records)
        return combined

    def get_basic_stats(self, seq, name=""):
        """Calculate basic statistics for a sequence"""
        length = len(seq)
        gc = (seq.count('G') + seq.count('C')) / length * 100 if length > 0 else 0
        a = seq.count('A') / length * 100 if length > 0 else 0
        t = seq.count('T') / length * 100 if length > 0 else 0
        c = seq.count('C') / length * 100 if length > 0 else 0
        g = seq.count('G') / length * 100 if length > 0 else 0
        n_count = seq.count('N')
        return {
            'name': name,
            'length': length,
            'gc_content': round(gc, 2),
            'A_pct': round(a, 2),
            'T_pct': round(t, 2),
            'C_pct': round(c, 2),
            'G_pct': round(g, 2),
            'N_count': n_count,
        }

    def kmer_similarity(self, seq1, seq2, k=31):
        """
        Fast similarity calculation based on k-mers.
        Returns the Jaccard similarity coefficient and shared k-mer ratio.
        This method is insensitive to sequence rearrangements, making it suitable
        for genome-level comparisons.
        """
        def get_kmers(seq, k):
            kmers = set()
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                if 'N' not in kmer:
                    kmers.add(kmer)
            return kmers

        kmers1 = get_kmers(seq1, k)
        kmers2 = get_kmers(seq2, k)

        if not kmers1 or not kmers2:
            return {'jaccard': 0.0, 'shared_ratio': 0.0, 'total1': 0, 'total2': 0}

        intersection = kmers1 & kmers2
        union = kmers1 | kmers2
        jaccard = len(intersection) / len(union) if union else 0.0
        shared_ratio = len(intersection) / len(kmers1) if kmers1 else 0.0

        return {
            'jaccard': round(jaccard, 6),
            'shared_ratio': round(shared_ratio, 6),
            'total_kmers_ref': len(kmers1),
            'total_kmers_asm': len(kmers2),
            'shared_kmers': len(intersection),
        }

    def sliding_window_identity(self, ref_seq, asm_seq, window_size=5000, step=2000):
        """
        Sliding window alignment: divides the reference into windows, uses k-mer matching
        to find the best corresponding region in the assembly, and computes the approximate
        identity for each window.
        Returns a list of window positions and corresponding identity values.
        """
        self._report_progress("Running sliding window analysis...", 40)

        # Build k-mer position index for the assembly
        k = 21
        asm_kmer_positions = {}
        for i in range(len(asm_seq) - k + 1):
            kmer = asm_seq[i:i + k]
            if 'N' not in kmer:
                if kmer not in asm_kmer_positions:
                    asm_kmer_positions[kmer] = []
                asm_kmer_positions[kmer].append(i)

        windows = []
        identities = []

        total_windows = max(1, (len(ref_seq) - window_size) // step + 1)
        processed = 0

        for start in range(0, len(ref_seq) - window_size + 1, step):
            ref_window = ref_seq[start:start + window_size]
            ref_kmers = []
            for i in range(0, len(ref_window) - k + 1, k):
                kmer = ref_window[i:i + k]
                if 'N' not in kmer:
                    ref_kmers.append(kmer)

            if not ref_kmers:
                windows.append(start)
                identities.append(0.0)
                continue

            match_count = 0
            for kmer in ref_kmers:
                if kmer in asm_kmer_positions:
                    match_count += 1

            identity = match_count / len(ref_kmers) * 100 if ref_kmers else 0.0
            windows.append(start)
            identities.append(round(identity, 2))

            processed += 1
            if processed % 10 == 0:
                pct = 40 + int(15 * processed / total_windows)
                self._report_progress(f"Sliding window: {processed}/{total_windows}", pct)

        return windows, identities

    def _global_alignment(self, ref_seq, asm_seq):
        """Perform global pairwise alignment on small sequences"""
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2

        best = next(aligner.align(ref_seq, asm_seq), None)
        if best is None:
            return None

        aligned_ref = str(best).split('\n')[0]
        aligned_asm = str(best).split('\n')[2]

        matches = sum(1 for a, b in zip(aligned_ref, aligned_asm) if a == b)
        mismatches = sum(1 for a, b in zip(aligned_ref, aligned_asm)
                        if a != b and a != '-' and b != '-')
        gaps_ref = aligned_ref.count('-')
        gaps_asm = aligned_asm.count('-')
        total_cols = len(aligned_ref)
        identity = matches / total_cols * 100 if total_cols > 0 else 0

        return {
            'identity': round(identity, 2),
            'matches': matches,
            'mismatches': mismatches,
            'gaps_ref': gaps_ref,
            'gaps_asm': gaps_asm,
            'aligned_length': total_cols,
            'score': best.score,
        }

    def _chunked_alignment(self, ref_seq, asm_seq, chunk_size=5000, overlap=500):
        """
        Chunked alignment: divides the reference into chunks and performs local alignment
        by finding the best matching region in the assembly for each chunk.
        Suitable for large sequences (>50 kb) to avoid memory overflow from full-sequence
        global alignment.
        """
        self._report_progress("Running chunked alignment (large genome mode)...", 20)

        # Build k-mer index for the assembly
        k = 21
        asm_kmer_pos = {}
        for i in range(len(asm_seq) - k + 1):
            kmer = asm_seq[i:i + k]
            if 'N' not in kmer:
                if kmer not in asm_kmer_pos:
                    asm_kmer_pos[kmer] = []
                asm_kmer_pos[kmer].append(i)

        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2

        total_matches = 0
        total_mismatches = 0
        total_gaps_ref = 0
        total_gaps_asm = 0
        total_aligned = 0
        total_score = 0
        chunk_count = 0

        step = chunk_size - overlap
        num_chunks = max(1, (len(ref_seq) - chunk_size) // step + 1)

        for chunk_idx, start in enumerate(range(0, len(ref_seq) - chunk_size + 1, step)):
            if chunk_idx % 5 == 0:
                pct = 20 + int(15 * chunk_idx / num_chunks)
                self._report_progress(f"Chunked alignment: {chunk_idx+1}/{num_chunks}", pct)

            ref_chunk = ref_seq[start:start + chunk_size]

            # Use k-mers to find the best matching region in the assembly
            chunk_kmers = []
            for i in range(0, len(ref_chunk) - k + 1, k):
                kmer = ref_chunk[i:i + k]
                if 'N' not in kmer:
                    chunk_kmers.append(kmer)

            if not chunk_kmers:
                continue

            # Find the most likely corresponding region in the assembly
            positions = []
            for kmer in chunk_kmers:
                if kmer in asm_kmer_pos:
                    positions.extend(asm_kmer_pos[kmer])

            if not positions:
                continue

            # Take the median position and extract the corresponding assembly window
            best_pos = int(median(positions))
            asm_start = max(0, best_pos - chunk_size // 2 - overlap)
            asm_end = min(len(asm_seq), asm_start + chunk_size + overlap * 2)
            asm_chunk = asm_seq[asm_start:asm_end]

            if len(asm_chunk) < chunk_size // 2:
                continue

            try:
                best = next(aligner.align(ref_chunk, asm_chunk), None)
                if best is not None:
                    aligned_ref = str(best).split('\n')[0]
                    aligned_asm = str(best).split('\n')[2]

                    m = sum(1 for a, b in zip(aligned_ref, aligned_asm) if a == b)
                    mm = sum(1 for a, b in zip(aligned_ref, aligned_asm)
                            if a != b and a != '-' and b != '-')
                    gr = aligned_ref.count('-')
                    ga = aligned_asm.count('-')
                    cols = len(aligned_ref)

                    total_matches += m
                    total_mismatches += mm
                    total_gaps_ref += gr
                    total_gaps_asm += ga
                    total_aligned += cols
                    total_score += best.score
                    chunk_count += 1
            except MemoryError:
                continue

        if chunk_count == 0 or total_aligned == 0:
            return None

        identity = total_matches / total_aligned * 100 if total_aligned > 0 else 0

        return {
            'identity': round(identity, 2),
            'matches': total_matches,
            'mismatches': total_mismatches,
            'gaps_ref': total_gaps_ref,
            'gaps_asm': total_gaps_asm,
            'aligned_length': total_aligned,
            'score': total_score,
            'chunks_aligned': chunk_count,
        }

    def biopython_alignment(self, ref_seq, asm_seq):
        """Smart alignment: global alignment for small sequences, chunked alignment for large ones"""
        max_len = max(len(ref_seq), len(asm_seq))
        if max_len <= 50000:
            self._report_progress("Running global alignment...", 20)
            return self._global_alignment(ref_seq, asm_seq)
        else:
            return self._chunked_alignment(ref_seq, asm_seq)

    def compute_coverage(self, ref_seq, asm_seq, k=21):
        """Estimate coverage of the reference by the assembly via k-mer matching"""
        self._report_progress("Computing coverage...", 60)

        asm_kmers = set()
        for i in range(len(asm_seq) - k + 1):
            kmer = asm_seq[i:i + k]
            if 'N' not in kmer:
                asm_kmers.add(kmer)

        ref_kmers = []
        for i in range(0, len(ref_seq) - k + 1, k):
            kmer = ref_seq[i:i + k]
            if 'N' not in kmer:
                ref_kmers.append(kmer)

        if not ref_kmers:
            return 0.0

        covered = sum(1 for kmer in ref_kmers if kmer in asm_kmers)
        coverage = covered / len(ref_kmers) * 100
        return round(coverage, 2)

    def run_full_analysis(self, ref_path, asm1_path, asm2_path):
        """
        Execute the complete analysis workflow:
        1. Read sequences
        2. Basic statistics
        3. K-mer similarity
        4. Sliding window analysis
        5. Biopython global alignment
        6. Coverage calculation
        7. Composite scoring
        """
        results = {}

        # Step 1: Read sequences
        self._report_progress("Reading reference sequence...", 5)
        ref_seq = self.read_fasta(ref_path)
        ref_stats = self.get_basic_stats(ref_seq, "Reference")
        results['reference'] = ref_stats

        self._report_progress("Reading Assembly 1...", 8)
        asm1_seq = self.read_fasta(asm1_path)
        asm1_stats = self.get_basic_stats(asm1_seq, "Assembly 1")
        results['assembly1'] = asm1_stats

        self._report_progress("Reading Assembly 2...", 11)
        asm2_seq = self.read_fasta(asm2_path)
        asm2_stats = self.get_basic_stats(asm2_seq, "Assembly 2")
        results['assembly2'] = asm2_stats

        # Step 2: K-mer similarity analysis
        self._report_progress("Computing k-mer similarity...", 14)
        results['kmer_asm1'] = self.kmer_similarity(ref_seq, asm1_seq, k=31)
        results['kmer_asm2'] = self.kmer_similarity(ref_seq, asm2_seq, k=31)

        # Step 3: Biopython global alignment
        if BIOPYTHON_AVAILABLE:
            results['alignment_asm1'] = self.biopython_alignment(ref_seq, asm1_seq)
            self._report_progress("Assembly 1 alignment complete. Running Assembly 2 alignment...", 30)
            results['alignment_asm2'] = self.biopython_alignment(ref_seq, asm2_seq)
        else:
            results['alignment_asm1'] = None
            results['alignment_asm2'] = None

        # Step 4: Coverage calculation
        results['coverage_asm1'] = self.compute_coverage(ref_seq, asm1_seq)
        results['coverage_asm2'] = self.compute_coverage(ref_seq, asm2_seq)

        # Step 5: Sliding window analysis
        self._report_progress("Running sliding window analysis for Assembly 1...", 70)
        windows1, identities1 = self.sliding_window_identity(ref_seq, asm1_seq)

        self._report_progress("Running sliding window analysis for Assembly 2...", 85)
        windows2, identities2 = self.sliding_window_identity(ref_seq, asm2_seq)

        results['window_asm1'] = {'positions': windows1, 'identities': identities1}
        results['window_asm2'] = {'positions': windows2, 'identities': identities2}

        # Step 6: Composite scoring (base score first)
        self._report_progress("Computing base composite scores...", 95)
        results['score1_base'] = self._calculate_score(results, asm=1)
        results['score2_base'] = self._calculate_score(results, asm=2)

        # Step 7: Quadripartite structure analysis (detect swap/inversion/junction errors)
        self._report_progress("Running quadripartite structure analysis...", 16)
        if len(ref_seq) > 30000:  # Only perform structure analysis for sequences > 30 kb
            results['structure1'] = self.run_structure_analysis(
                ref_seq, asm1_seq, asm_name="Assembly 1")
            results['structure2'] = self.run_structure_analysis(
                ref_seq, asm2_seq, asm_name="Assembly 2")
        else:
            results['structure1'] = {'structure_score': 5.0, 'structure_details': 'Sequence too short, structure analysis skipped'}
            results['structure2'] = {'structure_score': 5.0, 'structure_details': 'Sequence too short, structure analysis skipped'}

        # Final composite score = base score + structural score
        results['score1'] = {
            'total': round(results['score1_base']['total'] + results['structure1']['structure_score'], 2),
            'length_score': results['score1_base']['length_score'],
            'gc_score': results['score1_base']['gc_score'],
            'coverage_score': results['score1_base']['coverage_score'],
            'identity_score': results['score1_base']['identity_score'],
            'structure_score': results['structure1']['structure_score'],
        }
        results['score2'] = {
            'total': round(results['score2_base']['total'] + results['structure2']['structure_score'], 2),
            'length_score': results['score2_base']['length_score'],
            'gc_score': results['score2_base']['gc_score'],
            'coverage_score': results['score2_base']['coverage_score'],
            'identity_score': results['score2_base']['identity_score'],
            'structure_score': results['structure2']['structure_score'],
        }

        self._report_progress("Analysis complete!", 100)
        return results

    def _calculate_score(self, results, asm):
        """Calculate composite score (0-100 points)"""
        prefix = f'assembly{asm}'
        kmer_key = f'kmer_asm{asm}'
        align_key = f'alignment_asm{asm}'
        cov_key = f'coverage_asm{asm}'

        ref_len = results['reference']['length']
        asm_len = results[prefix]['length']
        ref_gc = results['reference']['gc_content']
        asm_gc = results[prefix]['gc_content']

        # Length similarity score (15%)
        length_ratio = min(ref_len, asm_len) / max(ref_len, asm_len)
        length_score = length_ratio * 15

        # GC content similarity score (15%)
        gc_diff = abs(ref_gc - asm_gc)
        gc_score = max(0, 15 - gc_diff * 1.5)

        # Coverage score (30%)
        coverage = results[cov_key]
        coverage_score = coverage / 100 * 30

        # Identity score (40%)
        if BIOPYTHON_AVAILABLE and results[align_key] is not None:
            identity = results[align_key]['identity']
            identity_score = identity / 100 * 40
        else:
            # Use k-mer similarity as fallback
            identity_score = results[kmer_key]['jaccard'] * 40

        total = length_score + gc_score + coverage_score + identity_score

        return {
            'total': round(total, 2),
            'length_score': round(length_score, 2),
            'gc_score': round(gc_score, 2),
            'coverage_score': round(coverage_score, 2),
            'identity_score': round(identity_score, 2),
        }

    # ============================================================
    # Chloroplast Genome Quadripartite Structure Analysis
    # (Detect swap/inversion/junction errors)
    # ============================================================

    def _revcomp(self, seq):
        """Return the reverse-complement of a sequence"""
        trans = str.maketrans('ACGTacgt', 'TGCAtgca')
        return seq[::-1].translate(trans)

    def _segment_map(self, ref_seq, asm_seq, num_segments=20, k=31):
        """
        Divide the reference into equal segments and map each segment to the
        assembly via k-mer matching.
        Returns: [(ref_start, ref_end, asm_best_pos, orientation, identity), ...]
          orientation: 'forward' or 'reversed'
        """
        self._report_progress("Running reference segment position mapping...", 17)
        seq_len = len(ref_seq)
        segment_size = seq_len // num_segments
        if segment_size < 1000:
            segment_size = 1000
            num_segments = max(1, seq_len // segment_size)

        # Build k-mer index for the assembly
        asm_kmer_pos = defaultdict(list)
        for i in range(len(asm_seq) - k + 1):
            kmer = asm_seq[i:i + k]
            if 'N' not in kmer:
                asm_kmer_pos[kmer].append(i)

        complement = str.maketrans('ACGTacgt', 'TGCAtgca')
        mappings = []

        for seg_idx in range(num_segments):
            seg_start = seg_idx * segment_size
            seg_end = min(seg_start + segment_size, seq_len)
            if seg_end - seg_start < 500:
                break
            seg_seq = ref_seq[seg_start:seg_end]

            # Forward matching
            fwd_positions = []
            fwd_matched = 0
            fwd_total = 0
            for i in range(0, len(seg_seq) - k + 1, 100):
                kmer = seg_seq[i:i + k]
                if 'N' in kmer:
                    continue
                fwd_total += 1
                if kmer in asm_kmer_pos:
                    fwd_matched += 1
                    fwd_positions.extend(asm_kmer_pos[kmer])

            # Reverse-complement matching
            rc_positions = []
            rc_matched = 0
            rc_total = 0
            for i in range(0, len(seg_seq) - k + 1, 100):
                kmer = seg_seq[i:i + k]
                if 'N' in kmer:
                    continue
                rc_kmer = kmer[::-1].translate(complement)
                rc_total += 1
                if rc_kmer in asm_kmer_pos:
                    rc_matched += 1
                    rc_positions.extend(asm_kmer_pos[rc_kmer])

            fwd_ratio = fwd_matched / fwd_total if fwd_total > 0 else 0
            rc_ratio = rc_matched / rc_total if rc_total > 0 else 0

            if fwd_ratio > rc_ratio:
                orientation = 'forward'
                best_pos = int(median(fwd_positions)) if fwd_positions else -1
                identity = fwd_ratio
            else:
                orientation = 'reversed'
                best_pos = int(median(rc_positions)) if rc_positions else -1
                identity = rc_ratio

            mappings.append((seg_start, seg_end, best_pos, orientation, round(identity, 4)))

        return mappings

    def _analyze_mapping_structure(self, mappings):
        """
        Analyze segment mapping results to detect structural errors.
        Returns: {
            'is_linear': bool,
            'is_global_inversion': bool,  # Whether the entire sequence is reverse-complement (global orientation)
            'inverted_segments': [],
            'position_breaks': [],
            'issues': [str, ...],
        }
        """
        valid = [m for m in mappings if m[2] >= 0]
        if len(valid) < 3:
            return {
                'is_linear': False,
                'is_global_inversion': False,
                'inverted_segments': [],
                'position_breaks': [],
                'num_breaks': 0,
                'issues': ['Unable to obtain sufficient valid mappings; sequences may be highly divergent'],
            }

        issues = []
        inverted_segments = []
        position_breaks = []

        # Check orientation
        forward_count = sum(1 for m in valid if m[3] == 'forward')
        reversed_count = sum(1 for m in valid if m[3] == 'reversed')
        total_valid = len(valid)

        # Determine if this is a global orientation inversion (nearly all segments reversed)
        is_global_inversion = (reversed_count >= total_valid * 0.85)

        if is_global_inversion:
            # Global inversion: the entire assembly is the reverse-complement of the reference;
            # not considered a structural error
            # Still check linear position order
            issues.append("Detected global orientation inversion (assembly is reverse-complement of reference)")
        elif reversed_count > total_valid * 0.3:
            # Partial inversion: structural inversion present
            issues.append("Some segments match in reverse-complement orientation (local inversions present)")
            # Only mark inverted segments when not a global inversion
            for i, m in enumerate(mappings):
                if m[2] >= 0 and m[3] == 'reversed':
                    inverted_segments.append(i)

        # Check local inversions: in non-global-inversion cases, find consecutive locally reversed segments
        if not is_global_inversion:
            in_inversion = False
            inv_start = None
            for i, m in enumerate(mappings):
                if m[2] >= 0 and m[3] == 'reversed':
                    if not in_inversion:
                        in_inversion = True
                        inv_start = i
                else:
                    if in_inversion and inv_start is not None:
                        in_inversion = False
                        if i - inv_start >= 2:
                            inverted_segments.extend(range(inv_start, i))
            if in_inversion and inv_start is not None and len(mappings) - inv_start >= 2:
                inverted_segments.extend(range(inv_start, len(mappings)))
            inverted_segments = list(set(inverted_segments))

        # Check linear position order
        # For global inversion, expect decreasing positions; otherwise expect increasing
        positions = [m[2] for m in valid]
        breaks = 0
        if is_global_inversion:
            # Global inversion: expect decreasing
            for i in range(1, len(positions)):
                if positions[i] > positions[i - 1]:  # increase = break
                    breaks += 1
                    position_breaks.append(i)
        else:
            for i in range(1, len(positions)):
                if positions[i] < positions[i - 1]:  # decrease = break
                    breaks += 1
                    position_breaks.append(i)

        is_linear = (breaks == 0)
        if breaks > 0:
            if breaks <= 2:
                issues.append(f"Found {breaks} position break(s); possible local rearrangement")
            else:
                issues.append(f"Found {breaks} position breaks; possible severe rearrangement or IR swap")

        # Detect IR swap signature: fold-back pattern in positions
        if not is_global_inversion and len(positions) >= 6:
            mid = len(positions) // 2
            first_half_inc = all(positions[i] <= positions[i + 1] for i in range(mid - 1))
            if first_half_inc:
                second_start = positions[mid]
                if any(p < second_start for p in positions[mid + 1:]):
                    issues.append("Detected possible fold-back pattern (IR swap signature)")

        return {
            'is_linear': is_linear,
            'is_global_inversion': is_global_inversion,
            'inverted_segments': inverted_segments,
            'position_breaks': position_breaks,
            'num_breaks': breaks,
            'issues': issues,
        }

    def _check_junction_kmers(self, ref_seq, asm_seq, k=31, flank=200):
        """
        Collect junction k-mer fingerprints at equally spaced points on the reference
        and check their order of occurrence in the assembly.
        No prior knowledge of IR boundaries required.
        """
        seq_len = len(ref_seq)
        # Sample 8 equally spaced "checkpoints" on the reference
        num_checks = 8
        check_positions = {}
        for i in range(num_checks):
            pos = (seq_len * (i + 1)) // (num_checks + 1)
            check_positions[f'checkpoint_{i+1}'] = pos

        # Extract k-mer sets from each checkpoint
        check_kmers = {}
        for name, pos in check_positions.items():
            start = max(0, pos - flank)
            end = min(seq_len, pos + flank)
            region = ref_seq[start:end]
            kmers = set()
            for i in range(len(region) - k + 1):
                kmer = region[i:i + k]
                if 'N' not in kmer:
                    kmers.add(kmer)
            check_kmers[name] = kmers

        # Find the best position of each checkpoint k-mer set in the assembly
        asm_positions = {}
        for name, kmers_set in check_kmers.items():
            if not kmers_set:
                continue
            positions = []
            for i in range(len(asm_seq) - k + 1):
                kmer = asm_seq[i:i + k]
                if 'N' not in kmer and kmer in kmers_set:
                    positions.append(i)
            if positions:
                asm_positions[name] = int(median(positions))

        # Check order
        ordered_names = sorted(check_positions.keys(),
                              key=lambda n: check_positions[n])
        asm_order = [(n, asm_positions.get(n, -1)) for n in ordered_names]
        valid_order = [x for x in asm_order if x[1] >= 0]

        issues = []
        if len(valid_order) >= 3:
            for i in range(1, len(valid_order)):
                if valid_order[i][1] < valid_order[i - 1][1]:
                    issues.append(
                        f"Checkpoint {valid_order[i][0]} position in assembly ({valid_order[i][1]:,}) "
                        f"appears before {valid_order[i-1][0]} ({valid_order[i-1][1]:,})"
                    )

        # Also detect local inversions: check whether each checkpoint matches forward or reverse
        complement = str.maketrans('ACGTacgt', 'TGCAtgca')
        inversion_flags = []
        for name, pos in asm_positions.items():
            if pos < 0:
                continue
            ref_pos = check_positions[name]
            ref_region = ref_seq[max(0, ref_pos - flank):min(seq_len, ref_pos + flank)]
            asm_start = max(0, pos - flank)
            asm_end = min(len(asm_seq), pos + flank)
            asm_region = asm_seq[asm_start:asm_end]

            # Forward match rate
            fwd = sum(1 for i in range(0, min(len(ref_region), len(asm_region)) - k + 1, 50)
                      if ref_region[i:i + k] == asm_region[i:i + k]) if min(len(ref_region), len(asm_region)) >= k else 0
            # Reverse match rate
            rc_ref = ref_region[::-1].translate(complement)
            rev = sum(1 for i in range(0, min(len(rc_ref), len(asm_region)) - k + 1, 50)
                      if rc_ref[i:i + k] == asm_region[i:i + k]) if min(len(rc_ref), len(asm_region)) >= k else 0

            if rev > fwd * 2:
                inversion_flags.append(f"Local inversion signal near checkpoint {name}")

        return {
            'check_positions': asm_positions,
            'order_issues': issues,
            'inversion_flags': inversion_flags,
        }

    def run_structure_analysis(self, ref_seq, asm_seq, asm_name=""):
        """
        Execute complete structure analysis (based on reference fragment mapping method).
        Does not rely on IR boundary detection; universally detects:
        1. Sequence rearrangements (swap/inversion)
        2. Local inversions
        3. Junction/order errors
        
        Returns a structure analysis result dictionary.
        """
        self._report_progress(f"Analyzing structural integrity of {asm_name}...", 16)
        
        result = {
            'name': asm_name,
            'structure_score': 5.0,
            'structure_details': '',
            'is_linear': True,
            'has_rearrangement': False,
            'inverted_segments': [],
            'order_issues': [],
            'inversion_flags': [],
            'ir_swapped': False,
            'ir_inverted': False,
            'large_inversions': [],
            'valid_junctions': 0,
            'ira_orientation': 'n/a',
            'irb_orientation': 'n/a',
            'ir_regions_found': False,
        }

        if len(ref_seq) < 30000:
            result['structure_details'] = 'Sequence too short, structure analysis skipped'
            return result

        # Step 1: Reference fragment position mapping
        mappings = self._segment_map(ref_seq, asm_seq, num_segments=20)
        result['_mappings'] = mappings

        # Step 2: Analyze mapping structure
        self._report_progress(f"Analyzing fragment mapping structure for {asm_name}...", 18)
        map_analysis = self._analyze_mapping_structure(mappings)
        result['is_linear'] = map_analysis['is_linear']
        result['has_rearrangement'] = not map_analysis['is_linear']
        result['is_global_inversion'] = map_analysis.get('is_global_inversion', False)
        result['inverted_segments'] = map_analysis['inverted_segments']
        
        # Convert inverted segments to position ranges for large_inversions
        if map_analysis['inverted_segments'] and not map_analysis.get('is_global_inversion', False):
            seg_size = len(ref_seq) // 20
            grouped = []
            sorted_inv = sorted(map_analysis['inverted_segments'])
            current_group = [sorted_inv[0]]
            for s in sorted_inv[1:]:
                if s == current_group[-1] + 1:
                    current_group.append(s)
                else:
                    grouped.append(current_group)
                    current_group = [s]
            grouped.append(current_group)
            for g in grouped:
                start = g[0] * seg_size
                end = min((g[-1] + 1) * seg_size, len(ref_seq))
                result['large_inversions'].append((start, end))
            result['ir_inverted'] = len(map_analysis['inverted_segments']) >= 2

        # Step 3: Checkpoint k-mer order verification
        self._report_progress(f"Verifying sequence order for {asm_name}...", 19)
        junction_check = self._check_junction_kmers(ref_seq, asm_seq)
        result['order_issues'] = junction_check['order_issues']
        result['inversion_flags'] = junction_check['inversion_flags']
        if junction_check['check_positions']:
            result['valid_junctions'] = len(junction_check['check_positions'])

        # Step 4: Comprehensive assessment of possible IR swap
        # IR swap hallmark: multiple segment position fold-backs + multiple inversions
        # Global inversion should not be treated as IR swap
        if not map_analysis.get('is_global_inversion', False):
            if (not map_analysis['is_linear'] and
                map_analysis.get('num_breaks', 0) >= 2 and
                len(map_analysis['inverted_segments']) >= 2):
                result['ir_swapped'] = True
            if junction_check['order_issues'] and map_analysis.get('num_breaks', 0) >= 2:
                result['ir_swapped'] = True

        # Step 5: Calculate structure score
        result['structure_score'] = self._calculate_structure_score(result)
        result['structure_details'] = self._format_structure_details(result, mappings)

        return result

    def _calculate_structure_score(self, struct_result):
        """
        Calculate structural integrity score (max 5 points).
        Distinguishes between global inversion (orientation only) and true structural errors.
        """
        score = 5.0

        # Global inversion: assembly is the reverse-complement of the reference; not a structural error
        is_global_inv = struct_result.get('is_global_inversion', False)

        # Rearrangement penalty
        if struct_result.get('has_rearrangement', False):
            if is_global_inv:
                score -= 0.5  # Global inversion + rearrangement: minor penalty
            else:
                score -= 1.5  # Non-global rearrangement: standard penalty

        # Inversion penalty (not counted for global inversion)
        if not is_global_inv:
            inv_count = len(struct_result.get('inverted_segments', []))
            score -= min(inv_count * 0.3, 2.0)  # Max deduction of 2 points

        # Order issue penalty
        order_issues = struct_result.get('order_issues', [])
        score -= len(order_issues) * 0.8

        # IR swap signature penalty
        if struct_result.get('ir_swapped', False) and not is_global_inv:
            score -= 2.0

        return max(0, round(score, 2))

    def _format_structure_details(self, struct_result, mappings=None):
        """Format structure analysis details as human-readable text"""
        lines = []

        if mappings:
            # Display mapping overview
            lines.append("Reference Fragment Mapping Results (20 equal segments):")
            seg_lines = []
            for i, (s, e, pos, orient, ident) in enumerate(mappings):
                mark = ""
                if pos < 0:
                    mark = " [unmapped]"
                elif orient == 'reversed':
                    mark = " [inverted]"
                seg_lines.append(f"  Seg{i+1:>2}: Ref[{s:>7}-{e:>7}] "
                               f"-> Asm ~{pos:>7,}bp  {orient}{mark}")
            lines.extend(seg_lines[:20])  # Display at most 20

        # Structure status summary
        lines.append(f"\nStructural Integrity Analysis:")
        if struct_result.get('is_global_inversion', False):
            lines.append("  (i) Global orientation inversion (assembly is reverse-complement of reference; not a structural error)")
        if struct_result.get('is_linear', True):
            lines.append("  [OK] Fragment mapping maintains linear order")
        else:
            lines.append("  [WARNING] Fragment mapping shows position breaks; possible rearrangement")

        inv_segs = struct_result.get('inverted_segments', [])
        if inv_segs and not struct_result.get('is_global_inversion', False):
            lines.append(f"  [WARNING] Detected {len(inv_segs)} locally inverted segment(s): {sorted(inv_segs)[:10]}")
        elif struct_result.get('is_global_inversion', False):
            lines.append("  (i) Global orientation inversion (not a local inversion)")
        else:
            lines.append("  [OK] No local inversions detected")

        order_issues = struct_result.get('order_issues', [])
        if order_issues:
            lines.append(f"  [WARNING] Order issues ({len(order_issues)} item(s)):")
            for oi in order_issues[:5]:
                lines.append(f"    - {oi}")
        else:
            lines.append("  [OK] Checkpoint k-mer order normal")

        inversion_flags = struct_result.get('inversion_flags', [])
        if inversion_flags:
            for f in inversion_flags[:3]:
                lines.append(f"  [WARNING] {f}")

        if struct_result.get('ir_swapped', False) and not struct_result.get('is_global_inversion', False):
            lines.append("  [WARNING] IR swap signature detected (multiple fold-backs + local inversions)")
        elif struct_result.get('is_global_inversion', False):
            lines.append("  (i) Global inversion mode; IR swap detection not applicable")
        else:
            lines.append("  [OK] No IR swap signature detected")

        large_inv = struct_result.get('large_inversions', [])
        if large_inv:
            lines.append(f"  Large inversion regions ({len(large_inv)} region(s)):")
            for s, e in large_inv[:5]:
                lines.append(f"    [{s:,}-{e:,}] (~{(e-s)//1000} kb)")

        lines.append(f"\nStructural Integrity Score: {struct_result.get('structure_score', 0):.2f} / 5.0")

        return '\n'.join(lines)


# ============================================================
# Graphical User Interface
# ============================================================
class CpGenomeComparatorApp:
    """Chloroplast genome comparison tool GUI"""

    def __init__(self, root):
        self.root = root
        self.root.title("Chloroplast Genome Assembly Comparison Tool")
        self.root.geometry("1000x750")
        self.root.minsize(900, 600)

        # Setup style
        self._setup_style()
        # Create widgets
        self._create_widgets()
        # Analysis state
        self.analysis_running = False
        self.results = None

    def _setup_style(self):
        style = ttk.Style()
        style.theme_use('clam')

        # Custom color scheme
        self.colors = {
            'bg': '#f5f6fa',
            'fg': '#2d3436',
            'accent': '#0984e3',
            'success': '#00b894',
            'warning': '#fdcb6e',
            'danger': '#d63031',
            'light': '#dfe6e9',
            'card': '#ffffff',
            'border': '#b2bec3',
        }

        self.root.configure(bg=self.colors['bg'])

        style.configure('Title.TLabel',
                        font=('Microsoft YaHei', 16, 'bold'),
                        foreground=self.colors['accent'],
                        background=self.colors['bg'])

        style.configure('Section.TLabel',
                        font=('Microsoft YaHei', 11, 'bold'),
                        foreground=self.colors['fg'],
                        background=self.colors['bg'])

        style.configure('Card.TFrame',
                        background=self.colors['card'],
                        relief='solid',
                        borderwidth=1)

        style.configure('Accent.TButton',
                        font=('Microsoft YaHei', 11, 'bold'),
                        background=self.colors['accent'])

        style.configure('Success.TLabel',
                        font=('Microsoft YaHei', 12, 'bold'),
                        foreground=self.colors['success'],
                        background=self.colors['card'])

        style.configure('Danger.TLabel',
                        font=('Microsoft YaHei', 10),
                        foreground=self.colors['danger'],
                        background=self.colors['card'])

    def _create_widgets(self):
        # Main container
        main_frame = tk.Frame(self.root, bg=self.colors['bg'])
        main_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=10)

        # ----- Title -----
        title_frame = tk.Frame(main_frame, bg=self.colors['bg'])
        title_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(title_frame, text="🧬 Chloroplast Genome Assembly Comparison Tool",
                  style='Title.TLabel').pack(side=tk.LEFT)
        ttk.Label(title_frame,
                  text="Compare two assemblies against a reference to determine which is more accurate",
                  font=('Microsoft YaHei', 9),
                  foreground='#636e72',
                  background=self.colors['bg']).pack(side=tk.LEFT, padx=(15, 0))

        # ----- File Selection Area -----
        file_card = ttk.Frame(main_frame, style='Card.TFrame')
        file_card.pack(fill=tk.X, pady=(0, 10), ipady=8)

        ttk.Label(file_card, text="📂 Select Sequence Files (FASTA format)",
                  style='Section.TLabel').pack(anchor=tk.W, padx=12, pady=(8, 5))

        # Reference sequence
        self._create_file_row(file_card, "Reference Sequence:", 'ref', 0)
        # Assembly 1
        self._create_file_row(file_card, "Assembly 1:", 'asm1', 1)
        # Assembly 2
        self._create_file_row(file_card, "Assembly 2:", 'asm2', 2)

        # ----- Control Buttons -----
        btn_frame = tk.Frame(main_frame, bg=self.colors['bg'])
        btn_frame.pack(fill=tk.X, pady=(0, 10))

        self.analyze_btn = tk.Button(
            btn_frame, text="▶  Start Analysis",
            font=('Microsoft YaHei', 12, 'bold'),
            bg=self.colors['accent'], fg='white',
            activebackground='#076bb5', activeforeground='white',
            relief=tk.FLAT, padx=25, pady=6,
            cursor='hand2', borderwidth=0,
            command=self._start_analysis
        )
        self.analyze_btn.pack(side=tk.LEFT, padx=(0, 10))

        self.clear_btn = tk.Button(
            btn_frame, text="🗑 Clear Results",
            font=('Microsoft YaHei', 10),
            bg=self.colors['light'], fg=self.colors['fg'],
            activebackground='#c8d6e5', activeforeground=self.colors['fg'],
            relief=tk.FLAT, padx=18, pady=5,
            cursor='hand2', borderwidth=0,
            command=self._clear_results
        )
        self.clear_btn.pack(side=tk.LEFT, padx=(0, 10))

        self.export_btn = tk.Button(
            btn_frame, text="📋 Export Report",
            font=('Microsoft YaHei', 10),
            bg=self.colors['light'], fg=self.colors['fg'],
            activebackground='#c8d6e5', activeforeground=self.colors['fg'],
            relief=tk.FLAT, padx=18, pady=5,
            cursor='hand2', borderwidth=0,
            command=self._export_report,
            state=tk.DISABLED
        )
        self.export_btn.pack(side=tk.LEFT)

        # Progress bar
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(
            btn_frame, variable=self.progress_var,
            maximum=100, length=200
        )
        self.progress_bar.pack(side=tk.RIGHT, padx=(10, 0))

        self.progress_label = ttk.Label(
            btn_frame, text="",
            font=('Microsoft YaHei', 9),
            foreground='#636e72',
            background=self.colors['bg']
        )
        self.progress_label.pack(side=tk.RIGHT, padx=(5, 0))

        # ----- Results Area (Notebook tabs) -----
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # == Tab 1: Summary ==
        self.summary_frame = ttk.Frame(self.notebook, style='Card.TFrame')
        self.notebook.add(self.summary_frame, text="📊 Summary")

        self.summary_text = scrolledtext.ScrolledText(
            self.summary_frame,
            font=('Consolas', 10),
            wrap=tk.WORD,
            bg='#ffffff',
            fg='#2d3436',
            relief=tk.FLAT,
            borderwidth=0,
            padx=12,
            pady=10,
        )
        self.summary_text.pack(fill=tk.BOTH, expand=True)
        self.summary_text.insert(tk.END, "Please select three FASTA files and click 「Start Analysis」.\n\n")
        self.summary_text.insert(tk.END, "Supported formats: .fasta, .fa, .fna, .txt\n")
        self.summary_text.insert(tk.END, "Select one file each for the reference and two assemblies.\n")
        self.summary_text.configure(state=tk.DISABLED)

        # == Tab 2: Details ==
        self.detail_frame = ttk.Frame(self.notebook, style='Card.TFrame')
        self.notebook.add(self.detail_frame, text="📋 Details")

        self.detail_text = scrolledtext.ScrolledText(
            self.detail_frame,
            font=('Consolas', 10),
            wrap=tk.WORD,
            bg='#ffffff',
            fg='#2d3436',
            relief=tk.FLAT,
            borderwidth=0,
            padx=12,
            pady=10,
        )
        self.detail_text.pack(fill=tk.BOTH, expand=True)
        self.detail_text.insert(tk.END, "Detailed comparison data will be displayed after analysis.\n")
        self.detail_text.configure(state=tk.DISABLED)

        # == Tab 3: Identity Plot ==
        self.canvas_frame = ttk.Frame(self.notebook, style='Card.TFrame')
        self.notebook.add(self.canvas_frame, text="📈 Identity Plot")

        canvas_container = tk.Frame(self.canvas_frame, bg='white')
        canvas_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.canvas = tk.Canvas(
            canvas_container,
            bg='white',
            highlightthickness=0,
            height=300,
        )
        self.canvas.pack(fill=tk.BOTH, expand=True)

        # Bottom status bar
        status_frame = tk.Frame(main_frame, bg=self.colors['bg'])
        status_frame.pack(fill=tk.X, pady=(5, 0))
        self.status_label = ttk.Label(
            status_frame, text="Ready",
            font=('Microsoft YaHei', 9),
            foreground='#636e72',
            background=self.colors['bg']
        )
        self.status_label.pack(side=tk.LEFT)

    def _create_file_row(self, parent, label_text, key, row):
        """Create a file selection row"""
        row_frame = tk.Frame(parent, bg=self.colors['card'])
        row_frame.pack(fill=tk.X, padx=12, pady=3)

        ttk.Label(row_frame, text=label_text,
                  font=('Microsoft YaHei', 10),
                  background=self.colors['card']).pack(side=tk.LEFT, padx=(0, 8))

        entry = ttk.Entry(row_frame, font=('Microsoft YaHei', 10), width=55)
        entry.pack(side=tk.LEFT, padx=(0, 5))

        def browse():
            filepath = filedialog.askopenfilename(
                title=f"Select {label_text.rstrip(':')} File",
                filetypes=[
                    ("FASTA files", "*.fasta *.fa *.fna *.fas *.txt"),
                    ("All files", "*.*"),
                ]
            )
            if filepath:
                entry.delete(0, tk.END)
                entry.insert(0, filepath)

        browse_btn = tk.Button(
            row_frame, text="Browse...",
            font=('Microsoft YaHei', 9),
            bg=self.colors['accent'], fg='white',
            activebackground='#076bb5', activeforeground='white',
            relief=tk.FLAT, padx=10, pady=2,
            cursor='hand2', borderwidth=0,
            command=browse
        )
        browse_btn.pack(side=tk.LEFT)

        setattr(self, f'{key}_entry', entry)

    def _start_analysis(self):
        """Start analysis (runs in a thread)"""
        if self.analysis_running:
            messagebox.showwarning("Notice", "Analysis is in progress. Please wait for completion.")
            return

        # Get file paths
        ref_path = self.ref_entry.get().strip()
        asm1_path = self.asm1_entry.get().strip()
        asm2_path = self.asm2_entry.get().strip()

        # Validate
        errors = []
        if not ref_path or not os.path.exists(ref_path):
            errors.append("Please select a valid reference sequence file.")
        if not asm1_path or not os.path.exists(asm1_path):
            errors.append("Please select a valid Assembly 1 file.")
        if not asm2_path or not os.path.exists(asm2_path):
            errors.append("Please select a valid Assembly 2 file.")

        if errors:
            messagebox.showerror("Input Error", "\n".join(errors))
            return

        if not BIOPYTHON_AVAILABLE:
            messagebox.showerror(
                "Missing Dependency",
                "Biopython library not detected.\n\n"
                "Please run the following command to install:\n"
                "  pip install biopython\n\n"
                "Restart the program after installation."
            )
            return

        # Disable buttons
        self.analysis_running = True
        self.analyze_btn.configure(
            text="⏳ Analyzing...",
            state=tk.DISABLED,
            bg='#636e72'
        )
        self.export_btn.configure(state=tk.DISABLED)
        self.progress_var.set(0)
        self.progress_label.configure(text="Initializing...")
        self.status_label.configure(text="Analyzing...")

        # Run analysis in a thread
        thread = threading.Thread(
            target=self._run_analysis_thread,
            args=(ref_path, asm1_path, asm2_path),
            daemon=True
        )
        thread.start()

    def _run_analysis_thread(self, ref_path, asm1_path, asm2_path):
        """Background analysis thread"""
        try:
            def progress_callback(msg, pct):
                self.root.after(0, self._update_progress, msg, pct)

            analyzer = SequenceAnalyzer(progress_callback=progress_callback)
            results = analyzer.run_full_analysis(ref_path, asm1_path, asm2_path)
            self.results = results
            self.root.after(0, self._display_results, results)
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, self._handle_error, error_msg)

    def _update_progress(self, msg, pct):
        """Update progress bar"""
        self.progress_var.set(pct)
        self.progress_label.configure(text=msg)

    def _handle_error(self, error_msg):
        """Handle analysis errors"""
        self.analysis_running = False
        self.analyze_btn.configure(
            text="▶  Start Analysis",
            state=tk.NORMAL,
            bg=self.colors['accent']
        )
        self.progress_label.configure(text="Analysis error")
        self.status_label.configure(text="Error")
        messagebox.showerror("Analysis Error", f"An error occurred during analysis:\n\n{error_msg}")

    def _display_results(self, results):
        """Display analysis results"""
        self.analysis_running = False
        self.analyze_btn.configure(
            text="▶  Start Analysis",
            state=tk.NORMAL,
            bg=self.colors['accent']
        )
        self.export_btn.configure(state=tk.NORMAL)
        self.progress_label.configure(text="Complete")
        self.status_label.configure(text="Analysis complete")

        score1 = results['score1']['total']
        score2 = results['score2']['total']

        if score1 > score2:
            winner = "Assembly 1"
            winner_score = score1
            loser = "Assembly 2"
            loser_score = score2
        elif score2 > score1:
            winner = "Assembly 2"
            winner_score = score2
            loser = "Assembly 1"
            loser_score = score1
        else:
            winner = "Both assemblies have the same score"
            winner_score = score1
            loser = ""
            loser_score = score2

        # Tab 1: Summary
        self._display_summary(results, winner, winner_score, loser, loser_score)
        # Tab 2: Details
        self._display_details(results)
        # Tab 3: Plot
        self._draw_window_plot(results)

        # Switch to Summary tab
        self.notebook.select(0)

    def _display_summary(self, results, winner, winner_score, loser, loser_score):
        """Display summary results"""
        self.summary_text.configure(state=tk.NORMAL)
        self.summary_text.delete(1.0, tk.END)

        # Configure text tags
        self.summary_text.tag_configure("title", font=('Consolas', 14, 'bold'),
                                        foreground='#0984e3')
        self.summary_text.tag_configure("winner", font=('Consolas', 13, 'bold'),
                                        foreground='#00b894', background='#f0fff4')
        self.summary_text.tag_configure("header", font=('Consolas', 11, 'bold'),
                                        foreground='#2d3436')
        self.summary_text.tag_configure("good", foreground='#00b894')
        self.summary_text.tag_configure("bad", foreground='#d63031')

        ref = results['reference']
        asm1 = results['assembly1']
        asm2 = results['assembly2']
        s1 = results['score1']
        s2 = results['score2']

        self.summary_text.insert(tk.END, "═" * 60 + "\n", "header")
        self.summary_text.insert(tk.END, "  Chloroplast Genome Assembly Comparison - Summary\n", "title")
        self.summary_text.insert(tk.END, "═" * 60 + "\n\n", "header")

        # Winner result
        if "both" not in winner.lower():
            self.summary_text.insert(tk.END,
                f"  ★★★  {winner} is closer to the reference ★★★\n\n", "winner")
        else:
            self.summary_text.insert(tk.END,
                f"  {winner}\n\n", "header")

        # Scoring table
        self.summary_text.insert(tk.END, "┌" + "─" * 56 + "┐\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^15} │ {:^15} │\n".format("Dimension", "Assembly 1", "Assembly 2"))
        self.summary_text.insert(tk.END, "├" + "─" * 56 + "┤\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}  │ {:^13.2f}  │\n".format("Identity (40%)", s1['identity_score'], s2['identity_score']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}  │ {:^13.2f}  │\n".format("Coverage (30%)", s1['coverage_score'], s2['coverage_score']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}  │ {:^13.2f}  │\n".format("Length Sim. (15%)", s1['length_score'], s2['length_score']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}  │ {:^13.2f}  │\n".format("GC Sim. (15%)", s1['gc_score'], s2['gc_score']))
        # New structural integrity score
        if 'structure_score' in s1:
            self.summary_text.insert(tk.END,
                "│ {:^20} │ {:^13.2f}  │ {:^13.2f}  │\n".format("Structure (+5)", s1['structure_score'], s2['structure_score']))
        self.summary_text.insert(tk.END, "├" + "─" * 56 + "┤\n")

        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}  │ {:^13.2f}  │\n".format("Total Score ★", s1['total'], s2['total']))
        self.summary_text.insert(tk.END, "└" + "─" * 56 + "┘\n\n")

        # Basic statistics comparison
        self.summary_text.insert(tk.END, "Basic Statistics:\n", "header")
        self.summary_text.insert(tk.END,
            "┌" + "─" * 56 + "┐\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^15} │ {:^15} │\n".format("Statistic", "Reference", "Assembly 1", "Assembly 2"))
        self.summary_text.insert(tk.END, "├" + "─" * 56 + "┤\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13,}bp │ {:^13,}bp │\n".format("Length", ref['length'], asm1['length']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13} │ {:^13,}bp │\n".format("", "", asm2['length']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}%  │ {:^13.2f}%  │\n".format("GC Content", ref['gc_content'], asm1['gc_content']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13} │ {:^13.2f}%  │\n".format("", "", asm2['gc_content']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13}  │ {:^13}  │\n".format("N Bases", ref['N_count'], asm1['N_count']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13}  │ {:^13}  │\n".format("", "", asm2['N_count']))
        self.summary_text.insert(tk.END, "└" + "─" * 56 + "┘\n\n")

        # K-mer analysis summary
        km1 = results['kmer_asm1']
        km2 = results['kmer_asm2']
        self.summary_text.insert(tk.END, "K-mer Similarity Analysis (k=31):\n", "header")
        self.summary_text.insert(tk.END,
            f"  Assembly 1: Jaccard={km1['jaccard']:.4f},  Shared Ratio={km1['shared_ratio']:.4f}\n")
        self.summary_text.insert(tk.END,
            f"  Assembly 2: Jaccard={km2['jaccard']:.4f},  Shared Ratio={km2['shared_ratio']:.4f}\n\n")

        # Coverage
        self.summary_text.insert(tk.END, "Reference Coverage:\n", "header")
        self.summary_text.insert(tk.END,
            f"  Assembly 1: {results['coverage_asm1']:.2f}%\n")
        self.summary_text.insert(tk.END,
            f"  Assembly 2: {results['coverage_asm2']:.2f}%\n\n")

        # Quadripartite structure analysis summary
        for st_key, st_label in [('structure1', 'Assembly 1'), ('structure2', 'Assembly 2')]:
            st = results.get(st_key)
            if not st or not isinstance(st, dict):
                continue
            if st.get('ir_regions_found'):
                self.summary_text.insert(tk.END, f"{st_label} Quadripartite Structure Analysis:\n", "header")
                ref_reg = st.get('ref_regions', {})
                if ref_reg:
                    self.summary_text.insert(tk.END,
                        f"  Reference Regions: LSC={ref_reg.get('lsc_length','?'):,}bp, "
                        f"IRa={ref_reg.get('ira_length','?'):,}bp, "
                        f"SSC={ref_reg.get('ssc_length','?'):,}bp, "
                        f"IRb={ref_reg.get('irb_length','?'):,}bp\n")
                self.summary_text.insert(tk.END,
                    f"  IR Orientation: IRa={st.get('ira_orientation','?')}, "
                    f"IRb={st.get('irb_orientation','?')}\n")
                issues = []
                if st.get('ir_swapped'):
                    issues.append("IRa/IRb swap")
                if st.get('ir_inverted'):
                    issues.append("IR inversion")
                inversions = st.get('large_inversions', [])
                if inversions:
                    issues.append(f"Large inversions ({len(inversions)} regions)")
                for ji in st.get('junction_issues', []):
                    issues.append(ji.replace('⚠ ', ''))
                if issues:
                    self.summary_text.insert(tk.END,
                        f"  [WARNING] Structure issues: {'; '.join(issues)}\n")
                else:
                    self.summary_text.insert(tk.END,
                        "  [OK] Quadripartite structure normal; no swap/inversion/junction errors detected\n")
                self.summary_text.insert(tk.END,
                    f"  Structure Score: {st.get('structure_score', 0):.2f} / 5.0\n\n")

        # Alignment results
        if BIOPYTHON_AVAILABLE:
            al1 = results.get('alignment_asm1')
            al2 = results.get('alignment_asm2')
            if al1 and al2:
                self.summary_text.insert(tk.END, "Global Alignment Results:\n", "header")
                self.summary_text.insert(tk.END,
                    f"  Assembly 1: Identity={al1['identity']:.2f}%,  Aligned Length={al1['aligned_length']:,}\n")
                self.summary_text.insert(tk.END,
                    f"  Assembly 2: Identity={al2['identity']:.2f}%,  Aligned Length={al2['aligned_length']:,}\n\n")

        # Conclusion
        self.summary_text.insert(tk.END, "═" * 60 + "\n", "header")
        self.summary_text.insert(tk.END,
            "Conclusion: Use the higher-scoring assembly for downstream analysis.\n", "header")
        self.summary_text.insert(tk.END,
            "If scores are close, examine the sliding window identity plot\n")
        self.summary_text.insert(tk.END,
            "for large-scale rearrangements or inversions.\n")

        self.summary_text.configure(state=tk.DISABLED)

    def _display_details(self, results):
        """Display detailed comparison information"""
        self.detail_text.configure(state=tk.NORMAL)
        self.detail_text.delete(1.0, tk.END)

        ref = results['reference']
        asm1 = results['assembly1']
        asm2 = results['assembly2']

        self.detail_text.insert(tk.END, "=" * 70 + "\n")
        self.detail_text.insert(tk.END, "  Detailed Comparison Report\n")
        self.detail_text.insert(tk.END, "=" * 70 + "\n\n")

        # Basic sequence info
        self.detail_text.insert(tk.END, "[Basic Sequence Information]\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  Reference:    Length={ref['length']:,} bp,  GC={ref['gc_content']:.2f}%,  "
            f"A={ref['A_pct']:.1f}% T={ref['T_pct']:.1f}% C={ref['C_pct']:.1f}% G={ref['G_pct']:.1f}%\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 1:   Length={asm1['length']:,} bp,  GC={asm1['gc_content']:.2f}%,  "
            f"A={asm1['A_pct']:.1f}% T={asm1['T_pct']:.1f}% C={asm1['C_pct']:.1f}% G={asm1['G_pct']:.1f}%\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 2:   Length={asm2['length']:,} bp,  GC={asm2['gc_content']:.2f}%,  "
            f"A={asm2['A_pct']:.1f}% T={asm2['T_pct']:.1f}% C={asm2['C_pct']:.1f}% G={asm2['G_pct']:.1f}%\n\n")

        # Length differences
        len_diff1 = abs(ref['length'] - asm1['length'])
        len_diff_pct1 = len_diff1 / ref['length'] * 100 if ref['length'] > 0 else 0
        len_diff2 = abs(ref['length'] - asm2['length'])
        len_diff_pct2 = len_diff2 / ref['length'] * 100 if ref['length'] > 0 else 0
        self.detail_text.insert(tk.END, "[Length Differences]\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 1 vs Reference: {len_diff1:,} bp ({len_diff_pct1:.2f}%)\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 2 vs Reference: {len_diff2:,} bp ({len_diff_pct2:.2f}%)\n\n")

        # GC content differences
        gc_diff1 = abs(ref['gc_content'] - asm1['gc_content'])
        gc_diff2 = abs(ref['gc_content'] - asm2['gc_content'])
        self.detail_text.insert(tk.END, "[GC Content Differences]\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 1 vs Reference: {gc_diff1:.2f}%\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 2 vs Reference: {gc_diff2:.2f}%\n\n")

        # K-mer analysis
        km1 = results['kmer_asm1']
        km2 = results['kmer_asm2']
        self.detail_text.insert(tk.END, "[K-mer Similarity Analysis (k=31)]\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 1:\n"
            f"    Reference k-mers: {km1['total_kmers_ref']:,}\n"
            f"    Assembly k-mers:  {km1['total_kmers_asm']:,}\n"
            f"    Shared k-mers:    {km1['shared_kmers']:,}\n"
            f"    Jaccard coeff:    {km1['jaccard']:.6f}\n"
            f"    Shared ratio:     {km1['shared_ratio']:.6f}\n\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 2:\n"
            f"    Reference k-mers: {km2['total_kmers_ref']:,}\n"
            f"    Assembly k-mers:  {km2['total_kmers_asm']:,}\n"
            f"    Shared k-mers:    {km2['shared_kmers']:,}\n"
            f"    Jaccard coeff:    {km2['jaccard']:.6f}\n"
            f"    Shared ratio:     {km2['shared_ratio']:.6f}\n\n")

        # Alignment details
        if BIOPYTHON_AVAILABLE:
            al1 = results.get('alignment_asm1')
            al2 = results.get('alignment_asm2')
            if al1 and al2:
                self.detail_text.insert(tk.END, "[Global Alignment Details (Biopython)]\n")
                self.detail_text.insert(tk.END, "-" * 70 + "\n")
                self.detail_text.insert(tk.END,
                    f"  Assembly 1:\n"
                    f"    Identity:       {al1['identity']:.2f}%\n"
                    f"    Matches:        {al1['matches']:,}\n"
                    f"    Mismatches:     {al1['mismatches']:,}\n"
                    f"    Ref gaps:       {al1['gaps_ref']:,}\n"
                    f"    Asm gaps:       {al1['gaps_asm']:,}\n"
                    f"    Aligned length: {al1['aligned_length']:,}\n"
                    f"    Score:          {al1['score']:.0f}\n\n")
                self.detail_text.insert(tk.END,
                    f"  Assembly 2:\n"
                    f"    Identity:       {al2['identity']:.2f}%\n"
                    f"    Matches:        {al2['matches']:,}\n"
                    f"    Mismatches:     {al2['mismatches']:,}\n"
                    f"    Ref gaps:       {al2['gaps_ref']:,}\n"
                    f"    Asm gaps:       {al2['gaps_asm']:,}\n"
                    f"    Aligned length: {al2['aligned_length']:,}\n"
                    f"    Score:          {al2['score']:.0f}\n\n")

        # Coverage
        self.detail_text.insert(tk.END, "[Reference Coverage]\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 1: {results['coverage_asm1']:.2f}%\n")
        self.detail_text.insert(tk.END,
            f"  Assembly 2: {results['coverage_asm2']:.2f}%\n\n")

        # Score breakdown
        self.detail_text.insert(tk.END, "[Composite Score Breakdown]\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        s1 = results['score1']
        s2 = results['score2']
        self.detail_text.insert(tk.END,
            f"  {'Dimension':<20} {'Weight':>8} {'Asm 1':>10} {'Asm 2':>10}\n")
        self.detail_text.insert(tk.END, f"  {'-'*50}\n")
        self.detail_text.insert(tk.END,
            f"  {'Identity Score':<20} {'40%':>8} {s1['identity_score']:>8.2f}  {s2['identity_score']:>8.2f}\n")
        self.detail_text.insert(tk.END,
            f"  {'Coverage Score':<20} {'30%':>8} {s1['coverage_score']:>8.2f}  {s2['coverage_score']:>8.2f}\n")
        self.detail_text.insert(tk.END,
            f"  {'Length Similarity':<20} {'15%':>8} {s1['length_score']:>8.2f}  {s2['length_score']:>8.2f}\n")
        self.detail_text.insert(tk.END,
            f"  {'GC Similarity':<20} {'15%':>8} {s1['gc_score']:>8.2f}  {s2['gc_score']:>8.2f}\n")
        if 'structure_score' in s1:
            self.detail_text.insert(tk.END,
                f"  {'Structure Integrity':<20} {'+5%':>8} {s1['structure_score']:>8.2f}  {s2['structure_score']:>8.2f}\n")
        self.detail_text.insert(tk.END, f"  {'-'*50}\n")
        self.detail_text.insert(tk.END,
            f"  {'Total Score':<20} {'105%':>8} {s1['total']:>8.2f}  {s2['total']:>8.2f}\n\n")

        # Quadripartite structure analysis details
        for st_key, st_label in [('structure1', 'Assembly 1'), ('structure2', 'Assembly 2')]:
            st = results.get(st_key)
            if not st or not isinstance(st, dict):
                continue
            details = st.get('structure_details', '')
            if details and details != 'Sequence too short, structure analysis skipped':
                self.detail_text.insert(tk.END, f"[Quadripartite Structure Analysis — {st_label}]\n")
                self.detail_text.insert(tk.END, "-" * 70 + "\n")
                self.detail_text.insert(tk.END, f"  {details}\n\n")

        self.detail_text.configure(state=tk.DISABLED)

    def _draw_window_plot(self, results):
        """Draw sliding window identity plot on Canvas"""
        self.canvas.delete("all")

        w1 = results['window_asm1']
        w2 = results['window_asm2']
        ref_len = results['reference']['length']

        canvas_width = self.canvas.winfo_width()
        canvas_height = self.canvas.winfo_height()

        if canvas_width < 100:
            canvas_width = 900
        if canvas_height < 100:
            canvas_height = 300

        margin_left = 80
        margin_right = 30
        margin_top = 30
        margin_bottom = 50
        plot_width = canvas_width - margin_left - margin_right
        plot_height = canvas_height - margin_top - margin_bottom

        # Background
        self.canvas.create_rectangle(
            margin_left, margin_top,
            margin_left + plot_width, margin_top + plot_height,
            fill='#fafafa', outline='#e0e0e0'
        )

        # Y-axis and labels
        for pct in [0, 25, 50, 75, 100]:
            y = margin_top + plot_height - (pct / 100 * plot_height)
            self.canvas.create_line(margin_left - 5, y, margin_left + plot_width, y,
                                    fill='#e0e0e0', dash=(2, 4))
            self.canvas.create_text(margin_left - 10, y,
                                    text=f"{pct}%", anchor=tk.E,
                                    font=('Consolas', 8), fill='#636e72')

        # X-axis labels
        num_labels = 8
        for i in range(num_labels + 1):
            pos = int(ref_len * i / num_labels)
            x = margin_left + (plot_width * i / num_labels)
            self.canvas.create_text(x, margin_top + plot_height + 15,
                                    text=f"{pos:,}", anchor=tk.N,
                                    font=('Consolas', 7), fill='#636e72')

        # Axis titles
        self.canvas.create_text(margin_left + plot_width // 2,
                                margin_top + plot_height + 35,
                                text="Reference Position (bp)", anchor=tk.N,
                                font=('Microsoft YaHei', 9), fill='#2d3436')
        self.canvas.create_text(15, margin_top + plot_height // 2,
                                text="I\nd\ne\nn\nt\ni\nt\ny\n%", anchor=tk.CENTER,
                                font=('Microsoft YaHei', 9), fill='#2d3436')

        # Draw data lines
        def draw_line(windows, identities, color, label):
            if not windows or len(windows) < 2:
                return
            points = []
            for pos, ident in zip(windows, identities):
                x = margin_left + (pos / ref_len * plot_width) if ref_len > 0 else margin_left
                y = margin_top + plot_height - (ident / 100 * plot_height)
                points.extend([x, y])
            if points:
                self.canvas.create_line(*points, fill=color, width=2, smooth=True)
                # Legend
                legend_y = margin_top + (20 if label == "Assembly 1" else 42)
                self.canvas.create_line(
                    margin_left + plot_width - 160, legend_y,
                    margin_left + plot_width - 120, legend_y,
                    fill=color, width=3
                )
                self.canvas.create_text(
                    margin_left + plot_width - 115, legend_y,
                    text=label, anchor=tk.W,
                    font=('Microsoft YaHei', 9), fill=color
                )

        draw_line(w1['positions'], w1['identities'], '#0984e3', 'Assembly 1')
        draw_line(w2['positions'], w2['identities'], '#d63031', 'Assembly 2')

        # Title
        self.canvas.create_text(margin_left + plot_width // 2, 10,
                                text="Sliding Window Identity Plot (window=5000 bp, step=2000 bp)",
                                anchor=tk.N,
                                font=('Microsoft YaHei', 10, 'bold'),
                                fill='#2d3436')

    def _clear_results(self):
        """Clear all results"""
        if self.analysis_running:
            messagebox.showwarning("Notice", "Analysis is in progress. Please wait for completion before clearing.")
            return

        self.results = None
        self.summary_text.configure(state=tk.NORMAL)
        self.summary_text.delete(1.0, tk.END)
        self.summary_text.insert(tk.END, "Please select three FASTA files and click 「Start Analysis」.\n\n")
        self.summary_text.insert(tk.END, "Supported formats: .fasta, .fa, .fna, .txt\n")
        self.summary_text.insert(tk.END, "Select one file each for the reference and two assemblies.\n")
        self.summary_text.configure(state=tk.DISABLED)

        self.detail_text.configure(state=tk.NORMAL)
        self.detail_text.delete(1.0, tk.END)
        self.detail_text.insert(tk.END, "Detailed comparison data will be displayed after analysis.\n")
        self.detail_text.configure(state=tk.DISABLED)

        self.canvas.delete("all")
        self.progress_var.set(0)
        self.progress_label.configure(text="")
        self.export_btn.configure(state=tk.DISABLED)
        self.status_label.configure(text="Results cleared")

    def _export_report(self):
        """Export report to file"""
        if not self.results:
            messagebox.showwarning("Notice", "No results to export. Please run an analysis first.")
            return

        filepath = filedialog.asksaveasfilename(
            title="Export Report",
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ],
            initialfile="chloroplast_genome_comparison_report.txt"
        )
        if not filepath:
            return

        try:
            results = self.results
            ref = results['reference']
            asm1 = results['assembly1']
            asm2 = results['assembly2']
            s1 = results['score1']
            s2 = results['score2']

            with open(filepath, 'w', encoding='utf-8') as f:
                f.write("=" * 70 + "\n")
                f.write("  Chloroplast Genome Assembly Comparison Report\n")
                f.write("=" * 70 + "\n\n")

                if s1['total'] > s2['total']:
                    f.write(f"Conclusion: Assembly 1 is closer to the reference (score {s1['total']:.2f} vs {s2['total']:.2f})\n\n")
                elif s2['total'] > s1['total']:
                    f.write(f"Conclusion: Assembly 2 is closer to the reference (score {s2['total']:.2f} vs {s1['total']:.2f})\n\n")
                else:
                    f.write(f"Conclusion: Both assemblies have the same score ({s1['total']:.2f})\n\n")

                f.write(f"[Basic Sequence Information]\n")
                f.write(f"  Reference: {ref['length']:,} bp, GC={ref['gc_content']:.2f}%\n")
                f.write(f"  Assembly 1: {asm1['length']:,} bp, GC={asm1['gc_content']:.2f}%\n")
                f.write(f"  Assembly 2: {asm2['length']:,} bp, GC={asm2['gc_content']:.2f}%\n\n")

                km1 = results['kmer_asm1']
                km2 = results['kmer_asm2']
                f.write(f"[K-mer Similarity (k=31)]\n")
                f.write(f"  Assembly 1: Jaccard={km1['jaccard']:.6f}, Shared Ratio={km1['shared_ratio']:.6f}\n")
                f.write(f"  Assembly 2: Jaccard={km2['jaccard']:.6f}, Shared Ratio={km2['shared_ratio']:.6f}\n\n")

                f.write(f"[Reference Coverage]\n")
                f.write(f"  Assembly 1: {results['coverage_asm1']:.2f}%\n")
                f.write(f"  Assembly 2: {results['coverage_asm2']:.2f}%\n\n")

                # Quadripartite structure analysis
                for st_key, st_label in [('structure1', 'Assembly 1'), ('structure2', 'Assembly 2')]:
                    st = results.get(st_key)
                    if not st or not isinstance(st, dict):
                        continue
                    details = st.get('structure_details', '')
                    if details and details != 'Sequence too short, structure analysis skipped':
                        f.write(f"[Quadripartite Structure Analysis — {st_label}]\n")
                        f.write(f"{details}\n\n")

                f.write(f"[Composite Score]\n")
                f.write(f"  Assembly 1: {s1['total']:.2f}")
                if 'structure_score' in s1:
                    f.write(f" (Base {results.get('score1_base', s1)['total']:.2f} + Structure {s1['structure_score']:.2f})")
                f.write(f"\n")
                f.write(f"  Assembly 2: {s2['total']:.2f}")
                if 'structure_score' in s2:
                    f.write(f" (Base {results.get('score2_base', s2)['total']:.2f} + Structure {s2['structure_score']:.2f})")
                f.write(f"\n")

            messagebox.showinfo("Export Successful", f"Report saved to:\n{filepath}")
            self.status_label.configure(text=f"Report exported: {os.path.basename(filepath)}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"Error exporting report:\n{str(e)}")


# ============================================================
# Main Entry Point
# ============================================================
def main():
    root = tk.Tk()

    # Set window icon (if available)
    try:
        root.iconbitmap(default='')
    except:
        pass

    app = CpGenomeComparatorApp(root)

    # Center the window
    root.update_idletasks()
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    w = 1000
    h = 750
    x = (sw - w) // 2
    y = (sh - h) // 2
    root.geometry(f"{w}x{h}+{x}+{y}")

    root.mainloop()


if __name__ == "__main__":
    main()
