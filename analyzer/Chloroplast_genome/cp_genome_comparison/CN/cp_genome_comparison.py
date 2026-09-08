# -*- coding: utf-8 -*-
"""
叶绿体基因组组装结果比对工具
用于比较两个组装得到的叶绿体基因组序列与参考序列，判断哪个组装结果更正确。
支持检测叶绿体基因组四区结构（LSC, IRa, SSC, IRb）的互换/倒位/接头错误。

依赖: Python 3.x, Biopython, tkinter (内置)
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
# 序列分析引擎
# ============================================================
class SequenceAnalyzer:
    """序列分析核心类，执行所有比对和评分计算"""

    def __init__(self, progress_callback=None):
        self.progress_callback = progress_callback

    def _report_progress(self, msg, pct):
        if self.progress_callback:
            self.progress_callback(msg, pct)

    def read_fasta(self, filepath):
        """读取FASTA文件，合并所有记录为一条序列"""
        records = []
        # 先尝试UTF-8，失败则尝试其他编码
        for encoding in ['utf-8', 'utf-8-sig', 'gbk', 'gb18030', 'latin-1']:
            try:
                with open(filepath, 'r', encoding=encoding) as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                break
            except (UnicodeDecodeError, UnicodeError):
                continue
        if not records:
            raise ValueError(f"文件 {filepath} 中未找到有效序列，或序列格式不正确")
        if len(records) == 1:
            return str(records[0].seq).upper()
        combined = "".join(str(r.seq).upper() for r in records)
        return combined

    def get_basic_stats(self, seq, name=""):
        """计算序列基本统计信息"""
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
        基于k-mer的快速相似度计算。
        返回Jaccard相似系数和共享k-mer比例。
        这种方法对序列重排不敏感，适合基因组级别比较。
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
        滑动窗口比对：将参考序列分成窗口，用k-mer匹配找到组装序列中的最佳对应区域，
        计算每个窗口的近似一致性。
        返回窗口位置列表和对应的一致性值。
        """
        self._report_progress("正在执行滑动窗口分析...", 40)

        # 为组装序列构建k-mer位置索引
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
                self._report_progress(f"滑动窗口分析: {processed}/{total_windows}", pct)

        return windows, identities

    def _global_alignment(self, ref_seq, asm_seq):
        """对小序列进行全局双序列比对"""
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
        分块比对：将参考序列分块，每块在组装序列中找到最佳匹配区域做局部比对。
        适用于大序列（>50kb），避免全序列全局比对的内存溢出。
        """
        self._report_progress("正在执行分块比对（大基因组模式）...", 20)

        # 为组装序列构建k-mer索引
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
                self._report_progress(f"分块比对: {chunk_idx+1}/{num_chunks}", pct)

            ref_chunk = ref_seq[start:start + chunk_size]

            # 用k-mer找到组装序列中最佳匹配区域
            chunk_kmers = []
            for i in range(0, len(ref_chunk) - k + 1, k):
                kmer = ref_chunk[i:i + k]
                if 'N' not in kmer:
                    chunk_kmers.append(kmer)

            if not chunk_kmers:
                continue

            # 找到组装序列中最可能对应此块的区域
            positions = []
            for kmer in chunk_kmers:
                if kmer in asm_kmer_pos:
                    positions.extend(asm_kmer_pos[kmer])

            if not positions:
                continue

            # 取中位数位置，提取组装序列的对应窗口
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
        """智能比对：小序列用全局比对，大序列用分块比对"""
        max_len = max(len(ref_seq), len(asm_seq))
        if max_len <= 50000:
            self._report_progress("正在执行全局比对...", 20)
            return self._global_alignment(ref_seq, asm_seq)
        else:
            return self._chunked_alignment(ref_seq, asm_seq)

    def compute_coverage(self, ref_seq, asm_seq, k=21):
        """通过k-mer匹配估算组装序列对参考序列的覆盖度"""
        self._report_progress("正在计算覆盖度...", 60)

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
        执行完整分析流程：
        1. 读取序列
        2. 基本统计
        3. k-mer相似度
        4. 滑动窗口分析
        5. Biopython全局比对
        6. 覆盖度计算
        7. 综合评分
        """
        results = {}

        # 步骤1：读取序列
        self._report_progress("正在读取参考序列...", 5)
        ref_seq = self.read_fasta(ref_path)
        ref_stats = self.get_basic_stats(ref_seq, "参考序列")
        results['reference'] = ref_stats

        self._report_progress("正在读取组装序列1...", 8)
        asm1_seq = self.read_fasta(asm1_path)
        asm1_stats = self.get_basic_stats(asm1_seq, "组装序列1")
        results['assembly1'] = asm1_stats

        self._report_progress("正在读取组装序列2...", 11)
        asm2_seq = self.read_fasta(asm2_path)
        asm2_stats = self.get_basic_stats(asm2_seq, "组装序列2")
        results['assembly2'] = asm2_stats

        # 步骤2：k-mer相似度分析
        self._report_progress("正在计算k-mer相似度...", 14)
        results['kmer_asm1'] = self.kmer_similarity(ref_seq, asm1_seq, k=31)
        results['kmer_asm2'] = self.kmer_similarity(ref_seq, asm2_seq, k=31)

        # 步骤3：Biopython全局比对
        if BIOPYTHON_AVAILABLE:
            results['alignment_asm1'] = self.biopython_alignment(ref_seq, asm1_seq)
            self._report_progress("组装序列1比对完成，正进行组装序列2比对...", 30)
            results['alignment_asm2'] = self.biopython_alignment(ref_seq, asm2_seq)
        else:
            results['alignment_asm1'] = None
            results['alignment_asm2'] = None

        # 步骤4：覆盖度计算
        results['coverage_asm1'] = self.compute_coverage(ref_seq, asm1_seq)
        results['coverage_asm2'] = self.compute_coverage(ref_seq, asm2_seq)

        # 步骤5：滑动窗口分析
        self._report_progress("正在执行组装序列1的滑动窗口分析...", 70)
        windows1, identities1 = self.sliding_window_identity(ref_seq, asm1_seq)

        self._report_progress("正在执行组装序列2的滑动窗口分析...", 85)
        windows2, identities2 = self.sliding_window_identity(ref_seq, asm2_seq)

        results['window_asm1'] = {'positions': windows1, 'identities': identities1}
        results['window_asm2'] = {'positions': windows2, 'identities': identities2}

        # 步骤6：综合评分（先做基础评分）
        self._report_progress("正在计算基础综合评分...", 95)
        results['score1_base'] = self._calculate_score(results, asm=1)
        results['score2_base'] = self._calculate_score(results, asm=2)

        # 步骤7：叶绿体四区结构分析（检测互换/倒位/接头错误）
        self._report_progress("正在进行四区结构分析...", 16)
        if len(ref_seq) > 30000:  # 仅对大于30kb的序列做结构分析
            results['structure1'] = self.run_structure_analysis(
                ref_seq, asm1_seq, asm_name="组装序列1")
            results['structure2'] = self.run_structure_analysis(
                ref_seq, asm2_seq, asm_name="组装序列2")
        else:
            results['structure1'] = {'structure_score': 5.0, 'structure_details': '序列太短，跳过结构分析'}
            results['structure2'] = {'structure_score': 5.0, 'structure_details': '序列太短，跳过结构分析'}

        # 最终综合评分 = 基础分 + 结构分
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

        self._report_progress("分析完成！", 100)
        return results

    def _calculate_score(self, results, asm):
        """计算综合评分（0-100分）"""
        prefix = f'assembly{asm}'
        kmer_key = f'kmer_asm{asm}'
        align_key = f'alignment_asm{asm}'
        cov_key = f'coverage_asm{asm}'

        ref_len = results['reference']['length']
        asm_len = results[prefix]['length']
        ref_gc = results['reference']['gc_content']
        asm_gc = results[prefix]['gc_content']

        # 长度相似度得分 (15%)
        length_ratio = min(ref_len, asm_len) / max(ref_len, asm_len)
        length_score = length_ratio * 15

        # GC含量相似度得分 (15%)
        gc_diff = abs(ref_gc - asm_gc)
        gc_score = max(0, 15 - gc_diff * 1.5)

        # 覆盖度得分 (30%)
        coverage = results[cov_key]
        coverage_score = coverage / 100 * 30

        # 一致性得分 (40%)
        if BIOPYTHON_AVAILABLE and results[align_key] is not None:
            identity = results[align_key]['identity']
            identity_score = identity / 100 * 40
        else:
            # 使用k-mer相似度替代
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
    # 叶绿体基因组四区结构分析（检测互换/倒位/接头错误）
    # ============================================================


    def _revcomp(self, seq):
        """返回序列的反向互补"""
        trans = str.maketrans('ACGTacgt', 'TGCAtgca')
        return seq[::-1].translate(trans)

    def _segment_map(self, ref_seq, asm_seq, num_segments=20, k=31):
        """
        将参考序列等分为多个片段，每个片段通过k-mer匹配映射到组装序列。
        返回: [(ref_start, ref_end, asm_best_pos, orientation, identity), ...]
          orientation: 'forward' 或 'reversed'
        """
        self._report_progress("正在进行参考片段位置映射...", 17)
        seq_len = len(ref_seq)
        segment_size = seq_len // num_segments
        if segment_size < 1000:
            segment_size = 1000
            num_segments = max(1, seq_len // segment_size)

        # 构建组装序列的k-mer索引
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

            # 正向匹配
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

            # 反向互补匹配
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
        分析片段映射结果，检测结构性错误。
        返回: {
            'is_linear': bool,
            'is_global_inversion': bool,  # 整个序列是否是反向互补（全局方向）
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
                'issues': ['无法获得足够的有效映射，可能序列差异较大'],
            }

        issues = []
        inverted_segments = []
        position_breaks = []

        # 检查方向
        forward_count = sum(1 for m in valid if m[3] == 'forward')
        reversed_count = sum(1 for m in valid if m[3] == 'reversed')
        total_valid = len(valid)

        # 判断是否为全局方向翻转（几乎所有片段都反向）
        is_global_inversion = (reversed_count >= total_valid * 0.85)

        if is_global_inversion:
            # 全局翻转：组装序列整体是参考的反向互补，不视为结构性错误
            # 但仍检查位置线性顺序
            issues.append("检测到全局方向翻转（组装序列整体为参考的反向互补）")
        elif reversed_count > total_valid * 0.3:
            # 部分反向：存在结构性倒位
            issues.append("部分片段以反向互补方式匹配（存在局部倒位）")
            # 只标记非全局翻转的倒位片段
            for i, m in enumerate(mappings):
                if m[2] >= 0 and m[3] == 'reversed':
                    inverted_segments.append(i)

        # 检查局部倒位：在非全局翻转的情况下，找出连续的局部反向片段
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

        # 检查位置线性顺序
        # 全局翻转时，期望位置递减；否则期望递增
        positions = [m[2] for m in valid]
        breaks = 0
        if is_global_inversion:
            # 全局翻转：期望递减
            for i in range(1, len(positions)):
                if positions[i] > positions[i - 1]:  # 递增=跳变
                    breaks += 1
                    position_breaks.append(i)
        else:
            for i in range(1, len(positions)):
                if positions[i] < positions[i - 1]:  # 递减=跳变
                    breaks += 1
                    position_breaks.append(i)

        is_linear = (breaks == 0)
        if breaks > 0:
            if breaks <= 2:
                issues.append(f"发现 {breaks} 处位置跳变，可能存在局部重排")
            else:
                issues.append(f"发现 {breaks} 处位置跳变，可能存在严重重排或IR互换")

        # 检测IR互换的特征：位置出现"折返"
        if not is_global_inversion and len(positions) >= 6:
            mid = len(positions) // 2
            first_half_inc = all(positions[i] <= positions[i + 1] for i in range(mid - 1))
            if first_half_inc:
                second_start = positions[mid]
                if any(p < second_start for p in positions[mid + 1:]):
                    issues.append("检测到可能的位置折返模式（IR互换特征）")

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
        在参考序列上等距采集接头k-mer指纹，检查它们在组装中的出现顺序。
        不依赖IR边界的先验知识。
        """
        seq_len = len(ref_seq)
        # 在参考序列上采样8个等距"检查点"
        num_checks = 8
        check_positions = {}
        for i in range(num_checks):
            pos = (seq_len * (i + 1)) // (num_checks + 1)
            check_positions[f'checkpoint_{i+1}'] = pos

        # 提取每个检查点的k-mer集
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

        # 在组装中寻找每个检查点k-mer集的最佳位置
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

        # 检查顺序
        ordered_names = sorted(check_positions.keys(), 
                              key=lambda n: check_positions[n])
        asm_order = [(n, asm_positions.get(n, -1)) for n in ordered_names]
        valid_order = [x for x in asm_order if x[1] >= 0]

        issues = []
        if len(valid_order) >= 3:
            for i in range(1, len(valid_order)):
                if valid_order[i][1] < valid_order[i - 1][1]:
                    issues.append(
                        f"检查点 {valid_order[i][0]} 在组装中的位置 ({valid_order[i][1]:,}) "
                        f"出现在 {valid_order[i-1][0]} ({valid_order[i-1][1]:,}) 之前"
                    )

        # 也检测局部倒位：对每个检查点，检查在组装中是正向还是反向匹配
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

            # 正向匹配率
            fwd = sum(1 for i in range(0, min(len(ref_region), len(asm_region)) - k + 1, 50)
                      if ref_region[i:i + k] == asm_region[i:i + k]) if min(len(ref_region), len(asm_region)) >= k else 0
            # 反向匹配率
            rc_ref = ref_region[::-1].translate(complement)
            rev = sum(1 for i in range(0, min(len(rc_ref), len(asm_region)) - k + 1, 50)
                      if rc_ref[i:i + k] == asm_region[i:i + k]) if min(len(rc_ref), len(asm_region)) >= k else 0

            if rev > fwd * 2:
                inversion_flags.append(f"检查点 {name} 附近存在局部倒位信号")

        return {
            'check_positions': asm_positions,
            'order_issues': issues,
            'inversion_flags': inversion_flags,
        }

    def run_structure_analysis(self, ref_seq, asm_seq, asm_name=""):
        """
        执行完整的结构分析（基于参考片段映射法）。
        不依赖IR边界检测，通用检测：
        1. 序列重排（互换/倒位）
        2. 局部倒位
        3. 接头/顺序错误
        
        返回结构分析结果字典。
        """
        self._report_progress(f"正在分析{asm_name}的结构完整性...", 16)
        
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
            result['structure_details'] = '序列太短，跳过结构分析'
            return result

        # 步骤1：参考片段位置映射
        mappings = self._segment_map(ref_seq, asm_seq, num_segments=20)
        result['_mappings'] = mappings

        # 步骤2：分析映射结构
        self._report_progress(f"正在分析{asm_name}的片段映射结构...", 18)
        map_analysis = self._analyze_mapping_structure(mappings)
        result['is_linear'] = map_analysis['is_linear']
        result['has_rearrangement'] = not map_analysis['is_linear']
        result['is_global_inversion'] = map_analysis.get('is_global_inversion', False)
        result['inverted_segments'] = map_analysis['inverted_segments']
        
        # 将倒位片段转换为位置范围的large_inversions
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

        # 步骤3：检查点k-mer顺序验证
        self._report_progress(f"正在验证{asm_name}的序列顺序...", 19)
        junction_check = self._check_junction_kmers(ref_seq, asm_seq)
        result['order_issues'] = junction_check['order_issues']
        result['inversion_flags'] = junction_check['inversion_flags']
        if junction_check['check_positions']:
            result['valid_junctions'] = len(junction_check['check_positions'])

        # 步骤4：综合判断是否可能发生IR互换
        # IR互换典型特征：多片段位置折返 + 多处倒位
        # 全局翻转不应被视为IR互换
        if not map_analysis.get('is_global_inversion', False):
            if (not map_analysis['is_linear'] and 
                map_analysis.get('num_breaks', 0) >= 2 and
                len(map_analysis['inverted_segments']) >= 2):
                result['ir_swapped'] = True
            if junction_check['order_issues'] and map_analysis.get('num_breaks', 0) >= 2:
                result['ir_swapped'] = True

        # 步骤5：计算结构得分
        result['structure_score'] = self._calculate_structure_score(result)
        result['structure_details'] = self._format_structure_details(result, mappings)

        return result

    def _calculate_structure_score(self, struct_result):
        """
        计算结构完整性得分（满分5分）。
        区分全局翻转（仅方向不同）和真正的结构性错误。
        """
        score = 5.0

        # 全局翻转：组装整体是参考的反向互补，不是结构性错误
        is_global_inv = struct_result.get('is_global_inversion', False)

        # 重排扣分
        if struct_result.get('has_rearrangement', False):
            if is_global_inv:
                score -= 0.5  # 全局翻转+重排：轻微扣分
            else:
                score -= 1.5  # 非全局翻转+重排：正常扣分

        # 倒位扣分（全局翻转不计入）
        if not is_global_inv:
            inv_count = len(struct_result.get('inverted_segments', []))
            score -= min(inv_count * 0.3, 2.0)  # 最多扣2分

        # 顺序问题扣分
        order_issues = struct_result.get('order_issues', [])
        score -= len(order_issues) * 0.8

        # IR互换特征扣分
        if struct_result.get('ir_swapped', False) and not is_global_inv:
            score -= 2.0

        return max(0, round(score, 2))

    def _format_structure_details(self, struct_result, mappings=None):
        """格式化结构分析详情为可读文本"""
        lines = []

        if mappings:
            # 显示映射概览
            lines.append("参考片段映射结果 (20个等距片段):")
            seg_lines = []
            for i, (s, e, pos, orient, ident) in enumerate(mappings):
                mark = ""
                if pos < 0:
                    mark = " [未定位]"
                elif orient == 'reversed':
                    mark = " [倒位]"
                seg_lines.append(f"  片段{i+1:>2}: 参考[{s:>7}-{e:>7}] "
                               f"-> 组装约{pos:>7,}bp  {orient}{mark}")
            lines.extend(seg_lines[:20])  # 最多显示20个

        # 结构状态摘要
        lines.append(f"\n结构完整性分析:")
        if struct_result.get('is_global_inversion', False):
            lines.append("  ℹ 全局方向翻转（组装整体为参考反向互补，非结构性错误）")
        if struct_result.get('is_linear', True):
            lines.append("  ✓ 片段映射保持线性顺序")
        else:
            lines.append("  ⚠ 片段映射存在位置跳变，可能发生重排")

        inv_segs = struct_result.get('inverted_segments', [])
        if inv_segs and not struct_result.get('is_global_inversion', False):
            lines.append(f"  ⚠ 检测到 {len(inv_segs)} 个局部倒位片段: {sorted(inv_segs)[:10]}")
        elif struct_result.get('is_global_inversion', False):
            lines.append("  ℹ 全局方向翻转（非局部倒位）")
        else:
            lines.append("  ✓ 未检测到局部倒位")

        order_issues = struct_result.get('order_issues', [])
        if order_issues:
            lines.append(f"  ⚠ 顺序问题 ({len(order_issues)}项):")
            for oi in order_issues[:5]:
                lines.append(f"    - {oi}")
        else:
            lines.append("  ✓ 检查点k-mer顺序正常")

        inversion_flags = struct_result.get('inversion_flags', [])
        if inversion_flags:
            for f in inversion_flags[:3]:
                lines.append(f"  ⚠ {f}")

        if struct_result.get('ir_swapped', False) and not struct_result.get('is_global_inversion', False):
            lines.append("  ⚠ 检测到IR互换特征（多片段折返+局部倒位）")
        elif struct_result.get('is_global_inversion', False):
            lines.append("  ℹ 全局翻转模式，IR互换检测不适用")
        else:
            lines.append("  ✓ 未检测到IR互换特征")

        large_inv = struct_result.get('large_inversions', [])
        if large_inv:
            lines.append(f"  大片段倒位区域 ({len(large_inv)}处):")
            for s, e in large_inv[:5]:
                lines.append(f"    [{s:,}-{e:,}] (~{(e-s)//1000}kb)")

        lines.append(f"\n结构完整性得分: {struct_result.get('structure_score', 0):.2f} / 5.0")

        return '\n'.join(lines)
        return '\n'.join(lines)


# ============================================================
# 图形用户界面
# ============================================================
class CpGenomeComparatorApp:
    """叶绿体基因组比对工具GUI"""

    def __init__(self, root):
        self.root = root
        self.root.title("叶绿体基因组组装结果比对工具")
        self.root.geometry("1000x750")
        self.root.minsize(900, 600)

        # 设置样式
        self._setup_style()
        # 创建界面
        self._create_widgets()
        # 分析状态
        self.analysis_running = False
        self.results = None

    def _setup_style(self):
        style = ttk.Style()
        style.theme_use('clam')

        # 自定义颜色方案
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
        # 主容器
        main_frame = tk.Frame(self.root, bg=self.colors['bg'])
        main_frame.pack(fill=tk.BOTH, expand=True, padx=15, pady=10)

        # ----- 标题 -----
        title_frame = tk.Frame(main_frame, bg=self.colors['bg'])
        title_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(title_frame, text="🧬 叶绿体基因组组装结果比对工具",
                  style='Title.TLabel').pack(side=tk.LEFT)
        ttk.Label(title_frame,
                  text="比较两个组装序列与参考序列，自动判断哪个组装结果更正确",
                  font=('Microsoft YaHei', 9),
                  foreground='#636e72',
                  background=self.colors['bg']).pack(side=tk.LEFT, padx=(15, 0))

        # ----- 文件选择区 -----
        file_card = ttk.Frame(main_frame, style='Card.TFrame')
        file_card.pack(fill=tk.X, pady=(0, 10), ipady=8)

        ttk.Label(file_card, text="📂 选择序列文件（FASTA格式）",
                  style='Section.TLabel').pack(anchor=tk.W, padx=12, pady=(8, 5))

        # 参考序列
        self._create_file_row(file_card, "参考序列 (Reference):", 'ref', 0)
        # 组装序列1
        self._create_file_row(file_card, "组装序列1 (Assembly 1):", 'asm1', 1)
        # 组装序列2
        self._create_file_row(file_card, "组装序列2 (Assembly 2):", 'asm2', 2)

        # ----- 控制按钮 -----
        btn_frame = tk.Frame(main_frame, bg=self.colors['bg'])
        btn_frame.pack(fill=tk.X, pady=(0, 10))

        self.analyze_btn = tk.Button(
            btn_frame, text="▶  开始分析",
            font=('Microsoft YaHei', 12, 'bold'),
            bg=self.colors['accent'], fg='white',
            activebackground='#076bb5', activeforeground='white',
            relief=tk.FLAT, padx=25, pady=6,
            cursor='hand2', borderwidth=0,
            command=self._start_analysis
        )
        self.analyze_btn.pack(side=tk.LEFT, padx=(0, 10))

        self.clear_btn = tk.Button(
            btn_frame, text="🗑 清除结果",
            font=('Microsoft YaHei', 10),
            bg=self.colors['light'], fg=self.colors['fg'],
            activebackground='#c8d6e5', activeforeground=self.colors['fg'],
            relief=tk.FLAT, padx=18, pady=5,
            cursor='hand2', borderwidth=0,
            command=self._clear_results
        )
        self.clear_btn.pack(side=tk.LEFT, padx=(0, 10))

        self.export_btn = tk.Button(
            btn_frame, text="📋 导出报告",
            font=('Microsoft YaHei', 10),
            bg=self.colors['light'], fg=self.colors['fg'],
            activebackground='#c8d6e5', activeforeground=self.colors['fg'],
            relief=tk.FLAT, padx=18, pady=5,
            cursor='hand2', borderwidth=0,
            command=self._export_report,
            state=tk.DISABLED
        )
        self.export_btn.pack(side=tk.LEFT)

        # 进度条
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

        # ----- 结果展示区（使用Notebook分页） -----
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # == 第1页：综合结果 ==
        self.summary_frame = ttk.Frame(self.notebook, style='Card.TFrame')
        self.notebook.add(self.summary_frame, text="📊 综合结果")

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
        self.summary_text.insert(tk.END, "请选择三个FASTA文件后，点击「开始分析」按钮。\n\n")
        self.summary_text.insert(tk.END, "支持的文件格式：.fasta, .fa, .fna, .txt\n")
        self.summary_text.insert(tk.END, "参考序列和组装序列各选择一个文件。\n")
        self.summary_text.configure(state=tk.DISABLED)

        # == 第2页：详细比对 ==
        self.detail_frame = ttk.Frame(self.notebook, style='Card.TFrame')
        self.notebook.add(self.detail_frame, text="📋 详细比对")

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
        self.detail_text.insert(tk.END, "详细比对数据将在分析完成后显示。\n")
        self.detail_text.configure(state=tk.DISABLED)

        # == 第3页：窗口一致性图 ==
        self.canvas_frame = ttk.Frame(self.notebook, style='Card.TFrame')
        self.notebook.add(self.canvas_frame, text="📈 一致性分布")

        canvas_container = tk.Frame(self.canvas_frame, bg='white')
        canvas_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.canvas = tk.Canvas(
            canvas_container,
            bg='white',
            highlightthickness=0,
            height=300,
        )
        self.canvas.pack(fill=tk.BOTH, expand=True)

        # 底部状态栏
        status_frame = tk.Frame(main_frame, bg=self.colors['bg'])
        status_frame.pack(fill=tk.X, pady=(5, 0))
        self.status_label = ttk.Label(
            status_frame, text="就绪",
            font=('Microsoft YaHei', 9),
            foreground='#636e72',
            background=self.colors['bg']
        )
        self.status_label.pack(side=tk.LEFT)

    def _create_file_row(self, parent, label_text, key, row):
        """创建文件选择行"""
        row_frame = tk.Frame(parent, bg=self.colors['card'])
        row_frame.pack(fill=tk.X, padx=12, pady=3)

        ttk.Label(row_frame, text=label_text,
                  font=('Microsoft YaHei', 10),
                  background=self.colors['card']).pack(side=tk.LEFT, padx=(0, 8))

        entry = ttk.Entry(row_frame, font=('Microsoft YaHei', 10), width=55)
        entry.pack(side=tk.LEFT, padx=(0, 5))

        def browse():
            filepath = filedialog.askopenfilename(
                title=f"选择{label_text.split('(')[0].strip()}文件",
                filetypes=[
                    ("FASTA文件", "*.fasta *.fa *.fna *.fas *.txt"),
                    ("所有文件", "*.*"),
                ]
            )
            if filepath:
                entry.delete(0, tk.END)
                entry.insert(0, filepath)

        browse_btn = tk.Button(
            row_frame, text="浏览...",
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
        """启动分析（在线程中运行）"""
        if self.analysis_running:
            messagebox.showwarning("提示", "分析正在进行中，请等待完成。")
            return

        # 获取文件路径
        ref_path = self.ref_entry.get().strip()
        asm1_path = self.asm1_entry.get().strip()
        asm2_path = self.asm2_entry.get().strip()

        # 验证
        errors = []
        if not ref_path or not os.path.exists(ref_path):
            errors.append("请选择有效的参考序列文件。")
        if not asm1_path or not os.path.exists(asm1_path):
            errors.append("请选择有效的组装序列1文件。")
        if not asm2_path or not os.path.exists(asm2_path):
            errors.append("请选择有效的组装序列2文件。")

        if errors:
            messagebox.showerror("输入错误", "\n".join(errors))
            return

        if not BIOPYTHON_AVAILABLE:
            messagebox.showerror(
                "缺少依赖",
                "未检测到Biopython库。\n\n"
                "请运行以下命令安装：\n"
                "  pip install biopython\n\n"
                "安装后重新启动程序即可。"
            )
            return

        # 禁用按钮
        self.analysis_running = True
        self.analyze_btn.configure(
            text="⏳ 分析中...",
            state=tk.DISABLED,
            bg='#636e72'
        )
        self.export_btn.configure(state=tk.DISABLED)
        self.progress_var.set(0)
        self.progress_label.configure(text="初始化...")
        self.status_label.configure(text="正在分析...")

        # 在线程中运行分析
        thread = threading.Thread(
            target=self._run_analysis_thread,
            args=(ref_path, asm1_path, asm2_path),
            daemon=True
        )
        thread.start()

    def _run_analysis_thread(self, ref_path, asm1_path, asm2_path):
        """后台分析线程"""
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
        """更新进度条"""
        self.progress_var.set(pct)
        self.progress_label.configure(text=msg)

    def _handle_error(self, error_msg):
        """处理分析错误"""
        self.analysis_running = False
        self.analyze_btn.configure(
            text="▶  开始分析",
            state=tk.NORMAL,
            bg=self.colors['accent']
        )
        self.progress_label.configure(text="分析出错")
        self.status_label.configure(text="出错")
        messagebox.showerror("分析错误", f"分析过程中出现错误：\n\n{error_msg}")

    def _display_results(self, results):
        """显示分析结果"""
        self.analysis_running = False
        self.analyze_btn.configure(
            text="▶  开始分析",
            state=tk.NORMAL,
            bg=self.colors['accent']
        )
        self.export_btn.configure(state=tk.NORMAL)
        self.progress_label.configure(text="完成")
        self.status_label.configure(text="分析完成")

        score1 = results['score1']['total']
        score2 = results['score2']['total']

        if score1 > score2:
            winner = "组装序列1"
            winner_score = score1
            loser = "组装序列2"
            loser_score = score2
        elif score2 > score1:
            winner = "组装序列2"
            winner_score = score2
            loser = "组装序列1"
            loser_score = score1
        else:
            winner = "两个组装结果得分相同"
            winner_score = score1
            loser = ""
            loser_score = score2

        # 第1页：综合结果
        self._display_summary(results, winner, winner_score, loser, loser_score)
        # 第2页：详细比对
        self._display_details(results)
        # 第3页：图表
        self._draw_window_plot(results)

        # 切换到综合结果页
        self.notebook.select(0)

    def _display_summary(self, results, winner, winner_score, loser, loser_score):
        """显示综合结果"""
        self.summary_text.configure(state=tk.NORMAL)
        self.summary_text.delete(1.0, tk.END)

        # 为文本添加标签样式
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
        self.summary_text.insert(tk.END, "  叶绿体基因组组装比对 - 综合结果\n", "title")
        self.summary_text.insert(tk.END, "═" * 60 + "\n\n", "header")

        # 胜出结果
        if "两者" not in winner:
            self.summary_text.insert(tk.END,
                f"  ★★★  {winner} 更接近参考序列 ★★★\n\n", "winner")
        else:
            self.summary_text.insert(tk.END,
                f"  {winner}\n\n", "header")

        # 评分表格
        self.summary_text.insert(tk.END, "┌" + "─" * 56 + "┐\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^15} │ {:^15} │\n".format("评估维度", "组装序列1", "组装序列2"))
        self.summary_text.insert(tk.END, "├" + "─" * 56 + "┤\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}分 │ {:^13.2f}分 │\n".format("一致性得分(40%)", s1['identity_score'], s2['identity_score']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}分 │ {:^13.2f}分 │\n".format("覆盖度得分(30%)", s1['coverage_score'], s2['coverage_score']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}分 │ {:^13.2f}分 │\n".format("长度相似度(15%)", s1['length_score'], s2['length_score']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}分 │ {:^13.2f}分 │\n".format("GC相似度(15%)", s1['gc_score'], s2['gc_score']))
        # 新增结构完整性得分
        if 'structure_score' in s1:
            self.summary_text.insert(tk.END,
                "│ {:^20} │ {:^13.2f}分 │ {:^13.2f}分 │\n".format("结构完整性(+5%)", s1['structure_score'], s2['structure_score']))
        self.summary_text.insert(tk.END, "├" + "─" * 56 + "┤\n")

        tag1 = "good" if s1['total'] > s2['total'] else ""
        tag2 = "good" if s2['total'] > s1['total'] else ""

        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}分 │ {:^13.2f}分 │\n".format("综合得分 ★", s1['total'], s2['total']))
        self.summary_text.insert(tk.END, "└" + "─" * 56 + "┘\n\n")

        # 基本统计对比
        self.summary_text.insert(tk.END, "基本统计信息：\n", "header")
        self.summary_text.insert(tk.END,
            "┌" + "─" * 56 + "┐\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^15} │ {:^15} │\n".format("统计项", "参考序列", "组装序列1", "组装序列2"))
        self.summary_text.insert(tk.END, "├" + "─" * 56 + "┤\n")
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13,}bp │ {:^13,}bp │\n".format("序列长度", ref['length'], asm1['length']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13} │ {:^13,}bp │\n".format("", "", asm2['length']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13.2f}%  │ {:^13.2f}%  │\n".format("GC含量", ref['gc_content'], asm1['gc_content']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13} │ {:^13.2f}%  │\n".format("", "", asm2['gc_content']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13}  │ {:^13}  │\n".format("N碱基数", ref['N_count'], asm1['N_count']))
        self.summary_text.insert(tk.END,
            "│ {:^20} │ {:^13}  │ {:^13}  │\n".format("", "", asm2['N_count']))
        self.summary_text.insert(tk.END, "└" + "─" * 56 + "┘\n\n")

        # k-mer分析摘要
        km1 = results['kmer_asm1']
        km2 = results['kmer_asm2']
        self.summary_text.insert(tk.END, "K-mer相似度分析 (k=31)：\n", "header")
        self.summary_text.insert(tk.END,
            f"  组装序列1: Jaccard={km1['jaccard']:.4f},  共享比例={km1['shared_ratio']:.4f}\n")
        self.summary_text.insert(tk.END,
            f"  组装序列2: Jaccard={km2['jaccard']:.4f},  共享比例={km2['shared_ratio']:.4f}\n\n")

        # 覆盖度
        self.summary_text.insert(tk.END, "参考序列覆盖度：\n", "header")
        self.summary_text.insert(tk.END,
            f"  组装序列1: {results['coverage_asm1']:.2f}%\n")
        self.summary_text.insert(tk.END,
            f"  组装序列2: {results['coverage_asm2']:.2f}%\n\n")

        # 四区结构分析摘要
        for st_key, st_label in [('structure1', '组装序列1'), ('structure2', '组装序列2')]:
            st = results.get(st_key)
            if not st or not isinstance(st, dict):
                continue
            if st.get('ir_regions_found'):
                self.summary_text.insert(tk.END, f"{st_label} 四区结构分析：\n", "header")
                ref_reg = st.get('ref_regions', {})
                if ref_reg:
                    self.summary_text.insert(tk.END,
                        f"  参考四区: LSC={ref_reg.get('lsc_length','?'):,}bp, "
                        f"IRa={ref_reg.get('ira_length','?'):,}bp, "
                        f"SSC={ref_reg.get('ssc_length','?'):,}bp, "
                        f"IRb={ref_reg.get('irb_length','?'):,}bp\n")
                self.summary_text.insert(tk.END,
                    f"  IR方向: IRa={st.get('ira_orientation','?')}, "
                    f"IRb={st.get('irb_orientation','?')}\n")
                issues = []
                if st.get('ir_swapped'):
                    issues.append("IRa/IRb互换")
                if st.get('ir_inverted'):
                    issues.append("IR倒位")
                inversions = st.get('large_inversions', [])
                if inversions:
                    issues.append(f"大片段倒位({len(inversions)}处)")
                for ji in st.get('junction_issues', []):
                    issues.append(ji.replace('⚠ ', ''))
                if issues:
                    self.summary_text.insert(tk.END,
                        f"  ⚠ 结构问题: {'; '.join(issues)}\n")
                else:
                    self.summary_text.insert(tk.END,
                        "  ✓ 四区结构正常，未检测到互换/倒位/接头错误\n")
                self.summary_text.insert(tk.END,
                    f"  结构得分: {st.get('structure_score', 0):.2f} / 5.0\n\n")

        # 对齐结果
        if BIOPYTHON_AVAILABLE:
            al1 = results.get('alignment_asm1')
            al2 = results.get('alignment_asm2')
            if al1 and al2:
                self.summary_text.insert(tk.END, "全局比对结果：\n", "header")
                self.summary_text.insert(tk.END,
                    f"  组装序列1: 一致性={al1['identity']:.2f}%,  对齐长度={al1['aligned_length']:,}\n")
                self.summary_text.insert(tk.END,
                    f"  组装序列2: 一致性={al2['identity']:.2f}%,  对齐长度={al2['aligned_length']:,}\n\n")

        # 结论
        self.summary_text.insert(tk.END, "═" * 60 + "\n", "header")
        self.summary_text.insert(tk.END,
            "结论：建议使用得分更高的组装结果进行后续分析。\n", "header")
        self.summary_text.insert(tk.END,
            "若两个得分接近，建议检查滑动窗口一致性分布图表，\n")
        self.summary_text.insert(tk.END,
            "关注序列是否存在大片段重排或倒位。\n")

        self.summary_text.configure(state=tk.DISABLED)

    def _display_details(self, results):
        """显示详细比对信息"""
        self.detail_text.configure(state=tk.NORMAL)
        self.detail_text.delete(1.0, tk.END)

        ref = results['reference']
        asm1 = results['assembly1']
        asm2 = results['assembly2']

        self.detail_text.insert(tk.END, "=" * 70 + "\n")
        self.detail_text.insert(tk.END, "  详细比对报告\n")
        self.detail_text.insert(tk.END, "=" * 70 + "\n\n")

        # 序列基本信息
        self.detail_text.insert(tk.END, "【序列基本信息】\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  参考序列:   长度={ref['length']:,} bp,  GC={ref['gc_content']:.2f}%,  "
            f"A={ref['A_pct']:.1f}% T={ref['T_pct']:.1f}% C={ref['C_pct']:.1f}% G={ref['G_pct']:.1f}%\n")
        self.detail_text.insert(tk.END,
            f"  组装序列1:  长度={asm1['length']:,} bp,  GC={asm1['gc_content']:.2f}%,  "
            f"A={asm1['A_pct']:.1f}% T={asm1['T_pct']:.1f}% C={asm1['C_pct']:.1f}% G={asm1['G_pct']:.1f}%\n")
        self.detail_text.insert(tk.END,
            f"  组装序列2:  长度={asm2['length']:,} bp,  GC={asm2['gc_content']:.2f}%,  "
            f"A={asm2['A_pct']:.1f}% T={asm2['T_pct']:.1f}% C={asm2['C_pct']:.1f}% G={asm2['G_pct']:.1f}%\n\n")

        # 长度差异
        len_diff1 = abs(ref['length'] - asm1['length'])
        len_diff_pct1 = len_diff1 / ref['length'] * 100 if ref['length'] > 0 else 0
        len_diff2 = abs(ref['length'] - asm2['length'])
        len_diff_pct2 = len_diff2 / ref['length'] * 100 if ref['length'] > 0 else 0
        self.detail_text.insert(tk.END, "【长度差异】\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  组装序列1 vs 参考: {len_diff1:,} bp ({len_diff_pct1:.2f}%)\n")
        self.detail_text.insert(tk.END,
            f"  组装序列2 vs 参考: {len_diff2:,} bp ({len_diff_pct2:.2f}%)\n\n")

        # GC含量差异
        gc_diff1 = abs(ref['gc_content'] - asm1['gc_content'])
        gc_diff2 = abs(ref['gc_content'] - asm2['gc_content'])
        self.detail_text.insert(tk.END, "【GC含量差异】\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  组装序列1 vs 参考: {gc_diff1:.2f}%\n")
        self.detail_text.insert(tk.END,
            f"  组装序列2 vs 参考: {gc_diff2:.2f}%\n\n")

        # K-mer分析
        km1 = results['kmer_asm1']
        km2 = results['kmer_asm2']
        self.detail_text.insert(tk.END, "【K-mer相似度分析 (k=31)】\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  组装序列1:\n"
            f"    参考k-mer总数: {km1['total_kmers_ref']:,}\n"
            f"    组装k-mer总数: {km1['total_kmers_asm']:,}\n"
            f"    共享k-mer数:   {km1['shared_kmers']:,}\n"
            f"    Jaccard系数:   {km1['jaccard']:.6f}\n"
            f"    共享比例:      {km1['shared_ratio']:.6f}\n\n")
        self.detail_text.insert(tk.END,
            f"  组装序列2:\n"
            f"    参考k-mer总数: {km2['total_kmers_ref']:,}\n"
            f"    组装k-mer总数: {km2['total_kmers_asm']:,}\n"
            f"    共享k-mer数:   {km2['shared_kmers']:,}\n"
            f"    Jaccard系数:   {km2['jaccard']:.6f}\n"
            f"    共享比例:      {km2['shared_ratio']:.6f}\n\n")

        # 比对详情
        if BIOPYTHON_AVAILABLE:
            al1 = results.get('alignment_asm1')
            al2 = results.get('alignment_asm2')
            if al1 and al2:
                self.detail_text.insert(tk.END, "【全局比对详情 (Biopython)】\n")
                self.detail_text.insert(tk.END, "-" * 70 + "\n")
                self.detail_text.insert(tk.END,
                    f"  组装序列1:\n"
                    f"    一致性:       {al1['identity']:.2f}%\n"
                    f"    匹配碱基数:   {al1['matches']:,}\n"
                    f"    错配碱基数:   {al1['mismatches']:,}\n"
                    f"    参考序列gap:  {al1['gaps_ref']:,}\n"
                    f"    组装序列gap:  {al1['gaps_asm']:,}\n"
                    f"    对齐总长度:   {al1['aligned_length']:,}\n"
                    f"    比对得分:     {al1['score']:.0f}\n\n")
                self.detail_text.insert(tk.END,
                    f"  组装序列2:\n"
                    f"    一致性:       {al2['identity']:.2f}%\n"
                    f"    匹配碱基数:   {al2['matches']:,}\n"
                    f"    错配碱基数:   {al2['mismatches']:,}\n"
                    f"    参考序列gap:  {al2['gaps_ref']:,}\n"
                    f"    组装序列gap:  {al2['gaps_asm']:,}\n"
                    f"    对齐总长度:   {al2['aligned_length']:,}\n"
                    f"    比对得分:     {al2['score']:.0f}\n\n")

        # 覆盖度
        self.detail_text.insert(tk.END, "【参考序列覆盖度】\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        self.detail_text.insert(tk.END,
            f"  组装序列1: {results['coverage_asm1']:.2f}%\n")
        self.detail_text.insert(tk.END,
            f"  组装序列2: {results['coverage_asm2']:.2f}%\n\n")

        # 得分明细
        self.detail_text.insert(tk.END, "【综合评分明细】\n")
        self.detail_text.insert(tk.END, "-" * 70 + "\n")
        s1 = results['score1']
        s2 = results['score2']
        self.detail_text.insert(tk.END,
            f"  {'维度':<20} {'权重':>8} {'组装1':>10} {'组装2':>10}\n")
        self.detail_text.insert(tk.END, f"  {'-'*50}\n")
        self.detail_text.insert(tk.END,
            f"  {'一致性得分':<20} {'40%':>8} {s1['identity_score']:>8.2f}分 {s2['identity_score']:>8.2f}分\n")
        self.detail_text.insert(tk.END,
            f"  {'覆盖度得分':<20} {'30%':>8} {s1['coverage_score']:>8.2f}分 {s2['coverage_score']:>8.2f}分\n")
        self.detail_text.insert(tk.END,
            f"  {'长度相似度':<20} {'15%':>8} {s1['length_score']:>8.2f}分 {s2['length_score']:>8.2f}分\n")
        self.detail_text.insert(tk.END,
            f"  {'GC相似度':<20} {'15%':>8} {s1['gc_score']:>8.2f}分 {s2['gc_score']:>8.2f}分\n")
        if 'structure_score' in s1:
            self.detail_text.insert(tk.END,
                f"  {'结构完整性':<20} {'+5%':>8} {s1['structure_score']:>8.2f}分 {s2['structure_score']:>8.2f}分\n")
        self.detail_text.insert(tk.END, f"  {'-'*50}\n")
        self.detail_text.insert(tk.END,
            f"  {'综合得分':<20} {'105%':>8} {s1['total']:>8.2f}分 {s2['total']:>8.2f}分\n\n")

        # 四区结构分析详情
        for st_key, st_label in [('structure1', '组装序列1'), ('structure2', '组装序列2')]:
            st = results.get(st_key)
            if not st or not isinstance(st, dict):
                continue
            details = st.get('structure_details', '')
            if details and details != '序列太短，跳过结构分析':
                self.detail_text.insert(tk.END, f"【四区结构分析 — {st_label}】\n")
                self.detail_text.insert(tk.END, "-" * 70 + "\n")
                self.detail_text.insert(tk.END, f"  {details}\n\n")

        self.detail_text.configure(state=tk.DISABLED)

    def _draw_window_plot(self, results):
        """在Canvas上绘制滑动窗口一致性图表"""
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

        # 背景
        self.canvas.create_rectangle(
            margin_left, margin_top,
            margin_left + plot_width, margin_top + plot_height,
            fill='#fafafa', outline='#e0e0e0'
        )

        # Y轴和标签
        for pct in [0, 25, 50, 75, 100]:
            y = margin_top + plot_height - (pct / 100 * plot_height)
            self.canvas.create_line(margin_left - 5, y, margin_left + plot_width, y,
                                    fill='#e0e0e0', dash=(2, 4))
            self.canvas.create_text(margin_left - 10, y,
                                    text=f"{pct}%", anchor=tk.E,
                                    font=('Consolas', 8), fill='#636e72')

        # X轴标签
        num_labels = 8
        for i in range(num_labels + 1):
            pos = int(ref_len * i / num_labels)
            x = margin_left + (plot_width * i / num_labels)
            self.canvas.create_text(x, margin_top + plot_height + 15,
                                    text=f"{pos:,}", anchor=tk.N,
                                    font=('Consolas', 7), fill='#636e72')

        # 轴标题
        self.canvas.create_text(margin_left + plot_width // 2,
                                margin_top + plot_height + 35,
                                text="参考序列位置 (bp)", anchor=tk.N,
                                font=('Microsoft YaHei', 9), fill='#2d3436')
        self.canvas.create_text(15, margin_top + plot_height // 2,
                                text="一\n致\n性\n%", anchor=tk.CENTER,
                                font=('Microsoft YaHei', 9), fill='#2d3436')

        # 绘制数据线
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
                # 图例
                legend_y = margin_top + (20 if label == "组装序列1" else 42)
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

        draw_line(w1['positions'], w1['identities'], '#0984e3', '组装序列1')
        draw_line(w2['positions'], w2['identities'], '#d63031', '组装序列2')

        # 标题
        self.canvas.create_text(margin_left + plot_width // 2, 10,
                                text="滑动窗口一致性分布图（窗口=5000bp, 步长=2000bp）",
                                anchor=tk.N,
                                font=('Microsoft YaHei', 10, 'bold'),
                                fill='#2d3436')

    def _clear_results(self):
        """清除所有结果"""
        if self.analysis_running:
            messagebox.showwarning("提示", "分析正在进行中，请等待完成后再清除。")
            return

        self.results = None
        self.summary_text.configure(state=tk.NORMAL)
        self.summary_text.delete(1.0, tk.END)
        self.summary_text.insert(tk.END, "请选择三个FASTA文件后，点击「开始分析」按钮。\n\n")
        self.summary_text.insert(tk.END, "支持的文件格式：.fasta, .fa, .fna, .txt\n")
        self.summary_text.insert(tk.END, "参考序列和组装序列各选择一个文件。\n")
        self.summary_text.configure(state=tk.DISABLED)

        self.detail_text.configure(state=tk.NORMAL)
        self.detail_text.delete(1.0, tk.END)
        self.detail_text.insert(tk.END, "详细比对数据将在分析完成后显示。\n")
        self.detail_text.configure(state=tk.DISABLED)

        self.canvas.delete("all")
        self.progress_var.set(0)
        self.progress_label.configure(text="")
        self.export_btn.configure(state=tk.DISABLED)
        self.status_label.configure(text="已清除结果")

    def _export_report(self):
        """导出报告到文件"""
        if not self.results:
            messagebox.showwarning("提示", "没有可导出的结果，请先进行分析。")
            return

        filepath = filedialog.asksaveasfilename(
            title="导出报告",
            defaultextension=".txt",
            filetypes=[
                ("文本文件", "*.txt"),
                ("所有文件", "*.*"),
            ],
            initialfile="叶绿体基因组比对报告.txt"
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
                f.write("  叶绿体基因组组装比对报告\n")
                f.write("=" * 70 + "\n\n")

                if s1['total'] > s2['total']:
                    f.write(f"结论：组装序列1 更接近参考序列（得分 {s1['total']:.2f} vs {s2['total']:.2f}）\n\n")
                elif s2['total'] > s1['total']:
                    f.write(f"结论：组装序列2 更接近参考序列（得分 {s2['total']:.2f} vs {s1['total']:.2f}）\n\n")
                else:
                    f.write(f"结论：两个组装结果得分相同（{s1['total']:.2f}）\n\n")

                f.write(f"【序列基本信息】\n")
                f.write(f"  参考序列: {ref['length']:,} bp, GC={ref['gc_content']:.2f}%\n")
                f.write(f"  组装序列1: {asm1['length']:,} bp, GC={asm1['gc_content']:.2f}%\n")
                f.write(f"  组装序列2: {asm2['length']:,} bp, GC={asm2['gc_content']:.2f}%\n\n")

                km1 = results['kmer_asm1']
                km2 = results['kmer_asm2']
                f.write(f"【K-mer相似度 (k=31)】\n")
                f.write(f"  组装序列1: Jaccard={km1['jaccard']:.6f}, 共享比例={km1['shared_ratio']:.6f}\n")
                f.write(f"  组装序列2: Jaccard={km2['jaccard']:.6f}, 共享比例={km2['shared_ratio']:.6f}\n\n")

                f.write(f"【参考序列覆盖度】\n")
                f.write(f"  组装序列1: {results['coverage_asm1']:.2f}%\n")
                f.write(f"  组装序列2: {results['coverage_asm2']:.2f}%\n\n")

                # 四区结构分析
                for st_key, st_label in [('structure1', '组装序列1'), ('structure2', '组装序列2')]:
                    st = results.get(st_key)
                    if not st or not isinstance(st, dict):
                        continue
                    details = st.get('structure_details', '')
                    if details and details != '序列太短，跳过结构分析':
                        f.write(f"【四区结构分析 — {st_label}】\n")
                        f.write(f"{details}\n\n")

                f.write(f"【综合评分】\n")
                f.write(f"  组装序列1: {s1['total']:.2f}分")
                if 'structure_score' in s1:
                    f.write(f" (基础 {results.get('score1_base', s1)['total']:.2f} + 结构 {s1['structure_score']:.2f})")
                f.write(f"\n")
                f.write(f"  组装序列2: {s2['total']:.2f}分")
                if 'structure_score' in s2:
                    f.write(f" (基础 {results.get('score2_base', s2)['total']:.2f} + 结构 {s2['structure_score']:.2f})")
                f.write(f"\n")

            messagebox.showinfo("导出成功", f"报告已保存到：\n{filepath}")
            self.status_label.configure(text=f"报告已导出: {os.path.basename(filepath)}")
        except Exception as e:
            messagebox.showerror("导出失败", f"导出报告时出错：\n{str(e)}")


# ============================================================
# 主入口
# ============================================================
def main():
    root = tk.Tk()

    # 设置窗口图标（如果存在）
    try:
        root.iconbitmap(default='')
    except:
        pass

    app = CpGenomeComparatorApp(root)

    # 窗口居中
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
