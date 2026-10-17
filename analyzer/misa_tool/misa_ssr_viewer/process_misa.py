#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
处理 MISA 输出文件，生成图 B / C / D 的绘图准备文件（可直接用于
GraphPad Prism、Excel、TBtools、R、Origin 等软件绘图）。

图B 热图   : 各物种 SSR 基序的存在(1)/缺失(0)矩阵
             - 基序按重复单元长度(1-10+ bp)排序，附基序信息表（用于上方分组子图）
             - 物种按 Ward 层次聚类排序（附聚类顺序表）
图C 饼图   : 所有物种 SSR 类别(p1-p6/c/c*)计数与百分比
图D 堆叠条 : 各物种 SSR 类别计数矩阵（按总数降序），附长格式表
"""

import os
import re
from collections import defaultdict, OrderedDict

# ------------------------- 配置 -------------------------
MISA_FILE = r"D:\零零碎碎\代码\NEW\0915\01_20.fasta.misa"
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(MISA_FILE)), "SSR_plot_data")

# 提取括号内的大写基序及其重复次数（小写为间隔序列，c* 中的 <...> 变异不会被匹配）
MOTIF_RE = re.compile(r'\(([ACGTN]{1,60})\)(\d+)')

# MISA 类型 -> (类别名称, 重复单元长度)
TYPE_INFO = OrderedDict([
    ('p1', ('Mono-nucleotide', 1)),
    ('p2', ('Di-nucleotide',  2)),
    ('p3', ('Tri-nucleotide', 3)),
    ('p4', ('Tetra-nucleotide', 4)),
    ('p5', ('Penta-nucleotide', 5)),
    ('p6', ('Hexa-nucleotide', 6)),
    ('c',  ('Compound SSR (c)', 0)),
    ('c*', ('Compound SSR (c*)', 0)),
])


# ------------------------- 解析 -------------------------
def parse_misa(path):
    """解析 .misa 文件，返回记录列表。"""
    records = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        fh.readline()  # 跳过表头
        for line in fh:
            line = line.rstrip('\r\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            records.append({'id': cols[0], 'type': cols[2], 'ssr': cols[3]})
    return records


def extract_motifs(ssr_string):
    """从 SSR 字符串中提取所有重复单元基序（复合 SSR 取其全部组成基序）。"""
    return [m.upper() for m, n in MOTIF_RE.findall(ssr_string)]


# ------------------- Ward 层次聚类（纯 Python 实现） -------------------
def ward_leaf_order(matrix):
    """
    对二值矩阵做 Ward 最小方差层次聚类（距离 = 平方欧氏距离，
    即差异基序数目），返回叶子（物种）顺序索引列表。
    """
    n = len(matrix)
    m = len(matrix[0]) if n else 0
    if n <= 1:
        return list(range(n))

    # 初始两两距离（平方欧氏）
    dist = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = sum((matrix[i][k] - matrix[j][k]) ** 2 for k in range(m))

    clusters = [{'id': i, 'members': [i], 'tree': (i,), 'size': 1} for i in range(n)]
    next_id = n

    while len(clusters) > 1:
        # 找距离最近的簇对
        best_a, best_b, best_v = 0, 1, None
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                ca, cb = clusters[a]['id'], clusters[b]['id']
                v = dist[(min(ca, cb), max(ca, cb))]
                if best_v is None or v < best_v:
                    best_v, best_a, best_b = v, a, b

        A, B = clusters[best_a], clusters[best_b]
        na, nb = A['size'], B['size']

        # Lance-Williams 更新 Ward 距离
        for C in clusters:
            if C is A or C is B:
                continue
            nc = C['size']
            ka = (min(A['id'], C['id']), max(A['id'], C['id']))
            kb = (min(B['id'], C['id']), max(B['id'], C['id']))
            d_new = ((na + nc) / (na + nb + nc)) * dist[ka] \
                  + ((nb + nc) / (na + nb + nc)) * dist[kb] \
                  - (nc / (na + nb + nc)) * dist[(min(A['id'], B['id']), max(A['id'], B['id']))]
            dist[(min(next_id, C['id']), max(next_id, C['id']))] = d_new

        merged = {
            'id': next_id,
            'members': A['members'] + B['members'],
            'tree': (A['tree'], B['tree']),
            'size': na + nb,
        }
        clusters = [c for idx, c in enumerate(clusters) if idx not in (best_a, best_b)]
        clusters.append(merged)
        next_id += 1

    def flatten(tree):
        if isinstance(tree, int):
            return [tree]
        out = []
        for sub in tree:
            out.extend(flatten(sub))
        return out

    return flatten(clusters[0]['tree'])


# ------------------------- 输出工具 -------------------------
def write_table(path, header, rows):
    with open(path, 'w', encoding='utf-8-sig', newline='\n') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')


def fmt_pct(x):
    return f"{x:.2f}%"


# ------------------------- 主流程 -------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    records = parse_misa(MISA_FILE)
    print(f"共解析 SSR 记录: {len(records)} 条")

    species = sorted({r['id'] for r in records})  # 物种（按编号自然排序）
    print(f"物种/样本数: {len(species)}  ({species[0]} ... {species[-1]})")

    # ---------- 统计 ----------
    species_type_count = defaultdict(lambda: defaultdict(int))  # 物种 -> 类型 -> 计数
    species_motifs = defaultdict(set)                            # 物种 -> 基序集合
    for r in records:
        species_type_count[r['id']][r['type']] += 1
        for mo in extract_motifs(r['ssr']):
            species_motifs[r['id']].add(mo)

    total_ssr = len(records)
    type_total = defaultdict(int)
    for r in records:
        type_total[r['type']] += 1

    # 基序排序：按重复单元长度，再按字母
    all_motifs = sorted({mo for s in species for mo in species_motifs[s]},
                        key=lambda x: (len(x), x))
    print(f"不同 SSR 基序总数: {len(all_motifs)}")

    # ================= 图B: 存在/缺失矩阵 =================
    # 原始物种顺序矩阵
    mat_orig = [[1 if mo in species_motifs[s] else 0 for mo in all_motifs]
                for s in species]

    # Ward 聚类后的物种顺序
    order_idx = ward_leaf_order(mat_orig)
    ward_species = [species[i] for i in order_idx]

    write_table(
        os.path.join(OUT_DIR, 'B1_heatmap_presence_absence_matrix.txt'),
        ['Species\\Motif'] + all_motifs,
        [[s] + mat_orig[i] for i, s in zip(order_idx, ward_species)],
    )
    write_table(
        os.path.join(OUT_DIR, 'B1b_heatmap_matrix_original_species_order.txt'),
        ['Species\\Motif'] + all_motifs,
        [[s] + row for s, row in zip(species, mat_orig)],
    )

    # 基序信息表（用于热图上方按重复单元长度分组的子图）
    write_table(
        os.path.join(OUT_DIR, 'B2_motif_info.txt'),
        ['MotifOrder', 'Motif', 'RepeatUnitLength_bp', 'SpeciesCount', 'Frequency_pct'],
        [[i + 1, mo, len(mo),
          sum(1 for s in species if mo in species_motifs[s]),
          fmt_pct(100.0 * sum(1 for s in species if mo in species_motifs[s]) / len(species))]
         for i, mo in enumerate(all_motifs)],
    )

    # 按重复单元长度分组的基序数量统计（上方子图数据）
    len_groups = defaultdict(int)
    for mo in all_motifs:
        key = str(len(mo)) if len(mo) <= 10 else '10+'
        len_groups[key] += 1
    len_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '10+']
    write_table(
        os.path.join(OUT_DIR, 'B3_motif_length_group_count.txt'),
        ['RepeatUnitLength_bp', 'MotifCount'],
        [[k, len_groups[k]] for k in len_order if k in len_groups],
    )

    # Ward 聚类物种顺序表
    write_table(
        os.path.join(OUT_DIR, 'B4_species_ward_cluster_order.txt'),
        ['ClusterOrder', 'Species', 'Total_SSR', 'MotifCount'],
        [[i + 1, s, sum(species_type_count[s].values()), len(species_motifs[s])]
         for i, s in enumerate(ward_species)],
    )

    # ================= 图C: 类别百分比（饼图） =================
    rows_c = []
    for t, (name, _) in TYPE_INFO.items():
        cnt = type_total.get(t, 0)
        rows_c.append([t, name, cnt, fmt_pct(100.0 * cnt / total_ssr)])
    # 汇总口径：复合 SSR 合并
    comp = type_total.get('c', 0) + type_total.get('c*', 0)
    rows_c2 = []
    for t in ['p1', 'p2', 'p3', 'p4', 'p5', 'p6']:
        cnt = type_total.get(t, 0)
        rows_c2.append([TYPE_INFO[t][0], cnt, fmt_pct(100.0 * cnt / total_ssr)])
    rows_c2.append(['Compound SSR (c+c*)', comp, fmt_pct(100.0 * comp / total_ssr)])

    write_table(os.path.join(OUT_DIR, 'C_pie_SSR_category_percentage.txt'),
                ['MISA_Type', 'Category', 'SSR_Count', 'Percentage'], rows_c)
    write_table(os.path.join(OUT_DIR, 'C_pie_SSR_category_percentage_merged.txt'),
                ['Category', 'SSR_Count', 'Percentage'], rows_c2)

    # ================= 图D: 各物种类别堆叠条形图 =================
    # 物种按总 SSR 数降序
    species_by_total = sorted(species, key=lambda s: -sum(species_type_count[s].values()))
    cats = list(TYPE_INFO.keys())
    write_table(
        os.path.join(OUT_DIR, 'D1_stackedbar_species_category_counts.txt'),
        ['Species'] + [f"{t}({TYPE_INFO[t][0]})" for t in cats] + ['Total_SSR'],
        [[s] + [species_type_count[s].get(t, 0) for t in cats]
             + [sum(species_type_count[s].values())] for s in species_by_total],
    )
    # 长格式（适合 R/ggplot2、Origin 等）
    rows_long = []
    for s in species_by_total:
        for t in cats:
            rows_long.append([s, t, TYPE_INFO[t][0],
                              species_type_count[s].get(t, 0),
                              sum(species_type_count[s].values())])
    write_table(os.path.join(OUT_DIR, 'D2_stackedbar_long_format.txt'),
                ['Species', 'MISA_Type', 'Category', 'Count', 'SpeciesTotal'],
                rows_long)

    # ---------------- 汇总 ----------------
    print('\n=== 输出文件（目录: %s） ===' % OUT_DIR)
    for fn in sorted(os.listdir(OUT_DIR)):
        print('  -', fn)
    print('\n类别计数（图C）:')
    for r in rows_c:
        print(f"  {r[0]:>3} {r[1]:<22} {r[2]:>6}  {r[3]}")
    print('\n物种 SSR 总数（图D，降序前10）:')
    for s in species_by_total[:10]:
        print(f"  {s}: {sum(species_type_count[s].values())}")


if __name__ == '__main__':
    main()
