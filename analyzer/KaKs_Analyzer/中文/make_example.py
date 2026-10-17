# -*- coding: utf-8 -*-
"""从原始数据生成示例文件, 用于演示 KaKs 分析工具."""
import sys
import pandas as pd
from pathlib import Path

sys.path.insert(0, r'D:\零零碎碎\代码\NEW\0921')

from kaks_analyzer import KaksProcessor

OUT_DIR = Path(r'D:\零零碎碎\代码\NEW\0921')

# 加载原始数据
src = OUT_DIR / "merge.axt.kaks"
df = KaksProcessor.load_file(str(src))
print(f"原始数据: {len(df)} 行, {df['Gene'].nunique()} 个基因, "
      f"{df['Species2'].nunique()} 个物种")

# 取前 8 个基因作为示例 (更精简, 涵盖各种缺失率)
sample_genes = sorted(df['Gene'].unique())[:8]
sample = df[df['Gene'].isin(sample_genes)].copy()
print(f"示例基因 ({len(sample_genes)}): {sample_genes}")
print(f"示例数据: {len(sample)} 行")

# 简化 Species1/Species2: 去掉数字 ID 后缀, 保留 genus + species
def short_name(s):
    parts = s.split('_')
    # 只保留前两个 token (genus + species), 如 Asparagus_acutifolius
    return '_'.join(parts[:2]) if len(parts) >= 2 else s

sample['Species1'] = sample['Species1'].apply(short_name)
sample['Species2'] = sample['Species2'].apply(short_name)

# 整理为 6 列: Species1, Species2, Gene, Ka, Ks, Ka/Ks
out = sample[['Species1', 'Species2', 'Gene', 'Ka', 'Ks', 'Ka/Ks']].copy()
# 物种名按字母序排序, 便于查看
out = out.sort_values(['Gene', 'Species1', 'Species2']).reset_index(drop=True)

# 写出 xlsx
out_xlsx = OUT_DIR / "example_kaks.xlsx"
out.to_excel(out_xlsx, index=False)
print(f"已生成: {out_xlsx}")

# 同时导出原始格式 (.kaks) - 读取原始 22 列版本
out_kaks = OUT_DIR / "example_kaks.kaks"
raw_df = pd.read_csv(src, sep='\t')
# 只保留示例基因 (按 _mafft 后缀的原始 gene 名)
raw_sample_genes = sorted(raw_df['Sequence'].apply(
    lambda s: str(s).split('|')[-1].replace('_mafft', '').strip()
).unique())[:8]
sample_raw = raw_df[raw_df['Sequence'].apply(
    lambda s: str(s).split('|')[-1].replace('_mafft', '').strip()
).isin(raw_sample_genes)].copy()
sample_raw.to_csv(out_kaks, sep='\t', index=False)
print(f"已生成: {out_kaks} (原始 .kaks 格式)")

# 输出示例统计
print()
print("示例数据中各基因 Ka/Ks 列缺失比例分布:")
sample["__miss__"] = sample['Ka/Ks'].isna() | (sample['Ka/Ks'] == 0)
miss_rate = sample.groupby('Gene')['__miss__'].mean()
print(miss_rate.describe())
print()
print("示例数据 Ka/Ks=0 的数量:", (sample['Ka/Ks']==0).sum())
print("示例数据 Ka/Ks NaN 的数量:", sample['Ka/Ks'].isna().sum())