#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

# ------------------------
# 用户配置
# ------------------------
TSV_FILTERED_DIR = Path("foldseek_all_vs_all/tsv_filtered")  # 第一步输出目录
OUT_FILE = Path("foldseek_all_vs_all/complex_tmscore.tsv")   # 输出 complex 结果

# ------------------------
# 读取所有 chunk TSV
# ------------------------
all_files = sorted(TSV_FILTERED_DIR.glob("*.tsv"))
if len(all_files) == 0:
    raise RuntimeError(f"No TSV files found in {TSV_FILTERED_DIR}")

print(f"Found {len(all_files)} chunk TSV files")

dfs = []
for f in all_files:
    df = pd.read_csv(
        f,
        sep="\t",
        header=None,
        names=["query","target","evalue","bits","alnlen","fident","qtmscore","ttmscore"]
    )
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

# ------------------------
# 去掉链后缀，得到复合物 ID
# 假设链后缀是 "_E", "_A", "_B" 等
# ------------------------
df["query_complex"] = df["query"].str.replace(r'_[A-Z]$', '', regex=True)
df["target_complex"] = df["target"].str.replace(r'_[A-Z]$', '', regex=True)

# ------------------------
# 计算复合物整体相似度
# 对每个复合物 pair：
#   1. 按 query 链分组，选每条链的最大 qtmscore
#   2. 对这些最大值求平均
# ------------------------
complex_scores = df.groupby(["query_complex","target_complex"]).apply(
    lambda g: g.groupby("query")["qtmscore"].max().mean()
).reset_index(name="complex_qtmscore")

# ------------------------
# 保存结果
# ------------------------
complex_scores.to_csv(OUT_FILE, sep="\t", index=False)

print(f"✅ Done! Complex-level TM-scores saved to {OUT_FILE}")
print(complex_scores.head())