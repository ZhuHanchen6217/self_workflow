#!/usr/bin/env python3
import subprocess
from pathlib import Path

########################
# 用户配置
########################
STRUCT_DIR = Path("/home/zhc/dalab/ab_ag_structures/SAbDab/PPI_DB/structures")
THREADS = 8
CHUNK_SIZE = 500   # 16GB 内存建议 300~800

########################
# 自动路径
########################
SCRIPT_DIR = Path(__file__).resolve().parent
OUTDIR = SCRIPT_DIR / "foldseek_all_vs_all"
TMPDIR = OUTDIR / "tmp"

DB = OUTDIR / "full_db"
RESULT_DIR = OUTDIR / "results"
TSV_DIR = OUTDIR / "tsv"
TSV_FILTERED_DIR = OUTDIR / "tsv_filtered"

########################
def run(cmd, desc=None):
    if desc:
        print(f"\n▶ {desc}")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

def chunked(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def filter_tsv(input_tsv: Path, output_tsv: Path):
    """
    过滤规则：
    1) 去掉 query == target
    2) 只保留上三角：仅保留 query < target
    3) min(qtmscore, ttmscore) >= 0.5
    """
    with input_tsv.open("r") as fin, output_tsv.open("w") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            q, t = parts[0], parts[1]
            qtmscore = float(parts[6])
            ttmscore = float(parts[7])

            if q == t:
                continue
            if q >= t:
                continue
            if min(qtmscore, ttmscore) < 0.5:
                continue

            fout.write(line)

def main():
    OUTDIR.mkdir(exist_ok=True)
    TMPDIR.mkdir(exist_ok=True)
    RESULT_DIR.mkdir(exist_ok=True)
    TSV_DIR.mkdir(exist_ok=True)
    TSV_FILTERED_DIR.mkdir(exist_ok=True)

    # 1️⃣ 收集结构
    structures = sorted(STRUCT_DIR.glob("*.cif.gz"))
    if len(structures) == 0:
        raise RuntimeError("No .cif.gz files found")
    print(f"Total structures: {len(structures)}")

    # 2️⃣ 创建全库 DB
    if not DB.exists():
        run([
            "foldseek", "createdb",
            str(STRUCT_DIR),
            str(DB),
            "--threads", str(THREADS)
        ], "Creating Foldseek database")

    # 3️⃣ 创建索引（只做一次）
    if not (DB.parent / (DB.name + ".idx")).exists():
        run([
            "foldseek", "createindex",
            str(DB),
            str(TMPDIR),
            "--threads", str(THREADS)
        ], "Creating index")

    # 4️⃣ 分块 all-vs-all：chunk（query） vs full（target）
    chunks = list(chunked(structures, CHUNK_SIZE))
    total_chunks = len(chunks)

    for chunk_id, chunk in enumerate(chunks, start=1):
        print(f"\n=== Chunk {chunk_id}/{total_chunks}: {len(chunk)} structures ===")

        chunk_dir = OUTDIR / f"chunk_{chunk_id}_dir"
        chunk_db = OUTDIR / f"chunk_{chunk_id}_db"
        chunk_result = RESULT_DIR / f"chunk_{chunk_id}_vs_full"
        chunk_tsv = TSV_DIR / f"chunk_{chunk_id}_vs_full.tsv"
        chunk_tsv_filtered = TSV_FILTERED_DIR / f"chunk_{chunk_id}_vs_full.filtered.tsv"

        # ✅ 断点续跑：如果过滤后的 TSV 已存在，直接跳过
        if chunk_tsv_filtered.exists():
            print(f"Skip chunk {chunk_id}, filtered TSV exists: {chunk_tsv_filtered.name}")
            continue

        # 4.1 为本 chunk 创建临时目录（放软链接）
        chunk_dir.mkdir(exist_ok=True)

        # 清理旧链接
        for p in chunk_dir.glob("*"):
            p.unlink()

        # 创建软链接到原始结构
        for p in chunk:
            (chunk_dir / p.name).symlink_to(p)

        # 4.2 createdb（从临时目录读取）
        run([
            "foldseek", "createdb",
            str(chunk_dir),
            str(chunk_db),
            "--threads", str(THREADS)
        ], f"Creating chunk DB {chunk_id}")

        # 4.3 search（chunk vs full）
        run([
            "foldseek", "search",
            str(chunk_db),
            str(DB),
            str(chunk_result),
            str(TMPDIR),
            "--threads", str(THREADS),
            "-s", "9.5",
            "--alignment-type", "1",
            "-a",
            "--max-seqs", "300"
        ], f"Searching chunk {chunk_id} vs full")

        # 4.4 convert to TSV
        run([
            "foldseek", "convertalis",
            str(chunk_db),
            str(DB),
            str(chunk_result),
            str(chunk_tsv),
            "--format-output",
            "query,target,evalue,bits,alnlen,fident,qtmscore,ttmscore"
        ], f"Converting chunk {chunk_id} results")

        # 4.5 过滤：去掉自比对 + 上三角 + min(TMs) >= 0.5
        print(f"▶ Filtering {chunk_tsv.name} -> {chunk_tsv_filtered.name}")
        filter_tsv(chunk_tsv, chunk_tsv_filtered)

    print("\n✅ DONE")
    print(f"Filtered TSV files in: {TSV_FILTERED_DIR}")

if __name__ == "__main__":
    main()