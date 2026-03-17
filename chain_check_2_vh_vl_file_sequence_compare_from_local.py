#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
并行比较 CSV 中 VH/VL 序列 与 本地 vh_vl_dataset 目录下文件中的 VH/VL 链序列，并输出单页翻页对齐查看器（HTML）。

输入（同目录）：
  final_structures_sift_v1_251203.csv
必须列：
  - pdb_id          （对应 vh_vl_dataset 下面的文件名，不包括后缀）
  - Hchain, VH_sequence
  - Lchain, VL_sequence

本地文件：
  - 默认目录：vh_vl_dataset/
  - 文件名：{pdb_id}.*  （支持 .fasta/.fa/.faa/.fna/.txt/.pdb/.cif；若同名多后缀，按优先级选择）

取链策略（严格区分大小写）：
  - 对 FASTA：优先匹配 header 中的 "[auth X]"（X=链名）
  - 若没命中，再用 "Chain(s) X"/"chain(s) X" 匹配
  - 若仍没命中：
      - 如果 FASTA 只有 1 条非空 entry：直接取（single-entry fallback）
      - 否则失败：NO_FILE_CHAIN

对齐方式：
  - Needleman–Wunsch 全局对齐
  - identity_no_gaps = matches/(matches+mismatches)，忽略 gaps

输出：
  1) vh_vl_file_vs_local_chain_report_2.csv
  2) alignments/chain_check_2.html（仅收录 status=DIFFERENT 的条目；H/L 链独立条目）
  3) alignments/chain_check_2_focus.html（仅收录 DIFFERENT 且 identity_no_gaps != 1 的条目）

注意：
  - alignments/ 文件夹若已存在，不做覆盖/清理；仅覆盖写入本次生成的两个 html 文件（同名则覆盖）。
  - 链名严格区分大小写（例如 'a' 和 'A' 不视为同一条链）。

依赖：
  pip install tqdm
"""

from __future__ import annotations

import csv
import html as html_mod
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from tqdm import tqdm


# =========================
# 配置
# =========================

INPUT_CSV = "final_structures_sift_v1_251203.csv"
DATASET_DIR = "vh_vl_dataset"

OUTPUT_CSV = "vh_vl_file_vs_local_chain_report_2.csv"

ALIGN_DIR = "alignments"
INDEX_HTML = "chain_check_2.html"
FOCUS_HTML = "chain_check_2_focus.html"

THREADS = 10

# Needleman–Wunsch 全局对齐参数
MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2

WRAP = 80
WHITESPACE_RE = re.compile(r"\s+")

# header 匹配（严格区分大小写）
AUTH_TAG_RE_TPL = r"\[auth\s+{chain}\]"

CHAIN_PATTERNS = [
    re.compile(r"\bChain\s+([A-Za-z0-9])\b"),
    re.compile(r"\bChains\s+([A-Za-z0-9, ]+)\b"),
    re.compile(r"\bchain\s+([A-Za-z0-9])\b"),
    re.compile(r"\bchains\s+([A-Za-z0-9, ]+)\b"),
]

# 支持的后缀优先级（越靠前优先）
FILE_SUFFIX_PRIORITY = [
    ".fasta",
    ".fa",
    ".faa",
    ".fna",
    ".txt",
    ".pdb",
    ".cif",
]


# =========================
# 数据结构
# =========================

def normalize_seq(seq: str) -> str:
    seq = WHITESPACE_RE.sub("", (seq or "").upper())
    return seq.replace("*", "")


def looks_like_fasta(text: str) -> bool:
    if not text:
        return False
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        return s.startswith(">")
    return False


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq_parts: List[str] = []
    for line in (text or "").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append((header, normalize_seq("".join(seq_parts))))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        entries.append((header, normalize_seq("".join(seq_parts))))
    return entries


def header_mentions_auth_chain(header: str, auth_chain: str) -> bool:
    h = header or ""
    c = (auth_chain or "").strip()
    if not c:
        return False
    auth_re = re.compile(AUTH_TAG_RE_TPL.format(chain=re.escape(c)))  # no IGNORECASE
    return bool(auth_re.search(h))


def header_mentions_chain_letter(header: str, chain_letter: str) -> bool:
    h = header or ""
    c = (chain_letter or "").strip()
    if not c:
        return False
    for pat in CHAIN_PATTERNS:
        m = pat.search(h)
        if not m:
            continue
        val = m.group(1)
        tokens = [t.strip() for t in re.split(r"[,	]+", val) if t.strip()]
        if c in tokens:
            return True
    return False


def pick_by_auth(entries: List[Tuple[str, str]], chain: str) -> Optional[str]:
    for header, seq in entries:
        if seq and header_mentions_auth_chain(header, chain):
            return seq
    return None


def pick_by_chain_letter(entries: List[Tuple[str, str]], chain: str) -> Optional[str]:
    for header, seq in entries:
        if seq and header_mentions_chain_letter(header, chain):
            return seq
    return None


def select_dataset_file(dataset_dir: Path, pdb_id_raw: str) -> Optional[Path]:
    """Select file by pdb_id (no suffix). If multiple matches, choose by suffix priority."""
    base = (pdb_id_raw or "").strip()
    if not base:
        return None

    # Prefer known suffixes
    for suf in FILE_SUFFIX_PRIORITY:
        p = dataset_dir / f"{base}{suf}"
        if p.exists() and p.is_file():
            return p

    # Fallback: any file whose stem == base
    # (e.g. weird extension)
    matches = [p for p in dataset_dir.glob(f"{base}.*") if p.is_file()]
    if not matches:
        return None

    # Deterministic order
    matches.sort(key=lambda x: x.name)
    return matches[0]


def fetch_chain_sequence_from_file(dataset_dir: Path, pdb_id_raw: str, chain: str) -> ChainSeqResult:
    chain = (chain or "").strip()  # keep case
    pdb_id_raw = (pdb_id_raw or "").strip()

    fp = select_dataset_file(dataset_dir, pdb_id_raw)
    if fp is None:
        return ChainSeqResult(pdb_id=pdb_id_raw, chain_id=chain, found=False, reason="no_file")

    try:
        text = fp.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return ChainSeqResult(
            pdb_id=pdb_id_raw,
            chain_id=chain,
            found=False,
            reason=f"error:read:{type(e).__name__}:{e}",
            file_path=str(fp),
        )

    if not looks_like_fasta(text):
        return ChainSeqResult(
            pdb_id=pdb_id_raw,
            chain_id=chain,
            found=False,
            reason="not_fasta_file",
            file_path=str(fp),
        )

    entries = parse_fasta(text)

    seq_auth = pick_by_auth(entries, chain)
    if seq_auth:
        return ChainSeqResult(
            pdb_id=pdb_id_raw,
            chain_id=chain,
            found=True,
            sequence=seq_auth,
            reason="file_auth_matched",
            file_path=str(fp),
        )

    seq_label = pick_by_chain_letter(entries, chain)
    if seq_label:
        return ChainSeqResult(
            pdb_id=pdb_id_raw,
            chain_id=chain,
            found=True,
            sequence=seq_label,
            reason="file_chain_letter_matched",
            file_path=str(fp),
        )

    nonempty = [(h, s) for (h, s) in entries if s]
    if len(nonempty) == 1:
        return ChainSeqResult(
            pdb_id=pdb_id_raw,
            chain_id=chain,
            found=True,
            sequence=nonempty[0][1],
            reason="file_single_entry_fallback",
            file_path=str(fp),
        )

    return ChainSeqResult(
        pdb_id=pdb_id_raw,
        chain_id=chain,
        found=False,
        reason="no_file_chain_match",
        file_path=str(fp),
    )


# =========================
# 全局对齐（Needleman–Wunsch）
# =========================

def needleman_wunsch(a: str, b: str, match: int, mismatch: int, gap: int) -> Tuple[str, str]:
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[0] * (m + 1) for _ in range(n + 1)]  # 0 diag, 1 up, 2 left

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
        tb[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap
        tb[0][j] = 2

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            bj = b[j - 1]
            s_diag = score[i - 1][j - 1] + (match if ai == bj else mismatch)
            s_up = score[i - 1][j] + gap
            s_left = score[i][j - 1] + gap
            best = s_diag
            move = 0
            if s_up > best:
                best = s_up
                move = 1
            if s_left > best:
                best = s_left
                move = 2
            score[i][j] = best
            tb[i][j] = move

    i, j = n, m
    a_out: List[str] = []
    b_out: List[str] = []
    while i > 0 or j > 0:
        move = tb[i][j]
        if i > 0 and j > 0 and move == 0:
            a_out.append(a[i - 1])
            b_out.append(b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or move == 1):
            a_out.append(a[i - 1])
            b_out.append("-")
            i -= 1
        else:
            a_out.append("-")
            b_out.append(b[j - 1])
            j -= 1

    return "".join(reversed(a_out)), "".join(reversed(b_out))


# =========================
# HTML Viewer
# =========================

def chunk(s: str, width: int) -> List[str]:
    return [s[i : i + width] for i in range(0, len(s), width)]


def make_colored_lines(local_block: str, file_block: str) -> Tuple[str, str, str]:
    local_spans: List[str] = []
    file_spans: List[str] = []
    mid: List[str] = []

    for v, p in zip(local_block, file_block):
        if v == p and v != "-":
            local_spans.append(f'<span class="match">{html_mod.escape(v)}</span>')
            file_spans.append(f'<span class="match">{html_mod.escape(p)}</span>')
            mid.append("|")
        elif v == "-" or p == "-":
            vcls = "gap" if v == "-" else "insdel"
            pcls = "gap" if p == "-" else "insdel"
            local_spans.append(f'<span class="{vcls}">{html_mod.escape(v)}</span>')
            file_spans.append(f'<span class="{pcls}">{html_mod.escape(p)}</span>')
            mid.append(" ")
        else:
            local_spans.append(f'<span class="mismatch">{html_mod.escape(v)}</span>')
            file_spans.append(f'<span class="mismatch">{html_mod.escape(p)}</span>')
            mid.append(".")

    return "".join(local_spans), "".join(file_spans), html_mod.escape("".join(mid))


def render_alignment_html(local_aln: str, file_aln: str, wrap: int) -> str:
    local_chunks = chunk(local_aln, wrap)
    file_chunks = chunk(file_aln, wrap)
    blocks = []
    for i, (vblk, pblk) in enumerate(zip(local_chunks, file_chunks), start=1):
        local_line, file_line, mid_line = make_colored_lines(vblk, pblk)
        blocks.append(
            f"""
            <div class="block">
              <div class="label">Block {i}</div>
              <pre class="seq"><span class="name">LOCAL</span>{local_line}</pre>
              <pre class="seq"><span class="name">FILE </span>{file_line}</pre>
              <pre class="mid"><span class="name">     </span>{mid_line}</pre>
            </div>
            """
        )
    return "\n".join(blocks)


def make_index_html(items: List[dict], wrap: int, title_text: str) -> str:
    rendered = []
    for it in items:
        rendered.append(
            {
                "title": it["title"],
                "meta": it["meta"],
                "stats": it["stats"],
                "html": render_alignment_html(it["local_aln"], it["file_aln"], wrap=wrap),
            }
        )

    data_json = json.dumps(rendered, ensure_ascii=False)

    css = """
    body { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono","Courier New", monospace; margin: 16px; }
    h1 { font-size: 16px; margin: 0 0 8px 0; }
    .toolbar { display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin: 10px 0 12px 0; }
    button { padding: 6px 10px; }
    select { padding: 6px 10px; min-width: 520px; }
    .meta { color: #444; margin: 6px 0 10px 0; }
    .stats { color: #333; margin: 6px 0 10px 0; }
    .block { margin: 14px 0; padding: 10px; border: 1px solid #ddd; border-radius: 8px; }
    .label { color: #666; margin-bottom: 6px; }
    pre { margin: 0; white-space: pre; overflow-x: auto; }
    .seq .name, .mid .name { display: inline-block; width: 5ch; color: #666; }
    .match { color: #111; }
    .mismatch { color: #b00020; font-weight: 700; }
    .gap { color: #999; background: #f3f3f3; }
    .insdel { color: #0b57d0; font-weight: 700; }
    .mid { color: #666; }
    .hint { color: #666; font-size: 12px; }
    .errorbox { color: #b00020; white-space: pre-wrap; border: 1px solid #f1b; padding: 10px; border-radius: 8px; }
    """

    js = """
    let items = [];
    let idx = 0;

    function showError(err) {
      const box = document.getElementById('content');
      box.innerHTML = '';
      const pre = document.createElement('pre');
      pre.className = 'errorbox';
      pre.textContent = String(err && err.stack ? err.stack : err);
      box.appendChild(pre);
    }

    function render() {
      const it = items[idx];
      document.getElementById('title').textContent = it.title ?? '';
      document.getElementById('meta').textContent = it.meta ?? '';
      document.getElementById('stats').textContent = it.stats ?? '';
      document.getElementById('content').innerHTML = (typeof it.html === 'string') ? it.html : String(it.html);
      document.getElementById('counter').textContent = `${{idx+1}} / ${items.length}`;
      document.getElementById('sel').value = String(idx);
    }

    function next() { if (idx < items.length-1) { idx++; render(); } }
    function prev() { if (idx > 0) { idx--; render(); } }

    function init() {
      try {
        const raw = document.getElementById('data').textContent;
        items = JSON.parse(raw);

        const sel = document.getElementById('sel');
        items.forEach((it, i) => {
          const opt = document.createElement('option');
          opt.value = String(i);
          opt.textContent = `${{i+1}}. ${it.title ?? ''}`;
          sel.appendChild(opt);
        });
        sel.addEventListener('change', (e) => {
          idx = parseInt(e.target.value, 10);
          render();
        });

        document.getElementById('prev').addEventListener('click', prev);
        document.getElementById('next').addEventListener('click', next);

        document.addEventListener('keydown', (e) => {
          if (e.key === 'ArrowLeft') prev();
          if (e.key === 'ArrowRight') next();
        });

        if (items.length === 0) {
          document.getElementById('title').textContent = 'No entries';
          document.getElementById('meta').textContent = '';
          document.getElementById('stats').textContent = '';
          document.getElementById('content').innerHTML = '<p>没有符合条件的条目。</p>';
          return;
        }
        render();
      } catch (e) {
        showError(e);
      }
    }

    window.addEventListener('load', init);
    """

    return f"""<!doctype html>
<html>
<head>
<meta charset=\"utf-8\