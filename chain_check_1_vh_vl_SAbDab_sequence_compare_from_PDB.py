#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
并行比较 CSV 中 VH/VL 序列与 RCSB PDB FASTA 链序列，并输出单页翻页对齐查看器（HTML）。

输入（同目录）：
  final_structures_sift_v1_251203.csv
必须列：
  - pdb_id          （取前四个字符作为 PDB entry id）
  - Hchain, VH_sequence
  - Lchain, VL_sequence

处理方式：
  - 输出 CSV 为“一行一个链”（chain_type=H/L）
  - 当链名或本地序列为空时：status=SKIP，note=missing_local_chain_or_sequence

取链策略（按优先级，严格区分大小写）：
1) 全局 display FASTA：/fasta/entry/{pdb}/display
   - 先严格匹配 “[auth X]”（X = chainId）
2) 若没命中，再用 “Chain(s) X” 匹配（严格大小写）
3) 再不行用 display?chainId=X fallback：
   - 若仅 1 条非空 entry：直接取
   - 否则再用 “Chain(s) X” 匹配
4) 都失败：NO_PDB

输出：
  1) vh_vl_vs_pdb_chain_report_1.csv
  2) alignments/chain_check_1.html（仅收录 status=DIFFERENT 的条目；H/L 链独立条目）
  3) alignments/chain_check_1_focus.html（仅收录 DIFFERENT 且 identity_no_gaps != 1 的条目）

注意：
  - alignments/ 文件夹若已存在，不做覆盖/清理；仅覆盖写入本次生成的两个 html 文件（同名则覆盖）。
  - 本脚本对链名比较严格区分大小写（例如 'a' 和 'A' 不视为同一条链）。

依赖：
  pip install requests tqdm
"""

from __future__ import annotations

import csv
import html as html_mod
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import requests
from tqdm import tqdm


# =========================
# 配置
# =========================

INPUT_CSV = "final_structures_sift_v1_251203.csv"
OUTPUT_CSV = "vh_vl_vs_pdb_chain_report_1.csv"

ALIGN_DIR = "alignments"
INDEX_HTML = "chain_check_1.html"
FOCUS_HTML = "chain_check_1_focus.html"

FASTA_DISPLAY_URL_TMPL = "https://www.rcsb.org/fasta/entry/{pdb_id}/display"

REQUEST_TIMEOUT_S = 20
RETRIES = 2
RETRY_BACKOFF_S = 1.5
THREADS = 10

# Needleman–Wunsch 全局对齐参数
MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2

WRAP = 80
WHITESPACE_RE = re.compile(r"\s+")


# =========================
# 数据结构
# =========================

@dataclass
class ChainSeqResult:
    pdb_id: str
    chain_id: str  # value from CSV (case-sensitive)
    found: bool
    sequence: Optional[str] = None
    reason: str = "ok"


# =========================
# FASTA 获取与解析
# =========================

def normalize_seq(seq: str) -> str:
    seq = WHITESPACE_RE.sub("", (seq or "").upper())
    return seq.replace("*", "")


def http_get_with_retry(url: str, session: requests.Session) -> requests.Response:
    last_err = None
    for attempt in range(RETRIES + 1):
        try:
            return session.get(url, timeout=REQUEST_TIMEOUT_S, headers={"Accept": "text/plain,*/*"})
        except Exception as e:
            last_err = e
            if attempt < RETRIES:
                time.sleep(RETRY_BACKOFF_S * (attempt + 1))
            else:
                raise
    raise last_err  # pragma: no cover


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


# =========================
# header 匹配（严格区分大小写）
# =========================

AUTH_TAG_RE_TPL = r"\[auth\s+{chain}\]"

CHAIN_PATTERNS = [
    re.compile(r"\bChain\s+([A-Za-z0-9])\b"),
    re.compile(r"\bChains\s+([A-Za-z0-9, ]+)\b"),
    re.compile(r"\bchain\s+([A-Za-z0-9])\b"),
    re.compile(r"\bchains\s+([A-Za-z0-9, ]+)\b"),
]


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
        tokens = [t.strip() for t in re.split(r"[,\s]+", val) if t.strip()]
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


# =========================
# 取链序列（按要求先后顺序）
# =========================

def fetch_chain_sequence(pdb_id: str, chain: str, session: requests.Session) -> ChainSeqResult:
    pdb_id = (pdb_id or "").strip().upper()
    chain = (chain or "").strip()  # keep case

    if len(pdb_id) != 4:
        return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=False, reason="bad_pdb_id")

    try:
        url_full = FASTA_DISPLAY_URL_TMPL.format(pdb_id=pdb_id)
        r_full = http_get_with_retry(url_full, session=session)
        if r_full.status_code == 404:
            return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=False, reason="no_entry")
        r_full.raise_for_status()
        if not looks_like_fasta(r_full.text):
            return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=False, reason="not_fasta_response")

        entries_full = parse_fasta(r_full.text)

        seq_auth = pick_by_auth(entries_full, chain)
        if seq_auth:
            return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=True, sequence=seq_auth, reason="full_display_auth_matched")

        seq_label = pick_by_chain_letter(entries_full, chain)
        if seq_label:
            return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=True, sequence=seq_label, reason="full_display_chain_letter_matched")

        url_chain = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display?chainId={chain}"
        r_chain = http_get_with_retry(url_chain, session=session)
        if r_chain.status_code == 404:
            return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=False, reason="no_entry")

        if r_chain.ok and looks_like_fasta(r_chain.text):
            entries_chain = parse_fasta(r_chain.text)

            nonempty = [(h, s) for (h, s) in entries_chain if s]
            if len(nonempty) == 1:
                return ChainSeqResult(
                    pdb_id=pdb_id,
                    chain_id=chain,
                    found=True,
                    sequence=nonempty[0][1],
                    reason="fallback_chainId_single_entry",
                )

            seq_chain_label = pick_by_chain_letter(entries_chain, chain)
            if seq_chain_label:
                return ChainSeqResult(
                    pdb_id=pdb_id,
                    chain_id=chain,
                    found=True,
                    sequence=seq_chain_label,
                    reason="fallback_chainId_chain_letter_matched",
                )

        return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=False, reason="no_chain_match")

    except Exception as e:
        return ChainSeqResult(pdb_id=pdb_id, chain_id=chain, found=False, reason=f"error:{type(e).__name__}:{e}")


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


def make_colored_lines(local_block: str, pdb_block: str) -> Tuple[str, str, str]:
    local_spans: List[str] = []
    pdb_spans: List[str] = []
    mid: List[str] = []

    for v, p in zip(local_block, pdb_block):
        if v == p and v != "-":
            local_spans.append(f'<span class="match">{html_mod.escape(v)}</span>')
            pdb_spans.append(f'<span class="match">{html_mod.escape(p)}</span>')
            mid.append("|")
        elif v == "-" or p == "-":
            vcls = "gap" if v == "-" else "insdel"
            pcls = "gap" if p == "-" else "insdel"
            local_spans.append(f'<span class="{vcls}">{html_mod.escape(v)}</span>')
            pdb_spans.append(f'<span class="{pcls}">{html_mod.escape(p)}</span>')
            mid.append(" ")
        else:
            local_spans.append(f'<span class="mismatch">{html_mod.escape(v)}</span>')
            pdb_spans.append(f'<span class="mismatch">{html_mod.escape(p)}</span>')
            mid.append(".")

    return "".join(local_spans), "".join(pdb_spans), html_mod.escape("".join(mid))


def render_alignment_html(local_aln: str, pdb_aln: str, wrap: int) -> str:
    local_chunks = chunk(local_aln, wrap)
    pdb_chunks = chunk(pdb_aln, wrap)
    blocks = []
    for i, (vblk, pblk) in enumerate(zip(local_chunks, pdb_chunks), start=1):
        local_line, pdb_line, mid_line = make_colored_lines(vblk, pblk)
        blocks.append(
            f"""
            <div class="block">
              <div class="label">Block {i}</div>
              <pre class="seq"><span class="name">LOCAL</span>{local_line}</pre>
              <pre class="seq"><span class="name">PDB  </span>{pdb_line}</pre>
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
                "html": render_alignment_html(it["local_aln"], it["pdb_aln"], wrap=wrap),
            }
        )

    # 放到 application/json，避免 <script> 直接内联 JSON 引发的解析问题
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
      document.getElementById('counter').textContent = `${idx+1} / ${items.length}`;
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
          opt.textContent = `${i+1}. ${it.title ?? ''}`;
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
<meta charset="utf-8"/>
<title>{html_mod.escape(title_text)}</title>
<style>{css}</style>
</head>
<body>
<h1 id="title">{html_mod.escape(title_text)}</h1>
<div class="toolbar">
  <button id="prev">Prev</button>
  <button id="next">Next</button>
  <span id="counter"></span>
  <select id="sel"></select>
</div>
<div class="hint">
策略：先找 “[auth X]”（严格大小写）；找不到就找 “Chain(s) X”（严格大小写）；再不行用 chainId fallback。红色：替换；灰色 '-'：gap；蓝色：与 gap 对应另一侧字符。
</div>
<div id="meta" class="meta"></div>
<div id="stats" class="stats"></div>
<div id="content"></div>

<script type="application/json" id="data">{data_json}</script>
<script>{js}</script>
</body>
</html>
"""

# =========================
# 并行任务：处理一条链（H or L）
# =========================

def process_chain(
    row_index: int,
    pdb_id_raw: str,
    chain_type: str,
    chain_id: str,
    local_seq: str,
) -> Tuple[dict, Optional[dict]]:
    session = requests.Session()

    pdb_entry = (pdb_id_raw or "").strip()[:4].upper()
    chain_id = (chain_id or "").strip()  # keep case
    local_seq_norm = normalize_seq(local_seq or "")

    if not chain_id or not local_seq_norm:
        csv_row = {
            "row_index": row_index,
            "pdb_id": pdb_id_raw,
            "pdb_entry": pdb_entry,
            "chain_type": chain_type,
            "chain_id": chain_id,
            "local_sequence_len": len(local_seq_norm),
            "PDB_sequence_len": "",
            "status": "SKIP",
            "identity_no_gaps": "",
            "aligned_len": "",
            "matches": "",
            "mismatches": "",
            "gaps_in_local": "",
            "gaps_in_PDB": "",
            "note": "missing_local_chain_or_sequence",
        }
        return csv_row, None

    fetched = fetch_chain_sequence(pdb_entry, chain_id, session=session)
    if not fetched.found or not fetched.sequence:
        csv_row = {
            "row_index": row_index,
            "pdb_id": pdb_id_raw,
            "pdb_entry": pdb_entry,
            "chain_type": chain_type,
            "chain_id": chain_id,
            "local_sequence_len": len(local_seq_norm),
            "PDB_sequence_len": "",
            "status": "NO_PDB",
            "identity_no_gaps": "",
            "aligned_len": "",
            "matches": "",
            "mismatches": "",
            "gaps_in_local": "",
            "gaps_in_PDB": "",
            "note": fetched.reason,
        }
        return csv_row, None

    pdb_seq = fetched.sequence

    local_aln, pdb_aln = needleman_wunsch(local_seq_norm, pdb_seq, MATCH_SCORE, MISMATCH_SCORE, GAP_SCORE)

    matches = 0
    mismatches = 0
    for x, y in zip(local_aln, pdb_aln):
        if x == "-" or y == "-":
            continue
        if x == y:
            matches += 1
        else:
            mismatches += 1

    denom = matches + mismatches
    identity = (matches / denom) if denom else 1.0

    gaps_in_local = local_aln.count("-")
    gaps_in_pdb = pdb_aln.count("-")

    status = "EXACT" if (local_seq_norm == pdb_seq) else "DIFFERENT"

    csv_row = {
        "row_index": row_index,
        "pdb_id": pdb_id_raw,
        "pdb_entry": pdb_entry,
        "chain_type": chain_type,
        "chain_id": chain_id,
        "local_sequence_len": len(local_seq_norm),
        "PDB_sequence_len": len(pdb_seq),
        "status": status,
        "identity_no_gaps": f"{identity:.4f}",
        "aligned_len": len(local_aln),
        "matches": matches,
        "mismatches": mismatches,
        "gaps_in_local": gaps_in_local,
        "gaps_in_PDB": gaps_in_pdb,
        "note": fetched.reason,
    }

    if status != "DIFFERENT":
        return csv_row, None

    title = f"{pdb_entry} {chain_type} chain {chain_id} | row {row_index}"
    meta = f"pdb_id={pdb_id_raw}"
    stats = (
        f"identity(no-gaps)={identity:.4f}, matches={matches}, mismatches={mismatches}, "
        f"gaps_in_local={gaps_in_local}, gaps_in_PDB={gaps_in_pdb}, aligned_len={len(local_aln)}"
    )

    viewer_item = {
        "title": title,
        "meta": meta,
        "stats": stats,
        "local_aln": local_aln,
        "pdb_aln": pdb_aln,
        "identity_no_gaps": identity,
    }

    return csv_row, viewer_item


# =========================
# 主程序
# =========================

def main() -> int:
    base = Path(__file__).resolve().parent
    in_path = base / INPUT_CSV
    if not in_path.exists():
        print(f"[ERROR] 找不到输入CSV：{in_path}", file=sys.stderr)
        return 2

    align_dir = base / ALIGN_DIR
    align_dir.mkdir(parents=True, exist_ok=True)

    with in_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        required = {"pdb_id", "Hchain", "VH_sequence", "Lchain", "VL_sequence"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            print(f"[ERROR] CSV缺少列：{sorted(missing)}；当前表头：{reader.fieldnames}", file=sys.stderr)
            return 3
        rows = list(reader)

    csv_rows: List[dict] = []
    viewer_items: List[dict] = []

    futures = []
    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        for idx, row in enumerate(rows, start=1):
            pdb_id_raw = (row.get("pdb_id") or "").strip()

            futures.append(
                ex.submit(
                    process_chain,
                    idx,
                    pdb_id_raw,
                    "H",
                    (row.get("Hchain") or "").strip(),
                    row.get("VH_sequence") or "",
                )
            )
            futures.append(
                ex.submit(
                    process_chain,
                    idx,
                    pdb_id_raw,
                    "L",
                    (row.get("Lchain") or "").strip(),
                    row.get("VL_sequence") or "",
                )
            )

        for fut in tqdm(as_completed(futures), total=len(futures), desc=f"Processing (threads={THREADS})", unit="chain"):
            csv_row, viewer_item = fut.result()
            csv_rows.append(csv_row)
            if viewer_item is not None:
                viewer_items.append(viewer_item)

    def _sort_key(d: dict) -> Tuple[int, str]:
        return int(d.get("row_index") or 0), str(d.get("chain_type") or "")

    csv_rows.sort(key=_sort_key)

    def _viewer_rownum(it: dict) -> int:
        m = re.search(r"\|\s*row\s+(\d+)", it.get("title", ""))
        return int(m.group(1)) if m else 0

    viewer_items.sort(key=lambda it: (_viewer_rownum(it), it.get("title", "")))

    out_csv = base / OUTPUT_CSV
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "row_index",
            "pdb_id",
            "pdb_entry",
            "chain_type",
            "chain_id",
            "local_sequence_len",
            "PDB_sequence_len",
            "status",
            "identity_no_gaps",
            "aligned_len",
            "matches",
            "mismatches",
            "gaps_in_local",
            "gaps_in_PDB",
            "note",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in csv_rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

    # DIFFERENT only
    main_items = viewer_items
    (align_dir / INDEX_HTML).write_text(
        make_index_html(main_items, wrap=WRAP, title_text="chain_check_1 (DIFFERENT only)"),
        encoding="utf-8",
    )

    # focus: DIFFERENT & identity != 1
    focus_items = [it for it in viewer_items if abs(float(it["identity_no_gaps"]) - 1.0) > 1e-12]
    (align_dir / FOCUS_HTML).write_text(
        make_index_html(focus_items, wrap=WRAP, title_text="chain_check_1_focus (DIFFERENT & identity_no_gaps != 1)"),
        encoding="utf-8",
    )

    print(f"\n[DONE] 输出完成：{out_csv}")
    print(f"[DONE] DIFFERENT 查看器：{align_dir / INDEX_HTML}")
    print(f"[DONE] FOCUS 查看器：{align_dir / FOCUS_HTML}")
    print(f"[INFO] DIFFERENT 条目数：{len(main_items)}")
    print(f"[INFO] FOCUS 条目数：{len(focus_items)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
