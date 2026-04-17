#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run: python3 dedup_tmgraph_090_nocli_filesize.py

Input rows (tab-separated preferred; comma also supported):
  1A00-assembly1    2H8D-assembly1    0.9742

Dedup strategy:
  - Build undirected graph with edges where tmmax >= THRESHOLD.
  - Split into connected components.
  - In each component, pick representatives via greedy dominating-set / set-cover:
      choose node that covers most uncovered nodes (itself+neighbors),
      tie-break by larger structure file size (.cif.gz) under STRUCT_DIR.

Outputs (all in OUTDIR):
  - representatives.tsv
  - rep_to_members.tsv
  - removed.tsv
  - summary.json
"""

import glob
import json
import os
import sys
from collections import defaultdict, deque

# =========================
# CONFIG (edit here if needed)
# =========================
INPUT_GLOB = "final_result_0.6_260409_*.csv"
THRESHOLD = 0.9

# Base directory containing structure files like "1A00-assembly1.cif.gz"
STRUCT_DIR = "/home/zhc/dalab/ab_ag_structures/SAbDab/PPI_DB/structures"

OUTDIR = "dedup_tmmax_0.9_out"
# =========================


def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)


def structure_path(structure_id: str) -> str:
    return os.path.join(STRUCT_DIR, structure_id + ".cif.gz")


def file_size_bytes(structure_id: str) -> int:
    """
    Return file size in bytes for STRUCT_DIR/<id>.cif.gz.
    If missing/unreadable, return -1 (worst).
    """
    p = structure_path(structure_id)
    try:
        return os.path.getsize(p)
    except OSError:
        return -1


def read_edges_build_graph(files, threshold: float):
    """
    Build adjacency list for edges >= threshold.
    Returns: adj dict node->set(neighbors), nodes set, stats dict
    """
    adj = defaultdict(set)
    nodes = set()
    total_rows = 0
    kept_edges = 0

    for fp in files:
        with open(fp, "r", encoding="utf-8") as f:
            for line in f:
                total_rows += 1
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 3:
                    parts = line.split(",")
                    if len(parts) < 3:
                        continue

                a, b, s = parts[0].strip(), parts[1].strip(), parts[2].strip()
                try:
                    tm = float(s)
                except ValueError:
                    continue

                nodes.add(a)
                nodes.add(b)

                if tm >= threshold and a != b:
                    adj[a].add(b)
                    adj[b].add(a)
                    kept_edges += 1

    stats = {
        "total_rows": total_rows,
        "kept_edges": kept_edges,
        "nodes": len(nodes),
        "files": files,
        "threshold": threshold,
    }
    return adj, nodes, stats


def connected_components(adj: dict, nodes: set):
    """
    Components on an undirected graph with possible isolated nodes.
    Returns list of lists.
    """
    seen = set()
    comps = []

    for n in nodes:
        if n in seen:
            continue
        q = deque([n])
        seen.add(n)
        comp = []
        while q:
            u = q.popleft()
            comp.append(u)
            for v in adj.get(u, ()):
                if v not in seen:
                    seen.add(v)
                    q.append(v)
        comps.append(comp)

    return comps


def quality_key_filesize(node: str, size_cache: dict):
    """
    Tie-break key:
      - Prefer larger file size (descending), so we use negative size in key.
      - Secondary: node string stable.
    Smaller key = better.
    """
    sz = size_cache.get(node)
    if sz is None:
        sz = file_size_bytes(node)
        size_cache[node] = sz
    return (-sz, node)


def greedy_dominating_set(component_nodes, adj, size_cache):
    """
    Greedy dominating-set / set-cover:
      U = all nodes in component.
      Cover(i) = {i} ∪ neighbors(i) (restricted to component)

    Choose node with maximum gain = |Cover(i) ∩ U|.
    Tie-break: larger .cif.gz file size.
    """
    comp_set = set(component_nodes)
    uncovered = set(component_nodes)

    # Precompute closed neighborhoods restricted to component
    neigh = {}
    for u in component_nodes:
        s = set(adj.get(u, ()))
        s.intersection_update(comp_set)
        s.add(u)
        neigh[u] = s

    reps = []
    covered_by = {}  # node -> rep

    while uncovered:
        best_u = None
        best_gain = -1
        best_q = None

        # candidates: uncovered only (fast)
        for u in uncovered:
            gain = len(neigh[u] & uncovered)
            if gain > best_gain:
                best_u = u
                best_gain = gain
                best_q = quality_key_filesize(u, size_cache)
            elif gain == best_gain and gain > 0:
                q = quality_key_filesize(u, size_cache)
                if q < best_q:
                    best_u = u
                    best_q = q

        if best_u is None:
            best_u = min(uncovered, key=lambda x: quality_key_filesize(x, size_cache))

        reps.append(best_u)
        newly = neigh[best_u] & uncovered
        for n in newly:
            covered_by[n] = best_u
        uncovered -= newly

    rep_to_members = defaultdict(list)
    for n, r in covered_by.items():
        rep_to_members[r].append(n)
    for r in rep_to_members:
        rep_to_members[r].sort()

    reps = sorted(set(reps), key=lambda x: quality_key_filesize(x, size_cache))
    return reps, rep_to_members


def main():
    ensure_outdir(OUTDIR)

    files = sorted(glob.glob(INPUT_GLOB))
    if not files:
        print(f"[ERROR] No files matched: {INPUT_GLOB}", file=sys.stderr)
        sys.exit(2)

    if not os.path.isdir(STRUCT_DIR):
        print(f"[ERROR] STRUCT_DIR not found: {STRUCT_DIR}", file=sys.stderr)
        sys.exit(2)

    print(f"[INFO] Reading {len(files)} files, threshold={THRESHOLD} ...", file=sys.stderr)
    adj, nodes, stats = read_edges_build_graph(files, THRESHOLD)
    print(f"[INFO] Total rows={stats['total_rows']}, kept_edges={stats['kept_edges']}, nodes={stats['nodes']}",
          file=sys.stderr)

    print("[INFO] Computing connected components ...", file=sys.stderr)
    comps = connected_components(adj, nodes)
    comps.sort(key=len, reverse=True)
    largest = len(comps[0]) if comps else 0
    print(f"[INFO] Components={len(comps)} largest={largest}", file=sys.stderr)

    # Deduplicate
    representatives = []
    rep_to_members_all = {}
    size_cache = {}  # structure_id -> filesize

    print("[INFO] Selecting representatives per component (tie-break by file size) ...", file=sys.stderr)
    for idx, comp in enumerate(comps, start=1):
        reps, rep_to_members = greedy_dominating_set(comp, adj, size_cache)
        representatives.extend(reps)
        rep_to_members_all.update(rep_to_members)

        if idx <= 5:
            print(f"[INFO] Component #{idx}: size={len(comp)} reps={len(reps)}", file=sys.stderr)

    representatives = sorted(set(representatives), key=lambda x: quality_key_filesize(x, size_cache))

    # Write outputs
    rep_fp = os.path.join(OUTDIR, "representatives.tsv")
    with open(rep_fp, "w", encoding="utf-8") as f:
        f.write("structure_id\tfilesize_bytes\n")
        for r in representatives:
            sz = size_cache.get(r)
            if sz is None:
                sz = file_size_bytes(r)
                size_cache[r] = sz
            f.write(f"{r}\t{sz}\n")

    map_fp = os.path.join(OUTDIR, "rep_to_members.tsv")
    with open(map_fp, "w", encoding="utf-8") as f:
        f.write("representative\tmember\tmember_filesize_bytes\n")
        for rep in representatives:
            members = rep_to_members_all.get(rep, [])
            for m in members:
                sz = size_cache.get(m)
                if sz is None:
                    sz = file_size_bytes(m)
                    size_cache[m] = sz
                f.write(f"{rep}\t{m}\t{sz}\n")

    reps_set = set(representatives)
    removed = sorted([n for n in nodes if n not in reps_set])
    rm_fp = os.path.join(OUTDIR, "removed.tsv")
    with open(rm_fp, "w", encoding="utf-8") as f:
        f.write("structure_id\n")
        for n in removed:
            f.write(f"{n}\n")

    summary_fp = os.path.join(OUTDIR, "summary.json")
    summary = {
        "threshold": THRESHOLD,
        "input_glob": INPUT_GLOB,
        "input_files": files,
        "struct_dir": STRUCT_DIR,
        "total_rows": stats["total_rows"],
        "kept_edges": stats["kept_edges"],
        "node_count": len(nodes),
        "component_count": len(comps),
        "largest_component": largest,
        "representative_count": len(representatives),
        "removed_count": len(removed),
        "tie_break": "prefer_larger_filesize_bytes",
        "missing_file_size_value": -1
    }
    with open(summary_fp, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"[DONE] Output written to: {OUTDIR}", file=sys.stderr)
    print(f"[DONE] representatives: {rep_fp}", file=sys.stderr)
    print(f"[DONE] rep_to_members:  {map_fp}", file=sys.stderr)
    print(f"[DONE] removed:        {rm_fp}", file=sys.stderr)
    print(f"[DONE] summary:        {summary_fp}", file=sys.stderr)


if __name__ == "__main__":
    main()