#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import sys
from collections import defaultdict, deque

import networkx as nx
import matplotlib.pyplot as plt

# =========================
# CONFIG
# =========================
INPUT_GLOB = "final_result_0.6_260409_*.csv"
THRESHOLD = 0.9

OUTDIR = "dedup_tmmax_0.9_out"
REPRESENTATIVES_TSV = os.path.join(OUTDIR, "representatives.tsv")

PLOT_DIR = os.path.join(OUTDIR, "component_plots_full")
DPI = 200

# Layout tuning
SEED = 42
NODE_SIZE_REP = 80
NODE_SIZE_NONREP = 10
EDGE_WIDTH = 0.2
EDGE_ALPHA = 0.15

# For very small components, edges can be drawn a bit darker
EDGE_ALPHA_SMALL = 0.35
EDGE_WIDTH_SMALL = 0.5
SMALL_COMP_N = 30
# =========================


def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


def read_representatives(fp: str) -> set:
    reps = set()
    if not os.path.exists(fp):
        print(f"[ERROR] Missing {fp}. Run dedup first.", file=sys.stderr)
        sys.exit(2)
    with open(fp, "r", encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            reps.add(line.split("\t")[0])
    return reps


def read_edges(files, threshold: float):
    """
    Return edges list (a,b) for tm>=threshold and set of nodes.
    Also store a weight dict for optional edge coloring.
    """
    nodes = set()
    edges = {}  # (min(a,b), max(a,b)) -> max_tm

    total_rows = 0
    kept_rows = 0

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

                nodes.add(a); nodes.add(b)
                if tm >= threshold and a != b:
                    kept_rows += 1
                    x, y = (a, b) if a <= b else (b, a)
                    prev = edges.get((x, y))
                    if prev is None or tm > prev:
                        edges[(x, y)] = tm

    print(f"[INFO] total_rows={total_rows} kept_rows(tm>={threshold})={kept_rows} unique_edges={len(edges)}",
          file=sys.stderr)
    return nodes, edges


def connected_components_from_adj(nodes, edges_dict):
    adj = defaultdict(set)
    for (a, b) in edges_dict.keys():
        adj[a].add(b)
        adj[b].add(a)

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

    comps.sort(key=len, reverse=True)
    return comps, adj


def build_component_graph(comp_nodes, adj):
    G = nx.Graph()
    G.add_nodes_from(comp_nodes)
    for u in comp_nodes:
        for v in adj.get(u, ()):
            if v in G:
                if u < v:
                    G.add_edge(u, v)
    return G


def plot_component(G, reps_set, outpath: str, title: str):
    n = G.number_of_nodes()

    # Layout: spring_layout is ok for <=~600 nodes; deterministic with seed
    pos = nx.spring_layout(G, seed=SEED, k=None)

    rep_nodes = [u for u in G.nodes() if u in reps_set]
    nonrep_nodes = [u for u in G.nodes() if u not in reps_set]

    # Edge style
    if n <= SMALL_COMP_N:
        edge_alpha = EDGE_ALPHA_SMALL
        edge_w = EDGE_WIDTH_SMALL
    else:
        edge_alpha = EDGE_ALPHA
        edge_w = EDGE_WIDTH

    plt.figure(figsize=(8, 8))
    plt.axis("off")
    plt.title(title, fontsize=10)

    # Draw edges first
    nx.draw_networkx_edges(G, pos, width=edge_w, alpha=edge_alpha, edge_color="#444444")

    # Draw non-reps
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=nonrep_nodes,
        node_size=NODE_SIZE_NONREP,
        node_color="#bdbdbd",
        linewidths=0
    )

    # Draw reps
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=rep_nodes,
        node_size=NODE_SIZE_REP,
        node_color="#d62728",
        linewidths=0.6,
        edgecolors="#7f0000"
    )

    # Label only reps
    labels = {u: u for u in rep_nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, font_color="#111111")

    plt.tight_layout()
    plt.savefig(outpath, dpi=DPI)
    plt.close()


def main():
    ensure_dir(PLOT_DIR)

    reps = read_representatives(REPRESENTATIVES_TSV)

    files = sorted(glob.glob(INPUT_GLOB))
    if not files:
        print(f"[ERROR] No files matched: {INPUT_GLOB}", file=sys.stderr)
        sys.exit(2)

    nodes, edges_dict = read_edges(files, THRESHOLD)
    comps, adj = connected_components_from_adj(nodes, edges_dict)

    print(f"[INFO] components={len(comps)} largest={len(comps[0]) if comps else 0}", file=sys.stderr)
    print(f"[INFO] writing plots to: {PLOT_DIR}", file=sys.stderr)

    # Plot each component
    for i, comp in enumerate(comps, start=1):
        rep_count = sum(1 for u in comp if u in reps)
        if rep_count == 1:
            continue  # skip trivial components with single representative

        G = build_component_graph(comp, adj)
        outpng = os.path.join(PLOT_DIR, f"component_{i:04d}_n{len(comp)}_reps{rep_count}.png")
        title = f"Component #{i}  n={len(comp)}  reps={rep_count}  (tmmax>={THRESHOLD})"
        plot_component(G, reps, outpng, title)

        if i % 100 == 0:
            print(f"[INFO] plotted {i}/{len(comps)}", file=sys.stderr)

    print("[DONE] all component plots generated.", file=sys.stderr)


if __name__ == "__main__":
    main()