#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Parallel + resume + rich dashboard (low flicker) for extracting protein-protein chain pairs
from assembly mmCIF (.cif.gz) using _atom_site (gemmi 0.7.3 compatible).

Your requirements implemented:
- Contact definition:
    representative point = Cβ, Gly uses Cα, fallback to Cα if Cβ missing
    residue contact if distance < 8 Å
- Chain-chain interaction if:
    |contact_residues(A)| + |contact_residues(B)| > 5   (sum version)
- Output: pair CIF files named like 1A0H_1.cif, 1A0H_2.cif, ...
- Output CSV: protein_pairs.csv with columns:
    file_name, pair_chain, contact_residue
- Keep non-canonical amino acids / modified residues:
    export keeps ATOM always; keeps HETATM only if that residue has backbone atoms (N/CA/C/O) in that chain
    drops waters + nucleic acids
- Parallel processing:
    default 10 worker processes
- Resume:
    skip PDBIDs that already have outputs in out_dir (PDBID_*.cif)
    BUT the last PDBID in alphabetical order is ALWAYS reprocessed (delete old outputs + remove its CSV rows)
- UI:
    rich Live dashboard with:
      - worker table (one row per worker slot showing current structure)
      - overall progress bar
    with reduced flicker (fixed layout/panels, only update renderables)

Dependencies:
  pip install gemmi numpy scipy rich
"""

import os
import re
import csv
import glob
import argparse
from itertools import combinations
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

import numpy as np
import gemmi
from scipy.spatial import cKDTree

from rich.console import Console
from rich.table import Table
from rich.live import Live
from rich.layout import Layout
from rich.panel import Panel
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn


# ---------- constants ----------
WATER_RESN = {"HOH", "WAT", "DOD"}

NA_RESN = {
    "A", "U", "G", "C", "I",
    "DA", "DT", "DG", "DC", "DI",
    "RA", "RU", "RG", "RC", "RI",
    "PSU",
}

NA_ATOM_HINTS = {
    "P", "OP1", "OP2", "OP3",
    "O3'", "O4'", "O5'",
    "C1'", "C2'", "C3'", "C4'", "C5'",
    "O2'",
}

PROT_BACKBONE_ATOMS = {"N", "CA", "C", "O"}


# ---------- helpers ----------
def pdb_id_from_filename(path: str) -> str:
    base = os.path.basename(path)
    m = re.match(r"^([A-Za-z0-9]{4})-assembly\d+\.cif\.gz$", base)
    if m:
        return m.group(1).upper()
    return os.path.splitext(os.path.splitext(base)[0])[0].upper()


def norm_chain(x: str) -> str:
    return str(x).strip()


def norm_resn(x: str) -> str:
    return str(x).strip().upper()


def norm_atom(x: str) -> str:
    return str(x).strip().upper()


def is_model1(x) -> bool:
    s = str(x).strip()
    return s == "1" or s == "1.0"


def safe_float(x) -> float:
    s = str(x).strip()
    if s in (".", "?", ""):
        return float("nan")
    return float(s)


# ---------- mmCIF atom_site loading (gemmi 0.7.3 compatible) ----------
def load_atom_site_table_and_idx(cif_gz_path: str):
    doc = gemmi.cif.read(cif_gz_path)
    block = doc.sole_block()

    # Try the full column set first (matches your 1I5K header list)
    rich = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "Cartn_x_esd",
        "Cartn_y_esd",
        "Cartn_z_esd",
        "occupancy_esd",
        "B_iso_or_equiv_esd",
        "pdbx_formal_charge",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    ]
    table = block.find("_atom_site.", rich)
    if table is not None:
        return table, {name: i for i, name in enumerate(rich)}

    # Fallback minimal set required by our algorithm and export
    minimal = [
        "group_PDB",
        "type_symbol",
        "label_atom_id",
        "label_comp_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    ]
    table = block.find("_atom_site.", minimal)
    if table is None:
        raise RuntimeError("Failed to read _atom_site table (missing required columns)")
    return table, {name: i for i, name in enumerate(minimal)}


# ---------- chain classification ----------
def build_chain_signatures(table, idx):
    chain_total_atoms = defaultdict(int)
    chain_bb_hits = defaultdict(int)
    chain_na_atom_hits = defaultdict(int)
    chain_resn_counts = defaultdict(lambda: defaultdict(int))

    i_chain = idx["auth_asym_id"]
    i_resn = idx["auth_comp_id"]
    i_atom = idx["auth_atom_id"]
    i_model = idx["pdbx_PDB_model_num"]

    for row in table:
        if not is_model1(row[i_model]):
            continue

        chain = norm_chain(row[i_chain])
        resn = norm_resn(row[i_resn])
        atom = norm_atom(row[i_atom])

        chain_total_atoms[chain] += 1
        chain_resn_counts[chain][resn] += 1

        if atom in PROT_BACKBONE_ATOMS:
            chain_bb_hits[chain] += 1
        if atom in NA_ATOM_HINTS:
            chain_na_atom_hits[chain] += 1

    return chain_total_atoms, chain_bb_hits, chain_na_atom_hits, chain_resn_counts


def classify_chain_is_protein(chain_id: str,
                              chain_total_atoms,
                              chain_bb_hits,
                              chain_na_atom_hits,
                              chain_resn_counts,
                              min_backbone_ratio=0.03,
                              min_backbone_atoms=8):
    """
    Protein chain detection without peptide length threshold.
    """
    total = chain_total_atoms.get(chain_id, 0)
    if total == 0:
        return False

    bb = chain_bb_hits.get(chain_id, 0)
    na_atom = chain_na_atom_hits.get(chain_id, 0)
    resn_map = chain_resn_counts.get(chain_id, {})

    # reject nucleic-acid-like chains
    na_res_atoms = sum(resn_map.get(r, 0) for r in NA_RESN)
    if na_res_atoms / total > 0.2:
        return False
    if na_atom / total > 0.05:
        return False

    # backbone enrichment
    if bb < min_backbone_atoms:
        return False
    if (bb / total) < min_backbone_ratio:
        return False

    return True


# ---------- residue representative points (CB/CA) ----------
def build_residue_rep_points_and_backbone_presence(table, idx, chain_id: str):
    i_chain = idx["auth_asym_id"]
    i_seq = idx["auth_seq_id"]
    i_resn = idx["auth_comp_id"]
    i_atom = idx["auth_atom_id"]
    i_x = idx["Cartn_x"]
    i_y = idx["Cartn_y"]
    i_z = idx["Cartn_z"]
    i_model = idx["pdbx_PDB_model_num"]

    ca = {}
    cb = {}
    backbone_presence = defaultdict(bool)

    for row in table:
        if not is_model1(row[i_model]):
            continue
        chain = norm_chain(row[i_chain])
        if chain != chain_id:
            continue

        seq = str(row[i_seq]).strip()
        resn = norm_resn(row[i_resn])
        atom = norm_atom(row[i_atom])
        key = (seq, resn)

        if atom in PROT_BACKBONE_ATOMS:
            backbone_presence[key] = True

        x = safe_float(row[i_x]); y = safe_float(row[i_y]); z = safe_float(row[i_z])
        if not np.isfinite(x + y + z):
            continue

        if atom == "CA":
            ca[key] = (x, y, z)
        elif atom == "CB":
            cb[key] = (x, y, z)

    keys = set(ca.keys()) | set(cb.keys())
    res_uids = []
    coords = []

    for (seq, resn) in keys:
        if resn == "GLY":
            if (seq, resn) not in ca:
                continue
            pt = ca[(seq, resn)]
        else:
            if (seq, resn) in cb:
                pt = cb[(seq, resn)]
            elif (seq, resn) in ca:
                pt = ca[(seq, resn)]
            else:
                continue

        uid = f"{chain_id}:{seq}:{resn}"
        res_uids.append(uid)
        coords.append(pt)

    if not coords:
        return [], np.zeros((0, 3), dtype=np.float32), backbone_presence
    return res_uids, np.asarray(coords, dtype=np.float32), backbone_presence


def compute_contacts(res_uids_A, coords_A, res_uids_B, coords_B, cutoff=8.0):
    if coords_A.shape[0] == 0 or coords_B.shape[0] == 0:
        return set(), set()

    treeB = cKDTree(coords_B)
    neighbors = treeB.query_ball_point(coords_A, r=cutoff)

    contactA, contactB = set(), set()
    for i, nb in enumerate(neighbors):
        if not nb:
            continue
        contactA.add(res_uids_A[i])
        for j in nb:
            contactB.add(res_uids_B[j])
    return contactA, contactB


# ---------- export pair CIF ----------
def write_pair_cif_from_atom_site(table, idx, out_path: str, chainA: str, chainB: str, backbone_presence_by_chain):
    keep_chains = {chainA, chainB}

    # output headers in the same order as loaded table
    tags_in_order = [None] * len(idx)
    for tag, i in idx.items():
        tags_in_order[i] = tag
    headers = [f"_atom_site.{t}" for t in tags_in_order]

    i_group = idx["group_PDB"]
    i_chain = idx["auth_asym_id"]
    i_seq = idx["auth_seq_id"]
    i_resn = idx["auth_comp_id"]
    i_model = idx["pdbx_PDB_model_num"]

    kept_rows = 0
    with open(out_path, "wt", encoding="utf-8") as f:
        f.write("data_pair\n#\nloop_\n")
        for h in headers:
            f.write(h + "\n")

        for row in table:
            if not is_model1(row[i_model]):
                continue

            group = str(row[i_group]).strip().upper()
            chain = norm_chain(row[i_chain])
            if chain not in keep_chains:
                continue

            resn = norm_resn(row[i_resn])
            if resn in WATER_RESN or resn in NA_RESN:
                continue

            seq = str(row[i_seq]).strip()
            key = (seq, resn)

            if group == "HETATM":
                # keep only HETATM residues that look like protein residues (have backbone atoms)
                if not backbone_presence_by_chain[chain].get(key, False):
                    continue

            f.write(" ".join(str(x) for x in row) + "\n")
            kept_rows += 1

    if kept_rows == 0:
        raise RuntimeError("pair CIF export kept 0 rows (after filtering)")


# ---------- worker ----------
def process_one_file_worker(payload):
    """
    Worker: process one cif.gz and write pair CIFs.
    Returns: (pdb_id, records, err_string)
      records: list of (out_name, chainA, chainB, contact_sum)
    """
    fp, out_dir, cutoff, min_contact_res = payload
    pdb_id = pdb_id_from_filename(fp)

    try:
        table, idx = load_atom_site_table_and_idx(fp)
        chain_total, chain_bb, chain_na, chain_resn_counts = build_chain_signatures(table, idx)

        protein_chains = []
        for chain_id in chain_total.keys():
            if classify_chain_is_protein(chain_id, chain_total, chain_bb, chain_na, chain_resn_counts):
                protein_chains.append(chain_id)
        protein_chains.sort()

        rep = {}
        bb_presence = {}
        for ch in protein_chains:
            res_uids, coords, backbone_presence = build_residue_rep_points_and_backbone_presence(table, idx, ch)
            rep[ch] = (res_uids, coords)
            bb_presence[ch] = backbone_presence

        pair_idx = 0
        records = []

        for c1, c2 in combinations(protein_chains, 2):
            resA, coA = rep[c1]
            resB, coB = rep[c2]
            contactA, contactB = compute_contacts(resA, coA, resB, coB, cutoff=cutoff)
            contact_sum = len(contactA) + len(contactB)

            if contact_sum > min_contact_res:
                pair_idx += 1
                out_name = f"{pdb_id}_{pair_idx}"
                out_path = os.path.join(out_dir, f"{out_name}.cif")

                write_pair_cif_from_atom_site(table, idx, out_path, c1, c2, bb_presence)
                records.append((out_name, c1, c2, contact_sum))

        return pdb_id, records, ""
    except Exception as e:
        return pdb_id, [], str(e)


# ---------- resume ----------
def existing_pdbids_in_outdir(out_dir: str):
    ids = set()
    for fp in glob.glob(os.path.join(out_dir, "*_*.cif")):
        base = os.path.basename(fp)
        m = re.match(r"^([A-Za-z0-9]{4})_\d+\.cif$", base)
        if m:
            ids.add(m.group(1).upper())
    return ids


def delete_existing_outputs_for_pdbid(out_dir: str, pdb_id: str):
    for fp in glob.glob(os.path.join(out_dir, f"{pdb_id}_*.cif")):
        try:
            os.remove(fp)
        except OSError:
            pass


def rewrite_csv_without_pdbid(csv_path: str, pdb_id: str):
    """
    Remove rows whose file_name starts with f"{pdb_id}_".
    Used when force reprocessing last ID.
    """
    if not os.path.exists(csv_path):
        return
    tmp = csv_path + ".tmp"
    with open(csv_path, "rt", newline="", encoding="utf-8") as fin, \
         open(tmp, "wt", newline="", encoding="utf-8") as fout:
        r = csv.reader(fin)
        w = csv.writer(fout)
        header = next(r, None)
        if header:
            w.writerow(header)
        for row in r:
            if not row:
                continue
            if row[0].startswith(pdb_id + "_"):
                continue
            w.writerow(row)
    os.replace(tmp, csv_path)


# ---------- rich UI ----------
def render_worker_table(slots):
    t = Table(expand=True)
    t.add_column("Slot", justify="right", width=6)
    t.add_column("Status", width=10)
    t.add_column("PDBID", width=8)
    t.add_column("File", overflow="fold")

    for i, s in enumerate(slots):
        t.add_row(str(i), s.get("status", "IDLE"), s.get("pdbid", ""), s.get("file", ""))
    return t


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", default="/home/zhc/dalab/ab_ag_structures/SAbDab/PPI_DB/PPI_dataset_260416/assembly_structures")
    ap.add_argument("--out_dir", default="/home/zhc/dalab/ab_ag_structures/SAbDab/PPI_DB/PPI_dataset_260416/protein_pair_dataset")
    ap.add_argument("--csv_path", default=None, help="Default: <out_dir>/protein_pairs.csv")
    ap.add_argument("--cutoff", type=float, default=8.0)
    ap.add_argument("--min_contact_res", type=int, default=5)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--overwrite_csv", action="store_true")
    ap.add_argument("--resume", action="store_true", default=True)
    ap.add_argument("--refresh_hz", type=float, default=3.0, help="UI refresh rate (Hz). Lower = less flicker.")
    args = ap.parse_args()

    console = Console()
    os.makedirs(args.out_dir, exist_ok=True)
    csv_path = args.csv_path or os.path.join(args.out_dir, "protein_pairs.csv")

    cif_files = sorted(
        os.path.join(args.in_dir, fn)
        for fn in os.listdir(args.in_dir)
        if fn.endswith(".cif.gz")
    )
    if not cif_files:
        raise SystemExit(f"No *.cif.gz found in: {args.in_dir}")

    file_items = [(pdb_id_from_filename(fp), fp) for fp in cif_files]
    last_pdbid = file_items[-1][0]  # alphabetical by filename sort

    # resume: infer done IDs from out_dir (unless overwrite)
    done_ids = set()
    if args.resume and (not args.overwrite_csv):
        done_ids = existing_pdbids_in_outdir(args.out_dir)

    # force reprocess last id
    if last_pdbid in done_ids:
        delete_existing_outputs_for_pdbid(args.out_dir, last_pdbid)
        done_ids.discard(last_pdbid)

    if (not args.overwrite_csv) and os.path.exists(csv_path):
        rewrite_csv_without_pdbid(csv_path, last_pdbid)

    tasks = []
    for pdbid, fp in file_items:
        if args.resume and (pdbid in done_ids):
            continue
        tasks.append((fp, args.out_dir, args.cutoff, args.min_contact_res))

    if not tasks:
        console.print("Nothing to do.")
        console.print(f"Last ID (forced reprocess each run): {last_pdbid}")
        return

    # CSV open (main process writes, avoids concurrent write corruption)
    mode = "wt" if args.overwrite_csv or (not os.path.exists(csv_path)) else "at"
    need_header = (mode == "wt")

    progress = Progress(
        TextColumn("[bold]Total[/bold]"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    )
    total_task = progress.add_task("all", total=len(tasks))

    # slots for worker table
    slots = [{"status": "IDLE", "pdbid": "", "file": ""} for _ in range(args.workers)]

    # ---- Build a fixed layout/panels ONCE to avoid flickering ----
    layout = Layout()
    layout.split_column(
        Layout(name="top", ratio=3),
        Layout(name="bottom", ratio=1),
    )
    workers_panel = Panel(render_worker_table(slots), title="Workers (current file)", border_style="cyan")
    progress_panel = Panel(progress, title="Overall progress", border_style="green")
    layout["top"].update(workers_panel)
    layout["bottom"].update(progress_panel)

    total_pairs = 0

    with open(csv_path, mode, newline="", encoding="utf-8") as fcsv:
        w = csv.writer(fcsv)
        if need_header:
            w.writerow(["file_name", "pair_chain", "contact_residue"])
            fcsv.flush()

        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            task_iter = iter(tasks)
            pending = []  # list of (future, slot_id)

            # submit initial tasks to fill slots
            for slot_id in range(args.workers):
                try:
                    payload = next(task_iter)
                except StopIteration:
                    break
                fp = payload[0]
                pdbid = pdb_id_from_filename(fp)
                slots[slot_id] = {"status": "RUN", "pdbid": pdbid, "file": os.path.basename(fp)}
                fut = ex.submit(process_one_file_worker, payload)
                pending.append((fut, slot_id))

            # Live dashboard (reduced flicker: fixed panels/layout)
            with Live(layout, console=console, refresh_per_second=float(args.refresh_hz), transient=False) as live:
                while pending:
                    # update only the panel contents
                    workers_panel.renderable = render_worker_table(slots)
                    progress_panel.renderable = progress
                    live.refresh()

                    done, _ = wait([f for f, _ in pending], return_when=FIRST_COMPLETED)

                    new_pending = []
                    for fut, slot_id in pending:
                        if fut not in done:
                            new_pending.append((fut, slot_id))
                            continue

                        pdbid, records, err = fut.result()
                        if err:
                            console.print(f"[yellow][WARN][/yellow] {pdbid}: {err}")
                        else:
                            for out_name, c1, c2, contact_sum in records:
                                w.writerow([out_name, f"{c1},{c2}", contact_sum])
                                total_pairs += 1
                            fcsv.flush()

                        progress.advance(total_task, 1)

                        # slot becomes idle
                        slots[slot_id] = {"status": "IDLE", "pdbid": "", "file": ""}

                        # submit next task into freed slot
                        try:
                            payload = next(task_iter)
                            fp = payload[0]
                            npdbid = pdb_id_from_filename(fp)
                            slots[slot_id] = {"status": "RUN", "pdbid": npdbid, "file": os.path.basename(fp)}
                            nfut = ex.submit(process_one_file_worker, payload)
                            new_pending.append((nfut, slot_id))
                        except StopIteration:
                            pass

                    pending = new_pending

                    # refresh after state changes
                    workers_panel.renderable = render_worker_table(slots)
                    progress_panel.renderable = progress
                    live.refresh()

    console.print(f"\nDone. Extracted {total_pairs} interacting pairs.")
    console.print(f"Pairs CIF dir: {args.out_dir}")
    console.print(f"CSV: {csv_path}")
    console.print(f"Last ID (forced reprocess each run): {last_pdbid}")


if __name__ == "__main__":
    main()