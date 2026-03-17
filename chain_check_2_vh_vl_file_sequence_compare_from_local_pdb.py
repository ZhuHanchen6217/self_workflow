import pandas as pd
from Bio import SeqIO, Align
from Bio.Align import PairwiseAligner
import os

# Load final structures data
data = pd.read_csv('final_structures_sift_v1_251203.csv')

# Initialize aligner for Needleman-Wunsch
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.gap_score = -2

results = []

for index, row in data.iterrows():
    pdb_id = row['pdb_id']
    Hchain = row['Hchain']
    Lchain = row['Lchain']
    
    # Read PDB file
    pdb_file_path = f'vh_vl_dataset/{pdb_id}.pdb'
    
    if os.path.exists(pdb_file_path):
        sequence_h = None
        sequence_l = None
        
        # Parse PDB file to get sequences
        with open(pdb_file_path, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                # Prefer SEQRES
                if 'SEQRES' in record.annotations:
                    sequence_h = str(record.annotations['SEQRES'][0])
                    sequence_l = str(record.annotations['SEQRES'][1])
                # Fallback to ATOM
                if 'ATOM' in record.annotations:
                    atom_seq = str(record.annotations['ATOM'])
                    if sequence_h is None:
                        sequence_h = atom_seq
                    if sequence_l is None:
                        sequence_l = atom_seq
        
        # Align sequences
        align_h = aligner.align(sequence_h, Hchain)
        align_l = aligner.align(sequence_l, Lchain)
        
        # Append results for CSV
        results.append({
            'pdb_id': pdb_id,
            'Hchain_score': align_h[0].score,
            'Lchain_score': align_l[0].score,
            'Hchain_alignment': str(align_h[0]),
            'Lchain_alignment': str(align_l[0]),
        })

# Output results to CSV
results_df = pd.DataFrame(results)
results_df.to_csv('vh_vl_file_vs_local_chain_report_2.csv', index=False)

# Generate HTML outputs (implementation needed for visualization)
# Sample code for generating HTML would go here...

# Fix JS template string bug
# Affected JS Code snippet:
# let example = `${idx+1}`;  # Correct usage instead of ${{idx+1}}