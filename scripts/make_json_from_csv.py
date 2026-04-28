from argparse import ArgumentParser
import pandas as pd
import json
import os, sys
from pathlib import Path
from typing import List, Union

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-f', '--input_file', type=Path, required=True, help='input csv file')
    parser.add_argument('-o', '--output_file', type=Path, required=False, default = None, help='output json file')
    parser.add_argument('-t', '--target_sequences', type=str, required=False, default=None, nargs='+', help='Target sequence. By default, will try to find a column in the dataframe containing the target sequence')
    parser.add_argument('-s', '--seeds', type=int, nargs='+', help= 'Seeds to use. [ex: --seeds 0 10 20]')
    parser.add_argument('-r', '--rows', type=int, default=[0, 10], nargs=2,
                        metavar=('START', 'END'), help = 'Start/end rows to select from file [ex: --rows 0 100]')
    parser.add_argument('--precomputed_msa_dir', type=str, default=None)
    parser.add_argument('--pairing_db', type=str, default= None, help='ex: uniref100')
    return parser.parse_args()

def find_target_col(df, candidates=None):
    """
    Find the column in a DataFrame that likely contains target sequences.

    Args:
        df: pandas DataFrame
        candidates: optional list of column name patterns to search for.
                    Defaults to common target sequence column names.

    Returns:
        str: the matching column name

    Raises:
        ValueError: if no matching column is found
    """
    if candidates is None:
        candidates = [
            'target', 'target_seq', 'target_sequence',
            'sequence', 'seq', 'protein_seq', 'protein_sequence',
            'amino_acid', 'aa_seq', 'fasta'
        ]

    columns = df.columns.str.lower()

    # 1. Try exact match first
    for c in candidates:
        matches = [col for col, col_lower in zip(df.columns, columns) if col_lower == c]
        if matches:
            return matches[0]

    # 2. Try partial match (column contains candidate string)
    for c in candidates:
        matches = [col for col, col_lower in zip(df.columns, columns) if c in col_lower]
        if matches:
            return matches[0]

    raise ValueError(
        f"Could not find a target sequence column. "
        f"Columns available: {list(df.columns)}"
    )


def make_msa(precomputed_msa_dir=None, pairing_db='uniref100'):
    res = {'pairing_db': pairing_db}
    if precomputed_msa_dir:
        res['precomputed_msa_dir'] = precomputed_msa_dir
    return res

def make_protein_chain(sequence, id, count=1,
                       precomputed_msa_dir=None, pairing_db='uniref100'):
    return {
        "proteinChain": {
            "sequence": sequence,
            "id": [id],
            "count": count,
            "msa": make_msa(precomputed_msa_dir, pairing_db)
        }
    }

# Placeholder for if we add constraints to the folding
def make_contraints():
    pass

def make_entry(protein_chains, name, seeds=None, constraints=None):
    entry = {"name": name}
    if seeds:
        entry["modelSeeds"] = seeds
    entry["sequences"] = [protein_chains]

    if constraints:
        entry["constraints"] = constraints
    return entry


def make_entry_from_row(vh:str, target: Union[List, str], name: str, seeds=None,
                        precomputed_msa_dir=None, pairing_db=None, constraints=None):
    """
    Used to create a full entry based on iterating on vh sequences
    Args:
        vh:
        target:
        name:
        seeds:
        precomputed_msa_dir:
        pairing_db:
        constraints:

    Returns:
        entry
    """
    # Start with VH with no precomputed msa
    protein_chains = [make_protein_chain(vh, 'H', 1, precomputed_msa_dir=None, pairing_db=None)]
    if type(target)==list:
        for i, t in enumerate(target):
            protein_chains.append(make_protein_chain(t, f'T{i}', 1, precomputed_msa_dir, pairing_db))
    elif type(target)==str:
        protein_chains.append(make_protein_chain(target, 'T', 1, precomputed_msa_dir, pairing_db))
    return make_entry(protein_chains, name, seeds, constraints)


def main():
    args = parse_args()
    start, end = args.rows
    df = pd.read_csv(args.input_file).iloc[start:end]
    # If not target is provided will try to parse from the df
    if args.target_sequences is None:
        target_col = find_target_col(df)
        # assumes a single target for every sequence in the df
        target = df[target_col].unique()[0]
    else:
        target = args.target_sequences

    # If no output_file specified, will save in the same directory as input_file with filename replaced as json
    if not args.output_file:
        parent = args.input_file.parent
        basename = args.input_file.stem
        outfile = parent / (basename + '.json')
    else:
        outfile = args.out_file

    df['seq_id'] = [f'{args.input_file.stem}_id_{i:06}' for i in range(len(df))]
    entries = []
    for _, row in df[['vh', 'seq_id']].iterrows():
        entries.append(make_entry_from_row(row['vh'], target,
                                           name = row['seq_id'], seeds=args.seeds,
                                           precomputed_msa_dir=args.precomputed_msa_dir,
                                           pairing_db=args.pairing_db, constraints=None))

    with open(outfile, 'w') as f:
        json.dump(entries, f, indent=4)

if __name__=='__main__':
    main()
