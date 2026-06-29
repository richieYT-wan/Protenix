"""

This script reads a csv file and tries to parse the target, then writes a json file compatible with Protenix's input format
Can also be passed the target sequences in the CLI:
example: 
python3 scripts/make_json_from_csv.py -f data/02_intermediate/1XIW_evaluation_engine_inputs/rfab_1XIW_eval_input_sequences.csv -t QTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCENCM MKIPIEELEDRVFVNCNTSITWVEGTVGTLLSDITRLDLGKRILDPRGIYRCNESTVQVHYRMCQS -s 0 10 20 -r 0 30000

"""
import hashlib
from argparse import ArgumentParser
import pandas as pd
import json
import os, sys
from pathlib import Path
from typing import List, Union, Dict, Tuple
import subprocess
import tempfile
from joblib import Parallel, delayed
from multiprocessing import Pool
from functools import partial
from tqdm.auto import tqdm

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


def find_target_col(df: pd.DataFrame, candidates: List[str]= None):
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


def process_target_msa_template(target: Union[List, str], output_dir: Path):
    """
    Takes the target and builds a JSON to then run `protenix mt` and get its MSA and template
    Args:
        target:

    Returns:
    """

    protein_chains = []
    if type(target)==list:
        for i, t in enumerate(target):
            protein_chains.append(make_protein_chain(t, f'T{i}', 1, None, None))
    elif type(target)==str:
        protein_chains.append(make_protein_chain(target, 'T', 1, None, None))

    entry = [make_entry(protein_chains, 'target', None, None)]
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as tmp:
        json.dump(entry, tmp, indent=2)
        tmp_path = Path(tmp.name)

    _ = subprocess.run(
        [
            "protenix", "mt",
            "-i", str(tmp_path),
            "-o", str(output_dir),
        ],
        text=True,
        check=True
    )

    target_file = tmp_path.parent / (tmp_path.stem + '-final-updated.json')
    with open(target_file, 'r') as f:
        target_json = json.load(f)

    return target_json

def make_msa(pairing_path: Path=None, nonpairing_path: Path=None):
    # Do default old behaviour
    if not (pairing_path and nonpairing_path):
        return {'pairing_db': 'uniref100'}
    else:
        return {'pairedMsaPath': str(pairing_path),
                'unpairedMsaPath': str(nonpairing_path)}

def make_dummy_msa(sequence: str, name: str, msa_dir: Path) -> Tuple[Path, Path]:
    """
    Creates a dummy MSA for a given sequence
    msa_dir should be the containing all the msas, e.g. path/to/msa
    then in this dir there will be path/to/msa/name/pairing.a3m, msa/<name>/non_pairing.a3m etc
    """
    os.makedirs(msa_dir / name, exist_ok=True)
    pairing_path = msa_dir / name / 'pairing.a3m'
    nonpairing_path = msa_dir / name / 'non_pairing.a3m'
    with open(msa_dir / name / 'pairing.a3m', 'w') as f:
        f.writelines(f'>{name}\n{sequence}\n')
    with open(msa_dir / name / 'non_pairing.a3m', 'w') as f:
        f.writelines(f'>{name}\n{sequence}\n')

    return pairing_path, nonpairing_path


def make_protein_chain(sequence, id, count=1,
                       pairing_path: Union[str,Path]=None,
                       nonpairing_path: Union[str, Path]=None):
    res = {
        "proteinChain": {
            "sequence": sequence,
            "id": [id],
            "count": count,
        }
    }
    if (pairing_path and nonpairing_path):
        res['proteinChain']['pairedMsaPath'] = str(pairing_path)
        res['proteinChain']['unpairedMsaPath'] = str(nonpairing_path)
        
    return res
        

# Placeholder for if we add constraints to the folding
def make_constraints():
    pass


def make_entry(protein_chains, name, seeds=None, constraints=None):
    entry = {"name": name}
    if seeds:
        entry["modelSeeds"] = seeds
    entry["sequences"] = protein_chains

    if constraints:
        entry["constraints"] = constraints
    return entry


def make_entry_from_row(vh:str, target_json:Dict, name: str, msa_dir: Path,
                        seeds: List[int]=None, constraints: List[str]=None,
                        ):
    """
    Used to create a full entry based on iterating on vh sequences
    Args:
        vh:
        target: Dictionary as processed by process_target_msa_template
        name:
        seeds:
        precomputed_msa_dir:
        pairing_db:
        constraints:

    Returns:
        entry
    """
    # Start with VH with dummy MSA
    pairing_path, non_pairing_path = make_dummy_msa(vh, name, msa_dir)
    protein_chains = [make_protein_chain(vh, 'H', 1, pairing_path, non_pairing_path)]
    protein_chains.extend(target_json[0]['sequences'])
    return make_entry(protein_chains, name, seeds, constraints)


# For multiprocessing
def wrapper_make_entry(row_tuple, target_json, seeds, msa_dir, constraints):
    _, row = row_tuple
    return make_entry_from_row(row['vh'], target_json,
                               name=f"vh_id_{generate_fablab_hash(row['vh'])}",
                               seeds=seeds, msa_dir=msa_dir, constraints=constraints)


# Taken from FabLab's hashing
def generate_fablab_hash(sequence, length=12):
    """Generate a deterministic, uppercase alphanumeric hash for a given input sequence.

    This function uses the MD5 hash of the input string, converts it into an integer,
    and then encodes it in a custom base-36 alphabet (0-9, A-Z). The final hash is
    truncated or left-padded with zeros to match the desired length (in case of anormaly
    small sequences). The base-36 alphabet is choosen to be more human friendly than
    full ascii character set, but  having higher cardinality than a simple base-16.

    Parameters:
        sequence (str): The input aa string to hash.
        length (int): The length of the resulting hash string (default is 10).
        length of 10 ensure no collision in 5 million sequence set

    Returns:
        str: An uppercase alphanumeric hash string of the specified length.

    Raises:
        ValueError: If the input is not a non-empty string.
    """

    # define possible output character set
    alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    base = len(alphabet)

    # Generate the MD5 hash
    md5_hash = hashlib.md5(sequence.encode()).digest()

    # Convert the hash to an integer
    hash_int = int.from_bytes(md5_hash, byteorder="big")

    # Encode the integer to the custom base
    encoded = []
    while hash_int > 0:
        hash_int, remainder = divmod(hash_int, base)
        encoded.append(alphabet[remainder])

    # Pad or trim to exact length
    encoded_str = "".join(reversed(encoded)).rjust(length, "0")[:length]

    return encoded_str


def main():
    args = parse_args()
    # Read the df and parse rows to be selected + creates outfile
    start, end = args.rows
    df = pd.read_csv(args.input_file)
    df = df.iloc[max(start, 0):min(len(df), end)]
    
    # If no output_file specified, will save in the same directory as input_file with filename replaced as json
    if not args.output_file:
        parent = args.input_file.parent
        basename = args.input_file.stem
        outfile = parent / (basename + '.json')
    else:
        outfile = args.out_file


    ###################################################
    #    Processing target for MSA/template search    #
    ###################################################
    # If no target is provided will try to parse from the df
    if args.target_sequences is None:
        target_col = find_target_col(df)
        # assumes a single target for every sequence in the df
        target = df[target_col].unique()[0]
    else:
        target = args.target_sequences

    target_json = process_target_msa_template(target, outfile.parent)

    ###################################################
    #      Processing VHH + target for final JSON     #
    ###################################################
    df['seq_id'] = [f'{args.input_file.stem}_id_{i:06}' for i in range(len(df))]
    entries = []
    
    # Multiprocessing 
    wrapper = partial(wrapper_make_entry, 
                      target_json=target_json, 
                      seeds=args.seeds, 
                      msa_dir=outfile.parent / 'msa', 
                      constraints=None)

    rows = list(df[['vh', 'seq_id']].iterrows())

    with Pool(16) as p:
        entries = list(tqdm(p.imap(wrapper, rows), total=len(rows)))


                                   
    # for _, row in df[['vh', 'seq_id']].iterrows():
    #     entries.append(make_entry_from_row(row['vh'], target_json,
    #                                        name = row['seq_id'], seeds=args.seeds,
    #                                        msa_dir = outfile.parent / 'msa', constraints=None))

    with open(outfile, 'w') as file:
        json.dump(entries, file, indent=2)


if __name__=='__main__':
    main()
