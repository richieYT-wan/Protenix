"""
Generates Protenix-compatible JSONs from a CSV of VH sequences against a fixed target.
Target MSA is computed once. VHH MSAs are built from a local SAbDAb nanobody database.
Local database is built from ~/scripts/setup_vhh_db.sh

Example:
    python3 scripts/make_json_from_csv.py \
        -f data/01_raw/<experiment_name>/<input_file>.csv \
        -o data/02_intermediate/<experiment_name> \
        --target_sequences <target_domain_0> <target_domain_1> [...] <target_domain_n> \
        --oas-db <database_path>/sabdab_nano_db \
        -s 0 13 30 42 1213 -r 0 100
"""
import hashlib
import json
import subprocess
import tempfile

from pathlib import Path
from argparse import ArgumentParser
from functools import partial
from multiprocessing import Pool

from typing import Dict, List, Tuple, Union

import pandas as pd
from tqdm.auto import tqdm


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-f', '--input_file', type=Path, required=True)
    parser.add_argument('-o', '--output_dir', type=Path, required=True)
    parser.add_argument('-t', '--target_sequences', type=str, nargs='+', default=None)
    parser.add_argument('-s', '--seeds', type=int, nargs='+')
    parser.add_argument('-n', '--name', type = str, required=False, default=None, help='Custom sample name')
    parser.add_argument('-r', '--rows', type=int, default=[0, 10], nargs=2,
                        metavar=('START', 'END'))
    parser.add_argument('--oas-db', type=str,
                        default='/data/databases/oas_nano/oas_nano_db',
                        help='Path to MMseqs2 OAS nanobody database')
    parser.add_argument('--skip-vhh-msa', action='store_true',
                        help='Skip VHH MSA building (fallback to dummy)')
    parser.add_argument('--n-jobs', type=int, default=16)
    return parser.parse_args()


from typing import Tuple

def build_vhh_msa(
    sequence: str,
    name: str,
    output_dir: Path,
    oas_db: str = "/data/databases/oas_nano/oas_nano_db",
    sensitivity: float = 7.5,
    max_seqs: int = 5000,
) -> Tuple[Path, Path]:
    """
    Run MMseqs2 search of a VHH against OAS, produce paired+unpaired A3M files.

    Returns:
        (paired_msa_path, unpaired_msa_path)
    """
    output_dir = Path(output_dir) / name
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        query_fasta = tmpdir / "query.fasta"
        query_fasta.write_text(f">{name}\n{sequence}\n")

        # Build query DB
        query_db = tmpdir / "queryDB"
        subprocess.run(["mmseqs", "createdb", str(query_fasta), str(query_db)], check=True)

        # Search
        result_db = tmpdir / "resultDB"
        subprocess.run([
            "mmseqs", "search",
            str(query_db), oas_db, str(result_db), str(tmpdir / "search_tmp"),
            "-s", str(sensitivity),
            "--max-seqs", str(max_seqs),
            "-a",  # store alignments
        ], check=True)

        # Convert to a3m
        a3m_out = output_dir / "msa.a3m"
        subprocess.run([
            "mmseqs", "result2msa",
            str(query_db), oas_db, str(result_db), str(a3m_out),
            "--msa-format-mode", "5",  # a3m
        ], check=True)

    # Protenix expects paired & unpaired MSAs. For VHH, no pairing makes sense
    # (no co-evolution with antigen), so we use the same file for both,
    # or put only the query in the paired file.
    unpaired_path = output_dir / "non_pairing.a3m"
    paired_path = output_dir / "pairing.a3m"

    a3m_text = a3m_out.read_text()
    unpaired_path.write_text(a3m_text)
    # Paired MSA: just the query (no pairing partner exists)
    paired_path.write_text(f">{name}\n{sequence}\n")

    return paired_path, unpaired_path

def find_target_col(df: pd.DataFrame, candidates: List[str] = None):
    if candidates is None:
        candidates = ['target', 'target_seq', 'target_sequence', 'sequence',
                      'seq', 'protein_seq', 'protein_sequence', 'amino_acid',
                      'aa_seq', 'fasta']
    cols_lower = df.columns.str.lower()
    for c in candidates:
        matches = [col for col, cl in zip(df.columns, cols_lower) if cl == c]
        if matches:
            return matches[0]
    for c in candidates:
        matches = [col for col, cl in zip(df.columns, cols_lower) if c in cl]
        if matches:
            return matches[0]
    raise ValueError(f"No target column found in {list(df.columns)}")


def process_target_msa_template(target: Union[List, str], output_dir: Path) -> Dict:
    """Run protenix mt on the target to get its MSA + templates. Runs once."""
    protein_chains = []
    if isinstance(target, list):
        for i, t in enumerate(target):
            protein_chains.append(make_protein_chain(t, f'T{i}', 1))
    else:
        protein_chains.append(make_protein_chain(target, 'T', 1))

    entry = [make_entry(protein_chains, 'target')]
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as tmp:
        json.dump(entry, tmp, indent=2)
        tmp_path = Path(tmp.name)

    subprocess.run(
        ["protenix", "mt", "-i", str(tmp_path), "-o", str(output_dir)],
        text=True, check=True
    )

    target_file = tmp_path.parent / (tmp_path.stem + '-final-updated.json')
    with open(target_file) as f:
        return json.load(f)


def make_dummy_msa(sequence: str, name: str, msa_dir: Path) -> Tuple[Path, Path]:
    """Fallback: dummy single-sequence MSA."""
    (msa_dir / name).mkdir(parents=True, exist_ok=True)
    pairing = msa_dir / name / 'pairing.a3m'
    nonpairing = msa_dir / name / 'non_pairing.a3m'
    pairing.write_text(f'>{name}\n{sequence}\n')
    nonpairing.write_text(f'>{name}\n{sequence}\n')
    return pairing, nonpairing


def make_protein_chain(sequence, id, count=1,
                       pairing_path: Union[str, Path] = None,
                       nonpairing_path: Union[str, Path] = None):
    res = {"proteinChain": {"sequence": sequence, "id": [id], "count": count}}
    if pairing_path and nonpairing_path:
        res['proteinChain']['pairedMsaPath'] = str(pairing_path)
        res['proteinChain']['unpairedMsaPath'] = str(nonpairing_path)
    return res


def make_entry(protein_chains, name, seeds=None, constraints=None):
    entry = {"name": name}
    if seeds:
        entry["modelSeeds"] = seeds
    entry["sequences"] = protein_chains
    if constraints:
        entry["constraints"] = constraints
    return entry


def make_entry_from_row(vh: str, target_json: Dict, name: str, msa_dir: Path,
                        oas_db: str, use_real_msa: bool = True,
                        seeds: List[int] = None, constraints: List[str] = None):
    """Build a Protenix entry with VHH + target chains. Both with real MSAs."""
    if use_real_msa:
        try:
            pairing_path, nonpairing_path = build_vhh_msa(
                sequence=vh, name=name, output_dir=msa_dir, oas_db=oas_db
            )
        except Exception as e:
            print(f"[warn] Real MSA failed for {name} ({e}), falling back to dummy")
            pairing_path, nonpairing_path = make_dummy_msa(vh, name, msa_dir)
    else:
        pairing_path, nonpairing_path = make_dummy_msa(vh, name, msa_dir)

    protein_chains = [make_protein_chain(vh, 'H', 1, pairing_path, nonpairing_path)]
    protein_chains.extend(target_json[0]['sequences'])
    return make_entry(protein_chains, name, seeds, constraints)


def wrapper_make_entry(row_tuple, target_json, seeds, msa_dir, oas_db,
                       use_real_msa, constraints):
    _, row = row_tuple
    name = f"{row['seq_id']}_{generate_fablab_hash(row['vh'])}"
    return make_entry_from_row(
        vh=row['vh'], target_json=target_json, name=name,
        msa_dir=msa_dir, oas_db=oas_db, use_real_msa=use_real_msa,
        seeds=seeds, constraints=constraints,
    )


def generate_fablab_hash(sequence, length=12):
    alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    base = len(alphabet)
    md5_hash = hashlib.md5(sequence.encode()).digest()
    hash_int = int.from_bytes(md5_hash, byteorder="big")
    encoded = []
    while hash_int > 0:
        hash_int, remainder = divmod(hash_int, base)
        encoded.append(alphabet[remainder])
    return "".join(reversed(encoded)).rjust(length, "0")[:length]


def main():
    args = parse_args()
    start, end = args.rows
    df = pd.read_csv(args.input_file)
    df = df.iloc[max(start, 0):min(len(df), end)]

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / (Path(args.input_file).stem + '.json')
    vhh_msa_dir = out_dir / 'vh_msa'

    ###################################################
    #    Processing target for MSA/template search    #
    #    (runs exactly once)                          #
    ###################################################
    target_json_path = out_dir / 'target_msa_cached.json'
    if target_json_path.exists():
        print(f"[info] Loading cached target MSA from {target_json_path}")
        with open(target_json_path) as f:
            target_json = json.load(f)
    else:
        if args.target_sequences is None:
            target_col = find_target_col(df)
            target = df[target_col].unique()[0]
        else:
            target = args.target_sequences
        print(f"[info] Running target MSA/template search...")
        target_json = process_target_msa_template(target, out_dir)
        with open(target_json_path, 'w') as f:
            json.dump(target_json, f, indent=2)
        print(f"[info] Cached target MSA at {target_json_path}")

    ###################################################
    #      Processing VHH + target for final JSON     #
    ###################################################
    if args.name:
        df['seq_id'] = [f'{args.name}_id_{i:06}' for i in range(len(df))]
    else:
        df['seq_id'] = [f'{Path(args.input_file).stem}_id_{i:06}' for i in range(len(df))]

    use_real_msa = not args.skip_vhh_msa
    if use_real_msa and not Path(args.oas_db + ".dbtype").exists():
        print(f"[warn] OAS DB not found at {args.oas_db}, falling back to dummy MSAs")
        use_real_msa = False

    wrapper = partial(
        wrapper_make_entry,
        target_json=target_json,
        seeds=args.seeds,
        msa_dir=vhh_msa_dir,
        oas_db=args.oas_db,
        use_real_msa=use_real_msa,
        constraints=None,
    )

    rows = list(df[['vh', 'seq_id']].iterrows())

    # NB: MMseqs2 is multi-threaded internally; using many parallel processes
    # that each call mmseqs can over-subscribe CPUs. Lower n_jobs if so.
    n_jobs = args.n_jobs if not use_real_msa else max(1, args.n_jobs // 4)
    print(f"[info] Building {len(rows)} entries with {n_jobs} workers "
          f"(real MSAs: {use_real_msa})")

    with Pool(n_jobs) as p:
        entries = list(tqdm(p.imap(wrapper, rows), total=len(rows)))

    with open(out_file, 'w') as file:
        json.dump(entries, file, indent=2)
    print(f"[info] Wrote {len(entries)} entries to {out_file}")


if __name__ == '__main__':
    main()