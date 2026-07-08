"""
    Wrapper to run the full pipeline from reading a config yaml file.
        - reading an input CSV (with columns `vh`)
        - generating the intermediate JSON file,
        - running the predictions

    Config file must be built as
    'prepare': # args for script/make_json_from_csv.py
        <prepare_keys> : <prepare_values>
    'predict': # args for script/batch_inference.py // protenix pred
        <predict_keys> : <predict_values>
    'parse': # args for script/parse_outputs.py
        <parse_keys> : <parse_values>
    Some arguments are built into the config for traceability and others are built into the argument parser for
    ease of use (e.g. for snakemake pipelines such as sabdab-db which needs to be pulled from GCS))
"""

import os, sys
import subprocess
from argparse import ArgumentParser
from pathlib import Path
from yaml import load as load_yaml, Loader
from typing import List, Dict, Tuple


def parse_args():
    p = ArgumentParser(description="""Full Protenix pipeline wrapper from reading a CSV file containing vh sequences,  
                                   config yaml file, preparing sequences, generating predictions and parsing outputs.  
                                   Config file must be built as  
                                        prepare: # args for script/make_json_from_csv.py
                                            <prepare_keys> : <prepare_values>
                                        predict: # args for script/batch_inference.py // protenix pred 
                                            <predict_keys> : <predict_values>,
                                   parse only requires command-line arguments (n_jobs)""")

    p.add_argument('-c', '--config_file', metavar='CONFIG', required=True, type=Path, help='Path to config file')
    p.add_argument('-i', '--input_file', metavar='INPUT_FILE', required=True, type=Path, help='Path to input csv file containing vh sequences')
    p.add_argument('-o', '--output_dir', metavar='OUTPUT_DIR', required=True, type=Path, help='Output directory')
    p.add_argument('-r', '--rows',metavar= ('START', 'END'), type = int, default = [0, 10], nargs = 2)
    p.add_argument('--sabdab-db', type=str, default='~/search_database/sabdab_nano/sabdab_nano_db', help='Path to MMseqs2 OAS nanobody database')
    p.add_argument('--n_jobs', type=int, default=16, help='n_jobs for pre/post processing')
    return p.parse_args()


def build_prepare_cmd(config, args, script_dir):
    # Assumes these scripts to be in the same directory as run.py itself
    cmd = ['python', str(script_dir / 'make_json_from_csv.py'),
           '-i', str(args.input_file),
           '-o', str(args.output_dir),
           '--output_json', f'{Path(args.input_file).stem}_rows_{args.rows[0]:04}_{args.rows[1]:04}',
           '-r', str(args.rows[0]), str(args.rows[1]),
           '--n_jobs', str(args.n_jobs),
           '--sabdab_db', str(args.sabdab_db)]
    target_sequences = config['prepare'].get('target_sequences')
    if target_sequences:
        cmd.extend(['-t', target_sequences])
    seeds = config['prepare'].get("seeds", [0])
    cmd.extend(['-s', *[str(s) for s in seeds]])
    name = config['prepare'].get('name')
    if name:
        cmd.extend(['-n', name])
    print(f'Running input preparation with command: {" ".join(cmd)}')
    return cmd

def build_predict_cmd(config, args, script_dir):
    # Assumes these scripts to be in the same directory as run.py itself
    cmd = ['python', str(script_dir.parent / 'runner' / 'batch_inference.py'),
           '-i', str(Path(args.output_dir) / f'{Path(args.input_file).stem}_rows_{args.rows[0]:04}_{args.rows[1]:04}.json'),
           '-o', str(args.output_dir / f'{Path(args.input_file).stem}_rows_{args.rows[0]:04}_{args.rows[1]:04}_outputs'),
           '--use_seeds_in_json', 'True']
    model = config['predict'].get("model", 'protenix-v2')
    cmd.extend(['-n', model])
    print(f'Running batch_inference with command: {" ".join(cmd)}')
    return cmd

def build_parse_cmd(args, script_dir):
    # Assumes these scripts to be in the same directory as run.py itself
    cmd = ['python', str(script_dir / 'parse_outputs.py'),
           '-i', str(args.output_dir / f'{Path(args.input_file).stem}_rows_{args.rows[0]:04}_{args.rows[1]:04}_outputs'),
           '--n_jobs', str(args.n_jobs)]
    print(f'Running output parsing with command: {" ".join(cmd)}')
    return cmd

def main():
    args = parse_args()
    with open(args.config_file, 'r') as f:
        config: Dict = load_yaml(f, Loader)
    script_dir = Path(__file__).parent
    subprocess.run(build_prepare_cmd(config, args, script_dir), check=True)
    subprocess.run(build_predict_cmd(config, args, script_dir), check=True)
    subprocess.run(build_parse_cmd(args, script_dir), check=True)

if __name__ == '__main__':
    main()


