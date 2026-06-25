#!/bin/bash
set -euo pipefail

# Run this from <project_root>:
#   bash scripts/launcher.sh

N_SHARDS=4
INPUT_JSON="./data/02_intermediate/260429_rfab_1XIW_filtered/rfab_hyp03_1XIW_filtered.json"
DUMP_BASE="./output/260429_RFab_1XIW_cofolds_no_constraints"

# 1. Shard the input JSON (sharder lives in scripts/)
python scripts/shard.py "$INPUT_JSON" $N_SHARDS

INPUT_STEM=$(basename "$INPUT_JSON" .json)
INPUT_DIR=$(dirname "$INPUT_JSON")

mkdir -p logs

# 2. Launch N parallel processes on the same GPU
for i in $(seq 0 $((N_SHARDS-1))); do
    CUDA_VISIBLE_DEVICES=0 python runner/inference.py \
        --model_name 'protenix_base_20250630_v1.0.0' \
        --seeds 13,30,1213 \
        --dump_dir "${DUMP_BASE}/shard${i}" \
        --input_json_path "${INPUT_DIR}/${INPUT_STEM}_shard${i}.json" \
        --model.N_cycle 10 \
        --sample_diffusion.N_sample 5 \
        --sample_diffusion.N_step 200 \
        --triangle_attention "cuequivariance" \
        --triangle_multiplicative "cuequivariance" \
        > logs/shard${i}.log 2>&1 &
done

wait
echo "All shards complete."