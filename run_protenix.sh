#!/usr/bin/bash
# -*- coding: utf-8 -*-
# ─────────────────────────────────────────────────────────────────────────────
# Run Protenix pipeline in parallel across multiple GPUs.
#
# Splits the input CSV into contiguous row chunks, one per GPU, and launches
# them concurrently. The target MSA cache is pre-warmed serially first to
# avoid races between parallel workers.
#
# Example:
#         ./launch_protenix.sh \
#             -i input.csv -c config.yaml -o output/ \
#             --n-gpus 8 --n-jobs 4 \
#             --sabdab-db /shared/dbs/sabdab_nano_v2/sabdab_nano_db
# ./launch_protenix.sh -i input.csv -c config.yaml -o output/ --n-gpus 16 --dry-run
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd /home/JV11_DK2/structure_sandbox/repos/Protenix/

# ─── Defaults ────────────────────────────────────────────────────────────────
INPUT=""
CONFIG=""
OUTDIR_ROOT=""
N_GPUS=""
N_JOBS=5
SABDAB_DB="/home/JV11_DK2/search_database/sabdab_nano/sabdab_nano_db"
DRY_RUN=0
FORCE=0

POLL_INTERVAL=30          # seconds between GPU-busy checks
MAX_WAIT=0                # 0 = wait indefinitely; else seconds

# ─── Usage ───────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $0 -i INPUT -c CONFIG -o OUTPUT_DIR --n-gpus N [options]

Required:
  -i, --input-file PATH      Path to input CSV
  -c, --config-file PATH     Path to config YAML
  -o, --output-dir PATH      Root output directory
      --n-gpus N             Number of GPUs to use (0..N-1)

Optional:
      --n-jobs N             Threads per GPU job (default: ${N_JOBS})
      --sabdab-db PATH       Path to MMseqs2 SAbDAb nanobody DB
                             (default: ${SABDAB_DB})
  -n, --dry-run              Print commands without executing
      --force                Don't wait for busy GPUs
  -h, --help                 Show this help
EOF
}

# ─── CLI parsing ─────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-file)   INPUT=$2;       shift 2 ;;
        -c|--config-file)  CONFIG=$2;      shift 2 ;;
        -o|--output-dir)   OUTDIR_ROOT=$2; shift 2 ;;
        --n-gpus)          N_GPUS=$2;      shift 2 ;;
        --n-jobs)          N_JOBS=$2;      shift 2 ;;
        --sabdab-db)       SABDAB_DB=$2;   shift 2 ;;
        -n|--dry-run)      DRY_RUN=1;      shift ;;
        --force)           FORCE=1;        shift ;;
        -h|--help)         usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

# ─── Validation ──────────────────────────────────────────────────────────────
missing=()
[[ -z "$INPUT"       ]] && missing+=("--input-file")
[[ -z "$CONFIG"      ]] && missing+=("--config-file")
[[ -z "$OUTDIR_ROOT" ]] && missing+=("--output-dir")
[[ -z "$N_GPUS"      ]] && missing+=("--n-gpus")
if (( ${#missing[@]} > 0 )); then
    echo "ERROR: missing required flags: ${missing[*]}" >&2
    usage; exit 2
fi

[[ -f "$INPUT"  ]] || { echo "ERROR: input file not found: $INPUT"   >&2; exit 1; }
[[ -f "$CONFIG" ]] || { echo "ERROR: config file not found: $CONFIG" >&2; exit 1; }
if ! [[ "$N_GPUS" =~ ^[0-9]+$ ]] || (( N_GPUS < 1 )); then
    echo "ERROR: --n-gpus must be a positive integer (got: $N_GPUS)" >&2; exit 2
fi
if ! [[ "$N_JOBS" =~ ^[0-9]+$ ]] || (( N_JOBS < 1 )); then
    echo "ERROR: --n-jobs must be a positive integer (got: $N_JOBS)" >&2; exit 2
fi

mkdir -p "${OUTDIR_ROOT}"

# Derive tag from config filename (strip 'run_config_' prefix, '.yaml' suffix)
TAG=$(basename "$CONFIG" .yaml | sed 's/^run_config_//')
OUTDIR="${OUTDIR_ROOT}/${TAG}"

# Total rows = we let the Python side handle bounds; here we just need to
# compute chunk sizes for N_GPUS workers. We use the CSV row count minus header.
TOTAL_ROWS=$(( $(wc -l < "$INPUT") - 1 ))
if (( TOTAL_ROWS < 1 )); then
    echo "ERROR: input CSV appears empty: $INPUT" >&2; exit 1
fi

if (( DRY_RUN )); then
    echo "[dry-run] Commands that would be executed:"
    echo
fi

echo "Config summary:"
echo "  input        : $INPUT"
echo "  config       : $CONFIG"
echo "  output_root  : $OUTDIR_ROOT"
echo "  tag / outdir : $TAG  →  $OUTDIR"
echo "  total rows   : $TOTAL_ROWS"
echo "  n_gpus       : $N_GPUS"
echo "  n_jobs/gpu   : $N_JOBS"
echo "  sabdab_db    : $SABDAB_DB"
echo

# ─── GPU busy check ──────────────────────────────────────────────────────────
is_gpu_free() {
    local gpu=$1 uuid count
    uuid=$(nvidia-smi --id=${gpu} --query-gpu=uuid --format=csv,noheader 2>/dev/null | tr -d ' ')
    if [[ -z "$uuid" ]]; then return 0; fi
    count=$(nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader 2>/dev/null \
            | grep -c "^${uuid}" || true)
    (( count == 0 ))
}

wait_for_gpu() {
    local gpu=$1 tag=$2 waited=0
    if (( FORCE )); then return 0; fi
    while ! is_gpu_free "$gpu"; do
        if (( waited == 0 )); then
            echo "[${tag} | GPU ${gpu}] busy, waiting..."
        elif (( waited % 300 == 0 )); then
            echo "[${tag} | GPU ${gpu}] still waiting (${waited}s)"
        fi
        if (( MAX_WAIT > 0 && waited >= MAX_WAIT )); then
            echo "[${tag} | GPU ${gpu}] ERROR: exceeded MAX_WAIT=${MAX_WAIT}s" >&2
            return 1
        fi
        sleep "$POLL_INTERVAL"
        waited=$(( waited + POLL_INTERVAL ))
    done
    if (( waited > 0 )); then
        echo "[${tag} | GPU ${gpu}] free (waited ${waited}s), launching."
    fi
}

# ─── Pre-warm target MSA cache (serial) ──────────────────────────────────────
prewarm_target() {
    mkdir -p "${OUTDIR}"

    if [[ -f "${OUTDIR}/target_msa_cached.json" && \
          -f "${OUTDIR}/target_msa_cached.fingerprint" ]]; then
        echo "[prewarm | ${TAG}] target cache already present, skipping"
        return 0
    fi

    echo "[prewarm | ${TAG}] building target MSA cache..."
    local logfile="${OUTDIR}/prewarm.log"

    local cmd=(
        python scripts/run.py
        -i "${INPUT}"
        -o "${OUTDIR}"
        -c "${CONFIG}"
        -r 0 1
        --n_jobs 1
        --sabdab-db "${SABDAB_DB}"
        --prepare-only
    )

    if (( DRY_RUN )); then
        echo "[dry-run] CUDA_VISIBLE_DEVICES=\"\" ${cmd[*]} > ${logfile} 2>&1"
        return 0
    fi

    # protenix mt is CPU/IO-bound; no GPU needed for prewarm.
    CUDA_VISIBLE_DEVICES="" "${cmd[@]}" > "${logfile}" 2>&1

    if [[ ! -f "${OUTDIR}/target_msa_cached.json" ]]; then
        echo "[prewarm | ${TAG}] ERROR: cache not created; see ${logfile}" >&2
        return 1
    fi
    echo "[prewarm | ${TAG}] done"
}

echo "═══ Pre-warming target MSA cache (serial) ═══"
prewarm_target || { echo "Prewarm failed"; exit 1; }
echo "═══ Prewarm complete; launching parallel workers ═══"
echo

# ─── Worker: process one contiguous chunk on one GPU ─────────────────────────
run_one() {
    local gpu=$1 start=$2 end=$3
    mkdir -p "${OUTDIR}"

    wait_for_gpu "$gpu" "$TAG" || return 1

    local args=(
        python scripts/run.py
        -i "${INPUT}"
        -o "${OUTDIR}"
        -c "${CONFIG}"
        -r "${start}" "${end}"
        --n_jobs "${N_JOBS}"
        --sabdab-db "${SABDAB_DB}"
    )

    local logfile="${OUTDIR}/run_gpu${gpu}_${start}_${end}.log"

    if (( DRY_RUN )); then
        echo "[${TAG} | GPU ${gpu}] CUDA_VISIBLE_DEVICES=${gpu} ${args[*]} > ${logfile} 2>&1"
    else
        echo "[${TAG} | GPU ${gpu}] rows ${start}..${end} (size $((end-start))) -> ${logfile}"
        CUDA_VISIBLE_DEVICES=${gpu} "${args[@]}" > "${logfile}" 2>&1
    fi
}

# ─── Compute chunk assignments; launch all GPU workers ───────────────────────
base=$(( TOTAL_ROWS / N_GPUS ))
rem=$((  TOTAL_ROWS % N_GPUS ))

pids=()
row_cursor=0
for (( gpu=0; gpu < N_GPUS; gpu++ )); do
    size=$(( base + ( gpu < rem ? 1 : 0 ) ))
    s=$row_cursor
    e=$(( row_cursor + size ))
    row_cursor=$e

    if (( size == 0 )); then
        echo "[${TAG} | GPU ${gpu}] no rows assigned, skipping"
        continue
    fi

    if (( DRY_RUN )); then
        run_one "$gpu" "$s" "$e"
    else
        run_one "$gpu" "$s" "$e" &
        pids+=($!)
    fi
done

if (( ! DRY_RUN )); then
    echo
    echo "Launched ${#pids[@]} GPU workers: ${pids[*]}"
    wait "${pids[@]}"
    echo "All GPU workers complete."
fi