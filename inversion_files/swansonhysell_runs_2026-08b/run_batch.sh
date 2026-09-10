#!/bin/bash
# Sensitivity runs, August 2026 — batch b.
#
# Repeats configurations from batch 08a with TAI_MRM4_z3 and TAI_MRM4_z5
# removed from the input data, to test how much those two grains determine the
# recovered t-T history. Everything else — script, seed, pinned environment —
# is identical to 08a, so a paired 08a/08b comparison isolates the two grains.
#
# Run from the inversion_files/ directory, AFTER batch 08a has finished:
#
#     nohup bash swansonhysell_runs_2026-08b/run_batch.sh > swansonhysell_runs_2026-08b/batch.log 2>&1 &

set -u
cd "$(dirname "$0")/.."          # inversion_files/
BATCH=swansonhysell_runs_2026-08b
DATA=mrm_no_MRM4_z3_z5.csv
JULIA="julia +1.12.5"
MAXJOBS=5

# Add entries here to widen the test; each pair must have an 08a counterpart
# at the same seed for the comparison to mean anything.
CONFIGS="no_jacobsville"
SEEDS="20260710"

# --- preflight -------------------------------------------------------------
# Guard against the two ways this run could silently be wrong: the shared
# script not yet supporting the input-file argument (it would fall back to the
# full mrm.csv and quietly reproduce 08a), or batch 08a still holding the CPUs.
if ! grep -q 'ARGS\[4\]' tTinversion_mrm.jl; then
    echo "ERROR: tTinversion_mrm.jl does not accept an input-file argument."
    echo "       Without it this batch would silently run on the full mrm.csv."
    exit 1
fi
if [ ! -f "$DATA" ]; then
    echo "ERROR: $DATA not found in $(pwd)"; exit 1
fi
if pgrep -qf "tTinversion_mrm.jl"; then
    echo "ERROR: inversion chains are still running (batch 08a?). Wait for them"
    echo "       to finish, or raise MAXJOBS deliberately if you mean to share."
    exit 1
fi
echo "[$(date '+%F %T')] input $DATA: $(($(wc -l < $DATA) - 1)) analyses"

run_one() {
    local config=$1 seed=$2
    local log="$BATCH/${config}_seed${seed}.log"
    local prefix="$BATCH/MRM_${config}_seed${seed}_no_MRM4_z3_z5"

    echo "[$(date '+%F %T')] START  $config seed=$seed"
    $JULIA --project="$BATCH" tTinversion_mrm.jl "$config" "$seed" "$BATCH" "$DATA" > "$log.raw" 2>&1
    local status=$?
    tr '\r' '\n' < "$log.raw" | grep -v -E '^MCMC (burn-in|collection)' > "$log"
    rm -f "$log.raw"

    if [ "$status" -ne 0 ]; then
        echo "[$(date '+%F %T')] FAILED $config seed=$seed (exit $status), see $log"
        return 1
    fi
    for trace in lldist ndist acceptance; do
        [ -f "${prefix}_${trace}.pdf" ] && \
            magick -density 150 "${prefix}_${trace}.pdf" -units PixelsPerInch -density 150 "${prefix}_${trace}.png"
    done
    echo "[$(date '+%F %T')] DONE   $config seed=$seed"
}

for seed in $SEEDS; do
    for config in $CONFIGS; do
        while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 30; done
        run_one "$config" "$seed" &
    done
done
wait
echo "[$(date '+%F %T')] ALL RUNS COMPLETE"
