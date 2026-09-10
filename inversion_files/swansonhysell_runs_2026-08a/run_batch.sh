#!/bin/bash
# Batch runner for the August 2026 inversion runs (new mrm.csv with raw He dates).
#
# Runs the 3 constraint configurations at the 3 seeds used in the July 2026
# batches, so that these chains differ from those only in the input data.
# Run from the inversion_files/ directory:
#
#     nohup bash swansonhysell_runs_2026-08a/run_batch.sh > swansonhysell_runs_2026-08a/batch.log 2>&1 &
#
# MAXJOBS chains run concurrently; each takes roughly 20 hours.

set -u
cd "$(dirname "$0")/.."          # inversion_files/
BATCH=swansonhysell_runs_2026-08a
JULIA="julia +1.12.5"
MAXJOBS=5

CONFIGS="unconstrained no_jacobsville all_constraints"
SEEDS="20260710 12345 777"

run_one() {
    local config=$1 seed=$2
    local log="$BATCH/${config}_seed${seed}.log"
    local prefix="$BATCH/MRM_${config}_seed${seed}"

    echo "[$(date '+%F %T')] START  $config seed=$seed"
    # Write unfiltered output straight to disk so the run can be monitored live;
    # a pipeline here would block-buffer and leave the log empty until exit.
    $JULIA --project="$BATCH" tTinversion_mrm.jl "$config" "$seed" "$BATCH" > "$log.raw" 2>&1
    local status=$?

    # Strip the ProgressMeter carriage-return output from the archived log
    tr '\r' '\n' < "$log.raw" | grep -v -E '^MCMC (burn-in|collection)' > "$log"
    rm -f "$log.raw"

    if [ "$status" -ne 0 ]; then
        echo "[$(date '+%F %T')] FAILED $config seed=$seed (exit $status), see $log"
        return 1
    fi

    # Archive the MCMC diagnostic traces as 150 dpi rasterizations; the vector
    # originals plot all 10^6 steps and run to ~15 MB each (see .gitignore).
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
