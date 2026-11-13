#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
cd "$root/.."
sim="$root/simulations"
settings="$sim/settings"
out="$root/outputs"
workers=1
save_inputs=0
skip=""
for arg in "$@"; do
    case "$arg" in
        --quick) settings="$sim/settings/quick"; out="$root/outputs/quick" ;;
        --workers=*) workers="${arg#*=}" ;;
        --save-inputs) save_inputs=1 ;;
        --skip=*) skip="${arg#*=}" ;;
        *) echo "unknown option: $arg" >&2; exit 1 ;;
    esac
done

step() {
    local name=$1
    shift
    case ",$skip," in
        *",$name,"*) echo "== $name: skipped"; return 0 ;;
    esac
    echo "== $name"
    "$@"
}

w="--workers=$workers"
step matrix Rscript "$sim/scripts/run_matrix.R" --settings="$settings/matrix.csv" \
    --out="$out/matrix" "$w" --save_inputs="$save_inputs"
step tensor Rscript "$sim/scripts/run_tensor.R" --settings="$settings/tensor.csv" \
    --out="$out/tensor" "$w" --save_inputs="$save_inputs"
step cv Rscript "$sim/scripts/run_cv.R" --settings_dir="$settings" --out="$out/cv" "$w"
step stars Rscript "$sim/scripts/run_stars.R" --settings_dir="$settings" --out="$out/stars" "$w"
step units Rscript "$sim/scripts/run_units.R" --settings="$settings/matrix.csv" \
    --out="$out/units" "$w"
step refit_pilot Rscript "$sim/scripts/run_refit_pilot.R" --settings="$settings/refit_pilot.csv" \
    --out="$out/refit_pilot" "$w"
step refit_comparison Rscript "$sim/scripts/run_refit_comparison.R" \
    --settings="$settings/refit_comparison.csv" --out="$out/refit_comparison" "$w" --save_inputs="$save_inputs"
step aggregate Rscript "$sim/scripts/summarize_results.R" --outputs="$out" --settings_dir="$settings"
step figures python3 -m replication.simulations.figures.draw_figures --results="$out/results"
