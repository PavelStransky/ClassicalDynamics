#!/bin/bash
#SBATCH --job-name=bhnc-sparse
#SBATCH --partition=ffa-preempt        # a requeued slice restarts from scratch; finished slices are kept
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4              # UMFPACK's dense frontal matrices use BLAS threads
#SBATCH --mem=32G                      # LU of the N = 30 sector: about 12 GB of L + U, plus the Krylov basis
#SBATCH --output=/home/%u/results/bh/number-conserving/quantum/3/sparse/logs/%x_%A_%a.out
#SBATCH --error=/home/%u/results/bh/number-conserving/quantum/3/sparse/logs/%x_%A_%a.err
#SBATCH --mail-user=pavel.stransky@matfyz.cuni.cz
#SBATCH --mail-type=END,FAIL

# The slow region of the Liouvillian spectrum at LARGE N (items C3 and C6 of
# number-conserving-BH-paper-TODO.md) on Chimera, through
# BHNumberConservingLiouvillianSparse.jl.  Only the sizes the laptop cannot
# hold go here - N = 30 and beyond; N <= 24 runs on the laptop (see
# number-conserving-BH-compute-plan.md).
#
# One submission does everything, in three stages chained by SLURM itself:
#
#   ./BHNumberConservingLiouvillianSparseChimera.sh <run> <N> <m> <w> [plan options]
#
#   stage plan   (this script, MODE=plan, one task): writes plan.txt - one
#                calibration slice fixes the spacing of the shifts unless
#                --spacing is given - then submits
#   stage slice  (MODE=slice, an ARRAY with one task per shift): slice_<k>.txt
#   stage merge  (MODE=merge, after the array, whatever its exit codes):
#                spectrum.txt, the coverage report and, if the coverage has
#                holes, a refine step that appends shifts and resubmits the
#                missing slices and another merge.  At most REFINE_ROUNDS times.
#
# Examples (C3: sector m = 1 and the window Re λ >= -16; C6: the gap, both
# sectors and a narrow window):
#
#   ./BHNumberConservingLiouvillianSparseChimera.sh c3-N30 30 1 16
#   ./BHNumberConservingLiouvillianSparseChimera.sh c6-N30-m0 30 0 2 --howmany 60
#   ./BHNumberConservingLiouvillianSparseChimera.sh c6-N30-m1 30 1 2 --howmany 60
#   ./BHNumberConservingLiouvillianSparseChimera.sh c3-N30-g8 30 1 16 --g -8
#
# The results land in $HOME/results/bh/number-conserving/quantum/3/sparse/<run>/;
# read them with  julia BHNumberConservingLiouvillianSparse.jl statistics <run>.
#
# Every slice is written atomically and skipped when present, so a preempted
# or failed task is repaired by the next merge round (or by resubmitting the
# array by hand).  Before the first submission warm the Julia cache with the
# same JULIA_CPU_TARGET as below and make sure KrylovKit is installed:
#   export JULIA_CPU_TARGET=generic
#   julia -e 'using Pkg; Pkg.add("KrylovKit"); Pkg.precompile()'
#   julia BHNumberConservingLiouvillianSparse.jl checks
#
# Memory and time scale steeply with N (fill of L + U ~ n^1.5 with
# n = C(N+2,2)^2/3): about 3 GB and minutes per slice at N = 24, about 12 GB
# and ~1 h per slice at N = 30.  For N = 36 raise --mem to 96G and --time.

set -euo pipefail

export JULIA_CPU_TARGET=generic
export CD_NO_PLOTS=true
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"    # UMFPACK's frontal matrices
REFINE_ROUNDS=${REFINE_ROUNDS:-3}
THROTTLE=${THROTTLE:-100}

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SELF="$SCRIPT_DIR/BHNumberConservingLiouvillianSparseChimera.sh"
RESULTS="$HOME/results/bh/number-conserving/quantum/3/sparse"
mkdir -p "$RESULTS/logs"

MODE="${MODE:-submit}"

case "$MODE" in
    submit)
        # Called by hand on the login node: queue the plan stage.
        [ $# -ge 4 ] || { echo "usage: $0 <run> <N> <m> <w> [plan options]" >&2; exit 1; }
        run="$1"; shift
        sbatch --job-name="plan-$run" --cpus-per-task=4 --time=4:00:00 \
               --export=ALL,MODE=plan,RUN="$run",PLAN_ARGS="$*" "$SELF"
        ;;

    plan)
        cd "$SCRIPT_DIR"
        # shellcheck disable=SC2086
        julia BHNumberConservingLiouvillianSparse.jl plan "$RUN" $PLAN_ARGS
        shifts=$(grep '^shifts' "$RESULTS/$RUN/plan.txt" | tr ',' '\n' | wc -l)
        array=$(sbatch --parsable --job-name="slice-$RUN" --array="0-$((shifts - 1))%$THROTTLE" \
                       --export=ALL,MODE=slice,RUN="$RUN" "$SELF")
        sbatch --job-name="merge-$RUN" --dependency="afterany:$array" --time=1:00:00 --mem=8G \
               --cpus-per-task=1 --export=ALL,MODE=merge,RUN="$RUN",ROUND=1 "$SELF"
        echo "plan written: $shifts shifts; slice array $array queued"
        ;;

    slice)
        cd "$SCRIPT_DIR"
        julia BHNumberConservingLiouvillianSparse.jl slice "$RUN"
        ;;

    merge)
        cd "$SCRIPT_DIR"
        julia BHNumberConservingLiouvillianSparse.jl merge "$RUN" | tee "$RESULTS/$RUN/merge_round_${ROUND}.txt"
        if grep -q "coverage INCOMPLETE\|missing slices" "$RESULTS/$RUN/merge_round_${ROUND}.txt" \
                && [ "$ROUND" -lt "$REFINE_ROUNDS" ]; then
            before=$(grep '^shifts' "$RESULTS/$RUN/plan.txt" | tr ',' '\n' | wc -l)
            julia BHNumberConservingLiouvillianSparse.jl refine "$RUN"
            after=$(grep '^shifts' "$RESULTS/$RUN/plan.txt" | tr ',' '\n' | wc -l)
            # the whole range is resubmitted: slices already on disk are skipped at once, which
            # also repairs any task of the previous round that was preempted or failed
            array=$(sbatch --parsable --job-name="slice-$RUN" --array="0-$((after - 1))%$THROTTLE" \
                           --export=ALL,MODE=slice,RUN="$RUN" "$SELF")
            sbatch --job-name="merge-$RUN" --dependency="afterany:$array" --time=1:00:00 --mem=8G \
                   --cpus-per-task=1 --export=ALL,MODE=merge,RUN="$RUN",ROUND=$((ROUND + 1)) "$SELF"
            echo "round $ROUND: $before -> $after shifts, array $array resubmitted"
        fi
        ;;

    *)
        echo "unknown MODE $MODE" >&2
        exit 1
        ;;
esac
