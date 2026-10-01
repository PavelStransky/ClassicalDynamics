#!/bin/bash
# run_ipnp36_number_conserving.sh
#
# Runs the IPNP36 part of number-conserving-BH-paper-TODO.md: every memory-hungry or long QUANTUM
# calculation, at the larger N that the 640 GB of this machine allow.  The split between Chimera
# (the classical Lyapunov sweeps), IPNP36 (this script) and the laptop (the light classical and
# consistency stages, and the analysis) is in number-conserving-BH-compute-plan.md.
#
# Usage, from the repository root (the stages run one after another; every one is resumable - the
# drivers skip what is already on disk - so an interrupted run is simply started again).  Run it
# detached from the terminal, it takes about a day and a half:
#
#   nohup bash run_ipnp36_number_conserving.sh core > ~/nc-ipnp36.log 2>&1 &
#   bash run_ipnp36_number_conserving.sh reference cuts          # selected stages
#   bash run_ipnp36_number_conserving.sh list
#
# Stages (wall-clock estimates on the 24 cores; `core` runs them in this order):
#   setup         KrylovKit, precompilation, consistency checks, sparse vs dense       ~ 30 min
#   reference     C4 C6 K6 K8  the four reference points, N = 6 ... 20 (dense)         ~  5 h
#   cuts          C1 A3  eta = 3 and g = -20 cuts, N = 8 ... 16                         ~  2 h
#   plane         C1   (g, eta) plane at N = 12, 14641 points                           ~  3 h
#   eta0          C2   eta = 0 line with parity, N = 8 ... 16 at every g, N = 18 at every 2nd  ~  3 h
#   sparse        C3 C6  slow region and gap at N = 18, 20, 24, 30 (sparse slicing)      ~ 13 h
#   modes         C5   observable weights with eigenvectors, N = 8 ... 16               ~  6 h
#   variants      C9 C11 K12  kappa = 0, Gamma_- = Gamma_+/2, kappa = 0.1 and 0.5, N <= 16 ~ 30 min
#   dsff          C10  41 spectra around the chaotic point, N = 12                     ~  5 min
#   steady        C7 A10  steady state vs classical measure, N = 8 ... 40              ~  2 h
#   trajectories  C8   quantum trajectories at the multistable point, N = 10 ... 100    ~  3 h
# Extensions (`extensions` runs them in this order; about a day and a half more):
#   plane14       C1   the plane at N = 14 on every second point of the grid (61 x 61)  ~  8 h
#   cuts18        A3   the eta = 3 cut at N = 18                                        ~  5 h
#   sparse36      C3   slow region at N = 36 (about 40 GB per slice, 12 at a time)      ~ 24 h
# Fallback for the classical sweeps if the Chimera queue is slow (any preset of
# BHMapNumberConservingChimeraSubmit.sh, run here 24 tasks at a time):
#   classical <preset>                                          e.g. classical cut-eta0
#
# Output under ~/results/bh/number-conserving/ (the same layout as on the laptop and on Chimera, so
# the directories can be copied back and merged); a log of every stage in .../logs/<stage>.log.

set -uo pipefail
cd "$(dirname "$0")"

WORKERS=${WORKERS:-24}             # Distributed workers of the dense Liouvillian driver
THREADS=${THREADS:-24}             # Julia threads of the multithreaded scripts
RESULTS="$HOME/results/bh/number-conserving"
LOGS="$RESULTS/logs"
mkdir -p "$LOGS"
export CD_NO_PLOTS=true
# Every driver writes under ~/results/bh/number-conserving/<its own subdirectory>; a BH_RESULTS_DIR
# left in the environment would send them all to one place and make the sparse plans below
# invisible to this script, which would then plan them again.
unset BH_RESULTS_DIR

run() {                            # run <log> <command...>
    local log="$1"; shift
    echo "==== $(date -Is)  $*" | tee -a "$LOGS/$log.log"
    "$@" 2>&1 | tee -a "$LOGS/$log.log"
}

map() {                            # map <task> [options]: the dense Liouvillian driver
    run "$1" julia BHNumberConservingLiouvillianMap.jl "$@" --workers "$WORKERS"
}

# Sparse slicing of one run with `processes` slice processes side by side: plan (unless present),
# slices, merge, and up to three refinement rounds until the coverage is complete.  One BLAS thread
# per process - the parallelism is across the slices.
sparse_run() {                     # sparse_run <name> <processes> <plan arguments...>
    local name="$1" processes="$2"; shift 2
    [ -f "$RESULTS/quantum/3/sparse/$name/plan.txt" ] ||
        run sparse julia BHNumberConservingLiouvillianSparse.jl plan "$name" "$@"
    for round in 0 1 2 3; do
        for ((part = 0; part < processes; part++)); do
            OPENBLAS_NUM_THREADS=1 julia BHNumberConservingLiouvillianSparse.jl slices "$name" \
                --part "$part" --parts "$processes" > "$LOGS/sparse_${name}_part${part}.log" 2>&1 &
        done
        wait
        report=$(run sparse julia BHNumberConservingLiouvillianSparse.jl merge "$name")
        grep -q "INCOMPLETE\|missing slices" <<< "$report" || break
        run sparse julia BHNumberConservingLiouvillianSparse.jl refine "$name" > /dev/null
    done
    run sparse julia BHNumberConservingLiouvillianSparse.jl statistics "$name"
}

stage() {
    case "$1" in
        list)
            sed -n '/^# Stages/,/^# Output/p' "$0" | sed '$d' ;;
        setup)
            run setup julia -e 'using Pkg; Pkg.add("KrylovKit"); Pkg.precompile()'
            run setup julia BHNumberConservingChecks.jl
            run setup julia BHNumberConservingLiouvillianSparse.jl validate 12 --w 16 ;;
        reference)
            map reference
            map reference --N 17,18,19,20 ;;
        cuts)
            map cut-eta3
            map cut-g20
            map cut-eta3 --N 16
            map cut-g20 --N 16 ;;
        plane)
            map plane ;;
        eta0)
            map eta0                      # N = 8 ... 14 at every g, N = 16 at every second
            map eta0 --N 16               # the remaining g at N = 16 (the done ones are skipped)
            map eta0 --N 18 --stride 2 ;;
        sparse)
            # C3: slow region of the clean sector at the chaotic point (statistics up to w = 16);
            # about 3 GB per slice at N = 24 and 14 GB at N = 30
            for N in 18 20 24; do sparse_run "c3-N$N" 24 "$N" 1 21; done
            sparse_run c3-N30 20 30 1 21
            # C6: the gap at the chaotic and the multistable point, both sectors
            for point in "-20 3 g-20e3" "-20 1 g-20e1"; do
                set -- $point
                for N in 18 20 24 30; do
                    for m in 0 1; do
                        sparse_run "c6-$3-N$N-m$m" 20 "$N" "$m" 3 --g "$1" --eta "$2" --howmany 60
                    done
                done
            done ;;
        modes)
            run modes julia BHNumberConservingLiouvillianModes.jl --N 8,10,12,14,16 --blas "$THREADS" ;;
        variants)
            for task in kappa0 asymmetric kappa; do map "$task"; done ;;
        dsff)
            map dsff ;;
        steady)
            run steady julia -t "$THREADS" BHNumberConservingSteadyState.jl --N 8,10,12,14,16,20,24,30,36,40 ;;
        trajectories)
            run trajectories julia -t "$THREADS" BHNumberConservingQuantumDynamics.jl trajectories \
                --N 10,20,30,40,50,60,70,80,90,100 --trajectories 24 --time 2000 ;;
        plane14)
            map plane --N 14 --stride 2 ;;
        cuts18)
            map cut-eta3 --N 18 ;;
        sparse36)
            sparse_run c3-N36 12 36 1 21 ;;
        core)
            for s in setup reference cuts plane eta0 sparse modes variants dsff steady trajectories; do
                stage "$s"
            done ;;
        extensions)
            for s in plane14 cuts18 sparse36; do stage "$s"; done ;;
        *)
            echo "unknown stage '$1' - run '$0 list'" >&2
            return 1 ;;
    esac
}

[ $# -eq 0 ] && set -- list
if [ "$1" = "classical" ]; then
    bash ./BHMapNumberConservingChimeraSubmit.sh "${2:?preset}" --local "${3:-$THREADS}"
    exit
fi
for s in "$@"; do stage "$s"; done
