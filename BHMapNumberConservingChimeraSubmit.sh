#!/bin/bash
# Submits one of the classical Lyapunov sweeps of number-conserving-BH-paper-TODO.md to Chimera,
# through BHMapNumberConservingChimera.sh / BHMapNumberConservingChimera.jl.
#
# Usage, from the repository root on the Chimera login node:
#   ./BHMapNumberConservingChimeraSubmit.sh <preset> [--dry-run]
#   ./BHMapNumberConservingChimeraSubmit.sh list
#
# or, on a machine WITHOUT SLURM (IPNP36, or the laptop in Git Bash), to run the same array tasks
# locally, P at a time - every task is a separate `julia BHMapNumberConservingChimera.jl` with its
# SLURM_ARRAY_TASK_ID set by hand, so the results are identical to the cluster ones and a preset
# can even be split between the cluster and a local machine (finished cells are skipped):
#   ./BHMapNumberConservingChimeraSubmit.sh <preset> --local P
#
# Each preset exports the BH_* variables read by BHMapNumberConservingChimera.jl (model, grid,
# trajectories, jitter, results tag), computes the array size from the grid - no Julia start on the
# login node - and passes --array, --time, --job-name and the log paths on the sbatch command line,
# where they override the #SBATCH defaults of BHMapNumberConservingChimera.sh.  Presets can be
# submitted one after another and wait in the queue together: every one of them writes into its own
# results directory, and none depends on the constants in the .jl file being edited.
#
# The presets (cost in laptop-core-hours at ~3 s per trajectory; the generic-CPU build on Chimera
# is roughly twice slower):
#
#   map-k0.3   C12   (g, eta) map, L = 3, kappa = 0.3, 121 x 121 cells, 100 trajectories,
#                    jittered - the main classical figure.                         ~ 1200 h
#   map-k0     C12/C9  the same plane at kappa = 0 (volume preserving: no attractors, the map
#                    shows where Hamiltonian-like chaos lives), 121 x 61, 25 traj.   ~ 150 h
#   cut-eta3   C1/C15/A3/A5  eta = 3, g = -50 ... -3 in steps of 0.05, JITTER = 0 (sharp
#                    windows, bifurcation points), 100 trajectories.                 ~   80 h
#   cut-eta0   C2/C13  eta = 0 exactly (the jittered map row eta = 0 is not), g = -50 ... -3
#                    in steps of 0.25, JITTER = 0.                                   ~   20 h
#   cut-g20    C1    g = -20, eta = 0 ... 6 in steps of 0.05, JITTER = 0 - the classical side of
#                    the quantum eta cut.                                            ~    7 h
#   L4-eta3    C15/K10  L = 4, eta = 3, g = -50 ... -20 in steps of 0.05, JITTER = 0 - where the
#                    classification runs strange -> torus -> undetermined.           ~   60 h
#   ray-J      the EXPERIMENTAL path: g = -20, eta = 3, kappa = 0.3 held fixed and the hopping
#                    J swept geometrically from 0.02 to 6 (121 values, J = 1 is the
#                    chaotic reference point), as a lattice depth would sweep it.  The
#                    MI threshold is crossed twice, at J ~ 0.028 and J ~ 4.4.  Every cell
#                    is integrated in units of 1/J (see BHMapNumberConservingChimera.jl),
#                    so the cost grows like 1/J: one cell per task.                  ~   45 h
#
# BEFORE THE FIRST SUBMISSION warm the Julia cache exactly as described in
# BHMapNumberConservingChimera.sh (JULIA_CPU_TARGET=generic, Pkg.precompile, one real block).
#
# Results:  $HOME/results/bh/number-conserving/<L>/<scan>/J_..._k_..._e_..._m_...[_<tag>]/
#           (ray-J: .../3/J/J_1.000_k_0.300_e_3.000_m_0.000_ray/)
# Logs:     $HOME/results/bh/number-conserving/<L>/logs/<preset>/bhnumcons_<task>.out/.err
# Analysis: python analyse_map_number_conserving.py <results dir>     (the two maps)
#           python analyse_cut_number_conserving.py <results dir>     (the one-dimensional cuts)

set -euo pipefail

preset="${1:-}"
dryRun="${2:-}"
parallel="${3:-}"          # with --local: the number of tasks running at once

# Throttle: at most this many array tasks of one submission run at once.
THROTTLE=300

# Axis given as first:last:count[:log] - only the count matters for the array size.
Count() { echo "$1" | cut -d: -f3; }

# Every preset starts from a clean slate: sbatch runs with --export=ALL, so a BH_* variable left
# exported in the login shell would otherwise leak into the job (BH_RESULTS_DIR would even send
# every preset into the same directory).
unset BH_L BH_SCAN BH_KAPPA BH_ETA BH_MODULATION BH_G BH_Y BH_TRAJECTORIES BH_RELAXATION_TIME \
      BH_INTEGRATION_TIME BH_JITTER BH_CELLS_PER_TASK BH_TAG BH_RESULTS_DIR BH_LOG_DIR
scan=eta
export BH_ETA=3.0          # the fixed η of the scans in g and κ; swept by the η scans

case "$preset" in
    map-k0.3)
        L=3; export BH_KAPPA=0.3 BH_G="-50:-2:121" BH_Y="0:6:121"
        export BH_TRAJECTORIES=100 BH_JITTER=1.0 BH_CELLS_PER_TASK=5 BH_TAG=""
        time="2:00:00" ;;
    map-k0)
        L=3; export BH_KAPPA=0.0 BH_G="-50:-2:121" BH_Y="0:6:61"
        export BH_TRAJECTORIES=25 BH_JITTER=1.0 BH_CELLS_PER_TASK=10 BH_TAG=""
        time="2:00:00" ;;
    cut-eta3)
        L=3; export BH_KAPPA=0.3 BH_G="-50:-3:941" BH_Y="3:3:1"
        export BH_TRAJECTORIES=100 BH_JITTER=0.0 BH_CELLS_PER_TASK=5 BH_TAG="cut"
        time="2:00:00" ;;
    cut-eta0)
        L=3; export BH_KAPPA=0.3 BH_G="-50:-3:189" BH_Y="0:0:1"
        export BH_TRAJECTORIES=100 BH_JITTER=0.0 BH_CELLS_PER_TASK=5 BH_TAG="cut"
        time="2:00:00" ;;
    cut-g20)
        L=3; export BH_KAPPA=0.3 BH_G="-20:-20:1" BH_Y="0:6:121"
        export BH_TRAJECTORIES=100 BH_JITTER=0.0 BH_CELLS_PER_TASK=5 BH_TAG="cut"
        time="2:00:00" ;;
    L4-eta3)
        L=4; export BH_KAPPA=0.3 BH_G="-50:-20:601" BH_Y="3:3:1"
        export BH_TRAJECTORIES=100 BH_JITTER=0.0 BH_CELLS_PER_TASK=3 BH_TAG="cut"
        time="3:00:00" ;;
    ray-J)
        L=3; scan=J; export BH_KAPPA=0.3 BH_ETA=3.0 BH_G="-20:-20:1" BH_Y="0.02:6:121:log"
        export BH_TRAJECTORIES=100 BH_JITTER=0.0 BH_CELLS_PER_TASK=1 BH_TAG="ray"
        time="8:00:00" ;;
    list|"")
        sed -n '/^# The presets/,/^# BEFORE/p' "$0" | sed '$d'
        exit 0 ;;
    *)
        echo "unknown preset '$preset' - run '$0 list'" >&2
        exit 1 ;;
esac

export BH_L="$L" BH_SCAN="$scan"

cells=$(( $(Count "$BH_G") * $(Count "$BH_Y") ))
tasks=$(( (cells + BH_CELLS_PER_TASK - 1) / BH_CELLS_PER_TASK ))

logDir="$HOME/results/bh/number-conserving/$L/logs/$preset"
export BH_LOG_DIR="$logDir"
mkdir -p "$logDir"     # SLURM opens --output/--error before the job script runs

echo "preset $preset: L = $L, scan $scan, kappa = $BH_KAPPA, eta = $BH_ETA, g = $BH_G, y = $BH_Y"
echo "  $cells cells, $BH_CELLS_PER_TASK per task -> --array=0-$((tasks - 1))%$THROTTLE, --time=$time"
echo "  logs in $logDir"

command=(sbatch
    --job-name="bhnc-$preset"
    --array="0-$((tasks - 1))%$THROTTLE"
    --time="$time"
    --output="$logDir/bhnumcons_%a.out"
    --error="$logDir/bhnumcons_%a.err"
    --export=ALL
    BHMapNumberConservingChimera.sh)

if [ "$dryRun" = "--dry-run" ]; then
    echo "  dry run: ${command[*]}"
elif [ "$dryRun" = "--local" ]; then
    parallel="${parallel:-$(( $(nproc) - 1 ))}"
    echo "  running the $tasks tasks locally, $parallel at a time; logs in $logDir"
    export CD_NO_PLOTS=true
    seq 0 $((tasks - 1)) | xargs -P "$parallel" -I{} sh -c \
        'SLURM_ARRAY_TASK_ID={} julia BHMapNumberConservingChimera.jl > "$BH_LOG_DIR/bhnumcons_{}.out" 2> "$BH_LOG_DIR/bhnumcons_{}.err"'
    echo "  done"
else
    "${command[@]}"
fi
