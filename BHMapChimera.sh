#!/bin/bash
#SBATCH --job-name=bhmap
#SBATCH --partition=ffa-preempt        # preemptible partition; job resumes cleanly (see note below)
#SBATCH --time=06:00:00                # generous single-CPU budget for one (J, E) pair; tune after a test run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G               # generous margin for DifferentialEquations/Plots precompilation; see note below
#SBATCH --array=0-650%300              # one task per (J, E) pair: 21 J values x 31 E values
#SBATCH --output=/home/%u/results/bh/lyapunov/5/logs/bhmap_%a.out
#SBATCH --error=/home/%u/results/bh/lyapunov/5/logs/bhmap_%a.err
# #SBATCH --mail-user=you@example.com  # uncomment and fill in to get end/fail notifications
# #SBATCH --mail-type=END,FAIL

# Submit from the repository root with:
#   sbatch BHMapChimera.sh
#
# Each array task computes exactly one (J, E) pair on a single CPU (index
# SLURM_ARRAY_TASK_ID, 0-based, flattened over the J-major, E-minor grid).
# Results are written per (J, U, E) file, and completed trajectory counts
# are read back from disk on restart, so a requeued/preempted task (this
# partition uses REQUEUE preemption) picks up where it left off instead of
# recomputing everything.
#
# SLURM's MaxArraySize (check with `scontrol show config | grep
# MaxArraySize`) may be smaller than the pair count above -- currently 651
# (21 J values x 31 E values), so unlikely to be an issue, but if J_VALUES /
# ENERGY_VALUES grow again and sbatch rejects the --array range, split the
# sweep into chunks with ARRAY_OFFSET, e.g. from the repo root:
#
#   chunk=1000
#   total=651   # = length(J_VALUES) * length(ENERGY_VALUES) in BHMapChimera.jl
#   for ((offset=0; offset<total; offset+=chunk)); do
#       last=$(( offset + chunk - 1 )); (( last >= total )) && last=$(( total - 1 ))
#       sbatch --array=0-$(( last - offset ))%300 --export=ALL,ARRAY_OFFSET=$offset BHMapChimera.sh
#   done
#
# Results and logs go under $HOME/results/bh/lyapunov/5/ (on Chimera, home is
# always /home/<login>). #SBATCH directives are parsed by sbatch itself, not
# a shell, so $HOME can't be used there directly -- %u (SLURM's own filename
# pattern for the submitting user) is used instead and expands the same way.
# The log directory must exist before the *first* submission (SLURM opens
# --output/--error when the job starts, before the mkdir -p below runs), so
# create it once by hand: mkdir -p "$HOME/results/bh/lyapunov/5/logs"
#
# IMPORTANT -- warm the cache before your first sbatch: with hundreds of
# array tasks starting close together, each one that finds no compiled
# cache tries to precompile DifferentialEquations/Plots/etc. itself; they
# race on the same lock files under ~/.julia/compiled ("stale pidfile"
# warnings) and can each briefly need well over 1G just to compile. Two
# steps, from the repo root:
#   salloc -n1 --mem=4G -p ffa-preempt
#   srun julia -e 'using Pkg; Pkg.precompile()'   # let it run to completion
#   srun julia BHMapChimera.jl                    # SLURM_ARRAY_TASK_ID unset -> runs
#                                                  # the full sweep; Ctrl-C after a few
#                                                  # (J, E) pairs finish, not immediately --
#                                                  # some extensions (autodiff/sparse
#                                                  # Jacobian) only get triggered partway
#                                                  # through a real solve, and Pkg.precompile()
#                                                  # alone doesn't always catch those.
# Re-run `srun julia BHMapChimera.jl` once more (still no SLURM_ARRAY_TASK_ID)
# and confirm no more "Being precompiled" lines before submitting the array.
#
# Adjust --partition/--time/--array/--mem-per-cpu to taste; see Chimera.md
# for the full partition table. If Julia is managed via a module on your
# account, uncomment the module load line below.

set -euo pipefail

# module load julia

mkdir -p "$HOME/results/bh/lyapunov/5/logs"

cd "$SLURM_SUBMIT_DIR"

# Diagnostics: if the array still wants to precompile despite a warm cache,
# compare this against the same commands run interactively (salloc) -- a
# different `julia` (PATH), $HOME, or active project would each explain it,
# since any of those points to a different (cold) depot/cache.
echo "host: $(hostname)"
echo "HOME: $HOME"
echo "julia: $(command -v julia)"
julia --version
julia -e 'println.(Base.DEPOT_PATH); using Pkg; println(Base.active_project())'

srun julia BHMapChimera.jl
