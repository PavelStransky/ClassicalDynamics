#!/bin/bash
#SBATCH --job-name=bhmap7
#SBATCH --partition=ffa-preempt        # preemptible partition; job resumes cleanly (see note below)
#SBATCH --time=0:40:00                 # PAIRS_PER_TASK (50) pairs x up to ~5 min each, under Chimera's 12 h cap; tune after a test run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G               # generous margin for DifferentialEquations/Plots precompilation; see note below
#SBATCH --array=0-4060%400              # one task per block of PAIRS_PER_TASK pairs: cld(101*301, 100) = 305 tasks (BHMapChimera.jl prints the exact range)
#SBATCH --output=/home/%u/results/bh/lyapunov/7/logs/bhmap_%a.out
#SBATCH --error=/home/%u/results/bh/lyapunov/7/logs/bhmap_%a.err
#SBATCH --mail-user=pavel.stransky@matfyz.cuni.cz  # uncomment and fill in to get end/fail notifications
#SBATCH --mail-type=END,FAIL

# Submit from the repository root with:
#   sbatch BHMapChimera.sh
#
# Each array task computes a contiguous block of at most PAIRS_PER_TASK
# pairs (set in BHMapChimera.jl, currently 100) on a single CPU: task t
# takes flat indices t*PAIRS_PER_TASK .. (t+1)*PAIRS_PER_TASK-1 over the
# J-major, E-minor grid. Batching amortises the multi-minute Julia +
# DifferentialEquations startup over ~100 solves instead of one, which is
# what keeps CPU efficiency reasonable. Results are written per (J, U, E)
# file, and completed trajectory counts are read back from disk on restart,
# so a requeued/preempted task (this partition uses REQUEUE preemption)
# picks up where it left off instead of recomputing everything.
#
# The array size is cld(length(J_VALUES) * length(ENERGY_VALUES),
# PAIRS_PER_TASK) -- currently cld(101*301, 100) = 305, so --array=0-304.
# Run `julia BHMapChimera.jl` once with no SLURM_ARRAY_TASK_ID and it
# prints the exact range to use. 305 is far below SLURM's MaxArraySize
# (`scontrol show config | grep MaxArraySize`), but if a much finer grid
# ever pushes past it, split the sweep with ARRAY_OFFSET (counted in tasks,
# not pairs), e.g. from the repo root:
#
#   chunk=200
#   ntasks=305   # = cld(length(J_VALUES)*length(ENERGY_VALUES), PAIRS_PER_TASK), printed by BHMapChimera.jl
#   for ((offset=0; offset<ntasks; offset+=chunk)); do
#       last=$(( offset + chunk - 1 )); (( last >= ntasks )) && last=$(( ntasks - 1 ))
#       sbatch --array=0-$(( last - offset ))%300 --export=ALL,ARRAY_OFFSET=$offset BHMapChimera.sh
#   done
#
# Results and logs go under $HOME/results/bh/lyapunov/4/ (on Chimera, home is
# always /home/<login>). #SBATCH directives are parsed by sbatch itself, not
# a shell, so $HOME can't be used there directly -- %u (SLURM's own filename
# pattern for the submitting user) is used instead and expands the same way.
# The log directory must exist before the *first* submission (SLURM opens
# --output/--error when the job starts, before the mkdir -p below runs), so
# create it once by hand: mkdir -p "$HOME/results/bh/lyapunov/4/logs"
#
# IMPORTANT -- warm the cache before your first sbatch, with the SAME
# JULIA_CPU_TARGET=generic set below (a cache built without it targets one
# specific node's CPU and won't be portable -- see that note). With hundreds
# of array tasks starting close together, each one that finds no valid
# cache tries to precompile DifferentialEquations/Plots/etc. itself; they
# race on the same lock files under ~/.julia/compiled ("stale pidfile"
# warnings) and can each briefly need well over 1G just to compile. Steps,
# from the repo root:
#   export JULIA_CPU_TARGET=generic
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
# If some nodes still recompile after this, they likely have a CPU feature
# JULIA_CPU_TARGET=generic doesn't cover either -- consider adding
# `#SBATCH --constraint=...` to pin the array to one consistent CPU family
# instead (see the node feature table in Chimera.md).
#
# Adjust --partition/--time/--array/--mem-per-cpu to taste; see Chimera.md
# for the full partition table. If Julia is managed via a module on your
# account, uncomment the module load line below.

set -euo pipefail

# module load julia

# Chimera's nodes span very different CPU generations (some lack AVX2/FMA
# entirely, others have AVX-512 -- see the node table in Chimera.md). Julia
# precompiles native code for the exact CPU it ran on by default, so a
# cache warmed on one node's CPU is invalid on a node with different
# features and gets silently recompiled there. `generic` makes the cache
# portable across all of them, at some cost to per-trajectory speed. This
# MUST also be set (identically) when warming the cache interactively --
# see the note above.
export JULIA_CPU_TARGET=generic

mkdir -p "$HOME/results/bh/lyapunov/7/logs"

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
