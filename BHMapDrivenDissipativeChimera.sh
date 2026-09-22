#!/bin/bash
#SBATCH --job-name=bhdriven2
#SBATCH --partition=ffa-preempt        # preemptible partition; job resumes cleanly (see note below)
#SBATCH --time=1:00:00                 # CELLS_PER_TASK (20) cells x typically 0.5-3 min each, plus rare slow cells; tune after a test run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G               # generous margin for DifferentialEquations precompilation; see note below
#SBATCH --array=0-4018%300             # one task per block of CELLS_PER_TASK cells: cld(241*241, 20) = 2905 tasks (BHMapDrivenDissipativeChimera.jl prints the exact range)
#SBATCH --output=/home/%u/results/bh/driven/2/logs/bhdriven_%a.out
#SBATCH --error=/home/%u/results/bh/driven/2/logs/bhdriven_%a.err
#SBATCH --mail-user=pavel.stransky@matfyz.cuni.cz
#SBATCH --mail-type=END,FAIL

# Submit from the repository root with:
#   sbatch BHMapDrivenDissipativeChimera.sh
#
# Each array task computes a contiguous block of at most CELLS_PER_TASK
# (Δ, f) cells (set in BHMapDrivenDissipativeChimera.jl, currently 20) on a
# single CPU: task t takes flat indices t*CELLS_PER_TASK ..
# (t+1)*CELLS_PER_TASK-1 over the Δ-major, f-minor grid. Batching amortises
# the multi-minute Julia + DifferentialEquations startup over the whole
# block, which is what keeps CPU efficiency reasonable. Results are written
# per (Δ, f) file, and completed trajectory counts are read back from disk on
# restart, so a requeued/preempted task (this partition uses REQUEUE
# preemption) picks up where it left off, losing at most the cell it was in.
# The same holds for a task killed by --time: resubmitting the array skips
# every finished cell.
#
# Timing: a cell is 500 trajectories. The desktop run of
# BHMapDrivenDissipative.jl (i9-13900HX) needed 47 CPU-s per cell at
# Δ ≈ 2.5-2.8, about twice that at Δ = 4-6; older Chimera CPUs with
# JULIA_CPU_TARGET=generic may add another factor ~2, so expect roughly
# 1-4 min per cell. About 1% of the cells contain much slower trajectories
# and took several minutes even on 10 workers, so --time leaves a wide margin
# for a block that catches one of them. Check `seff <jobid>_<task>` on a few
# finished tasks and tighten --time.
#
# The array size is cld(length(DELTA_VALUES) * length(F_VALUES),
# CELLS_PER_TASK) -- currently cld(241*241, 20) = 2905, so --array=0-2904.
# Run `julia BHMapDrivenDissipativeChimera.jl` once with no
# SLURM_ARRAY_TASK_ID and it prints the exact range to use (and exits
# without computing anything). If a finer grid ever pushes the array past
# SLURM's MaxArraySize (`scontrol show config | grep MaxArraySize`), split
# the sweep with ARRAY_OFFSET (counted in tasks, not cells), e.g. from the
# repo root:
#
#   chunk=2000
#   ntasks=2905   # printed by BHMapDrivenDissipativeChimera.jl
#   for ((offset=0; offset<ntasks; offset+=chunk)); do
#       last=$(( offset + chunk - 1 )); (( last >= ntasks )) && last=$(( ntasks - 1 ))
#       sbatch --array=0-$(( last - offset ))%500 --export=ALL,ARRAY_OFFSET=$offset BHMapDrivenDissipativeChimera.sh
#   done
#
# (%a restarts at 0 in every submission, so the log files of the chunks
# overwrite each other; the results do not.)
#
# Results go to $HOME/results/bh/driven/2/J_-1.000_g_2.000_k_1.000/ (the
# directory is derived from L, J, g, κ in BHMapDrivenDissipativeChimera.jl;
# set BH_RESULTS_DIR to override it), logs to $HOME/results/bh/driven/2/logs/.
# Cells already computed elsewhere are skipped if their files are copied into
# the results directory before submitting. #SBATCH directives are parsed by
# sbatch itself, not a shell, so $HOME can't be used there directly -- %u
# (SLURM's own filename pattern for the submitting user) is used instead and
# expands the same way on Chimera, where home is always /home/<login>.
# The log directory must exist before the *first* submission (SLURM opens
# --output/--error when the job starts, before the mkdir -p below runs), so
# create it once by hand: mkdir -p "$HOME/results/bh/driven/2/logs"
#
# IMPORTANT -- warm the cache before your first sbatch, with the SAME
# JULIA_CPU_TARGET=generic set below (a cache built without it targets one
# specific node's CPU and won't be portable -- see that note). With hundreds
# of array tasks starting close together, each one that finds no valid
# cache tries to precompile DifferentialEquations etc. itself; they race on
# the same lock files under ~/.julia/compiled ("stale pidfile" warnings) and
# can each briefly need well over 1G just to compile. Steps, from the repo
# root:
#   export JULIA_CPU_TARGET=generic
#   salloc -n1 --mem=4G -p ffa-preempt
#   srun julia -e 'using Pkg; Pkg.precompile()'   # let it run to completion
#   srun bash -c 'SLURM_ARRAY_TASK_ID=0 julia BHMapDrivenDissipativeChimera.jl'
#                                                  # computes the first block for real
#                                                  # (Δ = 0, the fast undriven edge); some
#                                                  # extensions only get triggered partway
#                                                  # through a real solve, and
#                                                  # Pkg.precompile() alone doesn't always
#                                                  # catch those. Ctrl-C after a few cells
#                                                  # have finished is fine.
# Re-run the last command once more and confirm no more "Being precompiled"
# lines before submitting the array. If some nodes still recompile after
# this, they likely have a CPU feature JULIA_CPU_TARGET=generic doesn't cover
# either -- consider adding `#SBATCH --constraint=...` to pin the array to
# one consistent CPU family instead (see the node feature table in
# Chimera.md).
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

mkdir -p "$HOME/results/bh/driven/2/logs"

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

srun julia BHMapDrivenDissipativeChimera.jl
