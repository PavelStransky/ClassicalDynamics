#!/bin/bash
#SBATCH --job-name=bhnumcons3
#SBATCH --partition=ffa-preempt        # preemptible partition; a requeued task resumes cleanly (see below)
#SBATCH --time=1:00:00                 # CELLS_PER_TASK (20) cells x roughly 2-5 min each; tune after a test run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G               # generous margin for DifferentialEquations precompilation; see below
#SBATCH --array=0-2928%300             # cld(121*121, 5) = 2928 tasks (BHMapNumberConservingChimera.jl prints the exact range)
#SBATCH --output=/home/%u/results/bh/number-conserving/3/logs/bhnumcons_%a.out
#SBATCH --error=/home/%u/results/bh/number-conserving/3/logs/bhnumcons_%a.err
#SBATCH --mail-user=pavel.stransky@matfyz.cuni.cz
#SBATCH --mail-type=END,FAIL

# Submit from the repository root with:
#   sbatch BHMapNumberConservingChimera.sh
#
# The (g, eta) map of the NUMBER-CONSERVING dissipative Bose-Hubbard model of
# number-conserving-BH.md; see BHMapNumberConserving.jl for the physics and the
# output format, and section 8 of the note for what the scan decides.
#
# Each array task computes a contiguous block of at most CELLS_PER_TASK cells
# (set in BHMapNumberConservingChimera.jl, currently 20) on a single CPU: task
# t takes flat indices t*CELLS_PER_TASK .. (t+1)*CELLS_PER_TASK-1 over the
# g-major, eta-minor grid. Batching amortises the multi-minute Julia +
# DifferentialEquations startup over the whole block, which is what keeps CPU
# efficiency reasonable. Results are written per (g, eta) file and the number
# of completed trajectories is read back from disk, so a requeued or preempted
# task (this partition uses REQUEUE preemption) picks up where it left off,
# losing at most the cell it was inside. The same holds for a task killed by
# --time: resubmitting the array skips every finished cell. Both properties
# were checked against the Distributed version - the two produce byte-identical
# result files, and a file truncated mid-cell is completed to exactly the file
# the uninterrupted run would have written.
#
# Timing: a cell is TRAJECTORIES (100) full Lyapunov spectra, each an ODE of
# dimension L(2L+1)+1 integrated over INTEGRATION_TIME (5000) with a QR
# reorthonormalisation every unit of time. On a desktop i9-13900HX one
# trajectory costs about 1 CPU-s at L = 3, so a cell is roughly 2 min and a
# block of 20 about 40 min; older Chimera CPUs with JULIA_CPU_TARGET=generic
# may add a factor ~2, hence the 2 h limit. There are no pathologically slow
# cells here the way there are in the driven map - the norm is conserved
# exactly, so no trajectory can run away - but a strongly chaotic cell needs
# smaller steps and can take several times the average. Check
# `seff <jobid>_<task>` on a few finished tasks and tighten --time. L = 4
# roughly doubles the cost per trajectory.
#
# The array size is cld(length(G_VALUES) * length(Y_VALUES), CELLS_PER_TASK) --
# currently cld(241*121, 20) = 1459, so --array=0-1458. Run
# `julia BHMapNumberConservingChimera.jl` once with no SLURM_ARRAY_TASK_ID and
# it prints the exact range to use (and exits without computing anything). If a
# finer grid ever pushes the array past SLURM's MaxArraySize
# (`scontrol show config | grep MaxArraySize`), split the sweep with
# ARRAY_OFFSET (counted in tasks, not cells), e.g. from the repo root:
#
#   chunk=1000
#   ntasks=1459   # printed by BHMapNumberConservingChimera.jl
#   for ((offset=0; offset<ntasks; offset+=chunk)); do
#       last=$(( offset + chunk - 1 )); (( last >= ntasks )) && last=$(( ntasks - 1 ))
#       sbatch --array=0-$(( last - offset ))%300 --export=ALL,ARRAY_OFFSET=$offset BHMapNumberConservingChimera.sh
#   done
#
# (%a restarts at 0 in every submission, so the log files of the chunks
# overwrite each other; the results do not.)
#
# Results go to
# $HOME/results/bh/number-conserving/3/eta/J_1.000_k_0.300_e_3.000_m_0.000/
# (the directory is derived from L, SCAN, J, κ, η and MODULATION in
# BHMapNumberConservingChimera.jl; set BH_RESULTS_DIR to override it), logs to
# $HOME/results/bh/number-conserving/3/logs/. Cells computed elsewhere are
# skipped if their files are copied into the results directory before
# submitting. The tasks also write parameters.txt (grid, constants and the
# section 5 threshold curve, which analyse_map_number_conserving.py reads);
# that write goes through a per-task temporary and an atomic rename, so the
# hundreds of tasks doing it at once cannot leave a truncated file behind.
#
# CHANGING THE PLANE. SCAN in BHMapNumberConservingChimera.jl selects (g, eta)
# at fixed kappa, (g, kappa) at fixed eta, or (g, modulation) at kappa = 0 -
# the last is question 2 of section 8. Y_VALUES and therefore the array size
# change with it (:modulation has 101 rows, so 1218 tasks), and PATH changes
# too, so the three sweeps do not collide. Re-run the script once without
# SLURM_ARRAY_TASK_ID after switching to get the new --array range, and update
# --output/--error if you change L.
#
# #SBATCH directives are parsed by sbatch itself, not a shell, so $HOME cannot
# be used there directly -- %u (SLURM's own filename pattern for the submitting
# user) is used instead and expands the same way on Chimera, where home is
# always /home/<login>. The log directory must exist before the *first*
# submission (SLURM opens --output/--error when the job starts, before the
# mkdir -p below runs), so create it once by hand:
#   mkdir -p "$HOME/results/bh/number-conserving/3/logs"
#
# IMPORTANT -- warm the cache before your first sbatch, with the SAME
# JULIA_CPU_TARGET=generic set below (a cache built without it targets one
# specific node's CPU and will not be portable -- see that note). With hundreds
# of array tasks starting close together, each one that finds no valid cache
# tries to precompile DifferentialEquations etc. itself; they race on the same
# lock files under ~/.julia/compiled ("stale pidfile" warnings) and can each
# briefly need well over 1G just to compile. Steps, from the repo root:
#   export JULIA_CPU_TARGET=generic
#   salloc -n1 --mem=4G -p ffa-preempt
#   srun julia -e 'using Pkg; Pkg.precompile()'   # let it run to completion
#   srun bash -c 'SLURM_ARRAY_TASK_ID=0 julia BHMapNumberConservingChimera.jl'
#                                                  # computes the first block for real
#                                                  # (g = -50, eta = 0, the regular corner);
#                                                  # some extensions only get triggered
#                                                  # partway through a real solve, and
#                                                  # Pkg.precompile() alone does not always
#                                                  # catch those. Ctrl-C after a couple of
#                                                  # cells have finished is fine.
# Re-run the last command once more and confirm no more "Being precompiled"
# lines before submitting the array. If some nodes still recompile after this,
# they likely have a CPU feature JULIA_CPU_TARGET=generic does not cover either
# -- consider adding `#SBATCH --constraint=...` to pin the array to one
# consistent CPU family instead (see the node feature table in Chimera.md).
#
# Adjust --partition/--time/--array/--mem-per-cpu to taste; see Chimera.md for
# the full partition table. If Julia is managed via a module on your account,
# uncomment the module load line below.

set -euo pipefail

# module load julia

# Chimera's nodes span very different CPU generations (some lack AVX2/FMA
# entirely, others have AVX-512 -- see the node table in Chimera.md). Julia
# precompiles native code for the exact CPU it ran on by default, so a cache
# warmed on one node's CPU is invalid on a node with different features and
# gets silently recompiled there. `generic` makes the cache portable across all
# of them, at some cost to per-trajectory speed. This MUST also be set
# (identically) when warming the cache interactively -- see the note above.
export JULIA_CPU_TARGET=generic

mkdir -p "$HOME/results/bh/number-conserving/3/logs"

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

srun julia BHMapNumberConservingChimera.jl
