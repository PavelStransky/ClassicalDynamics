#!/bin/bash
#SBATCH --job-name=bhmap
#SBATCH --partition=ffa-preempt        # preemptible partition; job resumes cleanly (see note below)
#SBATCH --time=06:00:00                # generous single-CPU budget for one (J, E) pair; tune after a test run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-15250                # one task per (J, E) pair: 101 J values x 151 E values, 300 running at once
#SBATCH --output=/home/strap8am/results/bh/lyapunov/5/logs/bhmap_%a.out
#SBATCH --error=/home/strap8am/results/bh/lyapunov/5/logs/bhmap_%a.err
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
# SLURM's MaxArraySize (check with
# `scontrol show config | grep MaxArraySize`) may be smaller than 15251. If
# sbatch rejects the --array range above, split the sweep into chunks with
# ARRAY_OFFSET, e.g. from the repo root:
#
#   chunk=1000
#   total=15251
#   for ((offset=0; offset<total; offset+=chunk)); do
#       last=$(( offset + chunk - 1 )); (( last >= total )) && last=$(( total - 1 ))
#       sbatch --array=0-$(( last - offset ))%300 --export=ALL,ARRAY_OFFSET=$offset BHMapChimera.sh
#   done
#
# Adjust --partition/--time/--array/--mem-per-cpu to taste; see Chimera.md
# for the full partition table. If Julia is managed via a module on your
# account, uncomment the module load line below.

set -euo pipefail

# module load julia

mkdir -p /home/strap8am/results/bh/lyapunov/5/logs

cd "$SLURM_SUBMIT_DIR"

srun julia BHMapChimera.jl
