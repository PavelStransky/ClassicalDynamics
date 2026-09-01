using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Printf

# Cluster-aware variant of BHMap.jl for the Chimera SLURM cluster.
#
# Unlike BHMap.jl (which parallelises trajectories within one job via
# Distributed/pmap), this version runs single-threaded and relies on SLURM
# to parallelise across (J, E) pairs instead: one array task computes a
# contiguous block of at most PAIRS_PER_TASK pairs sequentially, so each
# task needs only one CPU. Batching many pairs per task amortises the
# multi-minute Julia + DifferentialEquations startup over the whole block
# instead of paying it again for every single pair.
#
# The full (J, E) grid is flattened into a single 0-based index, J-major /
# E-minor. With SLURM_ARRAY_TASK_ID set, task t computes the pairs with
# flat index t*PAIRS_PER_TASK .. (t+1)*PAIRS_PER_TASK-1 (clamped to the
# grid); without it (e.g. running by hand), it falls back to the full
# sweep, same as BHMap.jl. The array therefore needs
# cld(length(J_VALUES) * length(ENERGY_VALUES), PAIRS_PER_TASK) tasks --
# the script prints the exact --array range on startup. ARRAY_OFFSET shifts
# the task id (counted in tasks, not pairs), so the sweep can be split
# across several sbatch submissions if it is larger than SLURM's
# MaxArraySize (see BHMapChimera.sh). This script never plots, so it skips
# loading Plots/ColorSchemes/PyPlot entirely -- see the CD_NO_PLOTS note in
# modules/ClassicalDynamics.jl.
ENV["CD_NO_PLOTS"] = "true"

include("models/BoseHubbardFull.jl")
include("modules/ClassicalDynamics.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))

# Random.seed!(1234)

# Constants and parameters
const TRAJECTORIES = 1000
const U = 1.0            # Float64 so modelParameters is a concrete NTuple -> type-stable EquationOfMotion!
const L = 4

const RESULTS_DIR = get(ENV, "BH_RESULTS_DIR", joinpath(homedir(), "results", "bh", "lyapunov", "$L"))

# Tangent-dynamics method for the Lyapunov exponent:
#   :vector -> single deviation vector, matrix-free J*v (fast; largest exponent only) -- default here
#   :matrix -> full 2f x 2f stability matrix + eigvals (slower; reproduces the original estimator exactly)
const TANGENT_DYNAMICS = :vector

function SingleTrajectory(energy, parameters, initialConditionEnergyTolerance, tangentDynamics)
    initialCondition = InitialCondition(energy, parameters, initialConditionEnergyTolerance)

    if initialCondition === nothing
        return -1
    end

    return TrajectoryLyapunov(initialCondition, parameters;
        sectionPlane=-1, maximumSectionPoints=-1, maximumIterations=1E6, tangentDynamics=tangentDynamics,
        manifoldProjection=BoseHubbardConservation!)[2]
end

function LyapunovMap(parameters, energy; initialConditionEnergyTolerance=0.0001, numTrajectories=100, tangentDynamics=TANGENT_DYNAMICS)
    time = @elapsed result = [SingleTrajectory(energy, parameters, initialConditionEnergyTolerance, tangentDynamics) for _ in 1:numTrajectories]

    L, J, U = parameters

    nonzero = filter(x -> x > 0, result)
    positive = filter(x -> x > 0.005, result)

    println("Finished J = $J, U = $U, E = $energy");
    println("Number of new trajectories: $(length(result)) ($(length(positive)) unstable)");

    if length(positive) > 0
        print("Λ = $(mean(positive))")
        if length(positive) > 1
            print(" ± $(var(positive))")
        end
        println()
    end

    if length(nonzero) > 0
        println("freg = ", 1 - length(positive) / length(nonzero))
    end

    println("Elapsed time: $time seconds")
    println()

    return result, positive
end

const J_VALUES = collect(LinRange(0.0, 0.5, 101))
const ENERGY_VALUES = collect(LinRange(-0.5, 0.0, 101))

# One SLURM array task computes a contiguous block of at most PAIRS_PER_TASK
# pairs from the flattened (J, E) grid. One pair takes up to ~5 min, so 100
# pairs stay well under Chimera's 12 h per-task limit.
const PAIRS_PER_TASK = 5

const N_ENERGY = length(ENERGY_VALUES)
const TOTAL_PAIRS = length(J_VALUES) * N_ENERGY
const N_TASKS = cld(TOTAL_PAIRS, PAIRS_PER_TASK)

# Flat 0-based index -> (J, E), J-major / E-minor.
pairAt(k) = (J_VALUES[k ÷ N_ENERGY + 1], ENERGY_VALUES[k % N_ENERGY + 1])

println("Grid: $(length(J_VALUES)) J x $N_ENERGY E = $TOTAL_PAIRS pairs; ",
        "$PAIRS_PER_TASK pairs/task -> submit with --array=0-$(N_TASKS - 1)")

# ARRAY_OFFSET shifts the task id (counted in tasks, not pairs), so the sweep
# can be split across several sbatch submissions (see BHMapChimera.sh).
offset = parse(Int, get(ENV, "ARRAY_OFFSET", "0"))

pairs = if haskey(ENV, "SLURM_ARRAY_TASK_ID")
    taskId = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) + offset # SLURM array indices are 0-based
    startIndex = taskId * PAIRS_PER_TASK

    if startIndex >= TOTAL_PAIRS
        @warn "Task index past the end of the grid; nothing to compute" taskId startIndex TOTAL_PAIRS
        Tuple{Float64, Float64}[]
    else
        stopIndex = min(startIndex + PAIRS_PER_TASK, TOTAL_PAIRS) - 1
        println("Task $taskId: flat pair indices $startIndex..$stopIndex ($(stopIndex - startIndex + 1) pairs)")
        [pairAt(k) for k in startIndex:stopIndex]
    end
else
    exit(1)
    [pairAt(k) for k in 0:(TOTAL_PAIRS - 1)]
end

mkpath(RESULTS_DIR)

for (j, energy) in pairs
    file = joinpath(RESULTS_DIR, @sprintf("%.3f_%.3f_%.3f", j, U, energy) * ".txt")
    trajectories = isfile(file) ? countlines(file) : 0

    println("Starting J = $j, U = $U, E = $energy");
    println("Computed trajectories: $trajectories, trajectories to compute: $(TRAJECTORIES - trajectories)");

    lyapunovs, positive = LyapunovMap((L, j, U), energy, numTrajectories=TRAJECTORIES - trajectories)

    if length(lyapunovs) == 0
        continue
    end

    open(file, "a") do io
        for lyapunov in lyapunovs
            println(io, lyapunov)
        end
    end
end
