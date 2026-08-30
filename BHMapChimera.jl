using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Printf

# Cluster-aware variant of BHMap.jl for the Chimera SLURM cluster.
#
# Unlike BHMap.jl (which parallelises trajectories within one job via
# Distributed/pmap), this version runs single-threaded and relies on SLURM
# to parallelise across (J, E) pairs instead: one array task computes one
# (J, E) pair sequentially, so each task needs only one CPU.
#
# The full (J, E) grid is flattened into a single index. With
# SLURM_ARRAY_TASK_ID set, this run only computes the pair at that index;
# without it (e.g. running by hand), it falls back to the full sweep, same
# as BHMap.jl. ARRAY_OFFSET shifts the index, so the sweep can be split
# across several sbatch submissions if the grid is larger than SLURM's
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
const TRAJECTORIES = 100
const U = 1.0            # Float64 so modelParameters is a concrete NTuple -> type-stable EquationOfMotion!
const L = 7

const RESULTS_DIR = get(ENV, "BH_RESULTS_DIR", ENV["HOME"] * "/results/bh/lyapunov/5/")

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

const J_VALUES = collect(LinRange(-0.5, 0.5, 21))
const ENERGY_VALUES = collect(LinRange(0.0, 1.5, 31))

offset = parse(Int, get(ENV, "ARRAY_OFFSET", "0"))

pairs = if haskey(ENV, "SLURM_ARRAY_TASK_ID")
    taskId = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) + offset # SLURM array indices are 0-based
    nEnergy = length(ENERGY_VALUES)
    jIndex = taskId ÷ nEnergy + 1
    eIndex = taskId % nEnergy + 1
    [(J_VALUES[jIndex], ENERGY_VALUES[eIndex])]
else
    [(j, energy) for j in J_VALUES for energy in ENERGY_VALUES]
end

for (j, energy) in pairs
    file = RESULTS_DIR * @sprintf("%.3f_%.3f_%.3f", j, U, energy) * ".txt"
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
