using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Distributed
using Printf

workers = 8

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

ENV["CD_NO_PLOTS"] = "true"

@everywhere using Logging
@everywhere global_logger(ConsoleLogger(stderr, Logging.Warn))

@everywhere include("models/BoseHubbardFull.jl")
@everywhere include("modules/ClassicalDynamics.jl")

# Random.seed!(1234)

# Constants and parameters
const TRAJECTORIES = 1000
const U = 1.0            # Float64 so modelParameters is a concrete NTuple -> type-stable EquationOfMotion!
const L = 3

const PATH = get(ENV, "BH_RESULTS_DIR", joinpath(homedir(), "results", "bh", "lyapunov", "$L"))

# Tangent-dynamics method for the Lyapunov exponent:
#   :vector -> single deviation vector, matrix-free J*v (fast; largest exponent only) -- default here
#   :matrix -> full 2f x 2f stability matrix + eigvals (slower; reproduces the original estimator exactly)
const TANGENT_DYNAMICS = :vector

function LyapunovMap(parameters, energy; initialConditionEnergyTolerance=0.0001, numTrajectories=100, tangentDynamics=TANGENT_DYNAMICS)
    # tangentDynamics is captured as a closure local (not a global) so it is serialised to the pmap workers
    function SingleTrajectory()
        initialCondition = InitialCondition(energy, parameters, initialConditionEnergyTolerance)

        if initialCondition === nothing
            return -1
        end

        lyapunov = TrajectoryLyapunov(initialCondition, parameters;
            sectionPlane=-1, maximumSectionPoints=-1, maximumIterations=1E6, tangentDynamics=tangentDynamics,
            manifoldProjection=BoseHubbardConservation!)[2]

        return lyapunov
    end

    input = 1:numTrajectories
    
    time = @elapsed result = pmap((args)->SingleTrajectory(), input)

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

for j in LinRange(-0.5, 0.5, 201)
    for energy in LinRange(0.0, 1.5, 301)
        mkpath(PATH)
        file = PATH * "/" * @sprintf("%.3f_%.3f_%.3f", j, U, energy) * ".txt"
        if isfile(file) 
            trajectories = countlines(file) 
        else 
            trajectories = 0
        end
        
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
end