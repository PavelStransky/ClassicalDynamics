using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Distributed
using Printf

workers = 10

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere using Logging
@everywhere global_logger(ConsoleLogger(stderr, Logging.Warn))

@everywhere include("models/BoseHubbardFull.jl")
@everywhere include("modules/ClassicalDynamics.jl")

Random.seed!(1234)

function LyapunovMap(parameters, energy; initialConditionTolerance = 0.0001, numTrajectories = 100)

    function SingleTrajectory()
        initialCondition = InitialCondition(energy, parameters, initialConditionTolerance)

        if initialCondition === nothing
            return -1
        end

        lyapunov = TrajectoryLyapunov(initialCondition, parameters; 
            sectionPlane=-1, maximumSectionPoints=-1, maximumIterations=1E6)[2]
        
        return lyapunov
    end

    input = 1:numTrajectories
    
    time = @elapsed result = pmap((args)->SingleTrajectory(), input)

    L, J, U = parameters

    result = filter(x -> x > 0, result)
    positive = filter(x -> x > 0.005, result)

    println("J = $J, E = $energy");
    println("Number of trajectories: $(length(result)) ($(length(positive)) unstable)");

    if length(positive) > 0
        print("Λ = $(mean(positive))")
        if length(positive) > 1
            print(" ± $(var(positive))")
        end
        println()
    end
    
    if length(result) > 0
        println("freg = ", 1 - length(positive) / length(result))
    end

    println("Elapsed time: $time seconds")
    println()

    return result, positive
end

for j in LinRange(-0.2, 0.2, 21)
    for energy in LinRange(0, 1.5, 76)
        lyapunovs, positive = LyapunovMap((3, j, 1), energy)

        if length(lyapunovs) == 0
            continue
        end

        file = "d:/results/bh/lyapunov/3/" * @sprintf("%.2f_%.2f", j, energy) * ".txt"

        open(file, "a") do io
            for lyapunov in lyapunovs
                println(io, lyapunov)
            end
        end
    end
end