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

@everywhere include("models/BoseHubbard.jl")
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

    println("Elapsed time: $time seconds")

    return result
end

for j in LinRange(-1, 1, 101)
    for energy in LinRange(-1, 1.5, 126)
        lyapunovs = LyapunovMap((2, j, 1), energy)

        file = "d:/results/bh/lyapunov/3/" * @sprintf("%.2f_%.2f", j, energy) * ".txt"

        open(file, "a") do io
            for lyapunov in lyapunovs
                println(io, lyapunov)
            end
        end
    end
end
