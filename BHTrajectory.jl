using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Plots

include("models/BoseHubbardFull.jl")
include("modules/ClassicalDynamics.jl")

Random.seed!(1234)

const TANGENT_DYNAMICS = :vector

bhParameters = (5, 0.5, 1)
energy = 1.2

# Trajectory(initialCondition, bhParameters; verbose=true, tolerance=1E-10)

for i in 1:10
    println("Step $i")
    initialCondition = InitialCondition(energy, bhParameters, 0.0001, maxInitialConditions=100000000)

    if initialCondition === nothing
        println("No initial condition found")
        continue
    else
        println("Energy: ", Energy(initialCondition, bhParameters))
    end

    calculationTime = @elapsed lyapunovs = TrajectoryLyapunov(initialCondition, bhParameters; 
        tangentDynamics=TANGENT_DYNAMICS,
        showFigures=true, sectionPlane=-1, maximumSectionPoints=-1, 
        regularThreshold=1e-3,
        timeInterval=(0, 1e6),
        relaxationTime=1000,
        historyLyapunovExponentLength=500,
        manifoldProjection=BoseHubbardConservation!,
        saveStep=2)

    println("Calculation time: ", calculationTime)
end

readline()
