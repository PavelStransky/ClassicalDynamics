using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Plots

include("models/BoseHubbardFull.jl")
include("modules/ClassicalDynamics.jl")

Random.seed!(1234)

const TANGENT_DYNAMICS = :vector

bhParameters = (3, 0.5, 1)
energy = 0.75
initialCondition = InitialCondition0(energy, bhParameters, 0.00001)

println("Initial condition: ", initialCondition)
println("Energy: ", Energy(initialCondition, bhParameters))

# Trajectory(initialCondition, bhParameters; verbose=true, tolerance=1E-10)

for i in 1:10
    println("Step $i")
lyapunovs = TrajectoryLyapunov(initialCondition, bhParameters; 
    tangentDynamics=TANGENT_DYNAMICS,
    showFigures=true, sectionPlane=-1, maximumSectionPoints=-1, 
    regularThreshold=1e-3,
    timeInterval=(0, 1e6),
    historyLyapunovExponentLength=1000,
    relativeFluctuationThreshold=1e-6,
    manifoldProjection=BoseHubbardConservation!,
    saveStep=2)

    readline()
end
