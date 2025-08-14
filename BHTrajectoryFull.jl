using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Plots
using Logging

include("models/BoseHubbardFull.jl")
include("modules/ClassicalDynamics.jl")

Random.seed!(1234)

bhParameters = (3,-0.2,1)
energy = 0.4
initialCondition = InitialCondition(energy, bhParameters, 0.0001)

println("Initial condition: ", initialCondition)
println("Energy: ", Energy(initialCondition, bhParameters))

Trajectory(initialCondition, bhParameters; verbose=true, tolerance=1E-10)

lyapunovs = TrajectoryLyapunov(initialCondition, bhParameters; 
 showFigures=true, sectionPlane=-1, maximumSectionPoints=-1, tolerance=1E-10, saveStep=1)