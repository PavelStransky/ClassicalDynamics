using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Plots

include("models/BoseHubbard.jl")
include("modules/ClassicalDynamics.jl")

bhParameters = [2,-1,1]
energy = 0.0
x0 = [0.0, 0.0, 0.0, 0.0]
initialCondition = InitialConditions(x0, energy, bhParameters, 1)

println("Initial condition: ", initialCondition)
println("Energy: ", Energy(initialCondition, bhParameters))

Trajectory(initialCondition, bhParameters; verbose=true)