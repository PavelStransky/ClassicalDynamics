using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Plots

include("models/Vibron.jl")
include("modules/ClassicalDynamics.jl")

# The following trajectory seems regular, but has an enormous LE
vibronParameters = [0.5,0.5,0.5]
initialCondition = [0.0, 1.0203574072397172, -0.9189955424295329, 0.23490092884445612]
Trajectory(initialCondition, vibronParameters; verbose=true)

initialCondition = [0.0, -0.7862068849266208, 0.7769513368418827, 0.7172308051665519]
TestTrajectory(initialCondition, vibronParameters)