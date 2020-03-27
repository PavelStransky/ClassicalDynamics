using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Plots

include("models/Vibron.jl")
include("modules/ClassicalDynamics.jl")

vibronParameters = [0.5,0.5,0.5]
energy = 0.5
numTrajectories = 10
PSPoints = 2000

initialCondition = [0.0, 0.1, missing, 0.1]
InitialCondition!(initialCondition, energy, vibronParameters)

initialCondition = [0.0, 1.0203574072397172, -0.9189955424295329, 0.23490092884445612]

pyplot(size = (1800,1000))
Trajectory(initialCondition, vibronParameters)