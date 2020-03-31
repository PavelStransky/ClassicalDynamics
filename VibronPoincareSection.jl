using Plots

include("models/Vibron.jl")
include("modules/ClassicalDynamics.jl")

vibronParameters = [0.1, -0.8, 0.7]
energy = 0.1
numTrajectories = 10
PSPoints = 1000

savePath = "d:\\results\\Vibron"

pyplot(size = (1200,1000))
fig = PoincareSection(energy, vibronParameters, numTrajectories; maxPSPoints=PSPoints, verbose=true, sectionPlane=1, min=-sqrt(2.0), max=sqrt(2.0), tolerance=1e-8)
#savefig(fig, savePath * "_PS_E=$(energy)_$vibronParameters.png")


# The following trajectory seems regular, but has an enormous LE
vibronParameters = [0.5,0.5,0.5]
energy = 0.5
initialCondition = [0.0, 1.0203574072397172, -0.9189955424295329, 0.23490092884445612]
#SolveTrajectory(initialCondition, vibronParameters; section=2, verbose=true, showFigures=true)

#initialCondition = [0.0, -0.7862068849266208, 0.7769513368418827, 0.7172308051665519]
#TestTrajectory(initialCondition, vibronParameters)