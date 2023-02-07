using Plots

include("models/Vibron.jl")
include("modules/ClassicalDynamics.jl")

vibronParameters = [0.1, -0.8, 0.7]
energy = 0.1
numTrajectories = 10
PSPoints = 1000

savePath = "d:\\results\\Vibron"

pyplot(size = (1200,1000))
fig = PoincareSection(energy, vibronParameters, numTrajectories; maximumSectionPoints=PSPoints, sectionCoordinateX=1, sectionCoordinateY=3, sectionPlane=4, tolerance=1e-8)
fig = PoincareSection(energy, vibronParameters, numTrajectories; maximumSectionPoints=PSPoints, sectionCoordinateX=2, sectionCoordinateY=4, sectionPlane=3, tolerance=1e-8)

#savefig(fig, savePath * "_PS_E=$(energy)_$vibronParameters.png")

