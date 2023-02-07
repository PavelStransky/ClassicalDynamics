using Plots

include("models/Dicke.jl")
include("modules/ClassicalDynamics.jl")

dickeParameters = [0.5,1,1,1]
energy = 0.5
numTrajectories = 10
PSPoints = 2000

savePath = "d:\\results\\Dicke"

pyplot(size = (1800,1000))
fig = PoincareSection(energy, dickeParameters, numTrajectories; maximumSectionPoints=PSPoints, minimumBound=-2.0, maximumBound=2.0, sectionPlane=2, sectionCoordinateX=1, sectionCoordinateY=3, maximumIterations=1E7)
fig = PoincareSection(energy, dickeParameters, numTrajectories; maximumSectionPoints=PSPoints, minimumBound=-8.0, maximumBound=8.0, sectionPlane=3, sectionCoordinateX=2, sectionCoordinateY=4, maximumIterations=1E7)

#savefig(fig, savePath * "_PS_E=$(energy)_$dickeParameters.png")