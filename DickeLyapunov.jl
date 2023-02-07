using Plots

include("models/Dicke.jl")
include("modules/ClassicalDynamics.jl")

dickeParameters = [2, 0.5, 1, 1]
energy = 1
numTrajectories = 20

savePath = "d:\\results\\Dicke"

pyplot(size = (1800, 1000))
fig = LyapunovExponents(energy, dickeParameters, numTrajectories; minimumBound=-2.0, maximumBound=2.0, sectionPlane=2, sectionCoordinateX=1, sectionCoordinateY=3)
savefig(fig, savePath * "_PS_E=$(energy)_$dickeParameters.png")