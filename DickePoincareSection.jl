using Plots

include("models/Dicke.jl")
include("modules/ClassicalDynamics.jl")

dickeParameters = [2,0.5,1,1]
energy = 50
numTrajectories = 5
PSPoints = 1000

savePath = "d:\\results\\Dicke"

pyplot(size = (1800,1000))
fig = PoincareSection(energy, dickeParameters, numTrajectories; sectionPlane=1, maxPSPoints=PSPoints, verbose=true, maxICNumber=10000, min=-sqrt(2.0*energy), max=sqrt(2.0*energy))
savefig(fig, savePath * "_PS_E=$(energy)_$dickeParameters.png")