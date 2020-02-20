using Plots

include("models\\Dicke.jl")
include("modules\\ClassicalDynamics.jl")

dickeParameters = [2,0.5,1,1]
energy = -1
numTrajectories = 10
PSPoints = 2000

savePath = "d:\\results\\Dicke"

pyplot(size = (1800,1000))
fig = PoincareSection(energy, dickeParameters, numTrajectories; maxPSPoints=PSPoints, verbose=true)
savefig(fig, savePath * "_PS_E=$(energy)_$dickeParameters.png")