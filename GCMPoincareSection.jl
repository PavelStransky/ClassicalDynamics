using Plots

include("models/GCM.jl")
include("modules/ClassicalDynamics.jl")

GCMParameters = [-1,0.6,1,1]
energy = 10
numTrajectories = 10
PSPoints = 2000

savePath = "d:\\results\\GCM"

pyplot(size = (1800,1000))
fig = PoincareSection(energy, GCMParameters, numTrajectories; maxPSPoints=PSPoints, verbose=true)
savefig(fig, savePath * "_PS_E=$(energy)_$GCMParameters.png")