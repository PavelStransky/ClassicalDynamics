using Plots

include("models/TwoFluidLipkin.jl")
include("modules/ClassicalDynamics.jl")

parameters = [0.5, 0.5, -0.5]
energy = 2/7
numTrajectories = 10
PSPoints = 1000

savePath = "d:\\results\\TwoFluid"

pyplot(size = (1200,1000))
fig = PoincareSection(energy, parameters, numTrajectories; maximumSectionPoints=PSPoints, sectionPlane=2, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), tolerance=1e-8, sectionCoordinateX=3, sectionCoordinateY=1)
#savefig(fig, savePath * "_PS_E=$(energy)_$vibronParameters.png")
