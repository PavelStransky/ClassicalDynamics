using Plots

include("models/TwoLipkin.jl")
include("modules/ClassicalDynamics.jl")

parameters = [1, 0.5, 0.5]
energy = -0.2
numTrajectories = 100
PSPoints = 1000

savePath = "d:\\results\\TwoLipkin"

pyplot(size = (1200,1000))
fig = PoincareSection(energy, parameters, numTrajectories; maximumSectionPoints=PSPoints, sectionPlane=2, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), tolerance=1e-8, sectionCoordinateX=3, sectionCoordinateY=1)
fig = PoincareSection(energy, parameters, numTrajectories; maximumSectionPoints=PSPoints, sectionPlane=2, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), tolerance=1e-8, sectionCoordinateX=4, sectionCoordinateY=1)
fig = PoincareSection(energy, parameters, numTrajectories; maximumSectionPoints=PSPoints, sectionPlane=4, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), tolerance=1e-8, sectionCoordinateX=3, sectionCoordinateY=1)
fig = PoincareSection(energy, parameters, numTrajectories; maximumSectionPoints=PSPoints, sectionPlane=3, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), tolerance=1e-8, sectionCoordinateX=4, sectionCoordinateY=2)
#savefig(fig, savePath * "_PS_E=$(energy)_$vibronParameters.png")
