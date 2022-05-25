using Plots

include("models/TwoFluidLipkin.jl")
include("modules/ClassicalDynamics.jl")

parameters = [0.5, 0.5, -0.5]
energy = 2/7
numTrajectories = 10
dimension = 100

savePath = "d:\\results\\TwoFluid\\"

pyplot(size = (1200,1000))

SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=1, sectionCoordinateX=2, sectionCoordinateY=3)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=1, sectionCoordinateX=2, sectionCoordinateY=4)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=1, sectionCoordinateX=3, sectionCoordinateY=4)

SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=2, sectionCoordinateX=1, sectionCoordinateY=3)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=2, sectionCoordinateX=1, sectionCoordinateY=4)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=2, sectionCoordinateX=3, sectionCoordinateY=4)

SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=3, sectionCoordinateX=1, sectionCoordinateY=2)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=3, sectionCoordinateX=1, sectionCoordinateY=4)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=3, sectionCoordinateX=2, sectionCoordinateY=4)

SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=4, sectionCoordinateX=1, sectionCoordinateY=2)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=4, sectionCoordinateX=1, sectionCoordinateY=3)
SolveEnergy(energy, parameters, dimension, savePath=savePath, timeout=10000, showFigures=true, randomize=true, minimumBound=-sqrt(2.0), maximumBound=sqrt(2.0), sectionPlane=4, sectionCoordinateX=2, sectionCoordinateY=3)
