using Distributed
using Plots

if nprocs() <= 4
    addprocs(5 - nprocs())
end

@everywhere include("models\\Dicke.jl")
@everywhere include("modules\\ClassicalDynamics.jl")

""" Calculates Poincaré sections  with Lyapunov exponents for various energies
    Parallel calculation, takes usually days to finish
"""
function DickeEnergies()
    savePath = "d:\\results\\Dicke"

    @sync @distributed for energy = -5:0.01:3
        parameters = [2.0, 1.0, 1.0, 1.0]
        dimension = 501

        sectionLyapunov, trajectories = SolveEnergy(energy, parameters, dimension, savePath=savePath)

        # Save results
        open(savePath * "_PS_E=$(energy)_$parameters.txt", "w") do io
            println(io, sectionLyapunov)
        end    
    end
end    
