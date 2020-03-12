using Distributed
using Plots

workers = 5

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere include("models/Dicke.jl")
@everywhere include("modules/ClassicalDynamics.jl")

""" Calculates Poincaré sections  with Lyapunov exponents for various energies
    Parallel calculation, takes usually days to finish
"""
function DickeEnergies()
    savePath = "d:\\results\\Dicke"

    @sync @distributed for energy = -1:0.1:10
        parameters = [1.0, 1.0, 1.0, 1.0]
        dimension = 101

        sectionLyapunov, trajectories = SolveEnergy(energy, parameters, dimension, savePath=savePath)

        # Save results
        open(savePath * "_PS_E=$(energy)_$parameters.txt", "w") do io
            println(io, sectionLyapunov)
        end    
    end
end    

savePath = "d:\\results\\Dicke"
parameters = [1.0, 1.0, 1.0, 1.0]
dimension = 101

c = Channel()
x = @async put!(remotecall_fetch(SolveEnergy, 2, 0, parameters, dimension, savePath="d:\\results\\Dicke"))
global i = 0.0
while !isready(c)
    global i
    i += 0.1
end
println(x)