using Distributed
using Plots

workers = 5

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere include("models\\Dicke.jl")
@everywhere include("modules\\ClassicalDynamics.jl")

""" Calculates PoincarDicke é sections  with Lyapunov exponents for various energies
    Parallel calculation, takes usually days to finish
"""

@everywhere function SolveItem(energy, λ)    
    parameters = [λ, 1.0, 1.0, 1.0]
    dimension = 101

    time = @elapsed sectionLyapunov, trajectories = SolveEnergy(energy, parameters, dimension)

    chaos = 0
    total = 0
    error = 0
    meanLyapunov = 0
    meanLyapunovChaos = 0
    maxLyapunov = 0

    for x in sectionLyapunov
        if x > 0.0
            total += 1
            meanLyapunov += x
        end
        if x == 0.0
            error += 1
        end
        if x > 0.01
            chaos += 1
            meanLyapunovChaos += x
        end        
        maxLyapunov = max(maxLyapunov, x)
    end

    result = [energy, λ, total, chaos, error, total > 0 ? meanLyapunov / total : 0, chaos > 0 ? meanLyapunovChaos / chaos : 0, maxLyapunov, myid(), time, trajectories]

    open("d:\\results\\Dicke_Map_$dimension.txt", "a") do io
        println(io, result)
    end    

    return result
end

function RunMap()
    input = [(energy, λ) for energy in 3:-0.1:-3, λ in 3:-0.1:0]
    result = pmap((args)->SolveItem(args...), input)
    println(result)
end

function ReadMap()
    open("d:\\results\\Dicke_tmp.txt", "r") do io
        println(read(io))
    end
end