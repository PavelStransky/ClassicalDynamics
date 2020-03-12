using Distributed
using Plots

workers = 5

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere include("models\\Dicke.jl")
@everywhere include("modules\\ClassicalDynamics.jl")

""" Calculates Poincaré sections with Lyapunov exponents for various energies
    Parallel calculation, takes usually days to finish
"""

@everywhere function SolveItem(energy, parameters, dimension; file="Dicke.txt", path="", alreadySolved=[])    
    for (E, λ) in alreadySolved
        if isapprox(E, energy) && isapprox(parameters[1], λ)
            println("Skipped $parameters, E=$energy, dim=$dimension (already calculated)")
            return
        end
    end

    time = @elapsed sectionLyapunov, trajectories = SolveEnergy(energy, parameters, dimension, savePath=path)

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

    result = [energy, parameters[1], total, chaos, error, total > 0 ? meanLyapunov / total : 0, chaos > 0 ? meanLyapunovChaos / chaos : 0, maxLyapunov, myid(), time, trajectories]

    open(path * file, "a") do io
        println(io, result)
    end

    return
end

function RunMap(; δ=1.0, ω=1.0, ω₀=1.0, path="", dimension=101, step=0.1)
    path *= "Dicke_"

    file = "Map_dim=$(dimension)_$([δ, ω, ω₀]).txt"
    alreadySolved = ReadMap(path * file)

    println("Already calculated $(length(alreadySolved)) points.")

    input = [(energy, [λ, δ, ω, ω₀], dimension) for energy in 3:-step:-3, λ in 3:-step:0]
    pmap((args)->SolveItem(args...; file=file, path=path, alreadySolved=alreadySolved), input)    

    return
end

function ReadMap(file="")
    Eλ = []

    for line in eachline(file)     
        line = replace(line, "[" => "")   
        line = replace(line, "]" => "")
        line = replace(line, "," => "")
        elements = split(line)

        append!(Eλ, [(parse(Float64, elements[1]), parse(Float64, elements[2]))])
    end

    return Eλ
end