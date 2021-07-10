using Distributed
using Random

workers = 8

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere include("models/Vibron.jl")
@everywhere include("modules/ClassicalDynamics.jl")

""" Calculates Poincaré sections with Lyapunov exponents for various energies
    Parallel calculation, takes usually days to finish
"""

@everywhere function SolveItem(energy, parameters, dimension; file="Vibron.txt", path="", alreadySolved=[], sectionPlane=2)    
    for (E, λ) in alreadySolved
        if isapprox(E, energy) && isapprox(parameters[1], λ)
            println("Skipped $parameters, E=$energy, dim=$dimension (already calculated)")
            return false
        end
    end

    time = @elapsed averageLyapunov, maximumLyapunov, freg, trajectories = SolveEnergy(energy, parameters, dimension, savePath=path, sectionPlane=sectionPlane, min=-sqrt(2.0), max=sqrt(2.0), timeout=7200, randomize=true)

    chaos = 0
    total = 0
    error = 0
    meanLyapunov = 0
    meanLyapunovChaos = 0

    for x in averageLyapunov
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
    end

    result = [energy, parameters[1], total, chaos, error, freg, total > 0 ? meanLyapunov / total : 0, chaos > 0 ? meanLyapunovChaos / chaos : 0, maximumLyapunov, myid(), time, trajectories]

    open(path * file, "a") do io
        println(io, result)
    end

    return true
end

function RunMapC(; C=0.2, path="", dimension=101, step=0.1, sectionPlane=1)
    path *= "Vibron_"

    file = "Map_dim=$(dimension)_$C.txt"
    alreadySolved = ReadMap(path * file)

    input = shuffle([(energy, [A, A-1, C], dimension) for energy in -1.2:step:1, A in 0:step:1])

    println("To be calculated $(length(input)).")
    println("Already calculated $(length(alreadySolved)) points.")

    pmap((args)->SolveItem(args...; file=file, path=path, alreadySolved=alreadySolved, sectionPlane=sectionPlane), input)

    return
end

function ReadMap(file="")
    Eλ = []

    try
        for line in eachline(file)     
            line = replace(line, "[" => "")   
            line = replace(line, "]" => "")
            line = replace(line, "," => "")
            elements = split(line)

            append!(Eλ, [(parse(Float64, elements[1]), parse(Float64, elements[2]))])
        end
    catch x
    end

    return Eλ
end