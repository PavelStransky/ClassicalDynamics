using Distributed
using Random
using Logging

workers = 24

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere include("models/TwoFluidLipkin.jl")
@everywhere include("modules/ClassicalDynamics.jl")

@everywhere disable_logging(Logging.Info)

""" Calculates Poincaré sections with Lyapunov exponents for various energies
    Parallel calculation, takes usually days to finish
"""

@everywhere function SolveItem(energy, parameters, dimension; file="2FL.txt", path="", alreadySolved=[])    
    for (E, λ) in alreadySolved
        if isapprox(E, energy) && isapprox(parameters[1], λ)
            println("Skipped $parameters, E=$energy, dim=$dimension (already calculated)")
            return false
        end
    end

    time = @elapsed averageLyapunov, freg, trajectories, lyapunovs = SolveEnergy(energy, parameters, dimension, savePath=path, timeout=10000, showFigures=false, randomize=true, tolerance=5e-9, sectionCoordinateX=1, sectionCoordinateY=3, sectionPlane=2)

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

    varλ = 0
    if length(lyapunovs) == 0
        maxλ = 0
        meanλ = 0
    else
        maxλ = maximum(lyapunovs)
        meanλ = mean(lyapunovs)
        if length(lyapunovs) > 1
            varλ = var(lyapunovs)
        end
    end

    result = [energy, parameters[1], total, chaos, error, freg, total > 0 ? meanLyapunov / total : 0, chaos > 0 ? meanLyapunovChaos / chaos : 0, length(lyapunovs), maxλ, meanλ, varλ, myid(), time, trajectories]

    open(path * file, "a") do io
        println(io, result)
    end

    return true
end

function RunMap(; η1=0.5, η2=-0.5, path="", dimension=101, step=0.1)
    path *= "2FT_"

    file = "Map_dim=$(dimension)_x1=$(η1)_x2=$(η2).txt"
    alreadySolved = ReadMap(path * file)

    input = shuffle([(energy, [ξ, η1, η2], dimension) for energy in -1:step:1, ξ in 0:step:1])

    println("Path: $path")
    println("File name: $file")
    println("To be calculated $(length(input)).")
    println("Already calculated $(length(alreadySolved)) points.")

    pmap((args)->SolveItem(args...; file=file, path=path, alreadySolved=alreadySolved), input)

    return
end

""" Calculates freg for one given set of parameters """
function RunEnergy(; ξ=0.5, η1=0.5, η2=-0.5, path="", dimension=101, step=0.01)
    path *= "2FL_"

    file = "Energy_dim=$(dimension)_$([ξ, η1, η2]).txt"
    alreadySolved = ReadMap(path * file)

    input = [(energy, [ξ, η1, η2], dimension) for energy in -1:step:1]

    println("To be calculated $(length(input)).")
    println("Already calculated $(length(alreadySolved)) points.")

    pmap((args)->SolveItem(args...; file=file, path=path, alreadySolved=alreadySolved), input)    

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