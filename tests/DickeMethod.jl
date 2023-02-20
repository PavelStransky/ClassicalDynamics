include("../models/Dicke.jl")
include("../modules/ClassicalDynamics.jl")

parameters = [1.0, 1.0, 1.0, 1.0]
energy = 30
dimension = 51

explicitRK = [DP5(), Tsit5(), TanYam7(), DP8(), TsitPap8()] #, Feagin10(), Feagin12(), Feagin14()]
explicitRKlazy = [BS5(), Vern6(), Vern7(), Vern8(), Vern9()]
multistepAB = [AB3(), AB4(), AB5(), ABM32(), ABM43(), ABM54()]
adaptiveA = [VCAB3(), VCAB4(), VCAB5(), VCABM3(), VCABM4(), VCABM5(), VCABM(), AN5()]

tolerance = 1e-8
path = "d:\\results\\Dicke\\test\\"

for solver in explicitRKlazy
    solverName = String(Symbol(solver))
    solverName = SubString(solverName, 1, findfirst('(', solverName) - 1)

    psPath = path * solverName * "_" * string(tolerance) * "_"
    println(psPath)
    
    time = @elapsed averageLyapunov, freg, trajectories, lyapunovs = SolveEnergy(energy, parameters, dimension, solver=solver, tolerance=tolerance, savePath=psPath, sectionPlane=2, sectionCoordinateX=3, sectionCoordinateY=1, minimumBound=-2, maximumBound=2)  

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

    result = [energy, parameters[1], total, chaos, error, freg, total > 0 ? meanLyapunov / total : 0, chaos > 0 ? meanLyapunovChaos / chaos : 0, length(lyapunovs), maxλ, meanλ, varλ, time, trajectories, tolerance, solverName]

    open(path * "Dicke.txt", "a") do io
        println(io, result)
    end

    println()
end