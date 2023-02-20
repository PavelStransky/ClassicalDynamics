using DifferentialEquations, DiffEqCallbacks
using LinearAlgebra
using Random
using Statistics
using Plots
using Printf
using Logging
using ColorSchemes

pyplot(size = (1200,1000))

mutable struct LyapunovIntegrationParameters
    dimension
    modelParameters
    energy
    relaxationTime
    relativeFluctuationThreshold
    regularThreshold
    maximumSectionPoints
    sectionPoints
    startTime
    timeout
    result
    lyapunovExponent
    historyLyapunovExponent
end    

struct SimpleIntegrationParameters
    modelParameters
    energy
end

""" 
    Rescales the Φ matrix and saves the Lyapunov exponent
    Called periodically (every "saveStep" time units) - SavingCallback
"""
function RescaleΦ!(u, t, integrator)
    integrationParameters = integrator.p
    dimension = integrationParameters.dimension

    Φ = reshape(integrator.u[(dimension + 1):end], dimension, dimension)
    
    maxEigenvalue = sqrt(maximum(eigvals(Φ * transpose(Φ))))
    integrator.u[(dimension + 1):end] /= maxEigenvalue

    lyapunov = 0

    relaxationTime = integrationParameters.relaxationTime
    if t > relaxationTime
        integrationParameters.lyapunovExponent += log(maxEigenvalue)
        lyapunov = integrationParameters.lyapunovExponent / (t - relaxationTime)

        lyapunovs = integrationParameters.historyLyapunovExponent[2:end]
        append!(lyapunovs, lyapunov)
        integrationParameters.historyLyapunovExponent = lyapunovs
    end

    return lyapunov
end

function TimeoutCondition(u, t, integrator)
    integrationParameters = integrator.p
    if(integrationParameters.timeout > 0 && time_ns() - integrationParameters.startTime > integrationParameters.timeout)
        integrationParameters.result = :Timeout
        return true
    end
    return false
end


""" Poincaré section conditions """
function SectionCondition(sectionPlane) 
    return (u, t, integrator) -> u[sectionPlane]
end


""" Hit Poincaré section - ContinuousCallback """
function Section!(integrator)
    integrationParameters = integrator.p
    maximumSectionPoints = integrationParameters.maximumSectionPoints

    integrationParameters.sectionPoints += 1
    if integrationParameters.sectionPoints >= maximumSectionPoints
        integrationParameters.result = :MaximumSectionPoints
        return terminate!(integrator)
    end

    lyapunovs = integrationParameters.historyLyapunovExponent
    lv = var(lyapunovs)
    lm = mean(lyapunovs)

    if (lm > 0) && ((lv / lm < integrationParameters.relativeFluctuationThreshold) || (lm < integrationParameters.regularThreshold))
        integrationParameters.result = :Converged
        return terminate!(integrator)
    end
end


""" Energy conservation - ManifoldProjection """
function EnergyConservation!(resid, u, integrationParameters, t)
    resid[1] = Energy(u, integrationParameters.modelParameters) - integrationParameters.energy
    resid[2:end] .= 0
end

""" The fourth coordinate completing the three given """
function FindMissingCoordinate(i, j, k)
    x = zeros(Float64, 4)
    x[i] = 1
    x[j] = 1
    x[k] = 1
    
    missingCoordinate = argmin(x)

    @info "Missing coordinate = $missingCoordinate"

    return missingCoordinate
end

""" Calculates Lyapunov exponent of an individual trajectory with given initial conditions.
    Works with systems with f=2 degrees of freedom and maximum 4 external parameters!
"""
function TrajectoryLyapunov(initialCondition, parameters; 
        solver=TsitPap8(), 
        sectionPlane=3, 
        savePath=nothing, 
        showFigures=false, 
        timeout=0,
        maximumIterations=2E6, 
        relaxationTime=100, 
        regularThreshold=1e-3, 
        relativeFluctuationThreshold=1e-5, 
        maximumSectionPoints=20000, 
        tolerance=1e-10, 
        saveStep=2, 
        lyapunovSeriesLength=500,
        timeInterval = (0, 1e5)
    )

    phaseSpaceDimension = length(initialCondition)

    # Initial condition, including the tangent dynamics
    x0 = zeros(phaseSpaceDimension * (phaseSpaceDimension + 1))
    x0[1:phaseSpaceDimension] = initialCondition
    x0[(phaseSpaceDimension + 1):end] = Matrix{Float64}(I, phaseSpaceDimension, phaseSpaceDimension)

    energy = Energy(x0, parameters)

    integrationParameters = LyapunovIntegrationParameters(phaseSpaceDimension, parameters, energy, relaxationTime, relativeFluctuationThreshold, regularThreshold, maximumSectionPoints, 0, time_ns(), 1E9 * timeout, :Start, 1, rand(lyapunovSeriesLength))

    sectionCondition = SectionCondition(sectionPlane) 
    lyapunovs = SavedValues(Float64, Float64)                                              # For a graph with the time evolution of Lyapunov exponents

    # Main part - calling the ODE solver
    fnc = ODEFunction(EquationOfMotion!)
    problem = ODEProblem(fnc, x0, timeInterval, integrationParameters)
    callback = CallbackSet(ManifoldProjection(EnergyConservation!, save=false), SavingCallback(RescaleΦ!, lyapunovs, saveat=saveStep:saveStep:1e6), ContinuousCallback(sectionCondition, Section!, nothing, save_positions=(false, true)), DiscreteCallback(TimeoutCondition, terminate!))
    time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, callback=callback, save_on=true, save_everystep=false, save_start=false, save_end=false, maxiters=maximumIterations, isoutofdomain=CheckDomain, verbose=true)

    # Print results
    lyapunov = mean(integrationParameters.historyLyapunovExponent)
    lv = var(integrationParameters.historyLyapunovExponent)

    if length(solution) > 0
        @info "Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(Energy(solution[end], parameters)), PS points = $(integrationParameters.sectionPoints) Λ = $lyapunov ± $lv"
    end

    @debug "retcode = $(solution.retcode), result = $(integrationParameters.result)"

    # Save all unstable or nonconvergent trajectories (for debug reasons)
    if !(solution.retcode == :Success || (solution.retcode == :Terminated && (integrationParameters.result == :Converged || integrationParameters.result == :MaximumSectionPoints)))
        if !isnothing(savePath)
            open(savePath * "Nonconvergent_Trajectories.txt", "a") do io
                println(io, "$parameters\t$energy\t$initialCondition\t$(solution.retcode)")
            end
        end

        @info "Nonconvergent trajectory with initialCondition = $initialCondition: retcode = $(solution.retcode), result = $(integrationParameters.result)"
        return [], 0, []    # Unstable trajectories can be distinguished later because their LE is exactly 0
    end

    if showFigures
        figure = plot(lyapunovs.t, lyapunovs.saveval, lw=2, title="Λ = $lyapunov ± $lv", label=nothing, xlabel="t", ylabel="Λ")
        display(plot(figure))
        if !isnothing(savePath)
            savefig(figure, savePath * "$initialCondition.png")
        end
    end

    return solution, lyapunov, zip(lyapunovs.t, lyapunovs.saveval)
end

""" Calculates Poincaré section for a given energy and number of trajectories """
function PoincareSection(energy, parameters, numTrajectories; 
        minimumBound=-sqrt(2.0), 
        maximumBound=sqrt(2.0), 
        sectionCoordinateX=2,
        sectionCoordinateY=4,
        sectionPlane=3,
        maximumInitialConditionsNumber=100,
        showFigures=true,
        kwargs...
    )
    
    figure = scatter()                  # We shall gradually fill the graph with chaotic trajectories

    missingCoordinate = FindMissingCoordinate(sectionCoordinateX, sectionCoordinateY, sectionPlane)

    trajectory = 0
    initialConditionsNumber = 0
    while trajectory < numTrajectories && initialConditionsNumber < maximumInitialConditionsNumber
        coordinateX = (maximumBound - minimumBound) * rand() + minimumBound
        coordinateY = (maximumBound - minimumBound) * rand() + minimumBound

        x = zeros(Float64, 4)
        x[sectionCoordinateX] = coordinateX
        x[sectionCoordinateY] = coordinateY

        initialConditions = InitialConditions(x, energy, parameters, missingCoordinate)
        initialConditionsNumber += 1

        if length(initialConditions) == 0
            continue                    # Failure of generating an initial condition
        end

        x[missingCoordinate] = rand(initialConditions)

        @info "$(length(initialConditions)) IC = $x, E = $(Energy(x, parameters))"

        solution, lyapunov = TrajectoryLyapunov(x, parameters; sectionPlane=sectionPlane, showFigures=false, relativeFluctuationThreshold=0, regularThreshold=0, kwargs...)
        if lyapunov <= 0
            continue                    # Unstable or nonconvergent trajectory
        end

        @info "Λ = $lyapunov, points = $(length(solution))"
        trajectory += 1

        if showFigures
            figure = scatter!(figure, solution[sectionCoordinateX,:], solution[sectionCoordinateY,:], label=nothing, markersize=3, markeralpha=0.8, markerstrokewidth=0)
            display(plot(figure))
        end
    end

    return figure
end


""" Solve all trajectories for a given energy at a lattice x, y = (1...dimension, 1...dimension) """
function SolveEnergy(energy, parameters, dimension; 
        minimumBound=-sqrt(2.0), 
        maximumBound=sqrt(2.0), 
        sectionCoordinateX=2,
        sectionCoordinateY=4,
        sectionPlane=3,
        regularThreshold=0.01, 
        maximumλ=0.25,
        showFigures=true,
        savePath=nothing,
        randomize=false, 
        timeout=0,
        kwargs...
    )

    missingCoordinate = FindMissingCoordinate(sectionCoordinateX, sectionCoordinateY, sectionPlane)

    averageLyapunov = zeros(Float64, dimension, dimension)
    countLyapunov = zeros(Int32, dimension, dimension)
    fregSection = fill(-1.0, (dimension, dimension))

    # All calculated Lyapunov exponents
    lyapunovs = zeros(Float64, 0)
    calculationTimes = zeros(Float64, 0)

    # Auxiliary variables just for info about the progress
    trajectories = 0
    infoInterval = 120E9
    startTime = time_ns()
    lastTime = time_ns()

    crossings = 0

    @info "Starting SolveEnergy: $parameters, E = $energy, dim = $dimension:"
    for ix = 1:dimension, iy = 1:dimension
        time = time_ns()

        if timeout > 0 && time - startTime > 1E9 * timeout
            print("TIMEOUT: ")
            break
        end

        if time - lastTime > infoInterval
            numNonzero = 0
            for i in averageLyapunov
                if i != 0
                    numNonzero += 1
                end
            end

            @printf("%.0f%% (%.0fs, %.0fs) %d trajectories (%.0f%% of the section covered)\n", 100 * (ix + iy / dimension) / dimension, (time - startTime) / 1E9, (time - lastTime) / 1E9, trajectories, 100 * numNonzero / (dimension * dimension))

            lastTime = time_ns()
        end

        if countLyapunov[ix, iy] > 0
            continue
        end
        
        # Randomize - IC selected randomly from a calculated cell, otherwise it is always taken from the cell centre
        x = (maximumBound - minimumBound) * (ix - (randomize ? rand() : 0.5)) / dimension  + minimumBound
        y = (maximumBound - minimumBound) * (iy - (randomize ? rand() : 0.5)) / dimension  + minimumBound
    
        ic = zeros(Float64, 4)
        ic[sectionCoordinateX] = x
        ic[sectionCoordinateY] = y

        m = InitialConditions(ic, energy, parameters, missingCoordinate)

        result = []
        lyapunov = -0.01                        # Unable to find initial condition - return small negative lyapunov exponent (in order to distinguish later kinematically inaccessible area)            

        for n in m
            ic[missingCoordinate] = n

            calculationTimeout = 0
            if length(calculationTimes) > 0
                calculationTimeout = 10 * mean(calculationTimes)
            end

            @info "Initial conditions = $ic, timeout = $calculationTimeout"

            calculationTime = @elapsed result, lyapunov = TrajectoryLyapunov(ic, parameters; timeout=calculationTimeout, savePath=savePath, showFigures=false, sectionPlane=sectionPlane, regularThreshold=regularThreshold, kwargs...)

            if lyapunov > 0 && length(result) > 0
                append!(calculationTimes, calculationTime)
                break
            else
                @info "Trajectory doesn't cross the section"
            end
        end
    
        averageLyapunov[ix, iy] = lyapunov
        append!(lyapunovs, lyapunov)

        if lyapunov > 0
            crossings += 1
            trajectories += 1                       # Number of convergent trajectories

            countLyapunov[ix, iy] = 1

            if lyapunov < regularThreshold
                fregSection[ix, iy] = 1
            else
                fregSection[ix, iy] = 0
            end
        else
            continue
        end

        currentTrajectory = zeros(Float64, dimension, dimension)
        currentTrajectory[ix, iy] = 1

        for i = 1:length(result)
            x = result.u[i][sectionCoordinateX]
            y = result.u[i][sectionCoordinateY]

            iix = convert(Int32, floor((x - minimumBound) * dimension / (maximumBound - minimumBound))) + 1
            iiy = convert(Int32, floor((y - minimumBound) * dimension / (maximumBound - minimumBound))) + 1

            if iix <= 0 || iix > dimension || iiy <= 0 || iiy > dimension
                continue
            end

            if averageLyapunov[iix, iiy] <= 0           # Nonzero Lyapunov exponent is always preferred, so any value replaces cells with nonpositive LE
                averageLyapunov[iix, iiy] = lyapunov
                countLyapunov[iix, iiy] = 1
                crossings += 1

                if fregSection[iix, iiy] < 0
                    fregSection[iix, iiy] = 0
                end
                if lyapunov > 0 && lyapunov < regularThreshold
                    fregSection[iix, iiy] += 1
                end

            elseif currentTrajectory[iix, iiy] == 0     # No visit yet to this point
                averageLyapunov[iix, iiy] += lyapunov
                countLyapunov[iix, iiy] += 1
                crossings += 1

                if fregSection[iix, iiy] < 0
                    fregSection[iix, iiy] = 0
                end
                if lyapunov > 0 && lyapunov < regularThreshold
                    fregSection[iix, iiy] += 1
                end
            end

            currentTrajectory[iix, iiy] = 1
        end
    end

    cl = replace(countLyapunov, 0=>1)                   # Missing trajectories - avoid division by zero
    averageLyapunov ./= cl
    fregSection ./= cl
    
    total = 0.0
    count = 0.0

    for x in fregSection
        if x < 0.0
            continue
        end

        total += 1.0
        count += x
    end

    freg = total > 0.0 ? count / total : 0.0        # Fraction of regularity (0 - totally chaotic, 1 - totally regular)

    lyapunovs = filter(x -> x > 0, lyapunovs)
    chaoticLyapunovs = filter(x -> x > regularThreshold, lyapunovs)

    xs = LinRange(minimumBound, maximumBound, dimension)

    if showFigures || !isnothing(savePath)
        al = deepcopy(averageLyapunov)
        al[al .>= maximumλ] .= maximumλ

        cl = replace(countLyapunov, 0=>-1)
        cl = cl  ./ maximum(cl)

        pannel1 = contourf(xs, xs, al, levels=LinRange(0, maximumλ, 100), c=:rainbow, clim=(0, maximumλ), title="Average λ [E = $(round(energy, digits=3)), trajectories = $trajectories]")
        pannel2 = contourf(xs, xs, fregSection, levels=LinRange(-0.1, 1, 100), c=:rainbow, clim=(-0.1, 1), title="freg = $(round(freg, digits=3)) (T = $(round((time_ns() - startTime) / 1E9, digits=0)))")
        pannel3 = contourf(xs, xs, cl, levels=LinRange(0, 1, 100), c=:rainbow, clim=(0, 1), title="Count (max = $(maximum(countLyapunov)))")
        if length(chaoticLyapunovs) > 0
            pannel4 = histogram(chaoticLyapunovs, bins=0:0.005:0.3, label=nothing, title="Λ = $(round(maximum(lyapunovs), sigdigits=3)), λ = $(round(mean(chaoticLyapunovs), sigdigits=3)) ± $(round(var(chaoticLyapunovs), sigdigits=2))")
        else
            pannel4 = plot()
        end

        figure = plot(pannel1, pannel2, pannel3, pannel4)

        if showFigures
            display(plot(figure))
        end
        
        if !isnothing(savePath)
            savefig(figure, savePath * "PS($(sectionCoordinateX)$(sectionCoordinateY)$(sectionPlane))_$(parameters)_E=$(round(energy, digits=3)).png")
        end
    end    
    
    if length(chaoticLyapunovs) > 0
        maximumLyapunov = maximum(chaoticLyapunovs)
    else
        maximumLyapunov = 0
    end

    @printf("%.0fs: Finished [λ, E] = [%.2f, %.2f], Λ = %.2f, freg = %.2f, %d trajectories, %d crossings\n", (time_ns() - startTime) / 1E9, parameters[1], energy, maximumLyapunov, freg, trajectories, crossings)

    return averageLyapunov, freg, trajectories, lyapunovs
end

""" Calculates Lyapunov exponents for a given energy and number of trajectories """
function LyapunovExponents(energy, parameters, numTrajectories; 
        minimumBound=-sqrt(2), 
        maximumBound=sqrt(2), 
        sectionCoordinateX=2,
        sectionCoordinateY=4,
        sectionPlane=3, 
        showFigures=true, 
        kwargs...
    )
    pannel1 = scatter()
    pannel2 = plot()
    figure = plot()

    missingCoordinate = FindMissingCoordinate(sectionCoordinateX, sectionCoordinateY, sectionPlane)
    @info "Missing coordinate = $missingCoordinate"

    trajectory = 0
    while trajectory < numTrajectories
        x = (maximumBound - minimumBound) * rand() + minimumBound
        y = (maximumBound - minimumBound) * rand() + minimumBound
   
        ic = zeros(Float64, 4)
        ic[sectionCoordinateX] = x
        ic[sectionCoordinateY] = y

        m = InitialConditions(ic, energy, parameters, missingCoordinate)

        if length(m) == 0
            continue                    # Failure of generating an initial condition
        end

        ic[missingCoordinate] = rand(m)

        @info "$(length(m)) IC = $ic, E = $(Energy(ic, parameters))"
        
        ps, lyapunov, lyapunovs = TrajectoryLyapunov(ic, parameters; showFigures=false, sectionPlane=sectionPlane, kwargs...)

        if lyapunov <= 0
            continue
        end

        trajectory += 1

        pannel1 = scatter!(pannel1, ps, label="Λ=$(round(lyapunov, digits=2))", markersize=3, markeralpha=0.8, markerstrokewidth=0)
        pannel2 = plot!(pannel2, collect(lyapunovs), lw=2, title="Lyapunov = $lyapunov")
        figure = plot(pannel1, pannel2, layout=2)

        if showFigures
            display(figure)
        end
    end

    return figure
end

""" Calculates an individual trajectory with given initial conditions"""
function Trajectory(initialCondition, parameters; 
        solver=TsitPap8(), 
        maxIterations=1E7, 
        timeInterval=(0.0, 100.0),
        verbose=true,
        tolerance=1e-10
    )

    x0 = zeros(20)
    x0[1:4] = initialCondition
    x0[5:end] = Matrix{Float64}(I, 4, 4)

    energy = Energy(x0, parameters)

    phaseSpaceDimension = length(initialCondition)
    integrationParameters = SimpleIntegrationParameters(parameters, energy)

    fnc = ODEFunction(EquationOfMotion!)
    problem = ODEProblem(fnc, x0, timeInterval, integrationParameters)
    time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, maxiters=maxIterations, isoutofdomain=CheckDomain, verbose=verbose)

    println("Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(Energy(solution[end], parameters))")

    # Show both pairs of canonically conjugated variables
    pannel1 = plot(solution, idxs=(1,3), lw=2)
    pannel2 = plot(solution, idxs=(2,4), lw=2)
    display(plot(pannel1, pannel2, layout=2))
end

""" Runs various methods to solve differential equations for the trajectory 
    For test purposes only!
""" 
function TestTrajectory(x, parameters; sectionPlane=3)
    pyplot(size = (2200,1000))
    display(plot())

    explicitRK = [DP5(), Tsit5(), TanYam7(), DP8(), TsitPap8(), Feagin10(), Feagin12(), Feagin14()]
    explicitRKlazy = [BS5(), Vern6(), Vern7(), Vern8(), Vern9()]
    multistepAB = [AB3(), AB4(), AB5(), ABM32(), ABM43(), ABM54()]
    adaptiveA = [VCAB3(), VCAB4(), VCAB5(), VCABM3(), VCABM4(), VCABM5(), VCABM(), AN5()]
    
    stiffdirk1 = [ImplicitEuler(), ImplicitMidpoint(), Trapezoid(), TRBDF2(), SDIRK2()]
    stiffdirk2 = [Kvaerno3(), Kvaerno4(), Kvaerno5(), Cash4(), KenCarp3(), KenCarp4(), KenCarp5(), Hairer4(), Hairer42()]
    rosenbrock = [ROS3P(), Rodas3(), RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(), Ros4LStab(), Rodas4(), Rodas42(), Rodas4P(), Rodas5()]
    rosenbrockW = [Rosenbrock23(), Rosenbrock32(), RosenbrockW6S4OS(), ROS34PW1a(), ROS34PW1b(), ROS34PW2(), ROS34PW3()]
    stabilizedRK = [ROCK2(), ROCK4(), RKC(), SERK2(), ESERK4(), ESERK5()]

    phaseSpaceDimension = length(initialCondition)
    energy = Energy(x, parameters)

    for method in explicitRK
        # Initial condition, including the tangent dynamics
        x0 = zeros(phaseSpaceDimension * (phaseSpaceDimension + 1))
        x0[1:phaseSpaceDimension] = initialCondition
        x0[(phaseSpaceDimension + 1):end] = Matrix{Float64}(I, phaseSpaceDimension, phaseSpaceDimension)

        integrationParameters = LyapunovIntegrationParameters(phaseSpaceDimension, parameters, energy, 0, 1E-5, 0.001, 10000, 0, time_ns(), 0, :Start, 1, rand(200))

        fnc = ODEFunction(EquationOfMotion!)
        problem = ODEProblem(fnc, x0, (0, 1e6), integrationParameters)

        lyapunovs = SavedValues(Float64, Float64)                                              # For a graph of Lyapunov exponents

        callback = CallbackSet(ManifoldProjection(EnergyConservation!, save=false), SavingCallback(RescaleΦ!, lyapunovs, saveat=2:2:1e6), DiscreteCallback(TimeoutCondition, terminate!))
        time = @elapsed solution = solve(problem, method, reltol=1e-8, abstol=1e-8, callback=callback, save_on=true, save_everystep=false, save_start=false, save_end=false, maxiters=1E8, isoutofdomain=CheckDomain, verbose=true)    

        lyapunov = mean(integrationParameters.historyLyapunovExponent)
        lv = var(integrationParameters.historyLyapunovExponent)
    
        println("Time = $time, Lyapunov = $lyapunov ± $lv")
        println(solution.retcode)

        p = plot!(lyapunovs.t, lyapunovs.saveval, lw=2, label="$method l=$(round(lyapunov, digits=3)) t=$(round(time, digits=1)) $(solution.retcode)", ylims=(0,1))
        display(p)
        #savefig(p, "d:\\results\\rosenbrockW.png")
    end
end
