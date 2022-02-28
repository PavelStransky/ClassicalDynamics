using DifferentialEquations, DiffEqCallbacks
using LinearAlgebra
using Random
using Statistics
using Plots
using Printf
using Logging

LogLevel(Logging.Debug)
pyplot(size = (1200,1000))

mutable struct IntegrationParameters
    dimension
    modelParameters
    energy
    relaxationTime
    relativeFluctuationThreshold
    regularThreshold
    maxPSPoints
    psPoints
    lyapunovExponent
    historyLyapunovExponent
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

""" Poincaré section conditions """
function SectionCondition(sectionPlane) 
    return (u, t, integrator) -> u[sectionPlane]
end

""" Hit Poincaré section - ContinuousCallback """
function Section!(integrator)
    integrationParameters = integrator.p
    maxPSPoints = integrationParameters.maxPSPoints

    integrationParameters.psPoints += 1
    if integrationParameters.psPoints >= maxPSPoints
        return terminate!(integrator)
    end

    lyapunovs = integrationParameters.historyLyapunovExponent
    lv = var(lyapunovs)
    lm = mean(lyapunovs)

    if (lm > 0) && ((lv / lm < integrationParameters.relativeFluctuationThreshold) || (lm < integrationParameters.regularThreshold))
        return terminate!(integrator)
    end
end

""" Energy conservation - ManifoldProjection """
function EnergyConservation!(resid, u, integrationParameters, t)
    resid[1] = Energy(u, integrationParameters.modelParameters) - integrationParameters.energy
    resid[2:end] .= 0
end

function FindMissingCoordinate(i, j, k)
    x = zeros(Float64, 4)
    x[i] = 1
    x[j] = 1
    x[k] = 1
    return argmin(x)
end

""" Calculates Poincaré section for a given energy and number of trajectories """
function PoincareSection(energy, parameters, numTrajectories; 
        minimumBound=-2.0, 
        maximumBound=2.0, 
        sectionCoordinateX=2,
        sectionCoordinateY=4,
        sectionPlane=3,
        maximumInitialConditionsNumber=1000,
        showFigures=true, 
        kwargs...
    )
    
    figure = [scatter()]                  # We shall gradually fill the graph with chaotic trajectories
    all = scatter()

    for i = 2:4
        push!(figure, scatter())
    end

    missingCoordinate = FindMissingCoordinate(sectionCoordinateX, sectionCoordinateY, sectionPlane)
    @info "Missing coordinate = $missingCoordinate"

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

        solution, lyapunov = TrajectoryLyapunov(x, parameters; maxSectionPoints=1000, sectionPlane=sectionPlane, showFigures=false, relativeFluctuationThreshold=0, regularThreshold=0, kwargs...)
        if lyapunov <= 0
            continue                    # Unstable or nonconvergent trajectory
        end

        @info "Λ = $lyapunov, points = $(length(solution))"

        trajectory += 1

        initialConditionIndex = zeros(Int64, length(solution))
        for i = 1:length(solution)
            initialConditions = InitialConditions(solution.u[i][1:4], energy, parameters, missingCoordinate)

            if length(initialConditions) == 0
                initialConditionIndex[i] = 1
                @info "No initial condition found for point $(solution.u[i][1:4]) with E = $(Energy(solution.u[i][1:4], parameters))"
                continue
            end

            distance = abs.(initialConditions .- solution.u[i][missingCoordinate])
            initialConditionIndex[i] = argmin(distance)

            if minimum(distance) > 1e-3
                @info "Minimum distance too big: $(solution.u[i][1:4]), $distance"
            end
        end

        @info "Initial condition index: $(maximum(initialConditionIndex))."

        resultX = Array{Array{Float64, 1}, 1}(undef, maximum(initialConditionIndex))
        resultY = Array{Array{Float64, 1}, 1}(undef, maximum(initialConditionIndex))
        for i = 1:maximum(initialConditionIndex)
            resultX[i] = zeros(Float64, 0)
            resultY[i] = zeros(Float64, 0)
        end

        for i = 1:length(solution)
            append!(resultX[initialConditionIndex[i]], solution.u[i][sectionCoordinateX])
            append!(resultY[initialConditionIndex[i]], solution.u[i][sectionCoordinateY])
        end

        for i = (length(figure) + 1):maximum(initialConditionIndex)
            push!(figure, scatter())
        end

        for i = 1:maximum(initialConditionIndex)
            figure[i] = scatter!(figure[i], resultX[i], resultY[i], label="Λ=$(round(lyapunov, digits=3))", markersize=3, markeralpha=0.8, markerstrokewidth=0)
            all = scatter!(all, resultX[i], resultY[i], label=nothing, markersize=3, markeralpha=0.8, markerstrokewidth=0)
        end

        if showFigures
            display(plot(figure[1], figure[2], figure[3], all))
        end
    end

    return figure
end


""" Calculates Lyapunov exponent of an individual trajectory with given initial conditions.
    Works with systems with f=2 degrees of freedom and maximum 4 external parameters!
"""
function TrajectoryLyapunov(initialCondition, parameters; 
        solver=TsitPap8(), 
        sectionPlane=3, 
        savePath=nothing, 
        showFigures=false, 
        maxIterations=2E6, 
        relaxationTime=10, 
        regularThreshold=1e-3, 
        relativeFluctuationThreshold=1e-5, 
        maxSectionPoints=10000, 
        tolerance=1e-10, 
        saveStep=2, 
        lyapunovSeriesLength=200,
        timeInterval = (0, 1e5)
    )

    phaseSpaceDimension = length(initialCondition)
    x0 = zeros(phaseSpaceDimension * (phaseSpaceDimension + 1))
    x0[1:phaseSpaceDimension] = initialCondition
    x0[(phaseSpaceDimension + 1):end] = Matrix{Float64}(I, phaseSpaceDimension, phaseSpaceDimension)

    energy = Energy(x0, parameters)

    integrationParameters = IntegrationParameters(phaseSpaceDimension, parameters, energy, relaxationTime, relativeFluctuationThreshold, regularThreshold, maxSectionPoints, 0, 1, rand(lyapunovSeriesLength))

    sectionCondition = SectionCondition(sectionPlane) 
    lyapunovs = SavedValues(Float64, Float64)                                              # For a graph with the time evolution of Lyapunov exponents

    # Main part - calling the ODE solver
    fnc = ODEFunction(EquationOfMotion!)
    problem = ODEProblem(fnc, x0, timeInterval, integrationParameters)
    callback = CallbackSet(ManifoldProjection(EnergyConservation!, save=false), SavingCallback(RescaleΦ!, lyapunovs, saveat=saveStep:saveStep:1e6), ContinuousCallback(sectionCondition, Section!, nothing, save_positions=(false, true)))
    time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, callback=callback, save_on=true, save_everystep=false, save_start=false, save_end=false, maxiters=maxIterations, isoutofdomain=CheckDomain, verbose=true)

    # Print results
    lyapunov = mean(integrationParameters.historyLyapunovExponent)
    lv = var(integrationParameters.historyLyapunovExponent)

    if length(solution) > 0
        @info "Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(Energy(solution[end], parameters)), PS points = $(integrationParameters.psPoints) Λ = $lyapunov ± $lv"
    end

    # Save all unstable or nonconvergent trajectories (for debug reasons)
    if !(solution.retcode == :Success || solution.retcode == :Terminated || (integrationParameters.psPoints >= maxPSPoints && relativeFluctuationThreshold >= 0))
        if !isnothing(savePath)
            open(savePath * "Nonconvergent_Trajectories.txt", "a") do io
                println(io, "$parameters\t$energy\t$initialCondition\t$(solution.retcode)")
            end
        end

        @warn "Nonconvergent trajectory with initialCondition = $initialCondition"
        return [], 0                # Unstable trajectories can be distinguished later because their LE is exactly 0
    end

    if showFigures
        figure = plot(lyapunovs.t, lyapunovs.saveval, lw=2, title="Lyapunov = $lyapunov ± $lv")
        display(figure)
    end

    return solution, lyapunov, zip(lyapunovs.t, lyapunovs.saveval)
end

#=
""" Solve all trajectories for a given energy at a lattice x, y = (1...dimension, 1...dimension) """
function SolveEnergy(energy, parameters, dimension; 
        minimum=-2.0, 
        maximum=2.0, 
        latticeCoordinate1=1,
        latticeCoordinate2=3,
        sectionPlane=2,
        missingCoordinate=4,
        regularLyapunov=0.01, 
        showFigures=false, 
        verbose=false, 
        randomize=false, 
        timeout=0,
        kwargs...
    )

    # Find the number of distinct IC
    maxNumIC = 0
    for ip = 1:dimension, iq = 1:dimension
        coordinate1 = (maximum - minimum) * (ip - 0.5) / dimension  + minimum
        coordinate2 = (maximum - minimum) * (iq - 0.5) / dimension  + minimum

        x0 = zeros(Float64, 4)
        x0[latticeCoordinate1] = coordinate1
        x0[latticeCoordinate2] = coordinate2

        maxNumIC = max(maxNumIC, length(InitialConditions(x0, energy, parameters, missingCoordinate)))
    end

    if maxNumIC == 0
        return nothing

    println("Number of different IC: $maxNumIC")

    averageLyapunov = zeros(Float64, maxNumIC, dimension, dimension)
    countLyapunov = zeros(Int32, maxNumIC, dimension, dimension)
    fregSection = zeros(Float64, maxNumIC, dimension, dimension)

    maximumLyapunov = 0.0

    # Auxiliary variables just for info about the progress
    trajectories = 0
    infoInterval = 120E9
    startTime = time_ns()
    lastTime = time_ns()

    crossings = 0

    println("Starting SolveEnergy: $parameters, E = $energy, dim = $dimension:")
    for iic = 1:maxNumIC, ip = 1:dimension, iq = 1:dimension
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

            @printf("%.0f%% (%.0fs, %.0fs) %d trajectories (%.0f%% of the section covered)\n", 100 * (ip + iq / dimension) / dimension, (time - startTime) / 1E9, (time - lastTime) / 1E9, trajectories, 100 * numNonzero / (dimension * dimension))

            lastTime = time_ns()
        end

        if countLyapunov[iic, ip, iq] > 0
            continue
        end
        
        # Randomize - IC selected randomly from a calculated cell, otherwise it is always taken from the cell centre
        coordinate1 = (maximum - minimum) * (ip - (randomize ? rand() : 0.5)) / dimension  + minimum
        coordinate2 = (maximum - minimum) * (iq - (randomize ? rand() : 0.5)) / dimension  + minimum
    
        x = zeros(Float64, 4)
        x[latticeCoordinate1] = coordinate1
        x[latticeCoordinate2] = coordinate2

        ic = InitialConditions(x, energy, parameters, missingCoordinate)
        if length(ic) >= iic
            x[missingCoordinate] = ic[iic]
            result, lyapunov = TrajectoryLyapunov(x, parameters; showFigures=showFigures, verbose=verbose, sectionPlane=sectionPlane, kwargs...)
        else
            result = []
            lyapunov = -0.01                        # Unable to find initial condition - return small negative lyapunov exponent (in order to distinguish later kinematically inaccessible area)
        end
    
        averageLyapunov[iic, ip, iq] = lyapunov
        maximumLyapunov = max(maximumLyapunov, lyapunov)
        countLyapunov[iic, ip, iq] = 1

        if lyapunov > 0
            crossings += 1
            trajectories += 1                       # Number of convergent trajectories

            if lyapunov < regularLyapunov
                fregSection[iic, ip, iq] = 1
            end
        else
            fregSection[iic, ip, iq] = -1
        end

        for p in ps
            ic = InitialConditions(p, energy, parameters, missingCoordinate)

            distance = abs.(ic - p[missingCoordinate])
            iiic = argmin(distance)

            if distance[iiic] > 1e-3 
                println("Index $iiic, distance $(distance[iiic]).")

            ii1 = convert(Int32, round((p[latticeCoordinate1] + 2) * (dimension - 1) / 4)) + 1
            ii2 = convert(Int32, round((p[latticeCoordinate2] + 2) * (dimension - 1) / 4)) + 1

            # Nonzero Lyapunov exponent is always preferred, so any value replaces cells with nonpositive LE
            if averageLyapunov[iiic, ii1, ii2] <= 0
                averageLyapunov[iiic, ii1, ii2] = lyapunov
                countLyapunov[iiic, ii1, ii2] = 1

                if lyapunov > 0 && lyapunov < regularLyapunov
                    fregSection[iiic, ii1, ii2] = 1
                end
            else
                averageLyapunov[iiic, ii1, ii2] += lyapunov
                countLyapunov[iiic, ii1, ii2] += 1

                if lyapunov > 0 && lyapunov < regularLyapunov
                    fregSection[iiic, ii1, ii2] += 1
                end
            end

            crossings += 1
        end
    end

    replace!(countLyapunov, 0=>1)                   # Missing trajectories - avoid division by zero
    averageLyapunov ./= countLyapunov
    fregSection ./= countLyapunov
    
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

    if showFigures
        display(contourf(averageLyapunov, title="E = $energy, trajectories = $trajectories"))
    end    
    
    @printf("%.0fs: Finished [λ, E] = [%.2f, %.2f], Λ = %.2f, freg = %.2f, %d trajectories, %d crossings\n", (time_ns() - startTime) / 1E9, parameters[1], energy, maximumLyapunov, freg, trajectories, crossings)

    return averageLyapunov, maximumLyapunov, freg, trajectories
end


""" Calculates Lyapunov exponents for a given energy and number of trajectories """
function LyapunovExponents(energy, parameters, numTrajectories; 
        min=-2.0, 
        max=2.0, 
        sectionPlane=2, 
        showFigures=true, 
        kwargs...
    )
    pannel1 = scatter()
    pannel2 = plot()
    figure = plot()

    trajectory = 0
    while trajectory < numTrajectories
        P = (max - min) * rand() + min
        Q = (max - min) * rand() + min

        x = SectionPlane(sectionPlane, P, Q)
        if !InitialCondition!(x, energy, parameters)
            continue
        end
        
        x = collect(skipmissing(x))

        ps, lyapunov, lyapunovs = TrajectoryLyapunov(x, parameters; showFigures=false, sectionPlane=sectionPlane, kwargs...)

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

    fnc = ODEFunction(EquationOfMotion!)
    problem = ODEProblem(fnc, x0, timeInterval, parameters)
    time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, maxiters=maxIterations, isoutofdomain=CheckDomain, verbose=verbose)

    println("Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(Energy(solution[end], parameters))")

    # Show both pairs of canonically conjugated variables
    pannel1 = plot(solution[1,:], solution[3,:], lw=2)
    pannel2 = plot(solution[2,:], solution[4,:], lw=2)
    display(plot(pannel1, pannel2, layout=2))
end

""" Runs various methods to solve differential equations for the trajectory 
    For test purposes only!
""" 
function TestTrajectory(x, parameters)
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

    for method in explicitRK
        x0 = zeros(20)
        x0[1:4] = x
        x0[5:end] = Matrix{Float64}(I, 4, 4)

        resize!(parameters, 4)
        params = [parameters..., Energy(x0, parameters), 0, 1e-5, 0.001, 10000, 0, 1]         
        # λ, δ, ω, ω₀, 5=energy, 6=relaxationTime, 7=relativeFluctuationThreshold, 8=regularThreshold, 9=maximum number of points in Poincare section, 
        #              10=num points in the Poincaré section, 11=current Lyapunov exponent 
        append!(params, rand(200))                                            

        fnc = ODEFunction(EquationOfMotion!)
        problem = ODEProblem(fnc, x0, (0, 1e6), params)

        lyapunovs = SavedValues(Float64, Float64)                                              # For a graph of Lyapunov exponents
        callback = CallbackSet(SavingCallback(rescale!, lyapunovs, saveat=2:2:1e6), ContinuousCallback(sectionCondition2, section!, nothing, save_positions=(false, true)), ManifoldProjection(energyConservation!, save=false))
        time = @elapsed solution = solve(problem, method, reltol=1e-8, abstol=1e-8, callback=callback, save_on=true, save_everystep=false, save_start=false, save_end=false, maxiters=1E8, isoutofdomain=CheckDomain, verbose=true)    

        lyapunov = mean(params[12:end])
        lv = var(params[12:end])

        println("Time = $time, Lyapunov = $lyapunov ± $lv")
        println(solution.retcode)

        p = plot!(lyapunovs.t, lyapunovs.saveval, lw=2, label="$method l=$(round(lyapunov, digits=3)) t=$(round(time, digits=1)) $(solution.retcode)", ylims=(0,1))
        display(p)
        #savefig(p, "d:\\results\\rosenbrockW.png")
    end
end
=#
