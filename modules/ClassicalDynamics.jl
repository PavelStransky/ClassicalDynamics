using DifferentialEquations, DiffEqCallbacks
using LinearAlgebra
using Random
using Statistics
using Plots
using Printf

""" 
    Rescales the Φ matrix and saves the Lyapunov exponent
    Called periodically (every "saveStep" time units) - SavingCallback
"""
function rescale!(u, t, integrator)
    Φ = reshape(integrator.u[5:end], 4, 4)
    maxEigenvalue = sqrt(maximum(eigvals(Φ * transpose(Φ))))

    integrator.u[5:end] /= maxEigenvalue

    lyapunov = 0

    relaxationTime = integrator.p[6]
    if t > relaxationTime
        integrator.p[11] += log(maxEigenvalue)
        lyapunov = integrator.p[11] / (t - relaxationTime)

        lyapunovs = integrator.p[13:end]
        append!(lyapunovs, lyapunov)
        integrator.p[12:end] = lyapunovs
    end

    return lyapunov
end

""" Poincaré section conditions """
sectionCondition1(u, t, integrator) = u[1]
sectionCondition2(u, t, integrator) = u[2]

""" Hit Poincaré section - ContinuousCallback """
function section!(integrator)
    maxPSPoints = integrator.p[9]

    integrator.p[10] += 1
    if integrator.p[10] >= maxPSPoints
        return terminate!(integrator)
    end

    lyapunovs = integrator.p[12:end]
    lv = var(lyapunovs)
    lm = mean(lyapunovs)

    relativeFluctuationThreshold = integrator.p[7]
    regularThreshold = integrator.p[8]
    if (lm > 0) && ((lv / lm < relativeFluctuationThreshold) || (lm < regularThreshold))
        return terminate!(integrator)
    end
end

""" Energy conservation - ManifoldProjection """
function energyConservation!(resid, u, params, t)
    resid[1] = Energy(u, params) - params[5]
    resid[2:end] .= 0
end

""" Calculates Lyapunov exponent of an individual trajectory with given initial conditions.
    Works with systems with f=2 degrees of freedom and maximum 4 external parameters!
"""
function TrajectoryLyapunov(initialCondition, parameters; 
        solver=TsitPap8(), 
        sectionPlane=2, 
        savePath=nothing, 
        verbose=false, 
        showFigures=false, 
        maxIterations=2E6, 
        relaxationTime=10, 
        regularThreshold=1e-3, 
        relativeFluctuationThreshold=1e-5, 
        maxPSPoints=10000, 
        tolerance=1e-10, 
        saveStep=2, 
        lyapunovSeriesLength=200,
        timeInterval = (0, 1e5)
    )

    x0 = zeros(20)
    x0[1:4] = initialCondition
    x0[5:end] = Matrix{Float64}(I, 4, 4)

    energy = Energy(x0, parameters)

    resize!(parameters, 4)          # !!! Works for 4 parameters at most !!!
    params = [parameters..., energy, relaxationTime, relativeFluctuationThreshold, regularThreshold, maxPSPoints, 0, 1]   
                                                                                            # 1...4=parameters, 5=energy, 6=relaxationTime, 7=relativeFluctuationThreshold, 8=regularThreshold, 9=maximum number of points in Poincare section, 
                                                                                            #                   10=num points in the Poincaré section, 11=current Lyapunov exponent 
    append!(params, rand(lyapunovSeriesLength))                                             # Queue with latest Lyapunov exponents (initialized as a random series)    

    sectionCondition = sectionPlane == 1 ? sectionCondition1 : sectionCondition2
    lyapunovs = SavedValues(Float64, Float64)                                              # For a graph with the time evolution of Lyapunov exponents

    # Main part - calling the ODE solver
    fnc = ODEFunction(EquationOfMotion!)
    problem = ODEProblem(fnc, x0, timeInterval, params)
    callback = CallbackSet(ManifoldProjection(energyConservation!, save=false), SavingCallback(rescale!, lyapunovs, saveat=saveStep:saveStep:1e6), ContinuousCallback(sectionCondition, section!, nothing, save_positions=(false, true)))
    time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, callback=callback, save_on=true, save_everystep=false, save_start=false, save_end=false, maxiters=maxIterations, isoutofdomain=CheckDomain, verbose=verbose)

    # Print results
    lyapunov = mean(params[12:end])
    lv = var(params[12:end])

    if verbose && length(solution) > 0
        println("Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(Energy(solution[end], parameters)), PS points = $(params[10]) Λ = $lyapunov ± $lv")
    end

    # Save all unstable or nonconvergent trajectories (for debug reasons)
    if !(solution.retcode == :Success || solution.retcode == :Terminated || (params[10] >= maxPSPoints && relativeFluctuationThreshold >= 0))
        if !isnothing(savePath)
            open(savePath * "Nonconvergent_Trajectories.txt", "a") do io
                println(io, "$parameters\t$energy\t$initialCondition\t$(solution.retcode)")
            end
        end

        return [], 0                # Unstable trajectories can be distinguished later because their LE is exactly 0
    end

    result = sectionPlane == 1 ? zip(solution[2,:], solution[4,:]) : zip(solution[1,:], solution[3,:])

    if showFigures
        pannel1 = scatter(result, title="P = $(x0[1]), Q = $(x0[3]) [$(length(solution))]")
        pannel2 = plot(lyapunovs.t, lyapunovs.saveval, lw=2, title="Lyapunov = $lyapunov ± $lv")
        display(plot(pannel1, pannel2, layout=2))
    end

    return result, lyapunov, zip(lyapunovs.t, lyapunovs.saveval)
end

""" Solve all trajectories for a given energy at a lattice P,Q = (1...dimension, 1...dimension) """
function SolveEnergy(energy, parameters, dimension; 
        min=-2.0, 
        max=2.0, 
        sectionPlane=2, 
        regularLyapunov=0.01, 
        showFigures=false, 
        verbose=false, 
        randomize=false, 
        kwargs...
    )

    averageLyapunov = zeros(Float64, dimension, dimension)
    countLyapunov = zeros(Int32, dimension, dimension)
    fregSection = zeros(Float64, dimension, dimension)

    maximumLyapunov = 0.0

    # Auxiliary variables just for info about the progress
    trajectories = 0
    infoStep = convert(Int32, round(dimension / 10))
    startTime = time_ns()
    lastTime = time_ns()

    crossings = 0

    println("Starting SolveEnergy: $parameters, E = $energy, dim = $dimension:")
    for ip = 1:dimension, iq = 1:dimension

        if ip % infoStep == 0 && iq == dimension
            numNonzero = 0
            for i in averageLyapunov
                if i != 0
                    numNonzero += 1
                end
            end

            @printf("%.0f%% (%.0fs, %.0fs) trajectories = %d (%.0f%% of the section covered)\n", 100 * ip / dimension, (time_ns() - startTime) / 1E9, (time_ns() - lastTime) / 1E9,trajectories, 100 * numNonzero / (dimension * dimension))

            lastTime = time_ns()
        end

        if countLyapunov[ip, iq] > 0
            continue
        end
        
        # Randomize - IC selected randomly from a calculated cell, otherwise it is always taken from the cell centre
        P = (max - min) * (ip - (randomize ? (1.5 + rand()) : 1.0)) / (dimension - 1.0) + min
        Q = (max - min) * (iq - (randomize ? (1.5 + rand()) : 1.0)) / (dimension - 1.0) + min
    
        x = sectionPlane == 1 ? [0.0, P, missing, Q] : [P, 0.0, Q, missing]
        if InitialCondition!(x, energy, parameters)
            x = collect(skipmissing(x))
            ps, lyapunov = TrajectoryLyapunov(x, parameters; showFigures=showFigures, verbose=verbose, sectionPlane=sectionPlane, kwargs...)
        else
            ps = []
            lyapunov = -0.01                        # Unable to find initial condition - return small negative lyapunov exponent (in order to distinguish later kinematically inaccessible area)
        end
    
        averageLyapunov[ip, iq] = lyapunov
        maximumLyapunov = maximum([maximumLyapunov, lyapunov])
        countLyapunov[ip, iq] = 1

        if lyapunov > 0
            crossings += 1
            trajectories += 1                       # Number of convergent trajectories

            if lyapunov < regularLyapunov
                fregSection[ip, iq] = 1
            end
        else
            fregSection[ip, iq] = -1
        end

        for p in ps
            iip = convert(Int32, round((p[1] + 2) * (dimension - 1) / 4)) + 1
            iiq = convert(Int32, round((p[2] + 2) * (dimension - 1) / 4)) + 1

            # Nonzero Lyapunov exponent is always preferred, so any value replaces cells with nonpositive LE
            if averageLyapunov[iip, iiq] <= 0
                averageLyapunov[iip, iiq] = lyapunov
                countLyapunov[iip, iiq] = 1

                if lyapunov > 0 && lyapunov < regularLyapunov
                    fregSection[iip, iiq] = 1
                end
            else
                averageLyapunov[iip, iiq] += lyapunov
                countLyapunov[iip, iiq] += 1

                if lyapunov > 0 && lyapunov < regularLyapunov
                    fregSection[iip, iiq] += 1
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
    
    @printf("Finished λ = %.3f, E = %.3f, trajectories = %d, crossings = %d, Λ = %.3f, freg = %.3f (%.0fs)", parameters[1], energy, trajectories, crossings, maximumLyapunov, freg, (time_ns() - startTime) / 1E9)

    return averageLyapunov, maximumLyapunov, freg, trajectories
end

""" Calculates Poincaré section for a given energy and number of trajectories """
function PoincareSection(energy, parameters, numTrajectories; 
        min=-2.0, 
        max=2.0, 
        sectionPlane=2, 
        maxICNumber=10000, 
        showFigures=true, 
        kwargs...
    )
    
    figure = scatter()                  # We shall gradually fill the graph with chaotic trajectories

    trajectory = 0
    while trajectory < numTrajectories && maxICNumber > 0
        P = (max - min) * rand() + min
        Q = (max - min) * rand() + min

        maxICNumber -= 1

        x = sectionPlane == 1 ? [0.0, P, missing, Q] : [P, 0.0, Q, missing]
        if !InitialCondition!(x, energy, parameters)
            continue                    # Failure of generating an initial condition
        end
        
        x = collect(skipmissing(x))

        println("IC = $x, E = $(Energy(x, parameters))")

        ps, lyapunov = TrajectoryLyapunov(x, parameters; sectionPlane=sectionPlane, showFigures=false, relativeFluctuationThreshold=0, regularThreshold=0, kwargs...)
        if lyapunov <= 0
            continue                    # Unstable or nonconvergent trajectory
        end

        trajectory += 1

        figure = scatter!(collect(ps), label="Λ=$(round(lyapunov, digits=3))", markersize=3, markeralpha=0.8, markerstrokewidth=0)

        if showFigures
            display(figure)
        end
    end

    return figure
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

        x = sectionPlane == 1 ? [0.0, P, missing, Q] : [P, 0.0, Q, missing]
        if !InitialCondition!(x, energy, parameters)
            continue
        end
        
        x = collect(skipmissing(x))

        ps, lyapunov, lyapunovs = TrajectoryLyapunov(x, parameters; showFigures=false, sectionPlane=sectionPlane, kwargs...)

        if lyapunov <= 0
            continue
        end

        trajectory += 1

        pannel1 = scatter!(pannel1, collect(ps), label="Λ=$(round(lyapunov, digits=2))", markersize=3, markeralpha=0.8, markerstrokewidth=0)
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
