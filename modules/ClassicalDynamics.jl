using DifferentialEquations, DiffEqCallbacks
using LinearAlgebra
using Random
using Statistics
using Plots

""" 
    Rescales the Φ matrix and saves the Lyapunov exponent
    Called periodically (every "saveStep" time units)
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

""" Hit Poincaré section """
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

""" ManifoldProjection method to conserve energy """
function energyConservation!(resid, u, params, t)
    resid[1] = Energy(u, params) - params[5]
    resid[2:end] .= 0
end

""" Solves individual trajectory with given initial conditions x"""
function SolveTrajectory(initialCondition, parameters; solver=TsitPap8(), section=2, savePath=nothing, verbose=false, showFigures=false, relaxationTime=10, regularThreshold=1e-3, relativeFluctuationThreshold=1e-5, maxPSPoints=10000, tolerance=1e-10, saveStep=2, lyapunovSeriesLength=200)
    x0 = zeros(20)
    x0[1:4] = initialCondition
    x0[5:end] = Matrix{Float64}(I, 4, 4)

    params = [parameters..., Energy(x0, parameters), relaxationTime, relativeFluctuationThreshold, regularThreshold, maxPSPoints, 0, 1]   
                                                                                            # λ, δ, ω, ω₀, 5=energy, 6=relaxationTime, 7=relativeFluctuationThreshold, 8=regularThreshold, 9=maximum number of points in Poincare section, 
                                                                                            #              10=num points in the Poincaré section, 11=current Lyapunov exponent 
    append!(params, rand(lyapunovSeriesLength))                                             # Queue with latest Lyapunov exponents (initialized as a random series)

    fnc = ODEFunction(EquationOfMotion!)
    problem = ODEProblem(fnc, x0, (0, 1e5), params)

    sectionCondition = section == 1 ? sectionCondition1 : sectionCondition2

    lyapunovs = SavedValues(Float64, Float64)                                              # For a graph with the time evolution of Lyapunov exponents
    callback = CallbackSet(ManifoldProjection(energyConservation!, save=false), SavingCallback(rescale!, lyapunovs, saveat=saveStep:saveStep:1e6), ContinuousCallback(sectionCondition, section!, nothing, save_positions=(false, true)))
    time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, callback=callback, save_on=true, save_everystep=false, save_start=false, save_end=false, maxiters=1E8, isoutofdomain=CheckDomain, verbose=verbose)

    # Save all unstable or nonconvergent trajectories
    if !(solution.retcode == :Success || solution.retcode == :Terminated || (params[10] >= maxPSPoints && relativeFluctuationThreshold >= 0))
        if !isnothing(savePath)
            open(savePath * "_nonconvergent_trajectories.txt", "a") do io
                println(io, "$(solution.retcode)\t$initialCondition")
            end
        end

        return [], 0    # Unstable trajectories can be distinguished later because their LE is exactly 0
    end

    lyapunov = mean(params[12:end])
    lv = var(params[12:end])

    if verbose
        println("Calculation time = $time, Trajectory time = $(solution.t[end]), PS points = $(params[10]) Λ = $lyapunov ± $lv")

#        m = reshape(solution.u[end][5:end], 4, 4)
#        println(eigvals(m * transpose(m)))

#        ev = 0.5*log.(eigvals(m * transpose(m))) 
#        ev .+= params[9]
#        ev ./= solution.t[end] 
    end

    if showFigures
        pannel1 = scatter(solution[1,:], solution[3,:], title="P = $P, Q = $Q [$(length(solution))]")
        pannel2 = plot(lyapunovs.t, lyapunovs.saveval, lw=2, title="Lyapunov = $lyapunov ± $lv")
        display(plot(pannel1, pannel2, layout=2))
    end

    result = section == 1 ? zip(solution[2,:], solution[4,:]) : zip(solution[1,:], solution[3,:])
    return result, lyapunov, zip(lyapunovs.t, lyapunovs.saveval)
end

""" Solve all trajectories for a given energy at a lattice P,Q = (1...dimension, 1...dimension) """
function SolveEnergy(energy, parameters, dimension; min=-2.0, max=2.0, section=2, showFigures=false, verbose=false, randomize=false, kwargs...)
    sectionLyapunov = zeros(Float64, dimension, dimension)
    countLyapunov = zeros(Int32, dimension, dimension)

    # Auxiliary variables just for info about the progress
    trajectories = 0
    infoStep = convert(Int32, round(dimension / 10))
    startTime = time_ns()

    crossings = 0

    println("Starting $parameters, E=$energy, dim=$dimension")
    for ip = 1:dimension, iq = 1:dimension

        if ip % infoStep == 0 && iq == dimension
            numNonzero = 0
            for i in sectionLyapunov
                if i != 0
                    numNonzero += 1
                end
            end

            println("Trajectories=$trajectories ($(round(100 * numNonzero / (dimension * dimension)))%)")
        end

        if countLyapunov[ip, iq] > 0
            continue
        end
        
        # Randomize - IC selected randomly from a calculated cell, otherwise it is always taken from the cell centre
        P = (max - min) * (ip - (randomize ? (1.5 + rand()) : 1.0)) / (dimension - 1.0) + min
        Q = (max - min) * (iq - (randomize ? (1.5 + rand()) : 1.0)) / (dimension - 1.0) + min
    
        x = section == 1 ? [0.0, P, missing, Q] : [P, 0.0, Q, missing]
        if InitialCondition!(x, energy, parameters)
            x = collect(skipmissing(x))
            ps, lyapunov = SolveTrajectory(x, parameters; showFigures=showFigures, verbose=verbose, section=section, kwargs...)
        else
            ps = []
            lyapunov = -0.01                        # Unable to find initial condition - return small negative lyapunov exponent
        end
    
        sectionLyapunov[ip, iq] = lyapunov
        countLyapunov[ip, iq] = 1

        if lyapunov > 0
            crossings += 1
        end

        if lyapunov > 0
            trajectories += 1                       # Number of convergent trajectories
        end

        for p in ps
            iip = convert(Int32, round((p[1] + 2) * (dimension - 1) / 4)) + 1
            iiq = convert(Int32, round((p[2] + 2) * (dimension - 1) / 4)) + 1

            # Nonzero Lyapunov exponent is always preferred, so any value replaces cells with nonpositive LE
            if sectionLyapunov[iip, iiq] <= 0
                sectionLyapunov[iip, iiq] = lyapunov
                countLyapunov[iip, iiq] = 1
            else
                sectionLyapunov[iip, iiq] += lyapunov
                countLyapunov[iip, iiq] += 1
            end

            crossings += 1
        end
    end

    replace!(countLyapunov, 0=>1)                   # Missing trajectories - avoid division by zero
    sectionLyapunov ./= countLyapunov

    if showFigures
        display(contourf(sectionLyapunov, title="E = $energy, trajectories = $trajectories", zlims=(-0.1, 1)))
    end

    println("Finished $parameters, E=$energy, trajectories=$trajectories, crossings=$crossings (Finished in $(round((time_ns()-startTime)/1e9))s)")

    return sectionLyapunov, trajectories
end

""" Calculates Poincaré section for a given energy and number of trajectories """
function PoincareSection(energy, parameters, numTrajectories; min=-2.0, max=2.0, section=2, maxICNumber=10000, showFigures=true, verbose=false, kwargs...)
    figure = scatter()

    trajectory = 0
    while trajectory < numTrajectories && maxICNumber > 0
        P = (max - min) * rand() + min
        Q = (max - min) * rand() + min

        maxICNumber -= 1

        x = section == 1 ? [0.0, P, missing, Q] : [P, 0.0, Q, missing]
        if !InitialCondition!(x, energy, parameters)
            continue
        end
        
        x = collect(skipmissing(x))

        println("IC = $x, E = $(Energy(x, parameters))")

        ps, lyapunov = SolveTrajectory(x, parameters; section=section, showFigures=false, verbose=verbose, relativeFluctuationThreshold=0, regularThreshold=0, kwargs...)
        trajectory += 1

        figure = scatter!(collect(ps), label="Λ=$(round(lyapunov, digits=3))", markersize=3, markeralpha=0.8, markerstrokewidth=0)

        if showFigures
            display(figure)
        end
    end

    return figure
end

""" Calculates Lyapunov exponents for a given energy and number of trajectories """
function LyapunovExponents(energy, parameters, numTrajectories; showFigures=true, verbose=false, kwargs...)
    pannel1 = scatter()
    pannel2 = plot()
    figure = plot()

    trajectory = 0
    while trajectory < numTrajectories
        P = 4.0 * rand() - 2.0
        Q = 4.0 * rand() - 2.0

        x = [P, 0.0, Q, missing]
        if !InitialCondition!(x, energy, parameters)
            continue
        end
        
        x = collect(skipmissing(x))
        ps, lyapunov, lyapunovs = SolveTrajectory(x, parameters; showFigures=false, verbose=verbose, kwargs...)

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

""" Runs various methods to solve differential equations for the trajectory """ 
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
