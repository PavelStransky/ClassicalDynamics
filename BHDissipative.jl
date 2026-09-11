# BHDissipative.jl
#
# Classical (N -> infinity) dynamics of the OPEN Bose-Hubbard model - pump, loss and
# dephasing - and the maximal Lyapunov exponent of its trajectories (Benettin algorithm).
#
# Notation, phase-space layout and normalisation are those of models/BoseHubbardFull.jl and
# BHTrajectory.jl:
#
#   parameters     (L, J, U), extended here to (L, J, U, κ, γ, σ) by DissipativeParameters
#   state          x = (p, q),  x[i] = p_i,  x[i + L] = q_i,  i = 1...L  (periodic chain)
#   deviation      x[2L + i] = δp_i,  x[2L + i + L] = δq_i   (the tangentDynamics = :vector layout)
#   normalisation  Σ_i (p_i² + q_i²) = 2,  i.e.  Σ_i I_i = 1  with  I_i = (p_i² + q_i²) / 2
#   amplitudes     ψ_i = (q_i + i p_i) / √2
#
# The Hamiltonian is exactly the one evaluated by Energy(x, parameters):
#
#   H = Σ_i [ -J (p_i p_{i+1} + q_i q_{i+1}) + (U/4) (p_i² + q_i²)² ]
#
# and the open-system equations of motion (Stratonovich noise, ∘ dW) read
#
#   dp_i = [ J (q_{i-1} + q_{i+1}) - U (p_i² + q_i²) q_i - (κ/2) p_i ] dt - √γ q_i ∘ dW_i + σ dV_i
#   dq_i = [-J (p_{i-1} + p_{i+1}) + U (p_i² + q_i²) p_i - (κ/2) q_i ] dt + √γ p_i ∘ dW_i + σ dY_i
#
#   κ = γ_loss - γ_pump   net damping,  Σ_i I_i(t) = e^{-κ t} Σ_i I_i(0)
#   γ                     dephasing rate (a random phase kick ψ_i -> e^{-i θ_i} ψ_i)
#   σ = √((γ_loss + γ_pump) / (2N))   finite-N (truncated Wigner) additive noise, expressed in the
#                         scaled variables above; σ = 0 is the N -> infinity limit of
#                         BHTrajectory.jl, and σ does NOT enter the tangent equations.
#
# With κ = γ = σ = 0 everything below reduces to the Hamiltonian problem of BHTrajectory.jl.
#
# Methods
#   :split  (default) Strang splitting. The deterministic part - hopping, interaction, damping and
#           the tangent dynamics - is an ODEProblem solved to `tolerance`; the dephasing is applied
#           EXACTLY, as a random rotation of every (p_i, q_i) and (δp_i, δq_i) plane, by a
#           PeriodicCallback every `noiseStep`. The norm law Σ I_i = e^{-κt} then holds to solver
#           tolerance.
#   :sde    SDEProblem with non-diagonal noise, solved with EulerHeun (Stratonovich). Simpler, but
#           the norm drifts at O(noiseStep). Useful as a cross-check.
#           Do NOT use Ito solvers (EM, SRIW1, SOSRI, ...) on these Stratonovich equations.

using DifferentialEquations
using StochasticDiffEq
using LinearAlgebra
using Random
using Statistics
using Printf

include("models/BoseHubbardFull.jl")
include("modules/ClassicalDynamics.jl")


""" Bose-Hubbard parameters (L, J, U) of models/BoseHubbardFull.jl extended by the rates of the
    three dissipative channels. The result destructures as `L, J, U = parameters`, so every
    function of the closed model (Energy, InitialCondition, EquationOfMotionTangentVector!, ...)
    accepts it unchanged, while κ, γ and σ stay reachable by name.

    κ = γ_loss - γ_pump (net damping), γ = dephasing rate, σ = additive (finite-N) noise. """
function DissipativeParameters(parameters; κ = 0.0, γ = 0.0, σ = 0.0)
    L, J, U = parameters
    return (L = L, J = float(J), U = float(U), κ = float(κ), γ = float(γ), σ = float(σ))
end


""" Equations of motion of the open Bose-Hubbard model together with a single deviation vector, in
    the layout of EquationOfMotionTangentVector! (models/BoseHubbardFull.jl).

    The conservative part - the trajectory and the matrix-free tangent dynamics - is taken over
    from the closed model; the only addition is the damping -κ/2, which the linearisation
    reproduces verbatim on the deviation vector. The dephasing and the additive noise are NOT here:
    they are applied either exactly by DephasingCallback (method = :split) or through NoiseTerm!
    (method = :sde). """
function EquationOfMotionDissipative!(dx, x, parameters, t)
    EquationOfMotionTangentVector!(dx, x, parameters, t)

    κ = parameters.modelParameters.κ
    if κ != 0
        halfκ = 0.5 * κ
        @inbounds @simd for i in eachindex(dx)
            dx[i] -= halfκ * x[i]
        end
    end

    return nothing
end


""" Exact dephasing step for method = :split - a PeriodicCallback firing every `noiseStep`.

    Dephasing rotates each (p_i, q_i) plane by a random angle √γ ΔW_i; the SAME rotation is applied
    to (δp_i, δq_i), so the trajectory and the deviation vector always see one and the same noise
    realisation. The rotation is norm-preserving, hence Σ I_i keeps decaying as e^{-κt} to solver
    tolerance. The additive noise σ acts on the trajectory only. """
function DephasingCallback(parameters, noiseStep, rng)
    L = parameters.L
    dimension = 2 * L
    dephasingAmplitude = sqrt(parameters.γ * noiseStep)
    additiveAmplitude = parameters.σ * sqrt(noiseStep)

    function Dephase!(integrator)
        x = integrator.u

        @inbounds for i = 1:L
            s, c = sincos(dephasingAmplitude * randn(rng))

            p, q = x[i], x[i + L]
            x[i] = c * p - s * q
            x[i + L] = s * p + c * q

            δp, δq = x[dimension + i], x[dimension + i + L]
            x[dimension + i] = c * δp - s * δq
            x[dimension + i + L] = s * δp + c * δq

            if additiveAmplitude > 0
                x[i] += additiveAmplitude * randn(rng)
                x[i + L] += additiveAmplitude * randn(rng)
            end
        end

        u_modified!(integrator, true)
    end

    return PeriodicCallback(Dephase!, float(noiseStep); save_positions=(false, false))
end


""" Noise term (the matrix G of dx = f dt + G ∘ dW) for method = :sde.
    Column i is the dephasing on site i, columns L + i and 2L + i the additive noise on p_i and q_i
    (present only when σ > 0). Entries that are never written stay zero. """
function NoiseTerm!(G, x, parameters, t)
    L, J, U, κ, γ, σ = parameters.modelParameters
    dimension = 2 * L
    dephasingAmplitude = sqrt(γ)

    @inbounds for i = 1:L
        G[i, i] = -dephasingAmplitude * x[i + L]
        G[i + L, i] = dephasingAmplitude * x[i]
        G[dimension + i, i] = -dephasingAmplitude * x[dimension + i + L]
        G[dimension + i + L, i] = dephasingAmplitude * x[dimension + i]

        if σ > 0
            G[i, L + i] = σ
            G[i + L, 2 * L + i] = σ
        end
    end

    return nothing
end


""" Non-mutating reader for the SavingCallback that records the energy and the total number of
    bosons Σ I_i - neither of which is conserved any more. """
function EnergyNorm(x, t, integrator)
    dimension = integrator.p.dimension
    return (Energy(x, integrator.p.modelParameters), 0.5 * sum(abs2, @view x[1:dimension]))
end


""" Calculates the largest Lyapunov exponent of an individual trajectory of the OPEN Bose-Hubbard
    model. The dissipative counterpart of TrajectoryLyapunov(...; tangentDynamics = :vector) from
    modules/ClassicalDynamics.jl, with which it shares the phase-space layout, the integration
    parameters and all the callbacks.

    Because neither the energy nor the norm is conserved, no ManifoldProjection and no Poincaré
    section are used, and the convergence test is switched off by default (regularThreshold =
    relativeFluctuationThreshold = 0), so the trajectory is always integrated over the whole
    `timeInterval`.

    method          :split (Strang splitting, exact dephasing - default) or :sde (EulerHeun)
    noiseStep       splitting interval (:split) or the fixed EulerHeun step (:sde)
    saveStep        how often the deviation vector is renormalised = resolution of the output
    relaxationTime  the exponent is accumulated only after this time (0 gives Λ = Σ log(growth) / t)
    seed            seed of the noise; identical seeds give identical realisations

    Returns - the solution, the Lyapunov exponent Λ [the damping contributes exactly -κ/2 to it, so
              Λ + κ/2 is the intrinsic exponent], the SavedValues of the running exponent and the
              SavedValues of (energy, norm). """
function TrajectoryLyapunovDissipative(initialCondition, parameters;
        method = :split,
        solver = DP8(),                     # deterministic solver of the :split method
        stochasticSolver = EulerHeun(),     # Stratonovich solver of the :sde method
        noiseStep = 1e-2,
        seed = 1,
        timeInterval = (0.0, 1e3),
        saveStep = 2,
        relaxationTime = 0,
        tolerance = 1e-10,
        maximumIterations = 2E6,
        timeout = 0,
        regularThreshold = 0,
        relativeFluctuationThreshold = 0,
        historyLyapunovExponentLength = 500,
        showFigures = false,
        savePath = nothing
    )

    L, J, U, κ, γ, σ = parameters
    phaseSpaceDimension = length(initialCondition)          # = 2L
    rng = Xoshiro(seed)

    # Initial condition + a single normalised deviation vector (Benettin method), i.e. exactly the
    # layout expected by EquationOfMotionTangentVector!
    x0 = zeros(2 * phaseSpaceDimension)
    x0[1:phaseSpaceDimension] = initialCondition
    deviation = @view x0[(phaseSpaceDimension + 1):end]
    deviation .= randn(rng, phaseSpaceDimension)
    deviation ./= sqrt(sum(abs2, deviation))

    energy = Energy(x0, parameters)
    norm0 = 0.5 * sum(abs2, @view x0[1:phaseSpaceDimension])

    noisy = γ > 0 || σ > 0
    if noisy
        maximumI = maximum(i -> 0.5 * (x0[i]^2 + x0[i + L]^2), 1:L)
        phaseAdvance = 2 * U * maximumI * noiseStep         # nonlinear phase advance per noise step
        if phaseAdvance > 0.2
            @warn "noiseStep = $noiseStep is too coarse: the nonlinear phase advance per step 2 U max(I) noiseStep = $phaseAdvance exceeds 0.2"
        end
    end

    # maximumSectionPoints = 0 -> no Poincaré section, the stopping decision is taken in AccumulateLyapunov!
    integrationParameters = LyapunovIntegrationParameters(phaseSpaceDimension, parameters, energy, relaxationTime,
        relativeFluctuationThreshold, regularThreshold, 0, 0, time_ns(), 1E9 * timeout, :Start,
        0.0, norm0, historyLyapunovExponentLength, Float64[])

    lyapunovs = SavedValues(Float64, Float64)               # The whole history of the immediate Lyapunov exponents (for a graph)
    observables = SavedValues(Float64, Tuple{Float64, Float64})     # Energy and norm, neither of them conserved

    rescale = PeriodicCallback(RescaleDeviationVector!, float(saveStep); save_positions=(false, false))
    record = SavingCallback(RunningLyapunov, lyapunovs, saveat=saveStep:saveStep:last(timeInterval))
    recordObservables = SavingCallback(EnergyNorm, observables, saveat=saveStep:saveStep:last(timeInterval))
    timeoutCallback = DiscreteCallback(TimeoutCondition, terminate!)

    if method === :split || !noisy
        callback = noisy ?
            CallbackSet(DephasingCallback(parameters, noiseStep, rng), rescale, record, recordObservables, timeoutCallback) :
            CallbackSet(rescale, record, recordObservables, timeoutCallback)

        problem = ODEProblem(ODEFunction(EquationOfMotionDissipative!), x0, timeInterval, integrationParameters)
        time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, callback=callback,
            save_on=true, save_everystep=false, save_start=true, save_end=true, maxiters=maximumIterations,
            isoutofdomain=CheckDomain, verbose=DEVerbosity())
    elseif method === :sde
        numberOfNoiseProcesses = σ > 0 ? 3 * L : L
        callback = CallbackSet(rescale, record, recordObservables, timeoutCallback)

        problem = SDEProblem(EquationOfMotionDissipative!, NoiseTerm!, x0, float.(timeInterval), integrationParameters;
            noise_rate_prototype = zeros(2 * phaseSpaceDimension, numberOfNoiseProcesses))
        time = @elapsed solution = solve(problem, stochasticSolver, dt=noiseStep, adaptive=false, seed=UInt64(seed),
            callback=callback, save_everystep=false, save_start=true, save_end=true, maxiters=maximumIterations)
    else
        error("method must be :split or :sde")
    end

    # Get results - average over the trailing window of running-exponent samples (empty -> exactly 0)
    history = integrationParameters.historyLyapunovExponent
    window = @view history[max(1, length(history) - historyLyapunovExponentLength + 1):end]
    lyapunov = isempty(window) ? 0.0 : mean(window)
    lv = length(window) > 1 ? var(window) : 0.0

    # Print result
    if length(solution.u) > 0
        finalState = solution.u[end]
        finalNorm = 0.5 * sum(abs2, @view finalState[1:phaseSpaceDimension])
        @info "Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(Energy(finalState, parameters)), Final norm = $finalNorm, Λ = $lyapunov ± $lv, Λ + κ/2 = $(lyapunov + 0.5 * κ)"
    end

    @debug "retcode = $(solution.retcode), result = $(integrationParameters.result)"

    # Save all unstable or nonconvergent trajectories (for debug reasons)
    if !(solution.retcode == DiffEqBase.ReturnCode.Success ||
         (solution.retcode == DiffEqBase.ReturnCode.Terminated && integrationParameters.result == :Converged))
        if !isnothing(savePath)
            open(savePath * "Nonconvergent_Trajectories.txt", "a") do io
                println(io, "$parameters\t$energy\t$initialCondition\t$(solution.retcode)")
            end
        end

        @info "Nonconvergent trajectory with initialCondition = $initialCondition: retcode = $(solution.retcode), result = $(integrationParameters.result)"
        return solution, 0.0, lyapunovs, observables
    end

    if showFigures
        panel1 = plot(lyapunovs.t, lyapunovs.saveval, lw=2, title="Λ = $lyapunov ± $lv", label=nothing, xlabel="t", ylabel="Λ")
        panel2 = plot(observables.t, last.(observables.saveval), lw=2, label="Σ I", xlabel="t")
        panel2 = plot!(panel2, observables.t, first.(observables.saveval), lw=2, label="E")
        figure = plot(panel1, panel2, layout=2)
        display(figure)

        if !isnothing(savePath)
            savefig(figure, savePath * "$initialCondition.png")
        end
    end

    return solution, lyapunov, lyapunovs, observables
end


""" Finite-time Lyapunov exponent on the window [t1, t2], reconstructed from the running exponents
    saved by TrajectoryLyapunovDissipative. For κ ≠ 0 the density decays, the trajectory becomes
    less and less nonlinear and the exponent drifts, so windowed exponents are more informative
    than the global average. """
function WindowLyapunov(lyapunovs, t1, t2; relaxationTime = 0)
    i1 = searchsortedfirst(lyapunovs.t, t1)
    i2 = searchsortedlast(lyapunovs.t, t2)

    if i1 >= i2 || i2 > length(lyapunovs.t)
        return 0.0
    end

    # saveval = Σ log(growth) / (t - relaxationTime), so the cumulative stretching is recovered back
    stretching1 = lyapunovs.saveval[i1] * (lyapunovs.t[i1] - relaxationTime)
    stretching2 = lyapunovs.saveval[i2] * (lyapunovs.t[i2] - relaxationTime)

    return (stretching2 - stretching1) / (lyapunovs.t[i2] - lyapunovs.t[i1])
end


""" Mean and standard error of the Lyapunov exponent over independent noise realisations (threaded).
    `initialCondition` is either a fixed phase-space point or a function seed -> point. """
function EnsembleLyapunov(initialCondition, parameters; seeds = 1:8, kwargs...)
    lyapunovs = zeros(Float64, length(seeds))

    Threads.@threads for i in eachindex(seeds)
        x = initialCondition isa Function ? initialCondition(seeds[i]) : initialCondition
        lyapunovs[i] = TrajectoryLyapunovDissipative(x, parameters; seed=seeds[i], kwargs...)[2]
    end

    return mean(lyapunovs), std(lyapunovs) / sqrt(length(lyapunovs))
end


""" Demonstration of all the dissipative channels, together with the consistency checks of the
    integration (energy drift, norm law, the exactly known exponent of the linear model). """
function Demo()
    Random.seed!(1234)

    bhParameters = (5, 0.5, 1)                  # (L, J, U), the model of BHTrajectory.jl
    energy = 0.4
    timeInterval = (0.0, 1000.0)

    initialCondition = InitialCondition(energy, bhParameters, 1e-3)

    if initialCondition === nothing
        println("No initial condition found")
        return
    end

    @printf("L = %d, J = %.2f, U = %.2f:  E = %+.4f, Σ I = %.4f\n", bhParameters[1], bhParameters[2], bhParameters[3],
        Energy(initialCondition, bhParameters), 0.5 * sum(abs2, initialCondition))

    energies(observables) = first.(observables.saveval)
    energyDrift(observables) = maximum(abs.(energies(observables) .- energies(observables)[1])) / abs(energies(observables)[1])
    normError(observables, κ) = maximum(abs.(last.(observables.saveval) ./ exp.(-κ .* observables.t) .- 1))

    # 1. Hamiltonian limit: κ = γ = σ = 0 reproduces BHTrajectory.jl (without the manifold projection)
    parameters = DissipativeParameters(bhParameters)
    _, lyapunov, _, observables = TrajectoryLyapunovDissipative(initialCondition, parameters; timeInterval=timeInterval)
    @printf("Hamiltonian:          Λ = %+.4f   rel. energy drift = %.1e   norm err = %.1e\n",
        lyapunov, energyDrift(observables), normError(observables, 0.0))

    # 2. Dephasing only - the norm is preserved and the two methods have to agree
    for method in (:split, :sde)
        parameters = DissipativeParameters(bhParameters; γ=0.1)
        _, lyapunov, _, observables = TrajectoryLyapunovDissipative(initialCondition, parameters;
            method=method, noiseStep=5e-3, timeInterval=timeInterval)
        @printf("γ = 0.1, %-6s:       Λ = %+.4f   norm err = %.1e\n", string(method), lyapunov, normError(observables, 0.0))
    end

    # 3. U = 0 - the dynamics is linear, the damping alone gives Λ = -κ/2 exactly
    parameters = DissipativeParameters((bhParameters[1], bhParameters[2], 0.0); κ=0.3, γ=0.5)
    _, lyapunov, _, _ = TrajectoryLyapunovDissipative(initialCondition, parameters; timeInterval=(0.0, 50.0))
    @printf("U = 0, κ = 0.3:       Λ = %+.6f   (exact: -κ/2 = -0.15)\n", lyapunov)

    # 4. Net damping - the density decays, so windowed exponents are used instead of the average
    κ = 0.02
    parameters = DissipativeParameters(bhParameters; κ=κ)
    _, _, lyapunovs, observables = TrajectoryLyapunovDissipative(initialCondition, parameters; timeInterval=(0.0, 1000.0))
    println("κ = $κ, finite-time exponents on time windows [t1, t2]:")
    for (t1, t2) in ((2, 50), (100, 150), (250, 300), (500, 1000))
        windowLyapunov = WindowLyapunov(lyapunovs, t1, t2)
        norms = last.(observables.saveval)[(observables.t .>= t1) .& (observables.t .<= t2)]
        @printf("   [%4d,%4d]:  Λ = %+.3f   Λ + κ/2 = %+.3f   ⟨Σ I⟩ = %.3f\n",
            t1, t2, windowLyapunov, windowLyapunov + κ / 2, mean(norms))
    end

    # 5. Average over independent noise realisations
    parameters = DissipativeParameters(bhParameters; γ=0.1)
    meanLyapunov, errorLyapunov = EnsembleLyapunov(initialCondition, parameters; seeds=1:4, timeInterval=(0.0, 300.0))
    @printf("ensemble, γ = 0.1:    Λ = %.3f ± %.3f\n", meanLyapunov, errorLyapunov)
end

# Runs only when the file is executed directly (julia BHDissipative.jl), so that
# include("BHDissipative.jl") from BHMapDissipative.jl and friends stays silent.
if abspath(PROGRAM_FILE) == @__FILE__
    Demo()
end
