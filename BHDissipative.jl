# BHDissipative.jl
#
# Classical (N -> infinity) dynamics of the OPEN Bose-Hubbard model - pump, loss and
# dephasing - and the maximal Lyapunov exponent of its trajectories (Benettin algorithm).
#
# Notation, phase-space layout and normalisation are those of models/BoseHubbardFull.jl and
# BHTrajectory.jl:
#
#   parameters     (L, J, U), extended here to (L, J, U, κ, γ, σ, Δ, f) by DissipativeParameters
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
#   dp_i = [ J (q_{i-1} + q_{i+1}) - U (p_i² + q_i²) q_i + Δ q_i - √2 f - (κ/2) p_i ] dt
#                                                                - √γ q_i ∘ dW_i + σ dV_i
#   dq_i = [-J (p_{i-1} + p_{i+1}) + U (p_i² + q_i²) p_i - Δ p_i        - (κ/2) q_i ] dt
#                                                                + √γ p_i ∘ dW_i + σ dY_i
#
#   κ = γ_loss - γ_pump   net damping
#   γ                     dephasing rate (a random phase kick ψ_i -> e^{-i θ_i} ψ_i)
#   σ = √((γ_loss + γ_pump) / (2N))   finite-N (truncated Wigner) additive noise, expressed in the
#                         scaled variables above; σ = 0 is the N -> infinity limit of
#                         BHTrajectory.jl, and σ does NOT enter the tangent equations.
#   Δ                     detuning, and f the coherent drive - see below.
#
# Coherent driving (bh_dissipation_driving.md). In the frame rotating with the pump the model gains
# a detuning term -Δ Σ_i n_i and a drive Σ_i (F b_i† + F* b_i); in the classical variables used here
# these are the two extra terms above, and the classical Hamiltonian becomes
#
#   H = Σ_i [ -J (p_i p_{i+1} + q_i q_{i+1}) + (U/4) (p_i² + q_i²)² - Δ I_i + √2 f q_i ]
#
# evaluated by DrivenEnergy below (Energy of models/BoseHubbardFull.jl is its f = Δ = 0 part).
#
#   Δ  is the note's Δ̃, the detuning from the BOTTOM OF THE BAND, i.e. from the k = 0 Bloch mode -
#      the only mode a uniform drive couples to, and the quantity every numerical result in the note
#      is quoted in. The bare detuning Δ_bare = Δ_pump - ω_0 that actually enters the equations of
#      motion is Δ_bare = Δ - J z = Δ - 2J for the periodic chain (z = 2); the code does that
#      subtraction internally, so pass Δ exactly as the note writes Δ̃.
#   f  is the drive amplitude per site, F/√N held fixed as N -> infinity (drive power ∝ N, matching
#      the O(N) loss). Real and positive without loss of generality.
#
# In the note's notation ψ_j = (q_j + i p_j)/√2 and g = 2U the equations above are the discrete
# Lugiato-Lefever equation  dψ_j = [iJ(ψ_{j+1}+ψ_{j-1}) + i Δ_bare ψ_j - i g|ψ_j|²ψ_j - i f
# - (κ/2)ψ_j] dt - i√γ ψ_j ∘ dW_j.
#
# What the drive changes:
#   * It breaks U(1), so Δ stops being a gauge parameter and becomes physical.
#   * The norm law is destroyed: d/dt Σ_i I_i = -κ Σ_i I_i - √2 f Σ_i p_i, so Σ I_i no longer decays
#     as e^{-κt}. Instead the flow has a compact ABSORBING BALL Σ_i I_i ≤ 4 f² L / κ², which is what
#     lets a nontrivial (and possibly strange) attractor exist at all.
#   * Trajectories therefore forget their initial condition and settle on an attractor. Measure the
#     exponent ON the attractor: give `relaxationTime` enough time for the transient to die out.
#   * The drive term is a CONSTANT, so it drops out of the linearisation and does not appear in the
#     tangent equations; the detuning is linear and does appear.
#
# With κ = γ = σ = Δ = f = 0 everything below reduces to the Hamiltonian problem of BHTrajectory.jl.
#
# Methods
#   :split  (default) Strang splitting. The deterministic part - hopping, interaction, damping and
#           the tangent dynamics - is an ODEProblem solved to `tolerance`; the dephasing is applied
#           EXACTLY, as a random rotation of every (p_i, q_i) and (δp_i, δq_i) plane, by a
#           PeriodicCallback every `noiseStep`. In the undriven model (f = 0) the norm law
#           Σ I_i = e^{-κt} then holds to solver tolerance.
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
    three dissipative channels and by the coherent drive. The result destructures as
    `L, J, U = parameters`, so every function of the closed model (Energy, InitialCondition,
    EquationOfMotionTangentVector!, ...) accepts it unchanged, while the rest stays reachable by name.

    κ = γ_loss - γ_pump (net damping), γ = dephasing rate, σ = additive (finite-N) noise,
    Δ = detuning from the bottom of the band (the note's Δ̃), f = drive amplitude F/√N. """
function DissipativeParameters(parameters; κ = 0.0, γ = 0.0, σ = 0.0, Δ = 0.0, f = 0.0)
    L, J, U = parameters
    return (L = L, J = float(J), U = float(U), κ = float(κ), γ = float(γ), σ = float(σ),
            Δ = float(Δ), f = float(f))
end


""" Classical Hamiltonian of the DRIVEN model,

        H = Σ_i [ -J (p_i p_{i+1} + q_i q_{i+1}) + (U/4)(p_i² + q_i²)² - Δ I_i + √2 f q_i ],

    i.e. Energy(x, parameters) of models/BoseHubbardFull.jl plus the detuning and the drive. It is
    conserved in the conservative limit κ = γ = σ = 0 whatever Δ and f are, which is what the check
    in Demo() uses; for Δ = f = 0 it coincides with Energy. """
function DrivenEnergy(x, parameters)
    L, J, U, κ, γ, σ, Δ, f = parameters

    energy = Energy(x, parameters)

    if Δ != 0 || f != 0
        bareDetuning = Δ - 2 * J             # Δ is the note's Δ̃, measured from the band bottom
        driveTerm = sqrt(2) * f
        @inbounds for i = 1:L
            energy += -bareDetuning * 0.5 * (x[i]^2 + x[i + L]^2) + driveTerm * x[i + L]
        end
    end

    return energy
end


""" Equations of motion of the open Bose-Hubbard model together with a single deviation vector, in
    the layout of EquationOfMotionTangentVector! (models/BoseHubbardFull.jl).

    The conservative part - the trajectory and the matrix-free tangent dynamics - is taken over
    from the closed model; added here are the damping -κ/2, the detuning Δ and the coherent drive f.
    The dephasing and the additive noise are NOT here: they are applied either exactly by
    DephasingCallback (method = :split) or through NoiseTerm! (method = :sde).

    Damping and detuning are linear, so the linearisation reproduces them verbatim on the deviation
    vector. The drive is a CONSTANT and therefore drops out of the tangent equations entirely - it
    shifts the attractor, not the stability around it. """
function EquationOfMotionDissipative!(dx, x, parameters, t)
    EquationOfMotionTangentVector!(dx, x, parameters, t)

    L, J, U, κ, γ, σ, Δ, f = parameters.modelParameters
    dimension = 2 * L

    if κ != 0
        halfκ = 0.5 * κ
        @inbounds @simd for i in eachindex(dx)
            dx[i] -= halfκ * x[i]
        end
    end

    # Δ ≠ 0 or f ≠ 0: the detuning is physical and has to be applied - note that the resonant driven
    # case is Δ = Δ̃ = 0, which still carries the bare detuning -2J. For Δ = f = 0 (the undriven
    # default) the term is skipped altogether, which is exact: without a drive Δ is pure gauge.
    if Δ != 0 || f != 0
        # Δ is the note's Δ̃, measured from the bottom of the band; the equations of motion take the
        # bare detuning, the hopping already supplying the J z = 2J offset for a uniform state.
        bareDetuning = Δ - 2 * J
        @inbounds for i = 1:L
            dx[i] += bareDetuning * x[i + L]                            # dp_i += Δ_bare q_i
            dx[i + L] -= bareDetuning * x[i]                            # dq_i -= Δ_bare p_i
            dx[dimension + i] += bareDetuning * x[dimension + i + L]    # the same on the deviation vector
            dx[dimension + i + L] -= bareDetuning * x[dimension + i]
        end
    end

    if f != 0
        driveTerm = sqrt(2) * f                                 # dp_i -= √2 f, constant: no tangent part
        @inbounds @simd for i = 1:L
            dx[i] -= driveTerm
        end
    end

    return nothing
end


""" Exact dephasing step for method = :split - a PeriodicCallback firing every `noiseStep`.

    Dephasing rotates each (p_i, q_i) plane by a random angle √γ ΔW_i; the SAME rotation is applied
    to (δp_i, δq_i), so the trajectory and the deviation vector always see one and the same noise
    realisation. The rotation is norm-preserving, so without a drive Σ I_i keeps decaying as e^{-κt}
    to solver tolerance. The additive noise σ acts on the trajectory only. """
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
    return (DrivenEnergy(x, integrator.p.modelParameters), 0.5 * sum(abs2, @view x[1:dimension]))
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
    relaxationTime  the exponent is accumulated only after this time (0 gives Λ = Σ log(growth) / t).
                    With a drive (f ≠ 0) the trajectory forgets its initial condition and settles on
                    an attractor, so make this long enough for the transient to die out - otherwise
                    the transient, not the attractor, is what gets measured.
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

    L, J, U, κ, γ, σ, Δ, f = parameters
    phaseSpaceDimension = length(initialCondition)          # = 2L
    rng = Xoshiro(seed)

    # Initial condition + a single normalised deviation vector (Benettin method), i.e. exactly the
    # layout expected by EquationOfMotionTangentVector!
    x0 = zeros(2 * phaseSpaceDimension)
    x0[1:phaseSpaceDimension] = initialCondition
    deviation = @view x0[(phaseSpaceDimension + 1):end]
    deviation .= randn(rng, phaseSpaceDimension)
    deviation ./= sqrt(sum(abs2, deviation))

    energy = DrivenEnergy(x0, parameters)
    norm0 = 0.5 * sum(abs2, @view x0[1:phaseSpaceDimension])

    noisy = γ > 0 || σ > 0
    if noisy
        # Without the drive Σ I_i never grows, so max(I) at t = 0 bounds it for all later times. The
        # drive pumps the ring instead, up to the absorbing ball Σ_i I_i ≤ 4 f² L / κ² of the note,
        # which is then the bound to respect. (For f > 0 and κ = 0 nothing bounds the norm at all.)
        maximumI = maximum(i -> 0.5 * (x0[i]^2 + x0[i + L]^2), 1:L)
        if f > 0
            maximumI = κ > 0 ? max(maximumI, 4 * f^2 * L / κ^2) : Inf
        end

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
        @info "Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(DrivenEnergy(finalState, parameters)), Final norm = $finalNorm, Λ = $lyapunov ± $lv, Λ + κ/2 = $(lyapunov + 0.5 * κ)"
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
    # Valid for f = 0 only: the drive destroys the norm law, see the header.
    normError(observables, κ) = maximum(abs.(last.(observables.saveval) ./ exp.(-κ .* observables.t) .- 1))

    # Largest Bogoliubov exponent of the uniform state, note sec. 5:
    #   λ_k = -κ/2 ± √(g²n² - A_k²),  A_k = ε̄_k - Δ + 2 g n,  ε̄_k = 2J(1 - cos k),  g = 2U
    function BogoliubovLyapunov(L, J, g, n, Δ, κ)
        largest = -Inf
        for m = 0:(L - 1)
            A = 2 * J * (1 - cos(2π * m / L)) - Δ + 2 * g * n
            radicand = g^2 * n^2 - A^2
            largest = max(largest, radicand > 0 ? -0.5κ + sqrt(radicand) : -0.5κ)
        end
        return largest
    end

    # Uniform fixed point of the driven model, note sec. 4: ψ = f / [(Δ - g n) + i κ/2] with
    # n [(Δ - g n)² + κ²/4] = f². Returned in the (p, q) layout used here.
    function UniformFixedPoint(L, g, n, Δ, κ)
        f = sqrt(n * ((Δ - g * n)^2 + κ^2 / 4))
        ψ = f / ((Δ - g * n) + im * κ / 2)
        return vcat(fill(sqrt(2) * imag(ψ), L), fill(sqrt(2) * real(ψ), L)), f
    end

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

    # ---------------------------------------------------------------- coherent driving
    println()

    # 6. Conservative driven limit - with κ = gamma = sigma = 0 the driven Hamiltonian
    #    DrivenEnergy is still exactly conserved, which tests detuning and drive together.
    parameters = DissipativeParameters(bhParameters; Δ=0.7, f=0.3)
    _, lyapunov, _, observables = TrajectoryLyapunovDissipative(initialCondition, parameters; timeInterval=timeInterval)
    @printf("driven, no dissipation:  Λ = %+.4f   rel. driven-energy drift = %.1e\n",
        lyapunov, energyDrift(observables))

    # 7. The uniform fixed point of note sec. 4, at the note's own parameters. It must be stationary
    #    to machine precision, and the exponent measured on it must reproduce the analytic
    #    Bogoliubov rate of note sec. 5.
    Lsites, Jhop, g, κ, Δ = 6, 1.0, -1.0, 1.0, -4.0
    println("uniform fixed point, L = 6, J = 1, g = -1, κ = 1, Δ = -4 (note sec. 4, 5, 7):")
    for n in (3.9684, 3.5, 2.5, 1.3650)              # 3.9684 and 1.3650 are the two folds
        x, f = UniformFixedPoint(Lsites, g, n, Δ, κ)
        parameters = DissipativeParameters((Lsites, Jhop, g / 2); κ=κ, Δ=Δ, f=f)

        _, lyapunov, _, observables = TrajectoryLyapunovDissipative(x, parameters;
            timeInterval=(0.0, 600.0), saveStep=1, relaxationTime=100, historyLyapunovExponentLength=200)

        @printf("   n = %6.4f, f = %6.4f:  Λ = %+.5f   Bogoliubov %+.5f   Σ I: %.4f -> %.4f\n",
            n, f, lyapunov, BogoliubovLyapunov(Lsites, Jhop, g, n, Δ, κ),
            n * Lsites, last(last.(observables.saveval)))
    end
    println("   (n = 3.0 is omitted: k = 0 is unstable there - the middle branch is a saddle - so the")
    println("    state slides along the uniform manifold to the upper branch before Lambda is measured.)")

    # 8. The chaotic attractor of note sec. 7 (L = 8, f = 3.2, λ_max = +0.74 there), and the
    #    absorbing ball Σ_i I_i <= 4 f² L / κ² that replaces the norm law once the drive is on.
    Lsites = 8
    f = 3.2
    parameters = DissipativeParameters((Lsites, Jhop, g / 2); κ=κ, Δ=Δ, f=f)
    Random.seed!(6)
    x = 0.5 .* randn(2 * Lsites)
    _, lyapunov, _, observables = TrajectoryLyapunovDissipative(x, parameters;
        timeInterval=(0.0, 4000.0), relaxationTime=1000)
    @printf("chaotic attractor, L = %d, f = %.1f:  Λ = %+.4f  (note: +0.74)   max Σ I = %.3f <= 4f^2 L/κ^2 = %.3f\n",
        Lsites, f, lyapunov, maximum(last.(observables.saveval)), 4 * f^2 * Lsites / κ^2)
    println("   (the note reports strong multistability here: other initial conditions land on other")
    println("    attractors, with Λ from negative - a fixed point - up to about +0.74.)")

end

# Runs only when the file is executed directly (julia BHDissipative.jl), so that
# include("BHDissipative.jl") from BHMapDissipative.jl and friends stays silent.
if abspath(PROGRAM_FILE) == @__FILE__
    Demo()
end
