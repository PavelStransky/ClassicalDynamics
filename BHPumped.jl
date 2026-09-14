# BHPumped.jl
#
# Classical (N -> infinity) dynamics of the INCOHERENTLY PUMPED Bose-Hubbard model with TWO-BODY
# LOSS - the discrete complex Ginzburg-Landau equation of pump_bh_cgle.py - and the Lyapunov
# SPECTRUM of its trajectories (Benettin algorithm with Gram-Schmidt reorthonormalisation).
#
# This is the UNDRIVEN counterpart of BHDissipative.jl. There the ring is filled by a coherent
# drive and emptied by a linear loss; here it is filled by an INCOHERENT pump and emptied by a
# TWO-BODY loss. The difference is not cosmetic:
#
#   * U(1) is UNBROKEN. The drive was the only term that fixed the global phase, so without it the
#     flow is equivariant under ψ_j -> e^{iθ} ψ_j. Every attractor therefore carries an exact ZERO
#     Lyapunov exponent (the phase mode), and Δ becomes a pure gauge parameter - see below. In
#     particular λ_max = 0 on a fixed point or a limit cycle here, not -κ/2 as in the driven model,
#     which is what the chaos threshold of BHMapPumped.jl has to separate chaos from.
#   * The gain is LINEAR and the loss NONLINEAR, the opposite of the driven model. The norm is not
#     bounded by an absorbing ball fed by the drive but by the saturation of the pump itself: the
#     origin repels at rate P/2 and the two-body loss stops the growth at |ψ|² = P/Γ.
#
# Lindblad model (the header of pump_bh_cgle.py):
#
#   H  = -J Σ_<ij> (b_i^† b_j + h.c.) - μ Σ_i n_i + (U_q/2) Σ_i n_i (n_i - 1)
#   L1 = √γ_p b_i^†  (incoherent pump)      L2 = √γ_l b_i          (one-body loss)
#   L3 = √Γ_2 b_i²   (two-body loss)        L4 = √Γ_b (b_i - b_j)  (bond-correlated loss)
#
# In the classical limit ψ_i = b_i/√N with U_q = g/N, Γ_2 = Γ/(2N) and γ_p, γ_l, Γ_b of order 1:
#
#   dψ_j/dt = (d + iJ) (ψ_{j+1} + ψ_{j-1} - 2ψ_j) + i Δ ψ_j + (P/2) ψ_j - (i g + Γ/2) |ψ_j|² ψ_j
#
#   P = γ_p - γ_l    net LINEAR GAIN. It is exactly the -κ of BHDissipative.jl (P = -κ): a pump
#                    that beats the one-body loss is a NEGATIVE damping. P > 0 is the interesting
#                    side, and it is what makes the vacuum unstable.
#   Γ = 2 N Γ_2      two-body loss - the nonlinear saturation, and the only thing that bounds |ψ|.
#   d = Γ_b / 2      real diffusion. It multiplies the discrete Laplacian, so unlike J it is not a
#                    Hamiltonian term: it damps every nonuniform mode, and its -2d on the diagonal
#                    is a genuine loss, not the gauge shift that the -2iJ of the same Laplacian is.
#   g = 2 U          the interaction of models/BoseHubbardFull.jl written as in the Python code and
#                    in bh_dissipation_driving.md (Kerr nonlinearity, the imaginary partner of Γ/2).
#   Δ                detuning, the note's Δ̃ (measured from the BOTTOM OF THE BAND, i.e. from the
#                    k = 0 Bloch mode), exactly as in BHDissipative.jl - and exactly the Dt of
#                    pump_bh_cgle.py, whose Laplacian absorbs the same 2J. Without a drive it is
#                    PURE GAUGE: ψ -> e^{iΔt} ψ removes it and no Lyapunov exponent depends on it.
#                    It is kept because it is what turns a rotating plane wave into a genuine FIXED
#                    POINT, which is how PlaneWaveGrowth is tested against the integrator in Demo().
#
# Together this is the discrete CGLE, i.e. a ring of Stuart-Landau (Hopf) oscillators coupled both
# conservatively (J) and diffusively (d).
#
# Phase-space layout, normalisation and amplitudes are those of models/BoseHubbardFull.jl,
# BHTrajectory.jl and BHDissipative.jl:
#
#   parameters     (L, J, U), extended here to (L, J, U, P, Γ, d, Δ) by PumpedParameters
#   state          x = (p, q),  x[i] = p_i,  x[i + L] = q_i,  i = 1...L  (periodic chain)
#   deviations     x[k·2L + i] = δp_i^(k),  x[k·2L + i + L] = δq_i^(k),  k = 1...deviations
#   amplitudes     ψ_i = (q_i + i p_i) / √2,  I_i = (p_i² + q_i²)/2 = |ψ_i|²,  S = Σ_i I_i
#
# NOT normalised to Σ_i I_i = 1: the pump sets the filling itself (S -> L P / Γ on the uniform
# state), so the initial condition only decides WHICH attractor is reached, never the density.
#
# In the (p, q) variables the equations of motion read
#
#   dp_j/dt =  J (q_{j-1} + q_{j+1}) + Δ_bare q_j - (g/2) r_j² q_j
#              + [P/2 - (Γ/4) r_j²] p_j + d (p_{j-1} + p_{j+1} - 2 p_j)
#   dq_j/dt = -J (p_{j-1} + p_{j+1}) - Δ_bare p_j + (g/2) r_j² p_j
#              + [P/2 - (Γ/4) r_j²] q_j + d (q_{j-1} + q_{j+1} - 2 q_j),      r_j² = p_j² + q_j²
#
# with Δ_bare = Δ - 2J, the same "bottom of the band" bookkeeping as in BHDissipative.jl (the
# hopping already supplies J z = 2J on a uniform state, z = 2). The first line of each equation is
# the Hamiltonian flow of
#
#   H = Σ_j [ -J (p_j p_{j+1} + q_j q_{j+1}) - (Δ_bare/2) r_j² + (g/8) r_j⁴ ]
#
# (PumpedEnergy below, = Energy of models/BoseHubbardFull.jl plus the detuning), the second is the
# Stuart-Landau gain/loss and the third the diffusion. With P = Γ = d = 0 everything reduces to the
# closed model of BHTrajectory.jl, and both H and S are then conserved.
#
# EXACT RESULTS, all checked in Demo() - they are what makes this model worth integrating:
#
#   plane waves      ψ_j = A e^{i(Qj - ωt)} exists iff |A|² = (P - 2 d D_Q)/Γ > 0, D_Q = 2(1-cos Q);
#                    its growth rate follows from 2x2 Bogoliubov blocks (PlaneWaveGrowth).
#   trapping region  for d = 0,  P/Γ ≤ S ≤ P L/Γ  is absorbing and the origin repels at rate P/2,
#                    so no trajectory escapes and none falls into the vacuum.
#   trace rule       Σ_k λ_k = (P - 2 d z) L - 2 Γ ⟨S⟩   (z = 2), an exact identity satisfied by the
#                    numerical spectrum - the sharpest test there is of the tangent dynamics.
#   Benjamin-Feir    in the continuum the uniform state is unstable iff Re[(d + iJ)(Γ/2 - i g)] < 0,
#                    i.e. J g < -d Γ/2. ON A LATTICE this is necessary but NOT sufficient: waves
#                    with |Q| > π/2 have an effective mass J cos Q of the opposite sign and stay
#                    stable at d = 0. A real d is what closes that escape (it makes |A|² fall with
#                    Q), and above d* = CriticalDiffusion(...) no plane wave is stable at all
#                    (d* ≈ 0.0939 at J = P = Γ = 1, g = -1; ≈ 0.073 on the L = 8 ring), so every
#                    trajectory has to end on something else - typically a strange attractor. That
#                    is not the only home of chaos: on the L = 8 map of BHMapPumped.jl strange
#                    attractors already coexist with the stable waves below d* for -1.7 ≤ g ≤ -0.4,
#                    and on a small ring a wave can turn stable again at larger d (CriticalDiffusion).
#
# Methods. The flow is deterministic - no dephasing and no truncated-Wigner noise, which is why
# there is no :split / :sde pair here; add them by reusing DephasingCallback and NoiseTerm! of
# BHDissipative.jl if they are ever needed. A plain ODEProblem of DifferentialEquations.jl is
# solved with DP8 and the tangent dynamics is carried along in the same state vector, exactly as in
# modules/ClassicalDynamics.jl. `deviations` selects how much of the spectrum is computed:
#
#   deviations = 1     largest exponent only (Benettin), the layout of EquationOfMotionTangentVector!
#   deviations = m     the m largest exponents; m = 2L gives the whole spectrum, hence n_pos, the
#                      trace rule and the Kaplan-Yorke dimension. The cost grows linearly with m.
#
# Reorthonormalisation is a modified Gram-Schmidt of the m deviation vectors, applied by a
# PeriodicCallback every `saveStep`: a QR done in place on views into u, so no matrix is ever formed
# and the growth factors are the diagonal of R.

using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Printf

include("models/BoseHubbardFull.jl")
include("modules/ClassicalDynamics.jl")


""" Bose-Hubbard parameters (L, J, U) of models/BoseHubbardFull.jl extended by the incoherent pump,
    the two-body loss and the diffusion. The result destructures as `L, J, U = parameters`, so every
    function of the closed model (Energy, InitialCondition, ...) accepts it unchanged, while the
    rest stays reachable by name.

    P = γ_pump - γ_loss (net linear gain, = -κ of BHDissipative.jl), Γ = two-body loss rate,
    d = diffusion (the real part of the coupling constant), Δ = detuning from the bottom of the
    band, pure gauge as long as there is no drive. The g of the Python code is 2U. """
function PumpedParameters(parameters; P = 0.0, Γ = 0.0, d = 0.0, Δ = 0.0)
    L, J, U = parameters
    return (L = L, J = float(J), U = float(U), P = float(P), Γ = float(Γ), d = float(d),
            Δ = float(Δ))
end


""" The same, straight from the constants of pump_bh_cgle.py, which works with g = 2U throughout and
    never mentions U. """
PumpedParametersG(L, J, g, P, Γ; d = 0.0, Δ = 0.0) =
    PumpedParameters((L, J, 0.5 * g); P = P, Γ = Γ, d = d, Δ = Δ)


""" Conservative part of the classical Hamiltonian,

        H = Σ_j [ -J (p_j p_{j+1} + q_j q_{j+1}) - (Δ_bare/2) r_j² + (g/8) r_j⁴ ],  Δ_bare = Δ - 2J,

    i.e. Energy(x, parameters) of models/BoseHubbardFull.jl plus the detuning (the `hamiltonian` of
    pump_bh_cgle.py). It is conserved in the conservative limit P = Γ = d = 0 only; with the pump on
    it is a diagnostic, saved alongside S by PumpedObservables. """
function PumpedEnergy(x, parameters)
    L, J, U, P, Γ, d, Δ = parameters

    energy = Energy(x, parameters)

    bareDetuning = Δ - 2 * J                    # Δ is measured from the bottom of the band
    if bareDetuning != 0
        @inbounds for i = 1:L
            energy -= 0.5 * bareDetuning * (x[i]^2 + x[i + L]^2)
        end
    end

    return energy
end


""" Total number of bosons S = Σ_i I_i = Σ_i |ψ_i|², the quantity the pump saturates. Conserved in
    the conservative limit only; on the uniform attractor it settles at L P / Γ. """
PumpedNorm(x, dimension) = 0.5 * sum(abs2, @view x[1:dimension])


""" Equations of motion of the incoherently pumped Bose-Hubbard model together with an arbitrary
    number of deviation vectors, in the layout of the header: u[1:2L] is the trajectory and
    u[k·2L + 1 : (k+1)·2L] the k-th deviation vector. Their number is read off the length of u, so
    one and the same right-hand side serves deviations = 1 (Benettin, the layout of
    EquationOfMotionTangentVector!) and deviations = 2L (the whole spectrum).

    Everything is matrix-free: the Jacobian entries of site j are evaluated once and applied to
    every deviation vector in turn, so the cost is linear in the number of deviation vectors and no
    2L x 2L matrix is ever formed. The linear terms (hopping, detuning, gain, diffusion) reproduce
    themselves on the deviation vectors; only the cubic terms have nontrivial derivatives

        ∂(dp)/∂p = P/2 - (Γ/4)(r² + 2p²) - g p q         ∂(dp)/∂q = -(g/2)(r² + 2q²) - (Γ/2) p q
        ∂(dq)/∂p =  (g/2)(r² + 2p²) - (Γ/2) p q          ∂(dq)/∂q = P/2 - (Γ/4)(r² + 2q²) + g p q

    whose trace P - Γ r_j² = P - 2 Γ I_j is what the trace rule of the header adds up. """
function EquationOfMotionPumped!(dx, x, parameters, t)
    L, J, U, P, Γ, d, Δ = parameters.modelParameters
    dimension = 2 * L

    g = 2 * U
    bareDetuning = Δ - 2 * J
    deviations = div(length(x), dimension) - 1

    @inbounds for i = 1:L
        p = x[i]
        q = x[i + L]
        r2 = p * p + q * q

        im1 = i == 1 ? L : i - 1
        ip1 = i == L ? 1 : i + 1

        kerr = 0.5 * g * r2                     # conservative phase rotation, = U r²
        relaxation = 0.5 * P - 0.25 * Γ * r2    # linear gain saturated by the two-body loss

        dx[i] = J * (x[im1 + L] + x[ip1 + L]) + (bareDetuning - kerr) * q + relaxation * p +
                d * (x[im1] + x[ip1] - 2 * p)
        dx[i + L] = -J * (x[im1] + x[ip1]) - (bareDetuning - kerr) * p + relaxation * q +
                    d * (x[im1 + L] + x[ip1 + L] - 2 * q)

        deviations == 0 && continue

        # Jacobian of the on-site (cubic + gain) terms, evaluated once per site for all deviations
        pq = p * q
        dpdp = relaxation - 0.5 * Γ * p * p - g * pq
        dpdq = -0.5 * g * (r2 + 2 * q * q) - 0.5 * Γ * pq
        dqdp = 0.5 * g * (r2 + 2 * p * p) - 0.5 * Γ * pq
        dqdq = relaxation - 0.5 * Γ * q * q + g * pq

        for k = 1:deviations
            offset = k * dimension

            δp = x[offset + i]
            δq = x[offset + i + L]

            dx[offset + i] = J * (x[offset + im1 + L] + x[offset + ip1 + L]) +
                             bareDetuning * δq + dpdp * δp + dpdq * δq +
                             d * (x[offset + im1] + x[offset + ip1] - 2 * δp)
            dx[offset + i + L] = -J * (x[offset + im1] + x[offset + ip1]) -
                                 bareDetuning * δp + dqdp * δp + dqdq * δq +
                                 d * (x[offset + im1 + L] + x[offset + ip1 + L] - 2 * δq)
        end
    end

    return nothing
end


""" Non-mutating reader for the SavingCallback that records the conservative energy and the total
    number of bosons S = Σ_i I_i. Neither is conserved once the pump is on; S is the interesting
    one, since it is the attractor that fixes it. """
function PumpedObservables(x, t, integrator)
    dimension = integrator.p.dimension
    return (PumpedEnergy(x, integrator.p.modelParameters), PumpedNorm(x, dimension))
end


""" Reorthonormalisation of the deviation vectors - the affect! of the PeriodicCallback that fires
    every `saveStep`, and the only place where the exponents are accumulated.

    A modified Gram-Schmidt applied in place to the views u[k·2L + 1 : (k+1)·2L], i.e. a QR without
    ever forming a matrix; log of the k-th growth factor (the k-th diagonal entry of R) is added to
    `spectrum[k]`, so after division by the elapsed time spectrum holds the `deviations` largest
    Lyapunov exponents in descending order.

    The largest one is additionally handed to AccumulateLyapunov! of modules/ClassicalDynamics.jl,
    which keeps the running-exponent history, the convergence test and the plot data shared with
    every other trajectory routine of the project. With deviations = 1 this reduces exactly to
    RescaleDeviationVector!.

    `accumulationWindow` records the measurement window actually covered: [1] is set to the time of
    the FIRST rescale after relaxationTime, whose own interval is discarded (it straddles
    relaxationTime and is the only one not entirely inside the window), and [2] is moved to every
    rescale counted afterwards. Dividing the sums by [2] - [1] rather than by the nominal
    (timeEnd - relaxationTime) is what makes the spectrum EXACT rather than short by whichever
    fraction of an interval the two ends of the run leave over - visible as a 0.1% miss in the trace
    rule at saveStep = 2 over a window of 2000.

    Discarding the straddling interval also sidesteps the division by t - relaxationTime inside
    AccumulateLyapunov!, which is a division by ~0 whenever a rescale happens to land on
    relaxationTime itself - and it lands there exactly when relaxationTime is a multiple of
    saveStep, i.e. in the most natural choice of both. """
function RescaleDeviations!(deviations, spectrum, accumulationWindow)
    function Rescale!(integrator)
        dimension = integrator.p.dimension
        growth = 1.0

        accumulate = integrator.t > integrator.p.relaxationTime
        if accumulate && isnan(accumulationWindow[1])
            accumulationWindow[1] = integrator.t
            accumulate = false
        end

        @inbounds for k = 1:deviations
            vk = view(integrator.u, (k * dimension + 1):((k + 1) * dimension))

            for l = 1:(k - 1)                   # orthogonalise against the already normalised ones
                vl = view(integrator.u, (l * dimension + 1):((l + 1) * dimension))
                axpy!(-dot(vl, vk), vl, vk)
            end

            growthK = sqrt(sum(abs2, vk))
            vk ./= growthK

            k == 1 && (growth = growthK)
            accumulate && (spectrum[k] += log(growthK))
        end

        u_modified!(integrator, true)

        if accumulate
            accumulationWindow[2] = integrator.t
            AccumulateLyapunov!(integrator, growth)
        end
    end

    return Rescale!
end


""" Random initial condition of the pumped model: a random direction in phase space scaled to a
    total density S = Σ_i I_i drawn uniformly in (0, L P / Γ].

    L P / Γ is the top of the trapping annulus P/Γ ≤ S ≤ L P/Γ of the header, i.e. the uniform
    state, and is the natural counterpart of the absorbing ball of the driven model: sampling it
    rather than starting near one particular state is what makes a chaotic fraction a basin measure.
    (With P ≤ 0 or Γ ≤ 0 there is no annulus and the fallback scale S ~ L is used.) """
function PumpedInitialCondition(parameters, rng = Random.default_rng())
    L, J, U, P, Γ, d, Δ = parameters

    maximumNorm = (P > 0 && Γ > 0) ? L * P / Γ : float(L)
    x = randn(rng, 2 * L)
    x .*= sqrt(2 * rand(rng) * maximumNorm / sum(abs2, x))

    return x
end


""" Largest Lyapunov exponent - or the `deviations` largest ones - of an individual trajectory of
    the incoherently pumped Bose-Hubbard model. The undriven, two-body-loss counterpart of
    TrajectoryLyapunovDissipative of BHDissipative.jl, with which it shares the phase-space layout,
    the integration parameters and all the callbacks.

    Nothing is conserved, so there is no ManifoldProjection and no Poincaré section, and the
    convergence test is switched off by default (regularThreshold = relativeFluctuationThreshold =
    0), i.e. the trajectory is always integrated over the whole `timeInterval`.

    deviations      1 = largest exponent only (cheapest); 2L = the whole spectrum, which is what the
                    trace rule and the Kaplan-Yorke dimension need
    saveStep        how often the deviation vectors are reorthonormalised = resolution of the output
    relaxationTime  the exponents are accumulated only after this time. The trajectory forgets its
                    initial condition and settles on an attractor, so make this long enough for the
                    transient to die out - otherwise the transient, not the attractor, is measured.

    Returns - the solution, the largest exponent Λ (trailing average of the running exponent), the
              whole spectrum (descending, length `deviations`, averaged over the accumulation
              window), the SavedValues of the running exponent and the SavedValues of (energy, S).

    A FAILED trajectory is flagged by Λ = NaN, not by Λ = 0 as in TrajectoryLyapunovDissipative:
    here zero is a perfectly physical answer. U(1) is unbroken, so every attractor carries the exact
    zero of the phase mode, and Λ = 0 is what a fixed point or a limit cycle actually gives. """
function TrajectoryLyapunovPumped(initialCondition, parameters;
        deviations = 1,
        solver = DP8(),
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

    L, J, U, P, Γ, d, Δ = parameters
    phaseSpaceDimension = length(initialCondition)          # = 2L
    rng = Xoshiro(seed)

    deviations >= 1 || error("at least one deviation vector is needed")
    deviations <= phaseSpaceDimension ||
        error("at most 2L = $phaseSpaceDimension deviation vectors span the tangent space")

    # Initial condition + `deviations` orthonormal deviation vectors; for deviations = 1 this is
    # exactly the layout of EquationOfMotionTangentVector! (tangentDynamics = :vector)
    x0 = zeros(phaseSpaceDimension * (1 + deviations))
    x0[1:phaseSpaceDimension] = initialCondition
    for k = 1:deviations
        vk = view(x0, (k * phaseSpaceDimension + 1):((k + 1) * phaseSpaceDimension))
        vk .= randn(rng, phaseSpaceDimension)
        for l = 1:(k - 1)
            vl = view(x0, (l * phaseSpaceDimension + 1):((l + 1) * phaseSpaceDimension))
            axpy!(-dot(vl, vk), vl, vk)
        end
        vk ./= sqrt(sum(abs2, vk))
    end

    energy = PumpedEnergy(x0, parameters)
    norm0 = PumpedNorm(x0, phaseSpaceDimension)

    # maximumSectionPoints = 0 -> no Poincaré section, the stopping decision is taken in AccumulateLyapunov!
    integrationParameters = LyapunovIntegrationParameters(phaseSpaceDimension, parameters, energy, relaxationTime,
        relativeFluctuationThreshold, regularThreshold, 0, 0, time_ns(), 1E9 * timeout, :Start,
        0.0, norm0, historyLyapunovExponentLength, Float64[])

    spectrum = zeros(deviations)                            # Σ log(growth) of every deviation vector
    accumulationWindow = [NaN, NaN]                         # [first, last] rescale actually counted

    lyapunovs = SavedValues(Float64, Float64)               # The whole history of the immediate Lyapunov exponents (for a graph)
    observables = SavedValues(Float64, Tuple{Float64, Float64})     # Energy and norm, neither of them conserved

    rescale = PeriodicCallback(RescaleDeviations!(deviations, spectrum, accumulationWindow), float(saveStep); save_positions=(false, false))
    record = SavingCallback(RunningLyapunov, lyapunovs, saveat=saveStep:saveStep:last(timeInterval))
    recordObservables = SavingCallback(PumpedObservables, observables, saveat=saveStep:saveStep:last(timeInterval))
    timeoutCallback = DiscreteCallback(TimeoutCondition, terminate!)

    callback = CallbackSet(rescale, record, recordObservables, timeoutCallback)

    problem = ODEProblem(ODEFunction(EquationOfMotionPumped!), x0, timeInterval, integrationParameters)
    time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, callback=callback,
        save_on=true, save_everystep=false, save_start=true, save_end=true, maxiters=maximumIterations,
        isoutofdomain=CheckDomain, verbose=DEVerbosity())

    # Get results - average over the trailing window of running-exponent samples (empty -> exactly 0)
    history = integrationParameters.historyLyapunovExponent
    window = @view history[max(1, length(history) - historyLyapunovExponentLength + 1):end]
    lyapunov = isempty(window) ? NaN : mean(window)     # NaN = nothing was measured, see the docstring
    lv = length(window) > 1 ? var(window) : 0.0

    # The spectrum is accumulated over the whole window from the first rescale past relaxationTime
    # to the end of the run, without any trailing average - a single exponent can be averaged over
    # the plateau of its running value, a whole spectrum cannot.
    accumulationTime = accumulationWindow[2] - accumulationWindow[1]
    spectrum = accumulationTime > 0 ? spectrum ./ accumulationTime : fill(NaN, deviations)

    if length(solution.u) > 0
        finalState = solution.u[end]
        @info "Calculation time = $time, Trajectory time = $(solution.t[end]), Final energy = $(PumpedEnergy(finalState, parameters)), Final norm = $(PumpedNorm(finalState, phaseSpaceDimension)), Λ = $lyapunov ± $lv"
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
        return solution, NaN, fill(NaN, deviations), lyapunovs, observables
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

    return solution, lyapunov, spectrum, lyapunovs, observables
end


""" Kaplan-Yorke (Lyapunov) dimension of a spectrum: the interpolated number of exponents whose
    running sum is still nonnegative. 0 when even the largest one is negative (a fixed point),
    length(spectrum) when the whole sum is - which on an attractor cannot happen. """
function KaplanYorkeDimension(spectrum)
    sorted = sort(collect(spectrum), rev=true)
    cumulative = cumsum(sorted)

    cumulative[1] < 0 && return 0.0

    k = findlast(>=(0), cumulative)
    return k == length(sorted) ? float(length(sorted)) : k + cumulative[k] / abs(sorted[k + 1])
end


# ---------------------------------------------------------------- plane waves and their stability

""" Squared amplitude of the plane wave ψ_j = A e^{i(Qj - ωt)}: |A|² = (P - 2 d D_Q)/Γ with
    D_Q = 2(1 - cos Q). Nonpositive means the wave does not exist - the diffusion damps that
    wavenumber faster than the pump feeds it. """
PlaneWaveAmplitude(P, Γ, d, Q) = (P - 2 * d * 2 * (1 - cos(Q))) / Γ


""" Detuning that makes the plane wave of wavenumber Q a genuine FIXED POINT instead of a rotating
    solution: Δ = J D_Q + g |A|², in the band-bottom convention of PumpedParameters. Used by Demo()
    to measure the Bogoliubov rate below directly with the integrator. """
PlaneWaveDetuning(J, g, P, Γ, d, Q) =
    J * 2 * (1 - cos(Q)) + g * PlaneWaveAmplitude(P, Γ, d, Q)


""" Largest growth rate max_k Re λ_k of the plane wave of wavenumber Q on a ring of L sites, from
    the analytic 2x2 Bogoliubov blocks (the pw_growth of pump_bh_cgle.py).

    The perturbation ψ_j = (A + u e^{ikj} + v* e^{-ikj}) e^{i(Qj - ωt)} mixes only the modes
    Q ± k, so the stability problem splits into 2x2 blocks, one per k = 2πm/L:

        A11 = W (D_Q - D_{Q+k}) - c,  A22 = conj(W) (D_Q - D_{Q-k}) - conj(c),
        W = d + i J,  c = (i g + Γ/2) |A|²,  D_q = 2(1 - cos q)

    with eigenvalues [tr ± √(tr² - 4 det)]/2, tr = A11 + A22, det = A11 A22 - |c|².
    Returns `nothing` when the wave does not exist. """
function PlaneWaveGrowth(L, J, g, P, Γ, d, Q)
    amplitude2 = PlaneWaveAmplitude(P, Γ, d, Q)
    amplitude2 > 0 || return nothing

    W = complex(d, J)
    c = (im * g + 0.5 * Γ) * amplitude2
    DQ = 2 * (1 - cos(Q))

    largest = -Inf
    for m = 0:(L - 1)
        k = 2π * m / L

        A11 = W * (DQ - 2 * (1 - cos(Q + k))) - c
        A22 = conj(W) * (DQ - 2 * (1 - cos(Q - k))) - conj(c)

        trace = A11 + A22
        determinant = A11 * A22 - c * conj(c)
        root = sqrt(trace * trace - 4 * determinant)

        largest = max(largest, real(0.5 * (trace + root)), real(0.5 * (trace - root)))
    end

    return largest
end


""" True if at least one of the L plane waves of the ring exists and is linearly stable at diffusion
    d. The U(1) phase mode makes the largest rate of a stable wave exactly 0, hence the tolerance
    instead of a strict inequality. """
function AnyStablePlaneWave(L, J, g, P, Γ, d)
    return any(0:(L - 1)) do m
        growth = PlaneWaveGrowth(L, J, g, P, Γ, d, 2π * m / L)
        growth !== nothing && growth < 1e-10
    end
end


""" Onset diffusion d*: the smallest d at which NONE of the L plane waves of the ring is linearly
    stable any more - the d* of pump_bh_cgle.py. Below it every trajectory has somewhere regular to
    land (d* ≈ 0.0939 at J = P = Γ = 1, g = -1 on the default L = 256, and ≈ 0.073 on the L = 8 ring
    of BHMapPumped.jl). When a stable wave survives all the way, as it does for J g > 0, d* = ∞ is
    reported as `maximumDiffusion`.

    Found by a SCAN in steps of maximumDiffusion/steps, refined by bisection inside the first
    bracket - not by bisection alone, as d_star of pump_bh_cgle.py does. Bisection assumes that
    stability, once lost, never comes back as d grows, and on a small ring that is false: at L = 8
    and -1.5 ≤ g ≤ -0.5 every wave is unstable from d ≈ 0.073 up to a second threshold (0.34 at
    g = -0.5, 0.92 at g = -1, 1.35 at g = -1.5), above which a wave is stable again. Bisection then
    lands in that upper window, reports ∞, and the onset disappears. A window without a stable wave
    that is narrower than one scan step can still be missed. """
function CriticalDiffusion(J, g, P, Γ; L = 256, maximumDiffusion = 3.0, steps = 600, iterations = 40)
    AnyStablePlaneWave(L, J, g, P, Γ, 0.0) || return 0.0

    stepSize = maximumDiffusion / steps

    for n = 1:steps
        high = n * stepSize
        AnyStablePlaneWave(L, J, g, P, Γ, high) && continue

        low = high - stepSize                   # the previous scan point, still stable
        for _ = 1:iterations
            middle = 0.5 * (low + high)
            low, high = AnyStablePlaneWave(L, J, g, P, Γ, middle) ? (middle, high) : (low, middle)
        end

        return high
    end

    return maximumDiffusion
end


""" Number of plane waves of the ring that exist, and of those that are linearly stable, at a given
    diffusion - section C of pump_bh_cgle.py. """
function PlaneWaveCount(L, J, g, P, Γ, d)
    growths = [PlaneWaveGrowth(L, J, g, P, Γ, d, 2π * m / L) for m = 0:(L - 1)]
    existing = filter(!isnothing, growths)

    return length(existing), count(<(1e-9), existing)
end


""" Demonstration of the pumped model, together with the consistency checks of the integration: the
    conservative limit, the analytic Bogoliubov rates of the plane waves, the lattice escape that
    the diffusion closes, and the trace rule and Kaplan-Yorke dimension of the full spectrum.
    Mirrors sections A-D of the `__main__` block of pump_bh_cgle.py, whose numbers it reproduces. """
function Demo()
    J, g, P, Γ = 1.0, -1.0, 1.0, 1.0            # the constants of pump_bh_cgle.py; J g < 0 is needed
    z = 2                                       # coordination number of the ring

    # 1. Conservative limit P = Γ = d = 0: the closed Bose-Hubbard model of BHTrajectory.jl, where
    #    both the Hamiltonian and the total density are conserved. Tests the hopping, the Kerr term
    #    and the detuning bookkeeping (Δ ≠ 0 must not break either conservation law).
    L = 8
    parameters = PumpedParametersG(L, J, g, 0.0, 0.0; Δ=0.7)
    Random.seed!(1234)
    x = randn(2 * L)

    _, lyapunov, _, _, observables = TrajectoryLyapunovPumped(x, parameters; timeInterval=(0.0, 1000.0))
    energies = first.(observables.saveval)
    norms = last.(observables.saveval)
    @printf("conservative limit:   Λ = %+.4f   rel. energy drift = %.1e   rel. norm drift = %.1e\n",
        lyapunov, maximum(abs.(energies .- energies[1])) / abs(energies[1]),
        maximum(abs.(norms .- norms[1])) / norms[1])

    # 2. Plane waves. With Δ = J D_Q + g|A|² the wave is a FIXED POINT, so the exponent measured by
    #    the integrator has to reproduce the analytic Bogoliubov rate of PlaneWaveGrowth - an
    #    end-to-end test of the equations of motion, of the tangent dynamics and of the conventions.
    #    All four waves below are UNSTABLE fixed points, so the measurement window has to end
    #    before the rounding error of the initial condition - 1e-16 of the amplitude - is amplified
    #    to something visible: after t the error is e^{growth t}, so a window of a fixed number of
    #    e-foldings, 20/growth, is both safe (1e-16 e^20 ≈ 5e-8) and equally accurate for every Q.
    #    That cap is also what limits the accuracy to a few per cent: a finite-time exponent
    #    measured over N e-foldings carries the alignment offset log(c₁)/N of the initial deviation
    #    vector, and N ≤ 20 here. The high-precision test of the same tangent dynamics is the trace
    #    rule of part 4 below, which needs no fixed point and holds to five digits.
    L, d = 16, 0.1
    println("plane waves, L = 16, d = 0.1 (measured vs the analytic Bogoliubov rate, few % apart):")
    for m = 0:2:6
        Q = 2π * m / L
        amplitude = sqrt(PlaneWaveAmplitude(P, Γ, d, Q))
        growth = PlaneWaveGrowth(L, J, g, P, Γ, d, Q)
        parameters = PumpedParametersG(L, J, g, P, Γ; d=d, Δ=PlaneWaveDetuning(J, g, P, Γ, d, Q))

        # ψ_j = A e^{iQj}, i.e. q_j = √2 A cos(Qj), p_j = √2 A sin(Qj)
        x = vcat([sqrt(2) * amplitude * sin(Q * (j - 1)) for j = 1:L],
                 [sqrt(2) * amplitude * cos(Q * (j - 1)) for j = 1:L])

        relaxationTime = 5 / growth             # time for the deviation vector to align
        timeEnd = 20 / growth
        saveStep = 0.05 * relaxationTime

        # ... and the running exponent is read on the last quarter of the window, where it has
        # settled - averaging it over the whole window would still be dragged down by the start.
        _, lyapunov, _, _, observables = TrajectoryLyapunovPumped(x, parameters;
            timeInterval=(0.0, timeEnd), saveStep=saveStep, relaxationTime=relaxationTime,
            historyLyapunovExponentLength=round(Int, 0.25 * (timeEnd - relaxationTime) / saveStep))

        @printf("   Q = %5.3f  |A|² = %+.4f   measured %+.6f   analytic %+.6f   S: %.4f -> %.4f\n",
            Q, amplitude^2, lyapunov, growth, L * amplitude^2, last(last.(observables.saveval)))
    end

    # 3. The lattice escape: the Benjamin-Feir criterion J g < -d Γ/2 is necessary but not
    #    sufficient, because short waves (|Q| > π/2) stay stable at d = 0. Pure analytics.
    println("stable plane waves of the ring (L = 64):")
    for d in (0.0, 0.05, 0.094, 0.1, 0.25)
        existing, stable = PlaneWaveCount(64, J, g, P, Γ, d)
        @printf("   d = %-6.3f exist %3d/64   stable %3d\n", d, existing, stable)
    end
    @printf("   d* (no stable plane wave above this) = %.5f\n", CriticalDiffusion(J, g, P, Γ))

    # 4. The whole Lyapunov spectrum, hence n_pos, the trace rule and the Kaplan-Yorke dimension.
    #    The trace rule Σ λ = (P - 2 d z) L - 2 Γ ⟨S⟩ is exact and is the sharpest test of the
    #    tangent dynamics there is; ⟨S⟩ is measured on the same run.
    L = 12
    println("Lyapunov spectra, L = $L (2L = $(2 * L) exponents):")
    for d in (0.0, 0.1, 0.25)
        parameters = PumpedParametersG(L, J, g, P, Γ; d=d)
        Random.seed!(5)
        x = PumpedInitialCondition(parameters)

        _, _, spectrum, _, observables = TrajectoryLyapunovPumped(x, parameters;
            deviations=2 * L, timeInterval=(0.0, 2500.0), relaxationTime=500.0)

        norms = last.(observables.saveval)[observables.t .>= 500.0]
        meanNorm = mean(norms)

        @printf("\n   d = %.2f   ⟨S⟩ = %.4f  (uniform value P L/Γ = %.1f)\n", d, meanNorm, P * L / Γ)
        @printf("     λ[1:6]  = %s\n", join((@sprintf("%+.4f", λ) for λ in spectrum[1:6]), "  "))
        @printf("     n_pos = %2d   Σ λ = %+.5f\n", count(>(1e-4), spectrum), sum(spectrum))
        @printf("     (P - 2 d z) L - 2 Γ ⟨S⟩ = %+.5f   <-- trace rule\n", (P - 2 * d * z) * L - 2 * Γ * meanNorm)
        @printf("     D_KY = %.4f  of %d\n", KaplanYorkeDimension(spectrum), 2 * L)
    end
end

# Runs only when the file is executed directly (julia BHPumped.jl), so that include("BHPumped.jl")
# from BHMapPumped.jl and friends stays silent.
if abspath(PROGRAM_FILE) == @__FILE__
    Demo()
end
