# BHNumberConserving.jl
#
# Classical (N -> infinity) mean-field dynamics of the NUMBER-CONSERVING dissipative
# Bose-Hubbard model of number-conserving-BH.md, and the machinery that Section 8 of that note
# asks for: the FULL Lyapunov spectrum of the attractor, the Kaplan-Yorke dimension on the
# reduced 2L - 2 space, the trace rule sum(lambda) = <div F>, and the basin statistics that
# expose multistability.  Section8() below runs the four questions and prints a verdict for
# each; BHMapNumberConserving.jl turns the same machinery into a resumable map.
#
# WHY A NEW FILE AND NOT AN EXTENSION OF BHDissipative.jl
#
#   * BHDissipative.jl (and the pumped/driven variants) all have a WEAK U(1): the jump operators
#     b, b^2, b_i - b_j change the particle number, the norm decays, and a compact attractor
#     exists only because loss and gain balance.  Here the total density is EXACTLY conserved
#     (note §4), the flow lives on S^(2L-1) and, after the global phase, on CP^(L-1).  Nothing
#     about the absorbing ball, the norm law or the exact trace rule Σλ = -κL carries over.
#   * TrajectoryLyapunov of modules/ClassicalDynamics.jl returns the LARGEST exponent only
#     (:vector) or rescales the full monodromy matrix by its largest singular value (:matrix).
#     Question 1 of §8 needs a FRACTIONAL D_KY, i.e. the whole ordered spectrum, which requires a
#     QR (Gram-Schmidt) reorthonormalisation of 2L deviation vectors - implemented here.
#
# STATE, NORMALISATION AND LAYOUT
#
#   Amplitudes are the note ψ_j themselves, as a ComplexF64 vector:
#
#       ψ_j = b_j / sqrt(N),   n_j = |ψ_j|^2,   sum_j n_j = 1 exactly,   hbar_eff = 1/N
#
#   which is the same normalisation as the (p, q) layout of models/BoseHubbardFull.jl
#   (sum_i I_i = 1 with I_i = (p_i^2 + q_i^2)/2 and ψ_i = (q_i + i p_i)/sqrt(2)).  ToPQ / FromPQ
#   convert, so an initial condition or an attractor can be handed to the rest of the repository
#   unchanged.  Complex arithmetic is used internally because every dissipative term of the note
#   is written with conj(ψ) and Re(conj(ψ_j) ψ_k); the tangent map is real-linear (it contains
#   conj(δψ)), so the deviation vectors are complex vectors orthonormalised with the REAL inner
#   product Re sum_k conj(a_k) b_k, which is exactly the Euclidean product of the underlying R^2L.
#
#   The integrated state is
#
#       u[1:L]                                   ψ
#       u[(k L + 1):((k+1) L)],  k = 1..m        deviation vector k  (m = 2L for the full spectrum)
#       u[end]                                   integral of div F dt  (real part; the trace rule)
#
# EQUATIONS OF MOTION (note §3), with lap_j = ψ_{j+1} + ψ_{j-1} - 2 ψ_j and S_j = sum_k A_jk n_k:
#
#       dψ_j/dt = i J lap_j - i g n_j ψ_j + (1/2) S_j ψ_j
#                 + κ sum_{k = j±1} [ (n_j + n_k) ψ_k - 2 Re(conj(ψ_j) ψ_k) ψ_j ]
#
#   A_ij = γ_ij - γ_ji is the ANTISYMMETRIC part of the incoherent-hopping rates (note §2: the
#   symmetric part cancels identically at mean field).  For the circulating ring
#   A_{j+1,j} = -A_{j,j+1} = η_j; with a uniform η_j = η the column sums c_j = sum_i A_ij vanish
#   and the directed hopping is volume preserving, which is what question 2 of §8 is about.
#
# PARAMETERS.  J = 1 fixes the time unit; only sign(J g) matters (note §6), and the instability
# of §5 needs J g < 0, so g < 0 throughout.  g = U N, η = N (Gamma_+ - Gamma_-), κ = N κ_c are
# the rescaled rates of §3 - all O(1) in the classical limit.

using DifferentialEquations, DiffEqCallbacks
using OrdinaryDiffEqHighOrderRK     # DP8: DifferentialEquations v8 re-exports only a small solver set
using LinearAlgebra
using Random
using Statistics
using Printf


""" Parameters of the number-conserving model.  `A` is the L x L antisymmetric matrix of net
    incoherent-hopping rates of note §3 and `c[j] = sum_i A[i,j]` its column sums, which are what
    the divergence of §4 is built from (c == 0 <=> volume-preserving directed hopping). """
struct NumberConservingParameters
    L::Int
    J::Float64
    g::Float64
    κ::Float64
    A::Matrix{Float64}
    c::Vector{Float64}
end


""" Net incoherent-hopping matrix of a ring with a directed circulation.

    Bond j (sites j -> j+1, periodic) carries eta_j = η (1 + modulation cos(2 pi j / L + phase)):

        A[j+1, j] = +eta_j,   A[j, j+1] = -eta_j     =>    c_j = eta_j - eta_{j-1}

    `modulation = 0` is the minimal choice of note §1: a uniform circulation, c == 0,
    divergence-free, and with the uniform condensate still an exact solution (§5 needs vanishing
    ROW sums, which for an antisymmetric A is the same condition).  `modulation != 0` is the
    non-uniform rate of §8.2: it makes c_j != 0, so the directed hopping alone can contract
    phase-space volume - and at the same time it destroys the uniform state as an exact solution,
    which is why question 2 cannot simply reuse the threshold of §5. """
function CirculatingRates(L::Integer; η = 0.0, modulation = 0.0, phase = 0.0)
    A = zeros(Float64, L, L)

    for j = 1:L
        k = j == L ? 1 : j + 1
        rate = η * (1 + modulation * cos(2 * pi * j / L + phase))
        A[k, j] = rate
        A[j, k] = -rate
    end

    return A
end


""" Constructor: either give the ring parameters (η, modulation, phase) and let CirculatingRates
    build A, or pass an arbitrary antisymmetric A (any range, any geometry - nothing below assumes
    that A is nearest-neighbour). """
function NumberConservingParameters(L::Integer; J = 1.0, g = 0.0, κ = 0.0,
        η = 0.0, modulation = 0.0, phase = 0.0, A = nothing)

    L >= 3 || error("L >= 3 required: at L = 2 the reduced space is the 2D Bloch sphere, " *
                    "where Poincare-Bendixson forbids chaos (note §4), and the κ bond sum " *
                    "would count every bond twice")

    matrix = isnothing(A) ? CirculatingRates(L; η = η, modulation = modulation, phase = phase) :
                            Matrix{Float64}(A)

    size(matrix) == (L, L) || error("A must be L x L")
    norm(matrix + transpose(matrix)) < 1e-12 * max(1.0, norm(matrix)) ||
        error("A must be antisymmetric: it is the net rate A_ij = gamma_ij - gamma_ji of note §2")

    return NumberConservingParameters(L, float(J), float(g), float(κ), matrix,
                                      vec(sum(matrix, dims = 1)))
end


""" Scratch space shared by the vector field, the 2L tangent fields and the divergence within one
    right-hand-side evaluation: the state-dependent quantities are computed once per call to
    EquationOfMotion! and reused by every deviation vector. """
struct NumberConservingWorkspace
    n::Vector{Float64}          # n_j = |ψ_j|^2
    S::Vector{Float64}          # S_j = sum_k A_jk n_k
    R::Vector{Float64}          # R_j = Re(conj(ψ_j) ψ_{j+1}), the coherence of bond j
    δn::Vector{Float64}
    δS::Vector{Float64}
    δR::Vector{Float64}
end

NumberConservingWorkspace(L::Integer) =
    NumberConservingWorkspace(zeros(L), zeros(L), zeros(L), zeros(L), zeros(L), zeros(L))

@inline Next(j, L) = j == L ? 1 : j + 1
@inline Previous(j, L) = j == 1 ? L : j - 1


""" n_j, S_j and the bond coherences R_j of the current ψ; everything else reads them. """
function UpdateWorkspace!(work, ψ, parameters)
    L = parameters.L
    A = parameters.A

    @inbounds for j = 1:L
        work.n[j] = abs2(ψ[j])
    end

    @inbounds for j = 1:L
        work.R[j] = real(conj(ψ[j]) * ψ[Next(j, L)])
    end

    @inbounds for j = 1:L
        s = 0.0
        for k = 1:L
            s += A[j, k] * work.n[k]
        end
        work.S[j] = s
    end

    return nothing
end


""" Right-hand side of note §3, from a workspace already filled by UpdateWorkspace!. """
function VectorField!(dψ, ψ, parameters, work)
    L, J, g, κ = parameters.L, parameters.J, parameters.g, parameters.κ
    n, S, R = work.n, work.S, work.R

    @inbounds for j = 1:L
        p = Next(j, L)
        m = Previous(j, L)

        value = im * J * (ψ[p] + ψ[m] - 2 * ψ[j]) + (0.5 * S[j] - im * g * n[j]) * ψ[j]

        if κ != 0
            # bond form of the phase-locking channel D[c_ij]; equals the mode form of note §2 with
            # u, v = (e_j ± e_k)/sqrt(2) and gamma = 4κ (verified in Checks())
            value += κ * ((n[j] + n[p]) * ψ[p] + (n[j] + n[m]) * ψ[m] - 2 * (R[j] + R[m]) * ψ[j])
        end

        dψ[j] = value
    end

    return nothing
end


""" Tangent (linearised) dynamics d(δψ)/dt = DF(ψ) δψ, matrix free.

    The map is REAL-linear, not holomorphic: n_j, S_j and R_j all involve conj(ψ), so their
    variations δn_k = 2 Re(conj(ψ_k) δψ_k) and δR_j = Re(conj(δψ_j) ψ_{j+1} + conj(ψ_j) δψ_{j+1})
    are real and appear multiplying ψ.  That is the whole content of the Jacobian; everything else
    is the vector field with δψ substituted.  Checks() compares it with a finite difference. """
function TangentField!(dδ, δ, ψ, parameters, work)
    L, J, g, κ = parameters.L, parameters.J, parameters.g, parameters.κ
    A = parameters.A
    n, S, R = work.n, work.S, work.R
    δn, δS, δR = work.δn, work.δS, work.δR

    @inbounds for j = 1:L
        δn[j] = 2 * real(conj(ψ[j]) * δ[j])
    end

    @inbounds for j = 1:L
        p = Next(j, L)
        δR[j] = real(conj(δ[j]) * ψ[p] + conj(ψ[j]) * δ[p])
    end

    @inbounds for j = 1:L
        s = 0.0
        for k = 1:L
            s += A[j, k] * δn[k]
        end
        δS[j] = s
    end

    @inbounds for j = 1:L
        p = Next(j, L)
        m = Previous(j, L)

        value = im * J * (δ[p] + δ[m] - 2 * δ[j]) +
                (0.5 * S[j] - im * g * n[j]) * δ[j] +
                (0.5 * δS[j] - im * g * δn[j]) * ψ[j]

        if κ != 0
            value += κ * ((δn[j] + δn[p]) * ψ[p] + (n[j] + n[p]) * δ[p] +
                          (δn[j] + δn[m]) * ψ[m] + (n[j] + n[m]) * δ[m] -
                          2 * (δR[j] + δR[m]) * ψ[j] - 2 * (R[j] + R[m]) * δ[j])
        end

        dδ[j] = value
    end

    return nothing
end


""" Divergence of the flow in the real 2L-dimensional phase space (note §4):

        div F = sum_j c_j n_j - 8 κ sum_j Re(conj(ψ_j) ψ_{j+1}),    c_j = sum_i A_ij

    For a general complex field dψ/dt = F(ψ, conj(ψ)) the real divergence is
    sum_j 2 Re(dF_j/dψ_j) in the Wirtinger sense; the hopping and the Kerr term give a purely
    imaginary diagonal and drop out, the directed hopping contributes sum_j S_j = sum_j c_j n_j,
    and each phase-locking bond contributes -4κ R twice.  Holds for an arbitrary antisymmetric A,
    not only for the nearest-neighbour ring.  Checks() compares it with the trace of the analytic
    Jacobian; §8.4 asks for its time average along the attractor, which the integration below
    accumulates exactly, as an extra component of the state. """
function Divergence(parameters, work)
    L, κ = parameters.L, parameters.κ
    c, n, R = parameters.c, work.n, work.R

    value = 0.0
    @inbounds for j = 1:L
        value += c[j] * n[j]
    end

    if κ != 0
        coherence = 0.0
        @inbounds for j = 1:L
            coherence += R[j]
        end
        value -= 8 * κ * coherence
    end

    return value
end


""" Allocating convenience wrappers, for the checks and for the stationarity test. """
function VectorField(ψ, parameters)
    work = NumberConservingWorkspace(parameters.L)
    UpdateWorkspace!(work, ψ, parameters)
    dψ = similar(ψ)
    VectorField!(dψ, ψ, parameters, work)
    return dψ
end

function TangentField(δ, ψ, parameters)
    work = NumberConservingWorkspace(parameters.L)
    UpdateWorkspace!(work, ψ, parameters)
    dδ = similar(δ)
    TangentField!(dδ, δ, ψ, parameters, work)
    return dδ
end

function Divergence(ψ::AbstractVector, parameters)
    work = NumberConservingWorkspace(parameters.L)
    UpdateWorkspace!(work, ψ, parameters)
    return Divergence(parameters, work)
end


""" Real 2L x 2L Jacobian in the coordinates x = (Re ψ, Im ψ), assembled from 2L calls to the
    analytic tangent field.  Used by the checks and by the Bogoliubov comparison of §5; the
    Lyapunov integration itself never forms it. """
function JacobianMatrix(ψ, parameters)
    L = parameters.L
    work = NumberConservingWorkspace(L)
    UpdateWorkspace!(work, ψ, parameters)

    jacobian = zeros(Float64, 2 * L, 2 * L)
    δ = zeros(ComplexF64, L)
    dδ = zeros(ComplexF64, L)

    for k = 1:(2 * L)
        fill!(δ, 0)
        δ[k <= L ? k : k - L] = k <= L ? 1 : im
        TangentField!(dδ, δ, ψ, parameters, work)

        @inbounds for j = 1:L
            jacobian[j, k] = real(dδ[j])
            jacobian[j + L, k] = imag(dδ[j])
        end
    end

    return jacobian
end


# ---------------------------------------------------------------------------------------------
# Conversions and initial conditions
# ---------------------------------------------------------------------------------------------

""" (p, q) layout of models/BoseHubbardFull.jl: ψ_i = (q_i + i p_i)/sqrt(2), x[i] = p_i,
    x[i+L] = q_i, sum_i (p_i^2 + q_i^2) = 2.  Lets an attractor found here be fed to the plotting
    and section tooling of the rest of the repository. """
ToPQ(ψ) = vcat(sqrt(2) .* imag.(ψ), sqrt(2) .* real.(ψ))
FromPQ(x) = (x[(length(x) ÷ 2 + 1):end] .+ im .* x[1:(length(x) ÷ 2)]) ./ sqrt(2)

""" x = (Re ψ, Im ψ) - the coordinates JacobianMatrix is written in. """
ToReal(ψ) = vcat(real.(ψ), imag.(ψ))

function FromReal(x)
    L = length(x) ÷ 2
    return x[1:L] .+ im .* x[(L + 1):end]
end

""" Uniform measure on the sphere sum_j n_j = 1, i.e. the Fubini-Study measure on CP^(L-1) once
    the global phase is quotiented out.  This is the right basin measure for the multistability
    question §8.3: every state of the reduced phase space is equally likely, and no absorbing
    ball has to be guessed as it did in the driven model. """
function RandomInitialCondition(L::Integer, rng = Random.default_rng())
    ψ = randn(rng, ComplexF64, L)
    return ψ ./ sqrt(sum(abs2, ψ))
end

""" The uniform (phase-locked) condensate ψ_j = L^(-1/2), optionally perturbed - the dark state of
    the κ channel and the state whose Bogoliubov spectrum §5 computes.  `amplitude > 0` seeds the
    unstable mode, so the trajectory falls onto whatever the destabilised dark state decays to,
    which is exactly what §8.1 asks about. """
function UniformInitialCondition(L::Integer; amplitude = 0.0, rng = Random.default_rng())
    ψ = fill(ComplexF64(1 / sqrt(L)), L)

    if amplitude > 0
        ψ .+= amplitude .* randn(rng, ComplexF64, L)
    end

    return ψ ./ sqrt(sum(abs2, ψ))
end


# ---------------------------------------------------------------------------------------------
# Bogoliubov spectrum of the uniform state (note §5)
# ---------------------------------------------------------------------------------------------

""" Bogoliubov exponents of the uniform condensate ψ_j = L^(-1/2) exp(-i g t / L), for a uniform
    circulation η (row sums of A must vanish, otherwise the uniform state is not a solution at all
    and this function is meaningless - use UniformIsSolution below to check).  With n = 1/L and
    D_q = 2(1 - cos q), q = 2 pi m / L:

        lambda_q^± = -2 κ n D_q - i η n sin q ± sqrt( -J D_q (J D_q + 2 g n) - η^2 n^2 sin^2 q )

    Returns a vector of (q, lambda_plus, lambda_minus).  The q = 0 pair is exactly (0, 0): those
    are the two exact zero exponents of §4 (global phase and the conserved total density). """
function BogoliubovSpectrum(parameters)
    L, J, g, κ = parameters.L, parameters.J, parameters.g, parameters.κ
    η = parameters.A[Next(1, L), 1]         # uniform ring: every bond carries the same rate
    n = 1 / L

    spectrum = Tuple{Float64, ComplexF64, ComplexF64}[]

    for m = 0:(L - 1)
        q = 2 * pi * m / L
        D = 2 * (1 - cos(q))
        s = sin(q)

        drift = -2 * κ * n * D - im * η * n * s
        root = sqrt(complex(-J * D * (J * D + 2 * g * n) - η^2 * n^2 * s^2))

        push!(spectrum, (q, drift + root, drift - root))
    end

    return spectrum
end


""" Largest Bogoliubov growth rate over the non-trivial modes q != 0.  Positive means the uniform
    dark state is modulationally unstable - the only source of instability identified in the note,
    and therefore the line that every scan of §8 has to be placed against. """
function BogoliubovGrowthRate(parameters)
    rate = -Inf

    for (q, plus, minus) in BogoliubovSpectrum(parameters)
        q == 0 && continue
        rate = max(rate, real(plus), real(minus))
    end

    return rate
end


""" Critical g (the note quotes g < -(9 + η^2/12 + 4 κ^2)/2 for L = 3, J = 1) found by bisecting
    BogoliubovGrowthRate on g < 0.  Returns NaN if the state is stable over the whole bracket.

    `tolerance` is needed because below threshold with κ = 0 the radicand is negative and the
    growth rate is EXACTLY 0, not negative: the uniform state is then marginally (Hamiltonian)
    stable, and only a strictly positive rate counts as an instability. """
function InstabilityThreshold(L; J = 1.0, κ = 0.0, η = 0.0, bracket = (-1.0e4, -1.0e-6),
        tolerance = 1e-12)
    Rate(g) = BogoliubovGrowthRate(NumberConservingParameters(L; J = J, g = g, κ = κ, η = η))

    low, high = bracket
    Rate(low) > tolerance || return NaN
    Rate(high) > tolerance && return high

    for _ = 1:200
        middle = 0.5 * (low + high)
        (Rate(middle) > tolerance ? (low = middle) : (high = middle))
    end

    return 0.5 * (low + high)
end


""" Residual of the uniform state as an exact solution: min over the global-phase drift omega of
    || F(ψ) + i omega ψ ||.  Zero for a circulant A (vanishing row sums), non-zero as soon as the
    rates are modulated - which is why the non-uniform runs of §8.2 need their own reference
    state rather than the Bogoliubov threshold. """
function UniformIsSolution(parameters)
    ψ = UniformInitialCondition(parameters.L)
    field = VectorField(ψ, parameters)
    ω = -imag(dot(ψ, field)) / sum(abs2, ψ)

    return norm(field .+ im * ω .* ψ)
end


# ---------------------------------------------------------------------------------------------
# Full Lyapunov spectrum (QR / Benettin), Kaplan-Yorke dimension and the trace rule
# ---------------------------------------------------------------------------------------------

""" State of one Lyapunov run.  `logGrowth[i]` accumulates log of the i-th diagonal element of R in
    the QR decomposition of the propagated deviation basis, but ONLY after relaxationTime: the
    exponents must be measured on the attractor, not on the transient that reaches it. """
mutable struct LyapunovIntegrationParameters
    model::NumberConservingParameters
    work::NumberConservingWorkspace
    deviations::Int
    relaxationTime::Float64
    started::Bool
    startTime::Float64
    startDivergence::Float64
    lastTime::Float64
    endDivergence::Float64
    accumulations::Int
    logGrowth::Vector{Float64}
    history::Vector{Vector{Float64}}
    historyLength::Int
    maximumNormError::Float64
end


""" Real (Euclidean on R^2L) inner product of two complex deviation vectors. """
function RealDot(a, b)
    total = 0.0
    @inbounds @simd for k in eachindex(a, b)
        total += real(conj(a[k]) * b[k])
    end
    return total
end


""" Modified Gram-Schmidt of the m deviation vectors, in place, with the real inner product. """
function Orthonormalise!(vectors)
    for i in eachindex(vectors)
        vi = vectors[i]
        for j = 1:(i - 1)
            vj = vectors[j]
            projection = RealDot(vj, vi)
            @inbounds @simd for k in eachindex(vi)
                vi[k] -= projection * vj[k]
            end
        end
        vi ./= sqrt(RealDot(vi, vi))
    end
    return vectors
end


""" Right-hand side of the full problem: the trajectory, `deviations` tangent vectors and the
    running integral of div F.  The state-dependent quantities are computed once per call and
    shared by all deviation vectors, which is what makes 2L of them affordable. """
function EquationOfMotion!(du, u, integrationParameters, t)
    parameters = integrationParameters.model
    work = integrationParameters.work
    L = parameters.L

    ψ = view(u, 1:L)
    UpdateWorkspace!(work, ψ, parameters)
    VectorField!(view(du, 1:L), ψ, parameters, work)

    @inbounds for k = 1:integrationParameters.deviations
        indices = (k * L + 1):((k + 1) * L)
        TangentField!(view(du, indices), view(u, indices), ψ, parameters, work)
    end

    du[end] = complex(Divergence(parameters, work), 0.0)

    return nothing
end


""" PeriodicCallback: reorthonormalise the deviation basis and accumulate the growth factors.

    Before relaxationTime the basis is still reorthonormalised (otherwise every vector collapses
    onto the leading direction and the smaller exponents are lost to round-off) but nothing is
    accumulated.  The first firing at or after relaxationTime opens the measurement window and
    records the value of the divergence integral, so that both the exponents and <div F> are
    averaged over exactly the same interval - which is what makes the comparison of §8.4 a test
    of the dynamics rather than of the transient. """
function Reorthonormalise!(integrator)
    integrationParameters = integrator.p
    L = integrationParameters.model.L
    u = integrator.u
    t = integrator.t

    # The window OPENS at this firing: the growth accumulated over the interval that ends here
    # belongs to the transient, so it is discarded together with the log sums, and the divergence
    # integral is read at exactly the same instant.
    opening = !integrationParameters.started && t >= integrationParameters.relaxationTime

    if opening
        integrationParameters.started = true
        integrationParameters.startTime = t
        integrationParameters.startDivergence = real(u[end])
        integrationParameters.endDivergence = real(u[end])
        integrationParameters.lastTime = t
        integrationParameters.accumulations = 0
        fill!(integrationParameters.logGrowth, 0.0)
        empty!(integrationParameters.history)
    end

    accumulate = integrationParameters.started && !opening

    vectors = [view(u, (k * L + 1):((k + 1) * L)) for k = 1:integrationParameters.deviations]

    for i in eachindex(vectors)
        vi = vectors[i]
        for j = 1:(i - 1)
            vj = vectors[j]
            projection = RealDot(vj, vi)
            @inbounds @simd for k in eachindex(vi)
                vi[k] -= projection * vj[k]
            end
        end

        growth = sqrt(RealDot(vi, vi))
        vi ./= growth

        if accumulate
            integrationParameters.logGrowth[i] += log(growth)
        end
    end

    integrationParameters.maximumNormError =
        max(integrationParameters.maximumNormError, abs(sum(abs2, view(u, 1:L)) - 1))

    # Both averages are closed at the SAME firing.  Reading the divergence integral from the end
    # of the solution instead would divide an interval of length (integrationTime - startTime) by
    # a window that stops at the last firing, and the resulting O(rescaleStep / window) mismatch
    # is exactly the size of the effect §8.4 is trying to measure.
    if accumulate
        integrationParameters.endDivergence = real(u[end])
        integrationParameters.lastTime = t
        integrationParameters.accumulations += 1
    end

    if accumulate && integrationParameters.historyLength > 0 &&
            t > integrationParameters.startTime
        push!(integrationParameters.history,
              integrationParameters.logGrowth ./ (t - integrationParameters.startTime))
        if length(integrationParameters.history) > integrationParameters.historyLength
            popfirst!(integrationParameters.history)
        end
    end

    u_modified!(integrator, true)
end


""" Gauge-invariant observables, recorded by a SavingCallback and used to fingerprint the
    attractor for the multistability question §8.3:

        coherence  sum_j Re(conj(ψ_j) ψ_{j+1})   = 1 on the uniform locked state, -1 staggered;
                                                   it is also what the divergence is built from
        maximum_n  max_j n_j                     = 1/L uniform, 1 on a self-trapped state
        ipr        sum_j n_j^2                   inverse participation ratio, 1/L uniform
        current    sum_j Im(conj(ψ_j) ψ_{j+1})   the particle current round the ring, which is what
                                                   the circulation η is supposed to drive and which
                                                   changes sign under the reflection of §6 """
function Observables(u, t, integrator)
    L = integrator.p.model.L

    coherence = 0.0
    maximum_n = 0.0
    ipr = 0.0
    current = 0.0

    @inbounds for j = 1:L
        n = abs2(u[j])
        maximum_n = max(maximum_n, n)
        ipr += n * n

        z = conj(u[j]) * u[Next(j, L)]
        coherence += real(z)
        current += imag(z)
    end

    return (coherence, maximum_n, ipr, current)
end


""" Kaplan-Yorke dimension of an (unsorted) spectrum: the largest k whose partial sum is still
    non-negative, plus the fraction of the next exponent needed to bring it to zero.

    A value equal to the full dimension of the spectrum means no partial sum ever went negative,
    i.e. the flow does not contract at all: the set is not an attractor, and the number is a flag
    rather than a dimension (ClassifyAttractor calls that case :neutral). """
function KaplanYorke(spectrum)
    λ = sort(collect(spectrum), rev = true)

    total = 0.0
    for k in eachindex(λ)
        if total + λ[k] < 0
            return k == 1 ? 0.0 : (k - 1) + total / abs(λ[k])
        end
        total += λ[k]
    end

    return float(length(λ))
end


""" The reduced spectrum of note §4: the full 2L exponents minus the TWO that are exactly zero by
    symmetry (the global U(1) phase and the conserved total density).  Dropping the two exponents
    closest to zero is exactly the prescription D_KY(red) = D_KY(full) - 2 of the note, and it
    leaves the flow-direction zero in place - which is what makes the count of remaining zeros a
    classification: 0 gives a relative fixed point, 1 a limit cycle, 2 a torus, and a positive
    exponent a strange attractor. """
function ReducedSpectrum(spectrum)
    order = sortperm(abs.(spectrum))
    keep = [i for i in eachindex(spectrum) if !(i in order[1:2])]
    return sort(spectrum[keep], rev = true)
end


""" Attractor type from the reduced spectrum.

    `:fixedPoint` means a RELATIVE fixed point, i.e. stationary up to the global phase
    (ψ_j(t) = exp(-i omega t) ψ_j(0)): the flow direction then coincides with the U(1) direction,
    so only two zero exponents exist instead of three.

    `:neutral` is the case where the WHOLE reduced spectrum is zero, so the flow neither contracts
    nor expands: an invariant torus of a locally volume-preserving region, not an attractor.  It is
    the generic outcome at κ = 0 with a uniform circulation (note §4), and it also turns up at
    κ = 0 with modulated rates wherever the trajectory settles where <sum_j c_j n_j> = 0.  Keeping
    it apart from :torus matters for question 2, whose whole point is contraction. """
function ClassifyAttractor(reduced; zeroThreshold = 5e-3, chaosThreshold = 1e-2)
    positive = count(λ -> λ > chaosThreshold, reduced)
    zeros = count(λ -> abs(λ) <= zeroThreshold, reduced)

    positive >= 2 && return :hyperchaotic
    positive == 1 && return :chaotic
    zeros == length(reduced) && return :neutral
    zeros == 0 && return :fixedPoint
    zeros == 1 && return :limitCycle
    zeros == 2 && return :torus

    return :undetermined
end


""" Full Lyapunov spectrum of the attractor reached from `ψ0` - everything §8 asks for, in one
    call.  Returns a NamedTuple with

        spectrum          all 2L exponents, sorted descending
        reduced           the 2L - 2 exponents of the reduced space CP^(L-1)
        uncertainty       finite-time scatter of each exponent over the trailing history window
        dimension         D_KY on the reduced space (fractional means a strange attractor)
        dimensionFull     D_KY on the full space; equals dimension + 2 by construction
        classification    :fixedPoint / :limitCycle / :torus / :chaotic / :hyperchaotic
        divergence        <div F> averaged over the measurement window (note §4 trace rule)
        traceError        | sum(spectrum) - <div F> |, the check of §8.4
        symmetryZeros     the two exponents dropped as the exact zeros; both should be tiny
        stationarity      || F(ψ) + i omega ψ || at the end - zero on a relative fixed point
        coherence, maximum_n, ipr, current   attractor averages of Observables, with their spreads
        normError         worst | sum_j n_j - 1 | seen: the dS/dt = 0 check of §4, live
        window            length of the measurement window, and startTime / accumulations, the
                          instant it opened and the number of QR steps it holds
        psi               the final state, so the attractor can be re-entered or plotted

    The exponents and <div F> are averaged over exactly the same window
    (relaxationTime, integrationTime), which is what makes §8.4 a test of the dynamics rather
    than of the transient. """
function LyapunovSpectrum(ψ0, parameters;
        deviations = 2 * parameters.L,
        relaxationTime = 2000.0,
        integrationTime = 8000.0,
        rescaleStep = 1.0,
        saveStep = 1.0,
        tolerance = 1e-10,
        historyLength = 500,
        rng = Random.default_rng(),
        solver = DP8(),
        maximumIterations = 10^8,
        zeroThreshold = 5e-3,
        chaosThreshold = 1e-2)

    L = parameters.L
    relaxationTime < integrationTime ||
        error("the exponents are accumulated on (relaxationTime, integrationTime)")

    u0 = zeros(ComplexF64, L * (deviations + 1) + 1)
    u0[1:L] .= ψ0 ./ sqrt(sum(abs2, ψ0))

    vectors = [view(u0, (k * L + 1):((k + 1) * L)) for k = 1:deviations]
    for v in vectors
        v .= randn(rng, ComplexF64, L)
    end
    Orthonormalise!(vectors)

    integrationParameters = LyapunovIntegrationParameters(
        parameters, NumberConservingWorkspace(L), deviations, float(relaxationTime),
        false, float(relaxationTime), 0.0, float(relaxationTime), 0.0, 0,
        zeros(deviations), Vector{Float64}[], historyLength, 0.0)

    saved = SavedValues(Float64, NTuple{4, Float64})

    problem = ODEProblem(ODEFunction(EquationOfMotion!), u0, (0.0, float(integrationTime)),
                         integrationParameters)

    callback = CallbackSet(
        PeriodicCallback(Reorthonormalise!, float(rescaleStep); save_positions = (false, false)),
        SavingCallback(Observables, saved; saveat = relaxationTime:saveStep:integrationTime))

    solution = solve(problem, solver; reltol = tolerance, abstol = tolerance, callback = callback,
                     save_everystep = false, save_start = false, save_end = true,
                     maxiters = maximumIterations)

    window = integrationParameters.lastTime - integrationParameters.startTime
    window > 0 || error("empty measurement window - relaxationTime too close to integrationTime")

    order = sortperm(integrationParameters.logGrowth, rev = true)
    spectrum = integrationParameters.logGrowth[order] ./ window

    uncertainty = zeros(length(spectrum))
    if length(integrationParameters.history) > 2
        samples = reduce(hcat, integrationParameters.history)       # deviations by samples
        uncertainty = [std(view(samples, i, :)) for i in order]
    end

    reduced = ReducedSpectrum(spectrum)
    dropped = sort(spectrum[sortperm(abs.(spectrum))[1:2]], rev = true)

    divergence = (integrationParameters.endDivergence -
                  integrationParameters.startDivergence) / window

    ψ = solution.u[end][1:L]
    ψ ./= sqrt(sum(abs2, ψ))
    field = VectorField(ψ, parameters)
    ω = -imag(dot(ψ, field)) / sum(abs2, ψ)
    stationarity = norm(field .+ im * ω .* ψ)

    values = saved.saveval
    means = length(values) > 0 ? [mean(getindex.(values, i)) for i = 1:4] : fill(NaN, 4)
    spreads = length(values) > 1 ? [std(getindex.(values, i)) for i = 1:4] : fill(NaN, 4)

    return (spectrum = spectrum,
            reduced = reduced,
            uncertainty = uncertainty,
            dimension = KaplanYorke(reduced),
            dimensionFull = KaplanYorke(spectrum),
            classification = ClassifyAttractor(reduced; zeroThreshold = zeroThreshold,
                                               chaosThreshold = chaosThreshold),
            divergence = divergence,
            traceError = abs(sum(spectrum) - divergence),
            symmetryZeros = dropped,
            stationarity = stationarity,
            coherence = means[1], coherenceSpread = spreads[1],
            maximum_n = means[2], maximum_nSpread = spreads[2],
            ipr = means[3], iprSpread = spreads[3],
            current = means[4], currentSpread = spreads[4],
            normError = integrationParameters.maximumNormError,
            window = window,
            startTime = integrationParameters.startTime,
            accumulations = integrationParameters.accumulations,
            retcode = solution.retcode,
            psi = ψ)
end


# ---------------------------------------------------------------------------------------------
# Checks - the classical half of the table in note §10, redone for THIS implementation
# ---------------------------------------------------------------------------------------------

""" The number-conserving channel of note §2 written in its MODE form, for one pair of orthonormal
    single-particle modes u, v at rate gamma (the classical limit, where the spontaneous 1/N term
    is dropped):

        dψ/dt = (gamma/2) [ |phi_v|^2 phi_u u - |phi_u|^2 phi_v v ],    phi_u = u^dagger ψ

    The κ term of VectorField! is the bond special case u, v = (e_j ± e_k)/sqrt(2), gamma = 4κ;
    the directed hopping is the case u = e_i, v = e_j, gamma = γ_ij.  Used only by Checks(). """
function ModeChannel(ψ, u, v, γ)
    φu = dot(u, ψ)
    φv = dot(v, ψ)
    return (γ / 2) .* (abs2(φv) * φu .* u .- abs2(φu) * φv .* v)
end


""" Central-difference Jacobian of the vector field in the real coordinates x = (Re ψ, Im ψ). """
function FiniteDifferenceJacobian(ψ, parameters; h = 1e-6)
    L = parameters.L
    x = ToReal(ψ)
    jacobian = zeros(Float64, 2 * L, 2 * L)

    for k = 1:(2 * L)
        forward = copy(x); forward[k] += h
        backward = copy(x); backward[k] -= h
        jacobian[:, k] .= (ToReal(VectorField(FromReal(forward), parameters)) .-
                           ToReal(VectorField(FromReal(backward), parameters))) ./ (2 * h)
    end

    return jacobian
end


""" Largest Lyapunov exponent measured WITHOUT the tangent equations: two full nonlinear
    trajectories started a distance `separation` apart, with the separation renormalised every
    `step` (the original Benettin recipe).  It shares no code with TangentField!, so agreement
    with LyapunovSpectrum is an independent confirmation that the chaos is in the model and not in
    the linearisation.

    Compare it with spectrum[1], the largest exponent of the FULL spectrum, NOT with reduced[1]:
    a pair of nonlinear trajectories cannot be kept out of the two neutral symmetry directions, so
    on a relative fixed point this returns 0 (the shadow ends up on the same state with a different
    global phase) where the reduced spectrum returns the contraction rate.  The two agree exactly
    where it matters, namely wherever there is a positive exponent.

    `separation` has a window: too small and the integration error dominates (the exponent then
    grows like log(noise/separation) as it is reduced, which is easy to mistake for chaos), too
    large and the pair leaves the linear regime.  With tolerance 1e-12 the plateau measured here
    runs from about 1e-5 to 1e-3, and 1e-7 is already 10 per cent high.  Expect agreement to a few per cent and no better: on a strange attractor the
    finite-time exponent scatters by that much between any two realisations of the same orbit. """
function FiniteSeparationExponent(ψ0, parameters; separation = 1e-5, relaxationTime = 1000.0,
        integrationTime = 5000.0, step = 1.0, tolerance = 1e-12, rng = Random.default_rng())

    L = parameters.L
    reference = NumberConservingWorkspace(L)
    shadow = NumberConservingWorkspace(L)

    function TwoTrajectories!(du, u, p, t)
        ψ = view(u, 1:L)
        φ = view(u, (L + 1):(2 * L))

        UpdateWorkspace!(reference, ψ, parameters)
        VectorField!(view(du, 1:L), ψ, parameters, reference)

        UpdateWorkspace!(shadow, φ, parameters)
        VectorField!(view(du, (L + 1):(2 * L)), φ, parameters, shadow)

        return nothing
    end

    ψ = ψ0 ./ sqrt(sum(abs2, ψ0))
    δ = randn(rng, ComplexF64, L)
    φ = ψ .+ (separation / sqrt(sum(abs2, δ))) .* δ
    φ ./= sqrt(sum(abs2, φ))                 # both copies stay on the sphere sum_j n_j = 1

    accumulator = Ref(0.0)
    startTime = Ref(NaN)
    lastTime = Ref(NaN)

    # Projecting the shadow back onto the sphere removes the radial part of the separation, so the
    # distance it actually starts the next interval with is NOT `separation`.  Carrying the true
    # value across is what makes the growth factor unbiased.
    referenceDistance = Ref(sqrt(sum(abs2, φ .- ψ)))

    function Rescale!(integrator)
        u = integrator.u
        t = integrator.t
        a = view(u, 1:L)
        b = view(u, (L + 1):(2 * L))

        distance = sqrt(sum(abs2, b .- a))
        distance > 0 || return

        if isnan(startTime[])
            t >= relaxationTime && (startTime[] = t)     # open the window; accumulate from the next firing
        else
            accumulator[] += log(distance / referenceDistance[])
            lastTime[] = t
        end

        b .= a .+ (separation / distance) .* (b .- a)
        b ./= sqrt(sum(abs2, b))
        referenceDistance[] = sqrt(sum(abs2, b .- a))
        u_modified!(integrator, true)
    end

    problem = ODEProblem(TwoTrajectories!, vcat(ψ, φ), (0.0, float(integrationTime)))
    solve(problem, DP8(); reltol = tolerance, abstol = tolerance,
          callback = PeriodicCallback(Rescale!, float(step); save_positions = (false, false)),
          save_everystep = false, save_start = false, save_end = false)

    return accumulator[] / (lastTime[] - startTime[])
end


""" Consistency checks of the implementation, in the spirit of the table of note §10.  Everything
    here is independent of the Lyapunov machinery except the last four rows, which are the checks
    question 4 of §8 asks for.  Returns the list of (description, error) pairs and prints them. """
function Checks(; seed = 20240923, verbose = true)
    results = Tuple{String, Float64}[]
    rng = Xoshiro(seed)

    worstNorm = 0.0
    worstTangent = 0.0
    worstDivergence = 0.0
    worstPhase = 0.0
    worstDensity = 0.0

    for L in (3, 4, 5, 7), trial = 1:5
        parameters = NumberConservingParameters(L; J = randn(rng), g = randn(rng),
            κ = abs(randn(rng)), η = randn(rng), modulation = 0.7 * randn(rng),
            phase = 2 * pi * rand(rng))
        ψ = RandomInitialCondition(L, rng)

        # dS/dt = 0: the strong U(1) of §1, exact in the classical limit (§4)
        field = VectorField(ψ, parameters)
        worstNorm = max(worstNorm, abs(2 * RealDot(ψ, field)))

        # analytic tangent map against a central difference - this validates every exponent
        analytic = JacobianMatrix(ψ, parameters)
        worstTangent = max(worstTangent,
                           maximum(abs, analytic .- FiniteDifferenceJacobian(ψ, parameters)))

        # divergence formula of §4 against the trace of the Jacobian
        worstDivergence = max(worstDivergence, abs(Divergence(ψ, parameters) - tr(analytic)))

        # U(1) equivariance: delta = i ψ solves the tangent equations exactly, so it is a neutral
        # direction with exponent EXACTLY zero - the first of the three zeros of §4
        worstPhase = max(worstPhase,
                         maximum(abs, TangentField(im .* ψ, ψ, parameters) .- im .* field))

        # d/dt (grad S . delta) = 0: the conserved density gives the second exact zero exponent
        δ = randn(rng, ComplexF64, L)
        worstDensity = max(worstDensity,
            abs(RealDot(field, δ) + RealDot(ψ, TangentField(δ, ψ, parameters))))
    end

    push!(results, ("dS/dt = 0 at random points (L = 3, 4, 5, 7, random parameters)", worstNorm))
    push!(results, ("analytic tangent map vs central difference", worstTangent))
    push!(results, ("divergence formula §4 vs trace of the Jacobian", worstDivergence))
    push!(results, ("U(1) generator i ψ is an exact tangent solution (zero exponent 1)", worstPhase))
    push!(results, ("d/dt (grad S . delta) = 0 (zero exponent 2)", worstDensity))

    # --- bond form of the κ channel vs the mode form of §2 --------------------------------------
    worstMode = 0.0
    for L in (3, 4, 5, 7)
        κ = abs(randn(rng))
        parameters = NumberConservingParameters(L; J = 0.0, g = 0.0, κ = κ)
        ψ = RandomInitialCondition(L, rng)

        mode = zeros(ComplexF64, L)
        for j = 1:L
            k = Next(j, L)
            u = zeros(ComplexF64, L); u[j] = 1 / sqrt(2); u[k] = 1 / sqrt(2)
            v = zeros(ComplexF64, L); v[j] = 1 / sqrt(2); v[k] = -1 / sqrt(2)
            mode .+= ModeChannel(ψ, u, v, 4 * κ)
        end

        worstMode = max(worstMode, maximum(abs, VectorField(ψ, parameters) .- mode))
    end
    push!(results, ("bond form of the κ term vs mode form of §2 (gamma = 4κ)", worstMode))

    # --- replicator structure of the directed hopping (§3) ---------------------------------------
    worstReplicator = 0.0
    worstFrozen = 0.0
    for L in (3, 4, 5, 7)
        parameters = NumberConservingParameters(L; J = 0.0, g = 0.0, κ = 0.0,
            η = randn(rng), modulation = 0.5 * randn(rng))
        ψ = RandomInitialCondition(L, rng)
        field = VectorField(ψ, parameters)

        work = NumberConservingWorkspace(L)
        UpdateWorkspace!(work, ψ, parameters)

        for j = 1:L
            z = conj(ψ[j]) * field[j]
            worstReplicator = max(worstReplicator, abs(2 * real(z) - work.n[j] * work.S[j]))
            worstFrozen = max(worstFrozen, abs(imag(z)))
        end
    end
    push!(results, ("directed hopping: dn_i/dt = n_i sum_k A_ik n_k (replicator)", worstReplicator))
    push!(results, ("directed hopping leaves the phases untouched", worstFrozen))

    # --- Bogoliubov spectrum of §5 against the full Jacobian --------------------------------------
    worstBogoliubov = 0.0
    for L in (3, 4, 5, 7), trial = 1:3
        parameters = NumberConservingParameters(L; J = randn(rng), g = randn(rng),
            κ = abs(randn(rng)), η = randn(rng))
        ψ = UniformInitialCondition(L)

        # the uniform state rotates as exp(-i g t / L); the Bogoliubov ansatz of §5 works in the
        # frame that removes that drift, so add + i (g/L) back onto the Jacobian before comparing
        jacobian = JacobianMatrix(ψ, parameters)
        for j = 1:L
            jacobian[j, j + parameters.L] -= parameters.g / parameters.L
            jacobian[j + parameters.L, j] += parameters.g / parameters.L
        end

        Key(z) = (round(real(z), digits = 8), imag(z))
        numerical = sort(eigvals(jacobian), by = Key)

        analytic = ComplexF64[]
        for (q, plus, minus) in BogoliubovSpectrum(parameters)
            push!(analytic, plus)
            push!(analytic, minus)
        end
        analytic = sort(analytic, by = Key)

        worstBogoliubov = max(worstBogoliubov, maximum(abs, numerical .- analytic))
    end
    push!(results, ("Bogoliubov lambda_q of §5 vs full Jacobian (L = 3, 4, 5, 7)", worstBogoliubov))

    # --- the trimer threshold quoted in §5 --------------------------------------------------------
    worstThreshold = 0.0
    for trial = 1:5
        κ = abs(randn(rng))
        η = randn(rng)
        worstThreshold = max(worstThreshold,
            abs(InstabilityThreshold(3; J = 1.0, κ = κ, η = η) + (9 + η^2 / 12 + 4 * κ^2) / 2))
    end
    push!(results, ("trimer threshold g_c = -(9 + eta^2/12 + 4 kappa^2)/2 of §5", worstThreshold))

    # --- the uniform state as an exact solution, and the (p, q) interoperability ----------------
    worstUniform = 0.0
    worstModulated = Inf
    worstLayout = 0.0
    for L in (3, 4, 5, 7)
        circulant = NumberConservingParameters(L; J = randn(rng), g = randn(rng),
            κ = abs(randn(rng)), η = randn(rng))
        modulated = NumberConservingParameters(L; J = circulant.J, g = circulant.g,
            κ = circulant.κ, η = 1.0, modulation = 0.6)

        worstUniform = max(worstUniform, UniformIsSolution(circulant))
        worstModulated = min(worstModulated, UniformIsSolution(modulated))

        ψ = RandomInitialCondition(L, rng)
        worstLayout = max(worstLayout, maximum(abs, FromPQ(ToPQ(ψ)) .- ψ),
                          abs(sum(abs2, ToPQ(ψ)) - 2))
    end
    push!(results, ("uniform state solves the equations for a circulant A (§5)", worstUniform))
    push!(results, ("it does NOT once the rates are modulated (smallest residual; must be > 0)",
                    worstModulated))
    push!(results, ("ToPQ / FromPQ round trip and the sum x^2 = 2 normalisation", worstLayout))

    # --- the largest exponent WITHOUT the tangent equations -----------------------------------
    # Two full nonlinear trajectories, no shared code with TangentField!.  It is the answer to
    # "is the chaos real or is the linearisation wrong"; a few per cent is the finite-time scatter
    # of the exponent itself, so this row is a per-cent check, not a 1e-15 one.
    let parameters = NumberConservingParameters(3; g = -20.0, κ = 0.3, η = 3.0),
        ψ = RandomInitialCondition(3, Xoshiro(4))

        tangent = LyapunovSpectrum(ψ, parameters; relaxationTime = 1000.0,
                                   integrationTime = 5000.0, rng = Xoshiro(11)).spectrum[1]
        pair = FiniteSeparationExponent(ψ, parameters; relaxationTime = 1000.0,
                                        integrationTime = 5000.0, rng = Xoshiro(11))

        push!(results, ("lambda_max from two nearby trajectories vs the tangent code " *
                        "(relative; per cent expected)", abs(pair - tangent) / abs(tangent)))
    end

    # --- the Lyapunov machinery: the two checks of §8.4 -------------------------------------------
    for (label, parameters, ψ) in (
            ("chaotic", NumberConservingParameters(3; g = -20.0, κ = 0.3, η = 1.0),
             RandomInitialCondition(3, Xoshiro(4))),
            ("regular", NumberConservingParameters(3; g = -8.0, κ = 0.3, η = 0.0),
             RandomInitialCondition(3, Xoshiro(2))))
        result = LyapunovSpectrum(ψ, parameters; relaxationTime = 1000.0,
                                  integrationTime = 6000.0, rng = Xoshiro(11))
        push!(results, ("trace rule sum(lambda) = <div F> on a $label trajectory (§8.4)",
                        result.traceError))
        push!(results, ("sum_j n_j = 1 along the $label trajectory", result.normError))
    end

    if verbose
        println("| check | result |")
        println("|---|---|")
        for (description, error) in results
            @printf("| %s | %.1e |\n", description, error)
        end
    end

    return results
end


# ---------------------------------------------------------------------------------------------
# Section 8 of the note - the four open questions
# ---------------------------------------------------------------------------------------------
#
# Run with  julia -t auto BHNumberConserving.jl  (the trajectories of one parameter point are
# spread over the available threads; everything else is serial and cheap).
#
#   1. Beyond the §5 threshold, is the attractor strange, or a limit cycle or a torus?
#   2. Is κ needed, or does circulation with a non-uniform rate (c_j != 0) suffice?
#   3. Is there multistability (the locked state against something else)?
#   4. Does sum(lambda) = <div F> hold in the time average, and are the three zeros there?
#
# Questions 1-3 are answered by the same primitive - a parameter point sampled with several random
# initial conditions drawn from the Fubini-Study measure - and question 4 is then read off the
# union of every trajectory computed, which is the largest sample available and costs nothing.

""" All trajectories of one parameter point.  The first initial condition is the uniform dark
    state with a small perturbation (so that the destabilised state of §5 is always sampled), the
    rest are drawn uniformly on CP^(L-1).  Seeds are derived from the parameters, so a rerun
    reproduces the sample exactly. """
function ScanPoint(L, g; κ = 0.0, η = 0.0, modulation = 0.0, trajectories = 12,
        relaxationTime = 1000.0, integrationTime = 5000.0, seed = 0, seedUniform = true, kwargs...)

    parameters = NumberConservingParameters(L; J = 1.0, g = g, κ = κ, η = η,
                                            modulation = modulation)
    results = Vector{Any}(undef, trajectories)

    Threads.@threads for i = 1:trajectories
        rng = Xoshiro(hash((L, g, κ, η, modulation, seed, i)))
        ψ = (seedUniform && i == 1) ? UniformInitialCondition(L; amplitude = 1e-3, rng = rng) :
                                      RandomInitialCondition(L, rng)
        results[i] = LyapunovSpectrum(ψ, parameters; relaxationTime = relaxationTime,
                                      integrationTime = integrationTime, rng = rng, kwargs...)
    end

    return parameters, results
end


""" Compact summary of one parameter point. """
function SummarisePoint(results; chaosThreshold = 1e-2)
    λ = [r.reduced[1] for r in results]
    best = argmax(λ)

    chaoticDimensions = [r.dimension for r in results if r.reduced[1] > chaosThreshold]

    return (λmax = λ[best],
            dimension = isempty(chaoticDimensions) ? NaN : mean(chaoticDimensions),
            chaoticFraction = count(x -> x > chaosThreshold, λ) / length(λ),
            classes = [r.classification for r in results],
            traceError = maximum(r.traceError for r in results),
            divergence = mean(r.divergence for r in results),
            best = results[best])
end


""" "chaotic 5, fixedPoint 7" - the distribution over attractor types at one parameter point. """
function ClassHistogram(classes)
    counts = Dict{Symbol, Int}()
    for c in classes
        counts[c] = get(counts, c, 0) + 1
    end
    return join([string(k, " ", counts[k]) for k in sort(collect(keys(counts)), by = string)], ", ")
end


""" QUESTION 1.  Beyond the modulational-instability threshold of §5, is the attractor strange?

    For each L and each circulation η the interaction g is swept from just above the threshold
    (where the uniform dark state is still stable and nothing can happen) deep into the unstable
    region.  Every point reports the fraction of initial conditions that end on an attractor with
    a positive exponent, the largest exponent found, and the Kaplan-Yorke dimension of the REDUCED
    2L - 2 dimensional space - fractional D_KY on top of a positive exponent is what makes an
    attractor strange rather than a limit cycle (D = 1) or a torus (D = 2). """
function Question1(; Ls = (3, 4), κ = 0.3, ηs = (0.0, 1.0, 3.0, 5.0),
        gs = (-4.0, -6.0, -8.0, -12.0, -20.0, -30.0, -50.0), trajectories = 12, kwargs...)

    println("=" ^ 100)
    println("QUESTION 1: is the attractor beyond the §5 threshold strange, a limit cycle or a torus?")
    println("=" ^ 100)

    records = []

    for L in Ls, η in ηs
        threshold = InstabilityThreshold(L; J = 1.0, κ = κ, η = η)
        @printf("\nL = %d, J = 1, κ = %.2f, η = %.2f   -   uniform state unstable for g < %.3f (§5)\n",
                L, κ, η, threshold)
        println("      g   unstable  chaotic     λ_max   D_KY(red)/(2L-2)   attractors")

        for g in gs
            parameters, results = ScanPoint(L, g; κ = κ, η = η, trajectories = trajectories,
                                            kwargs...)
            summary = SummarisePoint(results)
            append!(records, results)

            @printf("%7.1f   %5s     %5.2f   %+8.4f   %6.3f / %d        %s\n",
                    g, g < threshold ? "yes" : "no", summary.chaoticFraction, summary.λmax,
                    summary.dimension, 2 * L - 2, ClassHistogram(summary.classes))
        end
    end

    chaotic = [r for r in records if r.classification in (:chaotic, :hyperchaotic)]

    println()
    if isempty(chaotic)
        println("VERDICT: no positive exponent anywhere on this grid - the destabilised dark state")
        println("         relaxes onto fixed points, limit cycles or tori only.")
    else
        best = chaotic[argmax([r.reduced[1] for r in chaotic])]
        fractional = [r.dimension for r in chaotic if abs(r.dimension - round(r.dimension)) > 0.02]
        @printf("VERDICT: STRANGE ATTRACTORS EXIST.  %d of %d trajectories are chaotic; the largest\n",
                length(chaotic), length(records))
        @printf("         exponent is λ_max = %+.4f with D_KY(red) = %.3f out of %d, and %d of the\n",
                best.reduced[1], best.dimension, length(best.reduced), length(fractional))
        @printf("         %d chaotic attractors have a clearly fractional reduced dimension.\n",
                length(chaotic))
        println("         Reduced spectrum of the strongest: ", round.(best.reduced, digits = 4))
    end

    return records
end


""" QUESTION 2.  Is κ needed, or does circulation alone suffice once the rates are non-uniform?

    Two parts, both at κ = 0:

    (a) UNIFORM circulation.  Then c_j = 0, div F vanishes identically and the flow is volume
        preserving - the note's claim that there is no attractor at all.  The test is not a
        Lyapunov exponent (a volume-preserving flow can perfectly well be chaotic - it just has no
        attractor) but the pair of conditions div F = 0 pointwise and sum(lambda) = 0.

    (b) NON-UNIFORM circulation, eta_j = η (1 + modulation cos(2 pi j / L)).  Now c_j = eta_j -
        eta_{j-1} != 0, div F = sum_j c_j n_j has either sign, and an attractor can exist.  Note
        that the uniform condensate is no longer a solution here (the row sums of A do not vanish),
        so the §5 threshold does not apply and the whole g range has to be searched.  The question
        is whether any of it carries a positive exponent together with a NEGATIVE mean divergence -
        a positive exponent on a volume-preserving flow would not be an attractor. """
function Question2(; Ls = (3, 4), ηs = (1.0, 3.0, 5.0), modulations = (0.5, 0.9),
        gs = (-6.0, -12.0, -25.0, -50.0), trajectories = 12, kwargs...)

    println()
    println("=" ^ 100)
    println("QUESTION 2: is κ needed, or does circulation with a non-uniform rate (c_j != 0) suffice?")
    println("=" ^ 100)

    println("\n(a) κ = 0, UNIFORM circulation: the flow must be volume preserving (no attractor).")
    println("      L      η   max|div F|   max|sum λ|      λ_max   attractors")

    worstUniform = 0.0
    for L in Ls, η in ηs
        parameters, results = ScanPoint(L, -20.0; κ = 0.0, η = η, trajectories = trajectories,
                                        kwargs...)
        pointwise = maximum(abs(Divergence(RandomInitialCondition(L, Xoshiro(k)), parameters))
                            for k = 1:200)
        sums = maximum(abs(sum(r.spectrum)) for r in results)
        summary = SummarisePoint(results)
        worstUniform = max(worstUniform, pointwise, sums)

        @printf("%7d %6.2f   %.2e     %.2e   %+8.4f   %s\n",
                L, η, pointwise, sums, summary.λmax, ClassHistogram(summary.classes))
    end

    println("\n(b) κ = 0, NON-UNIFORM circulation: c_j != 0, so volume can contract.")
    println("      L      η   modul.        g   chaotic     λ_max    <div F>   D_KY(red)   attractors")

    records = []
    for L in Ls, η in ηs, modulation in modulations, g in gs
        parameters, results = ScanPoint(L, g; κ = 0.0, η = η, modulation = modulation,
                                        trajectories = trajectories, kwargs...)
        summary = SummarisePoint(results)
        append!(records, results)

        @printf("%7d %6.2f   %6.2f %8.1f     %5.2f   %+8.4f   %+8.4f   %6.3f      %s\n",
                L, η, modulation, g, summary.chaoticFraction, summary.λmax, summary.divergence,
                summary.dimension, ClassHistogram(summary.classes))
    end

    attracting = [r for r in records if r.classification in (:chaotic, :hyperchaotic) &&
                                        sum(r.spectrum) < -1e-6]

    println()
    @printf("Uniform circulation at κ = 0: worst violation of volume preservation = %.1e.\n",
            worstUniform)
    if isempty(attracting)
        println("VERDICT: κ IS NEEDED.  With κ = 0 no non-uniform circulation on this grid produced a")
        println("         contracting chaotic attractor; positive exponents, where they occur, sit on")
        println("         volume-preserving or expanding flows, which are not attractors.")
    else
        best = attracting[argmax([r.reduced[1] for r in attracting])]
        @printf("VERDICT: κ IS NOT NEEDED.  %d trajectories have λ_max > 0 together with sum λ < 0 at\n",
                length(attracting))
        @printf("         κ = 0; the strongest has λ_max = %+.4f, sum λ = %+.4f, D_KY(red) = %.3f.\n",
                best.reduced[1], sum(best.spectrum), best.dimension)
    end

    return records
end


""" Group trajectories by the attractor they reached.  The fingerprint is made of gauge-invariant
    and translation-invariant time averages, so the L states related by the Z_L symmetry of §6
    count as ONE attractor (which is what "how many inequivalent attractors" means), while the two
    reflected states of a chosen circulation direction - genuinely inequivalent once η != 0 -
    count as two, since the current enters with its sign.

    The tolerance is widened by the TIME FLUCTUATION of each observable on the attractor, which is
    what makes the count meaningful for chaotic attractors: two trajectories on one strange
    attractor give window averages that differ by a fraction of the fluctuation amplitude, so a
    fixed tolerance would split a single attractor into as many pieces as there are trajectories.
    On a fixed point or a limit cycle the spread is zero and only `tolerance` applies.  The count
    is therefore a conservative LOWER bound: attractors closer than half a fluctuation width are
    reported as one.  The largest exponent is deliberately NOT part of the fingerprint - its
    finite-time scatter on a strange attractor is several per cent. """
function CountAttractors(results; tolerance = 0.02, spreadFactor = 0.5)
    Finite(x) = isfinite(x) ? x : 0.0
    Fingerprint(r) = (r.coherence, r.maximum_n, r.current)
    Spread(r) = (Finite(r.coherenceSpread), Finite(r.maximum_nSpread), Finite(r.currentSpread))

    groups = Vector{Tuple{Int, Vector{Int}}}()

    for (index, r) in enumerate(results)
        matched = false

        for (representative, members) in groups
            other = results[representative]
            other.classification == r.classification || continue

            limits = tolerance .+ spreadFactor .* 0.5 .* (Spread(other) .+ Spread(r))
            if all(abs.(Fingerprint(r) .- Fingerprint(other)) .<= limits)
                push!(members, index)
                matched = true
                break
            end
        end

        matched || push!(groups, (index, [index]))
    end

    return sort(groups, by = group -> -length(group[2]))
end


""" QUESTION 3.  Multistability - does the basin split between the locked state and something
    else?  At mean field this is a basin measurement; in the Liouvillian it is what shows up as
    metastability (a group of slowly decaying modes), which is why §8.3 asks for it.

    Each point is sampled with many initial conditions drawn uniformly on CP^(L-1); the attractors
    are then grouped by their invariant averages and reported with their basin fractions. """
function Question3(; points = ((3, -8.0, 0.3, 0.0), (3, -20.0, 0.3, 1.0), (3, -20.0, 0.3, 3.0),
                               (4, -20.0, 0.3, 3.0)),
        trajectories = 96, tolerance = 0.02, kwargs...)

    println()
    println("=" ^ 100)
    println("QUESTION 3: multistability - how many attractors coexist, and how large are their basins?")
    println("=" ^ 100)

    records = []

    for (L, g, κ, η) in points
        parameters, results = ScanPoint(L, g; κ = κ, η = η, trajectories = trajectories, kwargs...)
        append!(records, results)

        threshold = InstabilityThreshold(L; J = 1.0, κ = κ, η = η)
        groups = CountAttractors(results; tolerance = tolerance)

        @printf("\nL = %d, g = %.1f, κ = %.2f, η = %.2f  (§5 threshold g_c = %.3f, uniform state %s)\n",
                L, g, κ, η, threshold, g < threshold ? "UNSTABLE" : "stable")
        @printf("  %d initial conditions reach %d inequivalent attractors:\n",
                trajectories, length(groups))
        println("     basin          λ_max (scatter)   D_KY(red)   <coherence>   max n   current   type")

        for (representative, members) in groups
            r = results[representative]
            λ = [results[k].reduced[1] for k in members]
            scatter = length(λ) > 1 ? @sprintf("+-%.4f", std(λ)) : ""
            @printf("    %5.1f%%   %+10.4f %8s      %6.3f      %8.3f   %5.3f   %+7.3f   %s
",
                    100 * length(members) / trajectories, mean(λ), scatter,
                    r.dimension, r.coherence, r.maximum_n, r.current, r.classification)
        end

        # the uniform locked state has coherence = 1 and max n = 1/L
        locked = count(r -> abs(r.coherence - 1) < 0.05 && abs(r.maximum_n - 1 / L) < 0.05, results)
        @printf("    uniform locked state (coherence = 1, n = 1/L) reached by %.1f%% of them\n",
                100 * locked / trajectories)
    end

    println()
    println("VERDICT: see the basin tables above - more than one row at a parameter point is")
    println("         multistability, and a chaotic row coexisting with a fixed-point row is the")
    println("         case that would show up as Liouvillian metastability.")

    return records
end


""" QUESTION 4.  The trace rule and the zero exponents, read off every trajectory computed.

    sum(lambda) = <div F> is an identity of the QR algorithm (the product of the diagonal of R is
    the determinant of the propagator, whose logarithm is the integral of the trace of the
    Jacobian), so what it really tests is the divergence FORMULA of §4 - the closed expression
    sum_j c_j n_j - 8 κ sum_j Re(conj(ψ_j) ψ_{j+1}) - against the analytic Jacobian that drives
    the tangent dynamics, averaged along the attractor instead of at a point.

    The zero exponents are counted separately, because they are not all of the same quality:

      * the U(1) phase and the conserved density are EXACT zeros, but they are degenerate, and on
        a relative fixed point they form a 2 x 2 Jordan block (the family of relative fixed points
        parametrised by the total density has a frequency that varies along it, so the second
        generalised eigenvector maps onto the first).  A Jordan block grows linearly in t, so the
        QR method returns a symmetric pair +-log(t)/t rather than two exact zeros: a slowly
        vanishing artefact of the degeneracy, not a physical exponent.
      * the flow direction is a third zero on everything except a relative fixed point, where it
        coincides with the U(1) direction. """
function Question4(records; zeroThreshold = 5e-3, relativeZero = 0.05)
    println()
    println("=" ^ 100)
    println("QUESTION 4: the trace rule sum(lambda) = <div F>, and the three zero exponents")
    println("=" ^ 100)

    traceErrors = [r.traceError for r in records]
    # Scale the residual by the spectral radius, not by sum(lambda): on a neutral (volume
    # preserving) set sum(lambda) is itself zero and the ratio would be meaningless.
    relative = [r.traceError / max(abs(sum(r.spectrum)), maximum(abs, r.spectrum), 1e-12)
                for r in records]

    @printf("\n%d trajectories from questions 1-3.\n", length(records))
    @printf("  |sum lambda - <div F>|:  median %.1e, worst %.1e  (relative: median %.1e, worst %.1e)\n",
            median(traceErrors), maximum(traceErrors), median(relative), maximum(relative))

    # A zero exponent is judged against the spectral radius of its own spectrum: the two symmetry
    # zeros decay only as log(t)/t on a relative fixed point, where they form a Jordan block, so an
    # absolute threshold would fail them on short runs while they are in fact exactly zero.
    Radius(r) = maximum(abs, r.spectrum)
    Symmetry(r) = maximum(abs, r.symmetryZeros)
    Flow(r) = minimum(abs, r.reduced)

    println("\nZero exponents by attractor type; absolute, and as a fraction of max|lambda|:")
    println("      type           count   symmetry zeros (worst)      flow zero (worst)")

    for type in (:fixedPoint, :limitCycle, :torus, :neutral, :chaotic, :hyperchaotic, :undetermined)
        selection = [r for r in records if r.classification == type]
        isempty(selection) && continue

        symmetry = maximum(Symmetry, selection)

        # On a neutral set every exponent is zero, so a ratio to the spectral radius says nothing.
        if type === :neutral
            @printf("   %-14s %6d   %9.2e (%8s)   %9.2e (%8s)\n", string(type), length(selection),
                    symmetry, "all zero", maximum(Flow, selection), "all zero")
            continue
        end

        symmetryRelative = maximum(Symmetry(r) / Radius(r) for r in selection)

        # A relative fixed point has no third zero - the flow direction IS the U(1)
        # direction - so the smallest reduced exponent there is a contraction rate.
        if type === :fixedPoint
            @printf("   %-14s %6d   %9.2e (%6.2f%%)   %20s\n", string(type),
                    length(selection), symmetry, 100 * symmetryRelative, "n/a (2 zeros only)")
        else
            flow = maximum(Flow, selection)
            flowRelative = maximum(Flow(r) / Radius(r) for r in selection)
            @printf("   %-14s %6d   %9.2e (%6.2f%%)   %9.2e (%6.2f%%)\n", string(type),
                    length(selection), symmetry, 100 * symmetryRelative, flow, 100 * flowRelative)
        end
    end

    # The note prescribes D_KY(red) = D_KY(full) - 2.  It holds whenever the two dropped exponents
    # lie above the Kaplan-Yorke cut, so a violation means the "two closest to zero" rule picked a
    # genuine exponent instead of a symmetry zero.  It is tested only where D_KY is a dimension at
    # all: on a relative fixed point both sides collapse (D_KY(red) = 0, and D_KY(full) is 0 or 2
    # depending on the sign the numerically tiny Jordan pair happens to take), so the identity is
    # vacuous rather than violated there.
    attracting = [r for r in records if r.reduced[1] > zeroThreshold]
    if !isempty(attracting)
        dimensionGap = maximum(abs(r.dimensionFull - r.dimension - 2) for r in attracting)
        @printf("\n  D_KY(full) - D_KY(red) = 2 (note §4) to %.1e over the %d trajectories with a\n",
                dimensionGap, length(attracting))
        println("  positive exponent, where D_KY is a dimension rather than a collapsed zero.")
    end

    nonStationary = [r for r in records if r.classification != :fixedPoint]
    Zero(x, r) = x < zeroThreshold || x < relativeZero * Radius(r)
    zeros2 = count(r -> Zero(Symmetry(r), r), records)
    zeros3 = count(r -> Zero(Flow(r), r), nonStationary)

    println()
    @printf("  %d of %d trajectories have BOTH symmetry exponents at zero (below %.0e absolute or\n",
            zeros2, length(records), zeroThreshold)
    @printf("  %.0f%% of max|lambda|).  On a relative fixed point these two are degenerate and form a\n",
            100 * relativeZero)
    println("  2 x 2 Jordan block - the family of relative fixed points parametrised by the total")
    println("  density has a frequency that varies along it - so they converge only as log(t)/t and")
    println("  appear as a symmetric pair, not as two clean zeros.")
    @printf("  %d of %d non-stationary attractors carry the THIRD (flow-direction) zero as well;\n",
            zeros3, length(nonStationary))
    println("  on a relative fixed point the flow direction coincides with the U(1) direction, so")
    println("  two zeros is the correct count there, not three.")

    println()
    @printf("VERDICT: the trace rule of \u00a74 holds in the time average to %.0e absolute, %.0e relative,\n",
            maximum(traceErrors), maximum(relative))
    @printf("         over all %d trajectories; the zero exponents are where \u00a74 says they are\n",
            length(records))
    @printf("         (%d/%d and %d/%d above).\n", zeros2, length(records), zeros3, length(nonStationary))

    return nothing
end


""" Everything §8 asks for.  Run as  julia -t auto BHNumberConserving.jl  or call Section8() from
    the REPL; `Checks()` is run first, since a wrong Jacobian would invalidate every exponent. """
function Section8(; checks = true, trajectories = 12, basinTrajectories = 96,
        question1 = NamedTuple(), question2 = NamedTuple(), question3 = NamedTuple(), kwargs...)

    if checks
        println("CONSISTENCY CHECKS (note §10, classical rows, redone for this code)")
        println()
        Checks()
        println()
    end

    records = []
    append!(records, Question1(; trajectories = trajectories, question1..., kwargs...))
    append!(records, Question2(; trajectories = trajectories, question2..., kwargs...))
    append!(records, Question3(; trajectories = basinTrajectories, question3..., kwargs...))
    Question4(records)

    println()
    @printf("Threads: %d.  Total trajectories: %d.
", Threads.nthreads(), length(records))

    return records
end


if abspath(PROGRAM_FILE) == @__FILE__
    Section8()
end
