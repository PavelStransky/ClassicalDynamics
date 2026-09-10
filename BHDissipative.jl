# dissipative_bh_lyapunov.jl
#
# Classical (N → ∞) Bose–Hubbard dynamics with pump, loss and dephasing, and the
# maximal Lyapunov exponent (Benettin algorithm), using DifferentialEquations.jl.
#
#   dq_i = [-J Σ_{j∼i} p_j + g I_i p_i - κ/2 q_i] dt + √γd p_i ∘ dW_i
#   dp_i = [ J Σ_{j∼i} q_j - g I_i q_i - κ/2 p_i] dt - √γd q_i ∘ dW_i     (Stratonovich)
#
#   I_i = (q_i² + p_i²)/2,   κ = γ_l - γ_p,   ψ_i = (q_i + i p_i)/√2
#
# State vector u = [q; p; δq; δp] (length 4L): the trajectory and one tangent vector are
# integrated together, so both automatically see the SAME noise realisation.
#
# Methods
#   :split (default) Strang splitting. The deterministic part (hopping, interaction,
#          damping and the tangent equations) is an ODEProblem solved to tolerance; the
#          dephasing is applied exactly, as a random rotation of every (q_i,p_i) and
#          (δq_i,δp_i) plane, at the midpoint of each interval Δt. The norm law
#          Σ I_i(t) = e^{-κt} Σ I_i(0) then holds to solver tolerance.
#   :sde   SDEProblem with non-diagonal noise, solved with EulerHeun (Stratonovich).
#          Simpler, but the norm drifts at O(Δt). Useful as a cross-check.
#          Do NOT use Itô solvers (EM, SRIW1, SOSRI, ...) on these Stratonovich equations.
#
# Optional σadd = √((γ_l + γ_p)/(2N)) adds the finite-N (truncated Wigner) additive noise.
# It does not enter the tangent equations. With κ = γd = σadd = 0 the code reduces to the
# Hamiltonian DNLS.

using DifferentialEquations
using StochasticDiffEq
using LinearAlgebra, Random, Statistics, Printf

# ------------------------------------------------------------------------- lattice
function chain_neighbors(L; periodic = true)
    nbrs = [Int[] for _ in 1:L]
    for i in 1:L-1
        push!(nbrs[i], i + 1)
        push!(nbrs[i+1], i)
    end
    if periodic && L > 2
        push!(nbrs[1], L)
        push!(nbrs[L], 1)
    end
    return nbrs
end

# ------------------------------------------------------- deterministic part + tangent
function drift!(du, u, prm, t)
    (; L, J, g, κ, nbrs) = prm
    @inbounds for i in 1:L
        q_i, p_i = u[i], u[L+i]
        dq_i, dp_i = u[2L+i], u[3L+i]
        sq = sp = sdq = sdp = 0.0
        for j in nbrs[i]
            sq += u[j];     sp += u[L+j]
            sdq += u[2L+j]; sdp += u[3L+j]
        end
        I = 0.5 * (q_i^2 + p_i^2)
        c = q_i * dq_i + p_i * dp_i                      # = δI_i
        du[i]    = -J * sp  + g * I * p_i  - 0.5κ * q_i
        du[L+i]  =  J * sq  - g * I * q_i  - 0.5κ * p_i
        du[2L+i] = -J * sdp + g * I * dp_i + g * p_i * c - 0.5κ * dq_i
        du[3L+i] =  J * sdq - g * I * dq_i - g * q_i * c - 0.5κ * dp_i
    end
    return nothing
end

# ------------------------------------------ noise matrix for the :sde method (4L × nW)
# Column i: dephasing on site i. Columns L+i and 2L+i: additive noise on q_i, p_i.
function noise!(G, u, prm, t)
    (; L, γd, σadd) = prm
    s = sqrt(γd)
    @inbounds for i in 1:L
        G[i, i]    =  s * u[L+i]
        G[L+i, i]  = -s * u[i]
        G[2L+i, i] =  s * u[3L+i]
        G[3L+i, i] = -s * u[2L+i]
        if σadd > 0
            G[i, L+i]    = σadd
            G[L+i, 2L+i] = σadd
        end
    end
    return nothing
end

# ------------------------------- exact dephasing step (+ additive noise) for :split
function dephase_kick!(u, L, sγ, sa, rng)
    @inbounds for i in 1:L
        s, c = sincos(sγ * randn(rng))                   # angle √γd ΔW_i
        q, p = u[i], u[L+i]
        u[i]    =  c * q + s * p
        u[L+i]  = -s * q + c * p
        dq, dp = u[2L+i], u[3L+i]                        # same rotation on tangent
        u[2L+i] =  c * dq + s * dp
        u[3L+i] = -s * dq + c * dp
        if sa > 0
            u[i]   += sa * randn(rng)
            u[L+i] += sa * randn(rng)
        end
    end
    return nothing
end

# ------------------------------------------------------------------- observables
function energy(u, prm)                                   # H_cl(q, p)
    (; L, J, g, nbrs) = prm
    E = 0.0
    @inbounds for i in 1:L
        I = 0.5 * (u[i]^2 + u[L+i]^2)
        E += 0.5g * I^2
        for j in nbrs[i]                                 # each bond visited twice
            E -= 0.5J * (u[i] * u[j] + u[L+i] * u[L+j])
        end
    end
    return E
end

function _energy(ψ::AbstractVector{<:Complex}, J, g, nbrs)   # same H_cl in terms of ψ
    E = 0.0
    @inbounds for i in eachindex(ψ)
        E += 0.5g * abs2(ψ[i])^2
        for j in nbrs[i]
            E -= J * real(conj(ψ[i]) * ψ[j])            # each bond visited twice
        end
    end
    return E
end
energy(ψ::AbstractVector{<:Complex}, J::Real, g::Real; periodic = true) =
    _energy(ψ, J, g, chain_neighbors(length(ψ); periodic))

mutable struct LyapRecord
    t::Vector{Float64}
    S::Vector{Float64}       # cumulative log-stretching
    dens::Vector{Float64}    # n(t) = Σ I_i / L
    E::Vector{Float64}
    Ssum::Float64
end
LyapRecord() = LyapRecord(Float64[], Float64[], Float64[], Float64[], 0.0)

# Benettin step. Returns false if the run should stop.
function renormalize!(integ, prm, rec, Δt)
    L = prm.L
    u = integ.u
    v = @view u[2L+1:4L]
    nv = norm(v)
    rec.Ssum += log(nv)
    v ./= nv
    u_modified!(integ, true)
    Imax = maximum(i -> 0.5 * (u[i]^2 + u[L+i]^2), 1:L)
    push!(rec.t, integ.t)
    push!(rec.S, rec.Ssum)
    push!(rec.dens, sum(abs2, @view u[1:2L]) / (2L))
    push!(rec.E, energy(u, prm))
    if !isfinite(nv) || prm.g * Imax * Δt > 0.2
        @warn "stopping at t = $(integ.t): non-finite tangent norm or g·max(I)·Δt > 0.2"
        return false
    end
    return true
end

# ------------------------------------------------------------------ main routine
"""
    lyapunov(ψ0, J, g; κ=0, γd=0, σadd=0, method=:split, Δt=1e-2, t_max=1e3,
             t_renorm=1.0, seed=1, alg=Tsit5(), abstol=1e-10, reltol=1e-10, periodic=true)

Maximal Lyapunov exponent for initial amplitudes `ψ0` (complex, ψ = (q+ip)/√2).
The density n = Σ|ψ0_i|²/L and the energy are set by `ψ0` (see `random_state`,
`state_with_energy`).

Time and precision:
  t_max          total integration time
  t_renorm       Benettin renormalisation interval = time resolution of the output
  abstol, reltol tolerances of the ODE solver (deterministic part, :split only)
  alg            ODE solver, e.g. Tsit5() or Vern9() for high precision
  Δt             noise/splitting interval (:split) or EulerHeun step (:sde); must divide
                 `t_renorm`. Without noise it is only used for the resolution check.
`seed` must be ≥ 1 (seed 0 means "random" in StochasticDiffEq).

Returns a NamedTuple with `t`, `S`, `λ = S/t`, `λ_shifted = λ + κ/2`, `dens`, `E`,
and the initial energy `E0` and density `n0`.
For κ ≠ 0 use windowed exponents, see `window_lambda`.
"""
function lyapunov(ψ0::AbstractVector{<:Complex}, J::Real, g::Real;
                  κ = 0.0, γd = 0.0, σadd = 0.0, method::Symbol = :split,
                  Δt = 1e-2, t_max = 1e3, t_renorm = 1.0, seed::Integer = 1,
                  alg = Tsit5(), abstol = 1e-10, reltol = 1e-10, periodic = true)
    L = length(ψ0)
    nsub = round(Int, t_renorm / Δt)
    nwin = round(Int, t_max / t_renorm)
    @assert nsub * Δt ≈ t_renorm "t_renorm must be a multiple of Δt"
    @assert seed ≥ 1
    prm = (L = L, J = float(J), g = float(g), κ = float(κ), γd = float(γd),
           σadd = float(σadd), nbrs = chain_neighbors(L; periodic))
    rng = Xoshiro(seed)
    v0 = randn(rng, 2L)
    v0 ./= norm(v0)
    u0 = vcat(sqrt(2) .* real.(ψ0), sqrt(2) .* imag.(ψ0), v0)
    tspan = (0.0, 2.0 * t_max)            # margin; the run is driven by step! below
    rec = LyapRecord()
    noisy = γd > 0 || σadd > 0

    if method === :sde && noisy
        nW = σadd > 0 ? 3L : L
        prob = SDEProblem(drift!, noise!, u0, tspan, prm;
                          noise_rate_prototype = zeros(4L, nW))
        integ = init(prob, EulerHeun(); dt = Δt, seed = UInt64(seed),
                     save_everystep = false, maxiters = 10^12)
        for _ in 1:nwin
            step!(integ, t_renorm, true)
            renormalize!(integ, prm, rec, Δt) || break
        end
    elseif method === :split || !noisy
        prob = ODEProblem(drift!, u0, tspan, prm)
        integ = init(prob, alg; abstol, reltol, save_everystep = false,
                     maxiters = 10^12)
        sγ, sa = sqrt(γd * Δt), σadd * sqrt(Δt)
        for _ in 1:nwin
            if noisy
                for _ in 1:nsub          # Strang: half flow, exact noise, half flow
                    step!(integ, Δt / 2, true)
                    dephase_kick!(integ.u, L, sγ, sa, rng)
                    u_modified!(integ, true)
                    step!(integ, Δt / 2, true)
                end
            else
                step!(integ, t_renorm, true)
            end
            renormalize!(integ, prm, rec, Δt) || break
        end
    else
        error("method must be :split or :sde")
    end
    λ = rec.S ./ rec.t
    return (t = rec.t, S = rec.S, λ = λ, λ_shifted = λ .+ κ / 2,
            dens = rec.dens, E = rec.E,
            E0 = energy(u0, prm), n0 = sum(abs2, ψ0) / L)
end

# Finite-time exponent on the window [t1, t2].
function window_lambda(res, t1, t2)
    i1 = searchsortedfirst(res.t, t1)
    i2 = searchsortedlast(res.t, t2)
    return (res.S[i2] - res.S[i1]) / (res.t[i2] - res.t[i1])
end

# Uniform on the norm shell Σ|ψ_i|² = nL (the long-time state produced by dephasing).
function random_state(L, n; rng = Random.default_rng())
    z = randn(rng, ComplexF64, L)
    return z .* (sqrt(n * L) / norm(z))
end

# Norm-preserving gradient flow of H_cl on the sphere Σ|ψ_i|² = Nrm², up or down in
# energy, until E_target is reached.
function _energy_flow(ψ0, E_target, J, g, nbrs, Nrm, tol, maxiter)
    ψ = ψ0 .* (Nrm / norm(ψ0))
    G = similar(ψ)
    E = _energy(ψ, J, g, nbrs)
    for _ in 1:maxiter
        abs(E - E_target) < tol * max(1, abs(E_target)) && return ψ
        for i in eachindex(ψ)                            # G = ∂H/∂ψ_i*
            G[i] = g * abs2(ψ[i]) * ψ[i]
            for j in nbrs[i]
                G[i] -= J * ψ[j]
            end
        end
        G .-= (real(dot(ψ, G)) / real(dot(ψ, ψ))) .* ψ    # tangent to the norm sphere
        G2 = real(dot(G, G))
        G2 < 1e-24 && error("gradient flow reached a stationary state at E = $E; " *
                            "E_target = $E_target is outside the reachable range")
        η = min(0.05 / (2J + abs(g) * maximum(abs2, ψ)), abs(E - E_target) / (2 * G2))
        ψ .+= sign(E_target - E) * η .* G
        ψ .*= Nrm / norm(ψ)
        E = _energy(ψ, J, g, nbrs)
    end
    error("energy flow did not converge")
end

"""
    state_with_energy(L, n, E_target, J, g; rng, tol=1e-10, periodic=true)

State with Σ|ψ_i|² = nL and H_cl = E_target (relative accuracy `tol`).
  * E(uniform) ≤ E_target ≤ E(random): bisection along ψ(s) ∝ (1-s)·uniform + s·random,
    i.e. the uniform (ground) state dressed with random fluctuations of all modes;
  * E_target > E(random): gradient ascent from the random state;
  * E_target < E(uniform): gradient descent from the uniform state (open chains, g < 0).
Errors if E_target lies outside the reachable energy range.
"""
function state_with_energy(L, n, E_target, J, g; rng = Random.default_rng(),
                           tol = 1e-10, maxiter = 10^6, periodic = true)
    nbrs = chain_neighbors(L; periodic)
    Nrm = sqrt(n * L)
    z = random_state(L, n; rng)
    u = fill(complex(sqrt(float(n))), L)
    Ez, Eu = _energy(z, J, g, nbrs), _energy(u, J, g, nbrs)
    if E_target >= Ez
        return _energy_flow(z, E_target, J, g, nbrs, Nrm, tol, maxiter)
    elseif E_target >= Eu
        mix = s -> begin
            ψ = (1 - s) .* u .+ s .* z
            ψ .* (Nrm / norm(ψ))
        end
        lo, hi = 0.0, 1.0                               # E(lo) ≤ E_target ≤ E(hi)
        for _ in 1:200
            mid = (lo + hi) / 2
            Em = _energy(mix(mid), J, g, nbrs)
            abs(Em - E_target) < tol * max(1, abs(E_target)) && return mix(mid)
            Em < E_target ? (lo = mid) : (hi = mid)
        end
        return mix((lo + hi) / 2)
    else
        return _energy_flow(u .+ 1e-3 .* z, E_target, J, g, nbrs, Nrm, tol, maxiter)
    end
end

# Mean and standard error of λ(t_max) over noise realisations (threaded).
function ensemble_lambda(make_ψ0, J, g; seeds = 1:8, kwargs...)
    λs = zeros(length(seeds))
    Threads.@threads for k in eachindex(seeds)
        res = lyapunov(make_ψ0(seeds[k]), J, g; seed = seeds[k], kwargs...)
        λs[k] = res.λ[end]
    end
    return mean(λs), std(λs) / sqrt(length(λs))
end

# ------------------------------------------------------------------ demo / checks
function demo()
    # n = Σ|ψ_i|²/L is the classical density; only the combination g·n/J matters.
    L, J, g, n = 4, 0.5, 2.0, 1.0
    ψ0 = random_state(L, n; rng = Xoshiro(1))
    n0 = sum(abs2, ψ0) / L
    normerr(r, κ) = maximum(abs.(r.dens ./ (n0 .* exp.(-κ .* r.t)) .- 1))

    @printf("E(ψ0) = %.3f   [ground state: %.3f,  ⟨E⟩ at T=∞: %.3f]\n",
            energy(ψ0, J, g), -2J * n * L + g * n^2 * L / 2, g * n^2 * L^2 / (L + 1))

    r = lyapunov(ψ0, J, g; t_max = 500)
    @printf("Hamiltonian:          λ = %+.4f   rel. energy drift = %.1e   norm err = %.1e\n",
            r.λ[end], maximum(abs.(r.E .- r.E[1])) / abs(r.E[1]), normerr(r, 0.0))

    ψlow = state_with_energy(L, n, 0.0, J, g; rng = Xoshiro(2))   # prescribed energy
    r = lyapunov(ψlow, J, g; t_max = 500)
    @printf("Hamiltonian, E0 = %.2f: λ = %+.4f\n", r.E0, r.λ[end])

    for m in (:split, :sde)
        r = lyapunov(ψ0, J, g; γd = 0.1, method = m, Δt = 5e-3, t_max = 1000)
        @printf("γd = 0.1, %-6s:       λ = %+.4f   norm err = %.1e\n", string(m), r.λ[end], normerr(r, 0.0))
    end

    r = lyapunov(ψ0, J, 0.0; κ = 0.3, γd = 0.5, t_max = 50)
    @printf("g = 0, κ = 0.3:       λ = %+.6f   (exact: -κ/2 = -0.15)\n", r.λ[end])

    κ = 0.02
    r = lyapunov(ψ0, J, g; κ, t_max = 1000)
    println("κ = $κ, finite-time exponents on time windows [t1, t2]:")
    for (t1, t2) in ((1, 50), (100, 150), (250, 300), (500, 1000))
        λw = window_lambda(r, t1, t2)
        nw = mean(r.dens[(r.t .>= t1) .& (r.t .<= t2)])
        @printf("   [%3d,%3d]:  λ = %+.3f   λ+κ/2 = %+.3f   ⟨n⟩ = %.3f\n", t1, t2, λw, λw + κ / 2, nw)
    end

    m, se = ensemble_lambda(s -> random_state(L, n; rng = Xoshiro(s)), J, g;
                            γd = 0.1, t_max = 300, seeds = 1:4)
    @printf("ensemble, γd = 0.1:   λ = %.3f ± %.3f\n", m, se)
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo()
end