# BHNumberConservingQuantumDynamics.jl
#
# Time evolution in the N-boson space of the number-conserving dissipative Bose-Hubbard ring - the
# quantum side of BHNumberConserving.jl in the TIME domain rather than in the spectrum.  Two items
# of number-conserving-BH-paper-TODO.md:
#
#   consistency   K1   the Lindblad master equation started from an SU(L) coherent state against
#                      the classical equations of motion started from the same point, for
#                      N = 10 ... 60.  The observables are the one-body density matrix
#                      D_ij = <b_i^dag b_j> / N, whose classical value is conj(ψ_i) ψ_j, and from it
#                      n_j, the bond coherence C and the current.  (<b_j> itself vanishes
#                      identically in a fixed-N space - the strong U(1) of §6 - so the one-body
#                      density matrix is the right object, not <b_j>.)  Up to the Ehrenfest time
#                      the deviation must fall like 1/N; a wrong rate convention (γ = N Γ,
#                      κ_c = κ / N, the factor 4 of c_j, the sign of η) shows up as an O(1)
#                      deviation that does NOT shrink with N.
#   trajectories  C8   quantum trajectories (Monte Carlo wave functions) at the multistable point
#                      (g, η) = (-20, 1), N = 10 ... 50: the state is followed through the
#                      single-trajectory expectation values, classified as SELF-TRAPPED (one site
#                      holding nearly all bosons, the classical fixed point with max n = 0.994) or
#                      NOT, and the dwell times give the escape/switching rates as functions of N,
#                      to be compared with the gap G_∞ = 0.165 (A8).  Half of the trajectories
#                      start from a random coherent state, half from the trapped state.
#
# Both propagate with Krylov exponentials (KrylovKit.exponentiate): the Kerr term makes the
# generator's norm grow like |g| N / 2, which a Runge-Kutta scheme would have to resolve step by
# step, while a Krylov exponential absorbs it in the subspace dimension.
#
# USAGE
#   julia -t 8  BHNumberConservingQuantumDynamics.jl consistency [--N 10,20,30,40]
#   julia -t 20 BHNumberConservingQuantumDynamics.jl trajectories [--N 10,20,30,40,50]
#                                                     [--trajectories 16] [--time 2000]
#
# OUTPUT in $BH_RESULTS_DIR (default ~/results/bh/number-conserving/quantum/3/dynamics).

ENV["CD_NO_PLOTS"] = "true"
isdefined(Main, :LiouvillianSector) || include(joinpath(@__DIR__, "BHNumberConservingLiouvillian.jl"))

# The classical model lives in its own module: both files define Checks() and a few other generic
# names, and nothing here should depend on which definition happens to win.
module Classical
    include(joinpath(@__DIR__, "BHNumberConserving.jl"))
end

using KrylovKit
using SparseArrays
using LinearAlgebra
using Random
using Statistics
using Printf

const DYNAMICS_RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "quantum", "3", "dynamics"))


# ---------------------------------------------------------------------------------------------
# Coherent states and one-body observables
# ---------------------------------------------------------------------------------------------

""" SU(L) coherent state |ψ; N> = (sum_j ψ_j b_j^dag)^N / sqrt(N!) |0> in the Fock basis:
    amplitude sqrt(N! / prod a_j!) prod ψ_j^a_j for the occupation vector a (ψ normalised). """
function CoherentState(basis::NumberBasis, ψ)
    ψ = ψ ./ norm(ψ)
    logFactorial = cumsum(vcat(0.0, log.(1:basis.N)))
    amplitudes = [exp(0.5 * (logFactorial[basis.N + 1] - sum(logFactorial[a .+ 1]))) *
                  prod(ψ[j]^a[j] for j in eachindex(a)) for a in basis.states]
    return amplitudes ./ norm(amplitudes)
end

""" One-body density matrix <b_i^dag b_j> / N of a density matrix or of a (normalised) pure state. """
function OneBody(hops, ρ::AbstractMatrix, N)
    L = size(hops, 1)
    return [tr(hops[i, j] * ρ) / N for i = 1:L, j = 1:L]
end

function OneBody(hops, ψ::AbstractVector, N)
    L = size(hops, 1)
    return [dot(ψ, hops[i, j] * ψ) / (N * real(dot(ψ, ψ))) for i = 1:L, j = 1:L]
end

""" b_i^dag b_j as sparse matrices; D_ij = <b_i^dag b_j> / N has the classical value conj(ψ_i) ψ_j. """
Hops(basis) = [Hop(basis, i, j) for i = 1:basis.L, j = 1:basis.L]

Densities(D) = real.(diag(D))
BondCoherence(D) = sum(real(D[j, mod1(j + 1, size(D, 1))]) for j = 1:size(D, 1))
Current(D) = sum(imag(D[j, mod1(j + 1, size(D, 1))]) for j = 1:size(D, 1))


# ---------------------------------------------------------------------------------------------
# K1: master equation against the classical equations of motion
# ---------------------------------------------------------------------------------------------

""" dρ/dt = K ρ + ρ K^dag + sum_k L_k ρ L_k^dag, on the vectorised density matrix. """
function LindbladAction(K, jumps, d)
    KAdjoint = sparse(adjoint(K))
    jumpAdjoints = [sparse(adjoint(Lk)) for Lk in jumps]
    return function (x)
        ρ = reshape(x, d, d)
        out = K * ρ + ρ * KAdjoint
        for (Lk, LkAdjoint) in zip(jumps, jumpAdjoints)
            out .+= Lk * (ρ * LkAdjoint)
        end
        return vec(out)
    end
end

""" One comparison: master equation and classical flow from the same point ψ0, sampled every
    `step` up to `time`.  Returns rows (t, max_ij |ΔD_ij|, n_1^q, n_1^cl, C^q, C^cl, I^q, I^cl,
    and max_ij |ΔD_ij| against each deliberately WRONG classical model in `controls`).

    `controls` are classical parameter sets that differ from the correct one by a convention
    error - the sign of η, a missing factor 4 in κ, a Kerr term U N instead of U (N - 1)... - and
    must stay O(1) away from the quantum curve however large N is.  They are what makes the test
    able to fail: a comparison that no model can lose proves nothing. """
function Consistency(N, g, η, κ; ψ0, time = 4.0, step = 0.1, L = 3, controls = ())
    basis = NumberBasis(L, N)
    d = length(basis)
    p = LiouvillianParameters(L, N; g = g, η = η, κ = κ)
    jumps = JumpOperators(basis, p)
    K = EffectiveGenerator(Hamiltonian(basis, p), jumps)
    action = LindbladAction(K, jumps, d)
    hops = Hops(basis)

    c = CoherentState(basis, ψ0)
    x = vec(c * c')
    times = 0.0:step:time

    function ClassicalTrajectory(gc, ηc, κc)
        classical = Classical.NumberConservingParameters(L; J = 1.0, g = gc, κ = κc, η = ηc)
        work = Classical.NumberConservingWorkspace(L)
        field!(du, u, q, t) = (Classical.UpdateWorkspace!(work, u, classical);
                               Classical.VectorField!(du, u, classical, work))
        return Classical.solve(Classical.ODEProblem(field!, ComplexF64.(ψ0 ./ norm(ψ0)), (0.0, time)),
                               Classical.DP8(); reltol = 1e-12, abstol = 1e-12, saveat = times).u
    end
    trajectory = ClassicalTrajectory(g, η, κ)
    wrong = [ClassicalTrajectory(control...) for control in controls]

    rows = []
    for (k, t) in enumerate(times)
        if k > 1
            x, info = exponentiate(action, step, x; tol = 1e-11, krylovdim = 60, ishermitian = false)
            info.converged == 0 && @warn "Krylov exponential not converged" N t
        end
        Dq = OneBody(hops, reshape(x, d, d), N)
        ψ = trajectory[k]
        Dc = conj.(ψ) * transpose(ψ)
        deviations = [maximum(abs, Dq .- conj.(w[k]) * transpose(w[k])) for w in wrong]
        push!(rows, (t, maximum(abs, Dq .- Dc), Densities(Dq)[1], Densities(Dc)[1],
                     BondCoherence(Dq), BondCoherence(Dc), Current(Dq), Current(Dc), deviations...))
    end
    return rows
end

""" K1.  Three tests, each at the times where the classical limit is supposed to hold.

    The Kerr term dephases a coherent state on the collapse time t_c ~ sqrt(N) / |g| (quantum
    phase diffusion), long before any chaotic Ehrenfest time: at g = -20 and N = 30 it is 0.27.
    The dissipative conventions are therefore tested where the interaction is weak and t_c long,
    and the Kerr convention at moderate g and short times:

      * (g, η, κ) = (-1, 3, 1)   - circulation and phase locking dominate; t <= 2
      * (-4, 3, 0.3)             - the Kerr term at moderate strength;       t <= 1
      * (-20, 3, 0.3)            - the chaotic point, only t <= 0.3 (t_c ~ 0.3 at N = 30)

    For each the columns give N max|ΔD| at a few times: constant in N means max|ΔD| ∝ 1/N, i.e.
    agreement.  The control columns use a classical model with η -> -η (sign of the circulation)
    and with κ -> κ / 4 (the factor 4 of c_j = (b_j^dag + b_{j+1}^dag)(b_j - b_{j+1})): their
    deviation must NOT shrink with N. """
function RunConsistency(Ns)
    directory = joinpath(DYNAMICS_RESULTS, "consistency")
    mkpath(directory)
    ψ0 = normalize(ComplexF64[0.8 + 0.1im, 0.3 - 0.4im, 0.2 + 0.25im])     # generic, fixed

    tests = ((-1.0, 3.0, 1.0, 2.0, (0.5, 1.0, 2.0)),
             (-4.0, 3.0, 0.3, 1.0, (0.25, 0.5, 1.0)),
             (-20.0, 3.0, 0.3, 0.3, (0.1, 0.2, 0.3)))

    for (g, η, κ, time, checkpoints) in tests
        controls = ((g, -η, κ), (g, η, κ / 4))
        @printf("
(g, η, κ) = (%.1f, %.1f, %.1f), t <= %.2f; collapse times sqrt(N)/|g| = %s
", g, η, κ,
                time, join([@sprintf("%.2f", sqrt(N) / abs(g)) for N in Ns], ", "))
        @printf("     N    %s    | control η -> -η at t = %.2f   κ -> κ/4
",
                join([@sprintf("N·max|ΔD|(t=%.2f)", t) for t in checkpoints], "  "), checkpoints[end])

        results = Vector{Any}(undef, length(Ns))
        Threads.@threads for i in eachindex(Ns)
            results[i] = Consistency(Ns[i], g, η, κ; ψ0 = ψ0, time = time, step = 0.05,
                                     controls = controls)
        end

        for (N, rows) in zip(Ns, results)
            at(t, column) = rows[argmin([abs(r[1] - t) for r in rows])][column]
            @printf("  %4d    %s    |   %.3f                          %.3f
", N,
                    join([@sprintf("%16.3f", N * at(t, 2)) for t in checkpoints], "  "),
                    at(checkpoints[end], 9), at(checkpoints[end], 10))
            file = joinpath(directory, @sprintf("N%03d_g%.2f_e%.2f_k%.2f.txt", N, g, η, κ))
            open(file, "w") do io
                println(io, "# t	max_abs_dD	n1_q	n1_cl	C_q	C_cl	I_q	I_cl	dD_wrong_eta_sign	dD_wrong_kappa_factor")
                for r in rows
                    println(io, join([@sprintf("%.10g", v) for v in r], "	"))
                end
            end
        end
    end
    println("
N·max|ΔD| constant in N: agreement (max|ΔD| ∝ 1/N).  The control columns are the")
    println("deviations max|ΔD| (not multiplied by N) from deliberately wrong classical models;")
    println("they must stay put as N grows - otherwise the test cannot tell right from wrong.")
end


# ---------------------------------------------------------------------------------------------
# C8: quantum trajectories and switching
# ---------------------------------------------------------------------------------------------

""" One Monte Carlo wave-function trajectory by the waiting-time method: the unnormalised state
    evolves with K (Krylov exponentials over `step`), a jump happens when |ψ|^2 drops to the
    random threshold r, at a time found by log-linear interpolation of |ψ|^2 inside the step and
    reached by one more exponential; the jump operator is drawn with probability ∝ |L_k ψ|^2.
    Every `record` time units the single-trajectory one-body density matrix is stored. """
function Trajectory(K, jumps, hops, ψ0, N, time, rng; step = 0.02, record = 0.5)
    ψ = ψ0 ./ norm(ψ0)
    r = rand(rng)
    t = 0.0
    nextRecord = 0.0
    records = NTuple{7, Float64}[]
    jumpsDone = 0
    generator = x -> K * x

    while t < time
        if t >= nextRecord
            D = OneBody(hops, ψ, N)
            n = Densities(D)
            push!(records, (t, n..., BondCoherence(D), Current(D), argmax(n)))
            nextRecord += record
        end

        next, _ = exponentiate(generator, step, ψ; tol = 1e-10, ishermitian = false)
        p0, p1 = real(dot(ψ, ψ)), real(dot(next, next))

        if p1 > r
            ψ = next
            t += step
            continue
        end

        # jump inside this step
        τ = step * log(p0 / r) / log(p0 / p1)
        ψ, _ = exponentiate(generator, τ, ψ; tol = 1e-10, ishermitian = false)
        t += τ
        weights = [real(dot(Lk * ψ, Lk * ψ)) for Lk in jumps]
        k = findfirst(cumsum(weights) .>= rand(rng) * sum(weights))
        ψ = jumps[k] * ψ
        ψ ./= norm(ψ)
        r = rand(rng)
        jumpsDone += 1
    end

    return records, jumpsDone
end

""" Trapped (1) when the most occupied site holds more than `threshold` of the bosons, averaged over
    a running window; otherwise 0.  The classical self-trapped fixed point has max n = 0.994, the
    chaotic set max n = 0.79 ± 0.11 (BHNumberConservingAttractors.jl basins). """
function TrappedSeries(records; threshold = 0.95, window = 10)
    maximum_n = [maximum(r[2:4]) for r in records]
    smoothed = [mean(maximum_n[max(1, i - window + 1):i]) for i in eachindex(maximum_n)]
    return [s > threshold ? 1 : 0 for s in smoothed]
end

""" Durations of the uninterrupted runs of each state (the first and last runs are censored and
    dropped) and the number of transitions. """
function Dwell(series, dt)
    runs = Dict(0 => Float64[], 1 => Float64[])
    start = 1
    for i = 2:length(series)
        if series[i] != series[i - 1]
            start > 1 && push!(runs[series[i - 1]], (i - start) * dt)
            start = i
        end
    end
    return runs
end

function RunTrajectories(Ns; trajectories = 16, time = 2000.0, g = -20.0, η = 1.0, κ = 0.3)
    directory = joinpath(DYNAMICS_RESULTS, "trajectories")
    mkpath(directory)

    # the classical self-trapped fixed point, as a coherent-state label
    classical = Classical.NumberConservingParameters(3; J = 1.0, g = g, κ = κ, η = η)
    trapped = nothing
    for seed = 1:200
        ψ = Classical.RandomInitialCondition(3, Xoshiro(seed))
        work = Classical.NumberConservingWorkspace(3)
        f!(du, u, q, t) = (Classical.UpdateWorkspace!(work, u, classical); Classical.VectorField!(du, u, classical, work))
        final = Classical.solve(Classical.ODEProblem(f!, ψ, (0.0, 3000.0)), Classical.DP8();
                                reltol = 1e-10, abstol = 1e-10, save_everystep = false, maxiters = 10^8).u[end]
        if maximum(abs2, final) > 0.99
            trapped = final
            break
        end
    end
    isnothing(trapped) && error("no initial condition reached the self-trapped state")
    @printf("self-trapped classical state: n = %s\n", string(round.(abs2.(trapped), digits = 4)))

    summary = []
    for N in Ns
        basis = NumberBasis(3, N)
        p = LiouvillianParameters(3, N; g = g, η = η, κ = κ)
        jumps = JumpOperators(basis, p)
        K = EffectiveGenerator(Hamiltonian(basis, p), jumps)
        hops = Hops(basis)
        trappedState = CoherentState(basis, trapped)

        results = Vector{Any}(undef, trajectories)
        elapsed = @elapsed Threads.@threads for i = 1:trajectories
            rng = Xoshiro(hash((N, i, "trajectory")))
            ψ0 = isodd(i) ? CoherentState(basis, Classical.RandomInitialCondition(3, rng)) : trappedState
            results[i] = Trajectory(K, jumps, hops, ψ0, N, time, rng)
        end

        escapes = Float64[]      # dwell times in the non-trapped state (ended by trapping)
        stays = Float64[]        # dwell times in the trapped state (ended by escape)
        hopsBetweenSites = 0
        trappedTime = 0.0
        for (i, (records, _)) in enumerate(results)
            open(joinpath(directory, @sprintf("N%03d_traj%02d.txt", N, i)), "w") do io
                println(io, "# t\tn1\tn2\tn3\tC\tcurrent\tmax_site")
                for r in records
                    println(io, join([@sprintf("%.6g", v) for v in r], "\t"))
                end
            end
            series = TrappedSeries(records)
            runs = Dwell(series, 0.5)
            append!(escapes, runs[0])
            append!(stays, runs[1])
            trappedTime += count(==(1), series) * 0.5
            sites = [Int(r[7]) for (r, s) in zip(records, series) if s == 1]
            hopsBetweenSites += count(i -> sites[i] != sites[i - 1], 2:length(sites))
        end

        totalTime = trajectories * time
        escapeRate = isempty(stays) ? 0.0 : length(stays) / trappedTime
        @printf("N = %3d (d_N = %4d): %.1f min; trapped %.1f%% of the time; %d escapes from the trapped state (rate %.2e), %d site changes while trapped; mean untrapped dwell %.1f\n",
                N, length(basis), elapsed / 60, 100 * trappedTime / totalTime, length(stays), escapeRate,
                hopsBetweenSites, isempty(escapes) ? NaN : mean(escapes))
        push!(summary, (N, trappedTime / totalTime, length(stays), escapeRate, hopsBetweenSites,
                        isempty(escapes) ? NaN : mean(escapes), mean(r[2] for r in results) / time))
    end

    open(joinpath(directory, "summary.txt"), "w") do io
        println(io, "# N\ttrapped_fraction\tescapes\tescape_rate\tsite_changes\tmean_untrapped_dwell\tjump_rate")
        for s in summary
            println(io, join([@sprintf("%.6g", v) for v in s], "\t"))
        end
    end
    println("compare the escape rates with the gap G_∞ = 0.165 at (g, η) = (-20, 1) (item A8)")
end


# ---------------------------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------------------------

function Options(arguments)
    options = Dict{String, String}()
    positional = String[]
    i = 1
    while i <= length(arguments)
        if startswith(arguments[i], "--")
            options[arguments[i][3:end]] = arguments[i + 1]
            i += 2
        else
            push!(positional, arguments[i])
            i += 1
        end
    end
    return positional, options
end

if abspath(PROGRAM_FILE) == @__FILE__
    positional, options = Options(ARGS)
    command = isempty(positional) ? "" : positional[1]
    Ns(default) = haskey(options, "N") ? parse.(Int, split(options["N"], ',')) : default
    @printf("Threads: %d\n", Threads.nthreads())
    BLAS.set_num_threads(1)

    if command == "consistency"
        RunConsistency(Ns([10, 20, 30, 40]))
    elseif command == "trajectories"
        RunTrajectories(Ns([10, 20, 30, 40, 50]);
                        trajectories = parse(Int, get(options, "trajectories", "16")),
                        time = parse(Float64, get(options, "time", "2000")))
    else
        println("usage: julia -t <threads> BHNumberConservingQuantumDynamics.jl consistency | trajectories [--N ...]")
    end
end
