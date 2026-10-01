# BHNumberConservingSteadyState.jl
#
# The steady state of the number-conserving dissipative Bose-Hubbard ring against the classical
# invariant measure on the attractor - item C7 of number-conserving-BH-paper-TODO.md (steady-state
# Husimi function on the population simplex) and the numbers A10 asks for (how close the two are,
# as a function of N).
#
# HUSIMI FUNCTION ON THE SIMPLEX.  With the SU(L) coherent states |ψ; N> the Husimi function is
# Q(ψ) = <ψ; N| ρ |ψ; N>.  Averaged over the L relative phases (the population simplex
# n_j = |ψ_j|^2 is what is plotted) every off-diagonal element of ρ in the Fock basis drops out,
# and what remains is exact and simple:
#
#     Q(n) = sum_a P(a) M(a; N, n),      P(a) = <a| ρ |a>,
#     M(a; N, n) = N! / prod_j a_j!  prod_j n_j^a_j            (the multinomial distribution)
#
# so only the POPULATIONS of the steady state are needed.  The classical counterpart is built with
# the SAME kernel: an ensemble of classical states m on the attractor (weighted by the basins,
# since the initial conditions are drawn uniformly on CP^(L-1)) is the mixture of the coherent
# states |m; N>, whose Fock populations are P_cl(a) = < M(a; N, m) >_m.  P and P_cl are then two
# distributions on the same d_N points and are compared directly:
#
#     total variation   TV = (1/2) sum_a |P(a) - P_cl(a)|
#     Hellinger         H  = sqrt(1 - sum_a sqrt(P(a) P_cl(a)))
#     coarse TV         TV_c, the same after binning n = a / N into cells of side 1/10
#
# TV and H compare the two distributions at the scale of the quantum fluctuations, 1/sqrt(N): they
# vanish only if the steady state IS the coherent-state-dressed invariant measure, including the
# shape of its O(1/sqrt(N)) spread, which the jump noise need not reproduce - so they may saturate.
# TV_c compares them at a FIXED resolution and must go to zero as N grows if the steady state
# concentrates on the classical invariant measure; it is the number for A10.  The mean occupation
# of the fullest site, <n_max>, is printed from P and from P_cl alike.
#
# THE STEADY STATE is the null vector of the q = 0 sector, found by inverse iteration on the sparse
# orbit-basis sector of BHNumberConservingLiouvillianSparse.jl, so N is not limited to the dense
# range: N = 24 needs a few GB, N = 30 about 12 GB.
#
# USAGE
#   julia -t 16 BHNumberConservingSteadyState.jl [--N 8,10,12,14,16,20,24] [--relaxation 8000]
#
# OUTPUT in $BH_RESULTS_DIR (default ~/results/bh/number-conserving/quantum/3/steady):
#   summary.txt                          point, N, TV, Hellinger, purity, quantum and classical <n_max>
#   husimi_g<g>_e<η>_N<N>.txt            the triangular grid: n1 n2 n3 Q_quantum Q_classical
#   populations_g<g>_e<η>_N<N>.txt       a1 a2 a3 P P_cl

include(joinpath(@__DIR__, "BHNumberConservingLiouvillianSparse.jl"))   # also the dense file

module Classical
    include(joinpath(@__DIR__, "BHNumberConserving.jl"))
end

using Printf
using Random
using Statistics

const STEADY_RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "quantum", "3", "steady"))

const POINTS = ((-20.0, 3.0), (-20.0, 1.0), (-4.0, 3.0))        # chaotic, multistable, regular
const GRID = 60                                                 # simplex grid: (GRID + 1)(GRID + 2)/2 points


""" Populations P(a) = <a|ρ_ss|a> of the steady state, from the null vector of the sparse q = 0
    sector.  An orbit O of length l carries the coefficient x_O N_O (L / l) on each of its pairs
    (BHNumberConservingLiouvillianSparse.jl), so the diagonal pair (a, a) gives ρ_aa directly.
    Returns P, the residual |M x| / |x| and the purity-free check that P is real and positive. """
function SteadyPopulations(p::LiouvillianParameters)
    basis = NumberBasis(p.L, p.N)
    orbits = OrbitBasis(basis, 0)
    M = SparseSector(p, 0; basis = basis, orbits = orbits)

    factor = lu(M - 1e-10 * I; control = FactorisationControl())
    x = randn(Xoshiro(1), ComplexF64, size(M, 1))
    for _ = 1:3                                 # inverse iteration: the gap separates λ = 0 by >~ 0.1
        x = factor \ x
        x ./= norm(x)
    end
    residual = norm(M * x)

    d = length(basis)
    P = zeros(ComplexF64, d)
    for a = 1:d
        linear = a + d * (a - 1)
        O = orbits.orbit[linear]
        l = round(Int, (orbits.normalisation[O] * p.L)^2)        # N_O = sqrt(l) / L
        P[a] = x[O] * orbits.normalisation[O] * p.L / l
    end
    P ./= sum(P)                                # fixes the arbitrary phase and scale of x

    imaginary = maximum(abs, imag.(P))
    negative = -min(0.0, minimum(real, P))
    return real.(P), basis, (residual = residual, imaginary = imaginary, negative = negative)
end


""" log of the multinomial M(a; N, n) for all Fock states at once. """
function LogMultinomial(basis::NumberBasis, n)
    logFactorial = cumsum(vcat(0.0, log.(1:basis.N)))
    return [logFactorial[basis.N + 1] - sum(logFactorial[a .+ 1]) +
            sum(a[j] == 0 ? 0.0 : a[j] * log(max(n[j], 1e-300)) for j in eachindex(a))
            for a in basis.states]
end

""" Classical Fock populations of the coherent-state mixture over the attractor(s): initial
    conditions uniform on CP^(L-1), relaxation, then samples every `step`. """
function ClassicalPopulations(basis, g, η; κ = 0.3, initial = 64, relaxation = 8000.0,
        window = 2000.0, step = 1.0)
    parameters = Classical.NumberConservingParameters(basis.L; J = 1.0, g = g, κ = κ, η = η)
    partial = [zeros(length(basis)) for _ = 1:initial]

    Threads.@threads for i = 1:initial
        work = Classical.NumberConservingWorkspace(basis.L)
        field!(du, u, q, t) = (Classical.UpdateWorkspace!(work, u, parameters);
                               Classical.VectorField!(du, u, parameters, work))
        ψ0 = Classical.RandomInitialCondition(basis.L, Xoshiro(hash((g, η, i, "steady"))))
        states = Classical.solve(Classical.ODEProblem(field!, ψ0, (0.0, relaxation + window)),
                                 Classical.DP8(); reltol = 1e-10, abstol = 1e-10,
                                 saveat = relaxation:step:(relaxation + window), maxiters = 10^8).u
        for ψ in states
            partial[i] .+= exp.(LogMultinomial(basis, abs2.(ψ) ./ sum(abs2, ψ)))
        end
        partial[i] ./= length(states)
    end

    Pcl = sum(partial) ./ initial
    return Pcl ./ sum(Pcl)
end

""" The Husimi function of a population vector on the triangular grid of the simplex. """
function HusimiGrid(basis, P)
    rows = []
    for i = 0:GRID, j = 0:(GRID - i)
        n = [i / GRID, j / GRID, (GRID - i - j) / GRID]
        push!(rows, (n..., sum(P .* exp.(LogMultinomial(basis, n)))))
    end
    return rows
end


function Run(Ns; relaxation = 8000.0)
    mkpath(STEADY_RESULTS)
    summary = joinpath(STEADY_RESULTS, "summary.txt")
    isfile(summary) || open(io -> println(io, "# g\teta\tN\tTV\tHellinger\tTV_coarse\tmax_n_quantum\tmax_n_classical\tresidual\tseconds"),
                            summary, "w")

    for (g, η) in POINTS, N in Ns
        time = @elapsed begin
            P, basis, check = SteadyPopulations(LiouvillianParameters(3, N; g = g, η = η, κ = 0.3))
            Pcl = ClassicalPopulations(basis, g, η; relaxation = relaxation)
        end
        (check.imaginary > 1e-8 || check.negative > 1e-8) &&
            @warn "steady-state populations not real and positive" g η N check

        tv = 0.5 * sum(abs, P .- Pcl)
        hellinger = sqrt(max(0.0, 1 - sum(sqrt.(P .* Pcl))))
        maximumQuantum = sum(P[a] * maximum(basis.states[a]) / N for a in eachindex(P))
        maximumClassical = sum(Pcl[a] * maximum(basis.states[a]) / N for a in eachindex(Pcl))

        # the same comparison at a FIXED resolution: n = a / N binned into cells of side 1/10
        coarse = Dict{Tuple{Int, Int}, Float64}()
        for (a, x, y) in zip(basis.states, P, Pcl)
            cell = (min(floor(Int, 10 * a[1] / N), 9), min(floor(Int, 10 * a[2] / N), 9))
            coarse[cell] = get(coarse, cell, 0.0) + x - y
        end
        tvCoarse = 0.5 * sum(abs, values(coarse))

        @printf("(g, η) = (%5.1f, %.1f), N = %2d: TV = %.4f, Hellinger = %.4f, coarse TV = %.4f, <n_max> quantum %.4f classical %.4f  (residual %.1e, %.1f s)\n",
                g, η, N, tv, hellinger, tvCoarse, maximumQuantum, maximumClassical, check.residual, time)
        open(io -> println(io, join([@sprintf("%.6g", v) for v in (g, η, N, tv, hellinger, tvCoarse,
                                                                   maximumQuantum, maximumClassical,
                                                                   check.residual, time)], "\t")),
             summary, "a")

        tag = @sprintf("g%.2f_e%.2f_N%02d", g, η, N)
        quantum = HusimiGrid(basis, P)
        classical = HusimiGrid(basis, Pcl)
        open(joinpath(STEADY_RESULTS, "husimi_$tag.txt"), "w") do io
            println(io, "# n1\tn2\tn3\tQ_quantum\tQ_classical")
            for (q, c) in zip(quantum, classical)
                @printf(io, "%.6f\t%.6f\t%.6f\t%.8e\t%.8e\n", q[1:3]..., q[4], c[4])
            end
        end
        open(joinpath(STEADY_RESULTS, "populations_$tag.txt"), "w") do io
            println(io, "# a1\ta2\ta3\tP\tP_cl")
            for (a, x, y) in zip(basis.states, P, Pcl)
                @printf(io, "%d\t%d\t%d\t%.10e\t%.10e\n", a..., x, y)
            end
        end
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    positional, options = SplitArguments(ARGS)
    Ns = haskey(options, "N") ? parse.(Int, split(options["N"], ',')) : [8, 10, 12, 14, 16, 20, 24]
    @printf("Threads: %d\n", Threads.nthreads())
    Run(Ns; relaxation = parse(Float64, get(options, "relaxation", "8000")))
end
