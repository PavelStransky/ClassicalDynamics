# BHNumberConservingLiouvillianModes.jl
#
# Which Liouvillian eigenmodes does an observable actually see?  Item C5 of
# number-conserving-BH-paper-TODO.md.  Any expectation value decomposes over the eigenmodes,
#
#     <O(t)> = sum_α exp(Λ_α t) Tr[O R_α] Tr[L_α^dag ρ(0)],        Tr[L_α^dag R_β] = δ_αβ,
#
# so a mode matters for O only if its OBSERVABLE WEIGHT Tr[O R_α] is appreciable.  Two numbers are
# written for every mode:
#
#   o_α = Tr[O R_α] / ||R_α||_1     the weight with R_α normalised in the TRACE norm.  Then
#                                   |o_α| <= ||O||_op, which is 1 for the three observables below,
#                                   whatever N - a threshold such as |o_α| > 0.05 means the same
#                                   thing at every N, which a Hilbert-Schmidt normalisation (whose
#                                   bound grows like sqrt(d_N)) would not.
#   c_α = Tr[O R_α] Tr[L_α^dag ρ0]  the full, normalisation-free amplitude of the mode in <O(t)>
#                                   after a quench from the coherent state ρ0 = |ψ0; N><ψ0; N| at a
#                                   generic point ψ0 (the one of the K1 test).  This needs the left
#                                   eigenvectors, obtained as the rows of the inverse of the matrix
#                                   of right eigenvectors (V \ r0).
#
# Observables (all with operator norm 1):  n̂_1 / N (it has a component in every momentum sector),
# the bond coherence Ĉ = (1/2N) sum_j (b_j^dag b_{j+1} + h.c.) and Σ_j n̂_j^2 / N^2 (translation
# invariant, so only the q = 0 sector carries them).
#
# USAGE   julia BHNumberConservingLiouvillianModes.jl [--N 8,10,12,14] [--blas 8]
#
# OUTPUT in $BH_RESULTS_DIR (default ~/results/bh/number-conserving/quantum/3/modes):
#   modes_g<g>_e<η>_N<N>_m<m>.txt   Re Λ, Im Λ, |o| and |c| for the three observables
#   summary.txt                     number of modes with |o_α| > 0.05, by point, N, sector and
#                                   observable - the count whose growth with N C5 asks for

ENV["CD_NO_PLOTS"] = "true"
isdefined(Main, :LiouvillianSector) || include(joinpath(@__DIR__, "BHNumberConservingLiouvillian.jl"))

using LinearAlgebra
using SparseArrays
using Printf

const MODES_RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "quantum", "3", "modes"))

const MODE_POINTS = ((-20.0, 3.0), (-4.0, 3.0), (-20.0, 0.0))       # chaotic and regular references
const THRESHOLD = 0.05
const QUENCH_STATE = normalize(ComplexF64[0.8 + 0.1im, 0.3 - 0.4im, 0.2 + 0.25im])


# ---------------------------------------------------------------------------------------------
# Between sector vectors and operators on the N-boson space
# ---------------------------------------------------------------------------------------------

""" Index ranges of the groups (k, k - m) of a sector, in the order LiouvillianSector uses. """
function SectorGroups(data::MomentumData, m)
    L = data.L
    groups = []
    offset = 0
    for k = 1:L
        rows = data.ranges[k]
        columns = data.ranges[mod1(k - m, L)]
        length(rows) * length(columns) == 0 && continue
        push!(groups, (rows = rows, columns = columns, range = (offset + 1):(offset + length(rows) * length(columns))))
        offset += length(rows) * length(columns)
    end
    return groups
end

""" The operator R on the N-boson space (Fock basis) of a vector of sector m.  Within a group the
    vector is the row-major flattening of the block X (the kron(A, conj(B)) convention of
    LiouvillianSector), and R = sum_k V_k X_k V_(k-m)^dag. """
function SectorOperator(data::MomentumData, m, vector)
    R = zeros(ComplexF64, size(data.V, 1), size(data.V, 1))
    for g in SectorGroups(data, m)
        X = permutedims(reshape(vector[g.range], length(g.columns), length(g.rows)))
        R .+= view(data.V, :, g.rows) * X * view(data.V, :, g.columns)'
    end
    return R
end

""" The component of an operator ρ in sector m, as a sector vector (inverse of SectorOperator). """
function SectorVector(data::MomentumData, m, ρ)
    ρMomentum = data.V' * ρ * data.V
    groups = SectorGroups(data, m)
    vector = zeros(ComplexF64, groups[end].range[end])
    for g in groups
        vector[g.range] = vec(permutedims(ρMomentum[g.rows, g.columns]))
    end
    return vector
end

""" Tr[O R] for a sparse O. """
function TraceProduct(O::SparseMatrixCSC, R)
    total = zero(ComplexF64)
    rows = rowvals(O)
    for column = 1:size(O, 2), index in nzrange(O, column)
        total += nonzeros(O)[index] * R[column, rows[index]]
    end
    return total
end

function Observables(basis::NumberBasis)
    N, L = basis.N, basis.L
    n1 = Hop(basis, 1, 1) ./ N
    coherence = sum(Hop(basis, j, mod1(j + 1, L)) + Hop(basis, mod1(j + 1, L), j) for j = 1:L) ./ (2 * N)
    squares = sum(Hop(basis, j, j)^2 for j = 1:L) ./ N^2
    return (n1 = n1, C = coherence, n2 = squares)
end

function CoherentProjector(basis::NumberBasis, ψ)
    logFactorial = cumsum(vcat(0.0, log.(1:basis.N)))
    c = [exp(0.5 * (logFactorial[basis.N + 1] - sum(logFactorial[a .+ 1]))) *
         prod(ψ[j]^a[j] for j in eachindex(a)) for a in basis.states]
    c ./= norm(c)
    return c * c'
end


# ---------------------------------------------------------------------------------------------
# The calculation
# ---------------------------------------------------------------------------------------------

""" Every mode of sector m at one point: eigenvalues, observable weights and quench amplitudes. """
function Modes(g, η, N, m; κ = 0.3)
    p = LiouvillianParameters(3, N; g = g, η = η, κ = κ)
    basis = NumberBasis(3, N)
    data = MomentumData(basis, p)
    observables = Observables(basis)

    decomposition = eigen(LiouvillianSector(data, m))
    r0 = SectorVector(data, m, CoherentProjector(basis, QUENCH_STATE))
    overlaps = decomposition.vectors \ r0                       # Tr[L_α^dag ρ0] for these R_α

    rows = []
    for α in eachindex(decomposition.values)
        R = SectorOperator(data, m, decomposition.vectors[:, α])
        traceNorm = sum(svdvals(R))
        weights = [TraceProduct(O, R) for O in (observables.n1, observables.C, observables.n2)]
        push!(rows, (real(decomposition.values[α]), imag(decomposition.values[α]),
                     (abs.(weights) ./ traceNorm)..., abs.(weights .* overlaps[α])...))
    end
    return rows
end

""" Consistency of the conversions at a small N: ρ0 reassembled from its sector components, and
    every R_α an eigen-operator of the full Liouvillian. """
function ModesChecks(; N = 5)
    p = LiouvillianParameters(3, N; g = -20.0, η = 3.0, κ = 0.3)
    basis = NumberBasis(3, N)
    data = MomentumData(basis, p)
    ρ0 = CoherentProjector(basis, QUENCH_STATE)
    rebuilt = sum(SectorOperator(data, m, SectorVector(data, m, ρ0)) for m = 0:2)
    @printf("ρ0 reassembled from its three sectors: %.1e\n", norm(rebuilt - ρ0))

    full = LiouvillianFull(basis, p)
    d = length(basis)
    worst = 0.0
    for m = 0:2
        decomposition = eigen(LiouvillianSector(data, m))
        for α in eachindex(decomposition.values)
            R = SectorOperator(data, m, decomposition.vectors[:, α])
            image = reshape(full * vec(R), d, d)
            worst = max(worst, norm(image - decomposition.values[α] * R) / norm(R))
        end
    end
    @printf("|L(R_α) - Λ_α R_α| / |R_α| over every mode of every sector: %.1e\n", worst)
end


function Run(Ns)
    mkpath(MODES_RESULTS)
    summary = joinpath(MODES_RESULTS, "summary.txt")
    isfile(summary) || open(io -> println(io, "# g\teta\tN\tm\tmodes\tcount_n1\tcount_C\tcount_n2\tslowest_heavy_n1\tseconds"),
                            summary, "w")

    for (g, η) in MODE_POINTS, N in Ns, m in (0, 1)
        file = joinpath(MODES_RESULTS, @sprintf("modes_g%.2f_e%.2f_N%02d_m%d.txt", g, η, N, m))
        isfile(file) && continue

        time = @elapsed rows = Modes(g, η, N, m)
        open(file, "w") do io
            println(io, "# re\tim\to_n1\to_C\to_n2\tc_n1\tc_C\tc_n2")
            for r in rows
                println(io, join([@sprintf("%.8e", v) for v in r], "\t"))
            end
        end

        counts = [count(r -> r[2 + k] > THRESHOLD, rows) for k = 1:3]
        heavy = [r[1] for r in rows if r[3] > THRESHOLD]
        slowest = isempty(heavy) ? NaN : minimum(abs, heavy)
        @printf("(g, η) = (%5.1f, %.1f), N = %2d, m = %d: %5d modes, |o| > %.2f for n1/N %4d, C %4d, Σn²/N² %4d  (%.1f s)\n",
                g, η, N, m, length(rows), THRESHOLD, counts..., time)
        open(io -> println(io, join([@sprintf("%.6g", v) for v in (g, η, N, m, length(rows), counts..., slowest, time)], "\t")),
             summary, "a")
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    options = Dict{String, String}()
    for i = 1:2:(length(ARGS) - 1)
        options[replace(ARGS[i], "--" => "")] = ARGS[i + 1]
    end
    BLAS.set_num_threads(parse(Int, get(options, "blas", "8")))
    ModesChecks()
    Run(haskey(options, "N") ? parse.(Int, split(options["N"], ',')) : [8, 10, 12, 14])
end
