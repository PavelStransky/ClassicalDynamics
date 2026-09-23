# BHNumberConservingLiouvillian.jl
#
# Liouvillian spectrum, spectral statistics and the nearest-neighbour spacing distribution of the
# NUMBER-CONSERVING dissipative Bose-Hubbard model of number-conserving-BH.md - the quantum side of
# BHNumberConserving.jl, which does the classical limit of the same model.
#
#     drho/dt = -i [H, rho] + sum_{i != j} Gamma_ij D[b_i^dag b_j] rho + kappa_c sum_j D[c_j] rho
#     H   = -J sum_j (b_j^dag b_{j+1} + h.c.) + (U/2) sum_j n_j (n_j - 1)
#     c_j = (b_j^dag + b_{j+1}^dag)(b_j - b_{j+1})
#
# PARAMETERS ARE THE CLASSICAL ONES.  The note fixes the semiclassical scaling in §3: g = U N,
# eta = N (Gamma_+ - Gamma_-), kappa = N kappa_c, with hbar_eff = 1/N.  This file takes
# (J, g, eta, kappa) and N and divides by N internally, so a sweep in N at fixed (J, g, eta, kappa)
# is the genuine hbar_eff -> 0 sequence of §7 - and the classical attractor that BHNumberConserving.jl
# finds at those very numbers is what the Liouvillian spectrum should be compared against.  There is
# no cutoff anywhere: the N-boson Hilbert space is finite and exact, d_N = C(N + L - 1, L - 1).
#
# SYMMETRIES - the whole point of the construction (note §6)
#
#   * STRONG U(1).  [N, H] = 0 and [N, L_k] = 0 for EVERY jump operator, so the Liouvillian does not
#     merely conserve the particle number on average: it is block diagonal in the pair (N, N') of
#     bra and ket particle numbers.  Working in the N-boson Hilbert space from the start therefore
#     already resolves this symmetry - the matrices built below ARE the (N, N) block, the one that
#     hosts the steady state.  Nothing is projected out by hand and nothing is truncated.
#     (The coherence blocks N != N' exist too and decay; they are not built here.)
#
#   * WEAK Z_L.  For uniform rates the translation T (T b_j^dag T^dag = b_{j+1}^dag) commutes with H
#     and permutes the jump operators among themselves, so the SUPEROPERATOR  rho -> T rho T^dag
#     commutes with the Liouvillian.  It is weak, not strong: individual jump operators are not
#     translation invariant, only the set is.  Its eigenvalues are exp(i q), q = 2 pi m / L, and a
#     state |k, a><k', b| sits in the sector q = k - k'.  The (N, N) block of dimension d_N^2 thus
#     splits into L sectors of about d_N^2 / L each - the last column of the table in §7.
#
#   * rho -> rho^dag is ANTIUNITARY and does NOT block diagonalise.  It maps sector q to sector -q,
#     so it constrains a sector only when q = -q mod 2 pi, i.e. q = 0, and q = pi for even L.  Those
#     SELF-CONJUGATE sectors are similar to a real matrix (split them into Hermitian and
#     anti-Hermitian parts), so a finite fraction of their eigenvalues is exactly real and their
#     spacing statistics are of the real-Ginibre kind.  Every other sector carries no antiunitary
#     constraint at all and is the clean place to look for Ginibre (class A) statistics.  This is
#     the one symmetry that is easy to satisfy oneself is "resolved" and is not: mixing a
#     self-conjugate sector into the histogram puts a spike of real eigenvalues into it.
#
#   * REFLECTION j -> -j exists only at eta = 0, where it reverses no directed hop because there is
#     none (note §6).  It then combines with Z_L into the dihedral group and, composed with the
#     adjoint, gives every sector an antiunitary symmetry, so the reference ensemble is no longer
#     Ginibre.  The code detects this and says so loudly rather than resolving it: eta = 0 is the
#     case where BHNumberConserving.jl finds no classical chaos either, so it is the control, not
#     the case of interest.  For eta != 0 the decomposition here is complete.
#
# WHY NEAREST-NEIGHBOUR SPACINGS NEED CARE HERE.  The spectrum is complex, so there is no ordering
# and no one-dimensional unfolding.  The distances are unfolded locally against the density
# estimated from the k-th neighbour, and only the densest part of the cloud is kept, because the
# edge of the spectrum has its own statistics.  The reference curves are the 2D Poisson (Rayleigh)
# law and the Grobe-Haake-Sommers law for the Ginibre ensemble, which has CUBIC level repulsion
# P(s) ~ s^3 - the dissipative counterpart of Wigner-Dyson.  The complex spacing RATIO of
# Sa, Ribeiro and Prosen is computed as well: it needs no unfolding at all and is the more reliable
# of the two, so the two should agree before anything is concluded.

using LinearAlgebra
using SparseArrays
using Random
using Statistics
using Printf

const BHL_LOAD_PLOTS = get(ENV, "CD_NO_PLOTS", "false") != "true"

if BHL_LOAD_PLOTS
    using Plots
    pyplot(size = (900, 700))
end


# ---------------------------------------------------------------------------------------------
# Fock basis of the N-boson sector - the strong U(1) of §6, resolved by construction
# ---------------------------------------------------------------------------------------------

""" Occupation-number basis of L sites at fixed total number N, d_N = C(N + L - 1, L - 1) states.
    `index` maps an occupation vector to its position, which is all the operators below need. """
struct NumberBasis
    L::Int
    N::Int
    states::Vector{Vector{Int}}
    index::Dict{Vector{Int}, Int}
end

function NumberBasis(L::Integer, N::Integer)
    states = Vector{Int}[]

    function Compositions!(prefix, remaining, sites)
        if sites == 1
            push!(states, vcat(prefix, remaining))
            return
        end
        for n = 0:remaining
            Compositions!(vcat(prefix, n), remaining - n, sites - 1)
        end
    end

    Compositions!(Int[], N, L)

    return NumberBasis(L, N, states, Dict(s => i for (i, s) in enumerate(states)))
end

Base.length(basis::NumberBasis) = length(basis.states)


""" b_i^dag b_j on the N-boson sector - the only operator anything below is built from, because
    every admissible jump operator of §1 is a combination of such pairs (note §2). """
function Hop(basis::NumberBasis, i::Integer, j::Integer)
    d = length(basis)
    rows = Int[]; cols = Int[]; values = ComplexF64[]

    for (column, state) in enumerate(basis.states)
        state[j] == 0 && continue

        if i == j
            push!(rows, column); push!(cols, column); push!(values, state[j])
            continue
        end

        target = copy(state)
        target[j] -= 1
        target[i] += 1

        push!(rows, basis.index[target]); push!(cols, column)
        push!(values, sqrt(state[j] * (state[i] + 1)))
    end

    return sparse(rows, cols, values, d, d)
end


""" Translation T: T b_j^dag T^dag = b_{j+1}^dag, i.e. the occupations are shifted round the ring.
    A permutation matrix; the weak Z_L symmetry of §6 is built from it. """
function TranslationOperator(basis::NumberBasis)
    d = length(basis)
    rows = [basis.index[circshift(state, 1)] for state in basis.states]

    return sparse(rows, collect(1:d), ones(ComplexF64, d), d, d)
end


# ---------------------------------------------------------------------------------------------
# Model: Hamiltonian, jump operators, effective non-Hermitian generator
# ---------------------------------------------------------------------------------------------

""" Parameters in the CLASSICAL normalisation of note §3, plus the boson number N that sets
    hbar_eff = 1/N.  The rates actually entering the Lindblad equation are U = g/N,
    Gamma_+ - Gamma_- = eta/N and kappa_c = kappa/N.

    `Γsym` is a symmetric incoherent hopping rate added to BOTH directions.  It cancels exactly in
    the mean-field limit (note §2: only the antisymmetric part A_ij survives), so it is invisible
    to BHNumberConserving.jl and is a purely quantum knob - at an unscaled rate it is plain
    decoherence.  `γd` is the dephasing D[n_j] of §1, also invisible to the deterministic mean
    field.  Both default to 0, which keeps the quantum model the exact counterpart of the
    classical one. """
struct LiouvillianParameters
    L::Int
    N::Int
    J::Float64
    g::Float64
    η::Float64
    κ::Float64
    Γsym::Float64
    γd::Float64
end

LiouvillianParameters(L, N; J = 1.0, g = -20.0, η = 3.0, κ = 0.3, Γsym = 0.0, γd = 0.0) =
    LiouvillianParameters(L, N, float(J), float(g), float(η), float(κ), float(Γsym), float(γd))


""" H = -J sum_j (b_j^dag b_{j+1} + h.c.) + (U/2) sum_j n_j (n_j - 1),  U = g/N. """
function Hamiltonian(basis::NumberBasis, p::LiouvillianParameters)
    d = length(basis)
    L = basis.L
    U = p.g / p.N

    H = spzeros(ComplexF64, d, d)

    for j = 1:L
        k = j == L ? 1 : j + 1
        H -= p.J * (Hop(basis, j, k) + Hop(basis, k, j))
    end

    for (column, state) in enumerate(basis.states)
        H[column, column] += 0.5 * U * sum(n * (n - 1) for n in state)
    end

    return H
end


""" The bond phase-locking operator c_jk = (b_j^dag + b_k^dag)(b_j - b_k)
    = n_j - n_k - b_j^dag b_k + b_k^dag b_j (Diehl et al. 2008; note §1).  Its dark state in each
    N sector is the uniform condensate. """
BondOperator(basis::NumberBasis, j, k) =
    Hop(basis, j, j) - Hop(basis, k, k) - Hop(basis, j, k) + Hop(basis, k, j)


""" All jump operators, already carrying their rates.  The SET is closed under translation, which
    is what makes Z_L a weak symmetry of the Liouvillian even though no single operator here is
    translation invariant. """
function JumpOperators(basis::NumberBasis, p::LiouvillianParameters)
    L = basis.L
    Γplus = max(p.η, 0.0) / p.N + p.Γsym          # hop j -> j+1
    Γminus = max(-p.η, 0.0) / p.N + p.Γsym        # hop j+1 -> j;  N (Γplus - Γminus) = η
    κc = p.κ / p.N

    operators = SparseMatrixCSC{ComplexF64, Int}[]

    for j = 1:L
        k = j == L ? 1 : j + 1
        Γplus > 0 && push!(operators, sqrt(Γplus) * Hop(basis, k, j))
        Γminus > 0 && push!(operators, sqrt(Γminus) * Hop(basis, j, k))
        κc > 0 && push!(operators, sqrt(κc) * BondOperator(basis, j, k))
    end

    if p.γd > 0
        for j = 1:L
            push!(operators, sqrt(p.γd) * Hop(basis, j, j))
        end
    end

    return operators
end


""" K = -i H - (1/2) sum_k L_k^dag L_k, the non-Hermitian generator.  The Liouvillian is then
    L[rho] = K rho + rho K^dag + sum_k L_k rho L_k^dag, which is the form every block below uses.
    K commutes with the translation even though the individual L_k do not, because sum_k L_k^dag L_k
    is translation invariant - that is exactly the statement that Z_L is a weak symmetry. """
function EffectiveGenerator(H, jumps)
    K = -im * H
    for Lk in jumps
        K -= 0.5 * (Lk' * Lk)
    end
    return K
end


# ---------------------------------------------------------------------------------------------
# Weak Z_L: momentum basis and the Liouvillian sectors (note §6)
# ---------------------------------------------------------------------------------------------

""" Everything the sector assembly needs, computed once: the momentum basis V (columns grouped by
    momentum k_m = 2 pi m / L), and K and the jump operators expressed in it.

    The basis comes from the projectors P_m = (1/L) sum_r exp(-i k_m r) T^r, whose range is the
    momentum-m eigenspace of T.  The subspace dimensions d_m sum to d_N and are the reason a sector
    is about d_N^2 / L and not d_N^2. """
struct MomentumData
    L::Int
    ranges::Vector{UnitRange{Int}}
    momenta::Vector{Float64}
    K::Matrix{ComplexF64}
    jumps::Vector{Matrix{ComplexF64}}
end

function MomentumData(basis::NumberBasis, p::LiouvillianParameters)
    d = length(basis)
    L = basis.L
    T = TranslationOperator(basis)

    powers = Vector{SparseMatrixCSC{ComplexF64, Int}}(undef, L)
    powers[1] = sparse(ComplexF64(1) * I, d, d)
    for r = 2:L
        powers[r] = powers[r - 1] * T
    end

    V = zeros(ComplexF64, d, d)
    ranges = UnitRange{Int}[]
    momenta = Float64[]
    offset = 0

    for m = 0:(L - 1)
        k = 2 * pi * m / L
        projector = zeros(ComplexF64, d, d)
        for r = 0:(L - 1)
            projector .+= exp(-im * k * r) .* powers[r + 1]
        end
        projector ./= L

        dimension = round(Int, real(tr(projector)))
        if dimension > 0
            V[:, (offset + 1):(offset + dimension)] = svd(projector).U[:, 1:dimension]
        end

        push!(ranges, (offset + 1):(offset + dimension))
        push!(momenta, k)
        offset += dimension
    end

    offset == d || error("momentum subspaces do not exhaust the Hilbert space: $offset of $d")
    norm(V' * V - I) < 1e-10 * d || error("momentum basis is not orthonormal")

    H = Hamiltonian(basis, p)
    jumps = JumpOperators(basis, p)
    K = EffectiveGenerator(H, jumps)

    return MomentumData(L, ranges, momenta, Matrix(V' * K * V),
                        [Matrix(V' * Lk * V) for Lk in jumps])
end


""" True when the sector q = 2 pi m / L is mapped to itself by rho -> rho^dag, i.e. when q = -q.
    Only m = 0 and, for even L, m = L/2.  See the header: those sectors are similar to a real
    matrix and their spacing statistics are NOT the Ginibre ones. """
SelfConjugateSector(L, m) = (mod(2 * m, L) == 0)


""" The reflection j -> -j of note §6 is a symmetry of the ring only when the circulation vanishes:
    it reverses every directed hop, so it maps the jump set onto itself iff Gamma_+ = Gamma_-, i.e.
    eta = 0.  (The bond operators c_j and the dephasing n_j are reflection covariant anyway, and a
    symmetric background Gamma_sym does not spoil it.)  When it holds, the Z_L decomposition below
    is no longer the whole story - see the warning printed by Analyse. """
ReflectionIsSymmetry(p::LiouvillianParameters) = p.η == 0


""" The Liouvillian restricted to the Z_L sector q = 2 pi m / L, as a dense matrix.

    A basis operator |k, a><k', b| has  T rho T^dag = exp(i (k - k')) rho, so the sector m collects
    the pairs whose row momentum exceeds the column momentum by m.  Within the sector the pairs are
    grouped by their ROW momentum k; the group (k, k - q) has dimension d_k d_{k-q}, and an
    operator acting as A on the row index and B on the column index contributes kron(A, conj(B)).
    With K momentum diagonal, K contributes only to the diagonal groups, while the jump operators -
    which are NOT translation invariant individually - connect different groups.

    The full (N, N) block is never formed: it would be d_N^2 x d_N^2, whereas this is about
    (d_N^2 / L)^2, the entry in the last column of the table in note §7. """
function LiouvillianSector(data::MomentumData, m::Integer)
    L = data.L
    rowRange(k) = data.ranges[k]
    columnRange(k) = data.ranges[mod1(k - m, L)]

    dimensions = [length(rowRange(k)) * length(columnRange(k)) for k = 1:L]
    offsets = cumsum(vcat(0, dimensions))
    total = offsets[end]

    block = zeros(ComplexF64, total, total)

    for k1 = 1:L, k2 = 1:L
        dimensions[k1] == 0 && continue
        dimensions[k2] == 0 && continue

        piece = zeros(ComplexF64, dimensions[k1], dimensions[k2])

        for Lk in data.jumps
            piece .+= kron(view(Lk, rowRange(k1), rowRange(k2)),
                           conj(view(Lk, columnRange(k1), columnRange(k2))))
        end

        if k1 == k2
            rows = length(rowRange(k1))
            columns = length(columnRange(k1))
            piece .+= kron(view(data.K, rowRange(k1), rowRange(k1)),
                           Matrix{ComplexF64}(I, columns, columns))
            piece .+= kron(Matrix{ComplexF64}(I, rows, rows),
                           conj(view(data.K, columnRange(k1), columnRange(k1))))
        end

        block[(offsets[k1] + 1):offsets[k1 + 1], (offsets[k2] + 1):offsets[k2 + 1]] = piece
    end

    return block
end


""" The whole (N, N) block in one piece, WITHOUT the Z_L decomposition: 
    L = I (x) K + conj(K) (x) I + sum_k conj(L_k) (x) L_k  in the column-stacking convention
    vec(A rho B) = kron(transpose(B), A) vec(rho).  Only for the checks and for tiny systems -
    it is d_N^2 x d_N^2 and the point of the sectors is not to build it. """
function LiouvillianFull(basis::NumberBasis, p::LiouvillianParameters)
    d = length(basis)
    identity = sparse(ComplexF64(1) * I, d, d)

    H = Hamiltonian(basis, p)
    jumps = JumpOperators(basis, p)
    K = EffectiveGenerator(H, jumps)

    superoperator = kron(identity, K) + kron(conj(K), identity)
    for Lk in jumps
        superoperator += kron(conj(Lk), Lk)
    end

    return superoperator
end


""" Eigenvalues of every Z_L sector.  Returns a vector of (m, q, eigenvalues, selfConjugate),
    one entry per sector, with the eigenvalues sorted by decreasing real part - so the first one of
    sector 0 is the steady state (exactly 0) and the next is the Liouvillian gap.

    With `useConjugation` (the default) only the sectors m = 0 ... L/2 are diagonalised and the rest
    are obtained by complex conjugation, since rho -> rho^dag maps sector q onto sector -q (verified
    to 1e-13 in Checks()).  That halves the work at L = 3 and is what makes N = 16 affordable;
    pass false to diagonalise every sector independently, which is the slower cross-check. """
function SectorSpectra(p::LiouvillianParameters; verbose = true, useConjugation = true)
    basis = NumberBasis(p.L, p.N)
    data = MomentumData(basis, p)

    spectra = Vector{Vector{ComplexF64}}(undef, p.L)
    computed = useConjugation ? (0:div(p.L, 2)) : (0:(p.L - 1))

    for m in computed
        time = @elapsed values = eigvals(LiouvillianSector(data, m))
        spectra[m + 1] = sort(values, by = z -> (-real(z), imag(z)))

        if verbose
            @printf("  sector m = %d (q = %.4f): %5d eigenvalues%s, %.1f s
",
                    m, 2 * pi * m / p.L, length(values),
                    SelfConjugateSector(p.L, m) ? ", SELF-CONJUGATE" : "", time)
        end
    end

    for m = 0:(p.L - 1)
        isassigned(spectra, m + 1) && continue
        spectra[m + 1] = sort(conj.(spectra[mod(p.L - m, p.L) + 1]), by = z -> (-real(z), imag(z)))
        verbose && @printf("  sector m = %d (q = %.4f): %5d eigenvalues, by conjugation of sector %d
",
                           m, 2 * pi * m / p.L, length(spectra[m + 1]), mod(p.L - m, p.L))
    end

    return [(m = m, q = 2 * pi * m / p.L, values = spectra[m + 1],
             selfConjugate = SelfConjugateSector(p.L, m)) for m = 0:(p.L - 1)]
end


# ---------------------------------------------------------------------------------------------
# Spectral statistics of a COMPLEX spectrum
# ---------------------------------------------------------------------------------------------

""" For every eigenvalue: the distances to its `neighbours` nearest partners, and the indices of
    the nearest and the next-to-nearest.  Computed row by row, so the full n x n distance matrix is
    never formed - the sectors reach tens of thousands of eigenvalues.  Both statistics below read
    this one pass. """
function NeighbourData(values, neighbours::Integer)
    n = length(values)
    neighbours < n || error("asking for $neighbours neighbours out of $n eigenvalues")

    distances = zeros(Float64, n, neighbours)
    nearest = zeros(Int, n)
    next = zeros(Int, n)
    buffer = zeros(Float64, n)
    order = zeros(Int, n)

    for i = 1:n
        @inbounds for j = 1:n
            buffer[j] = abs(values[i] - values[j])
        end
        buffer[i] = Inf                                  # exclude the eigenvalue itself

        @inbounds for j = 1:n
            order[j] = j
        end
        partialsort!(order, 1:neighbours; by = j -> buffer[j])

        @inbounds for k = 1:neighbours
            distances[i, k] = buffer[order[k]]
        end
        nearest[i] = order[1]
        next[i] = order[2]
    end

    return distances, nearest, next
end


""" Unfolded nearest-neighbour spacings of a complex spectrum.

    A complex spectrum cannot be ordered, so the one-dimensional unfolding by a smooth staircase is
    not available.  Instead the local density is estimated from the k-th neighbour,
    rho_i = k / (pi d_k(i)^2), and the spacing is measured in units of rho_i^{-1/2}.  Only the
    `keep` fraction with the HIGHEST local density is retained, because the edge of the eigenvalue
    cloud has its own statistics and would otherwise contaminate the histogram; `neighbours` has to
    be large enough to average over the fluctuations and small enough that the density is still
    local (20-50 works over the whole range used here).

    The steady state is removed first: it is an exact zero forced by trace preservation, not a
    member of the statistical ensemble. """
function NearestNeighbourSpacings(values; neighbours = 30, keep = 0.5, window = nothing,
        zeroTolerance = 1e-8)

    values = filter(z -> abs(z) > zeroTolerance, values)
    length(values) > 4 * neighbours || error("too few eigenvalues for these settings")

    distances, _, _ = NeighbourData(values, neighbours)
    scale = distances[:, neighbours] ./ sqrt(neighbours / pi)       # 1 / sqrt(local density)
    spacings = distances[:, 1] ./ scale

    selected = SelectReferences(spacings, values, scale, keep, window)
    return selected ./ mean(selected)
end


""" Which eigenvalues enter a statistic.  Two ways of choosing, and they answer different questions.

    `window = nothing` keeps the `keep` fraction with the highest local density: the BULK of the
    eigenvalue cloud, away from its edge, which is the standard choice when one wants the generic
    statistics of the operator.

    `window = -c` keeps instead every eigenvalue with Re lambda >= -c: the spectral region CLOSE TO
    THE STEADY STATE.  Corps and Relano (arXiv:2609.18464) show that this is the region that
    matters physically - observables decompose as <O(t)> = sum_alpha exp(lambda_alpha t) Tr[O R_alpha]
    Tr[L_alpha^dag rho(0)], and the weights Tr[O R_alpha] are negligible for very negative
    Re lambda, so fast-decaying modes carry no observable consequence.  In their model the full
    spectrum reaches Ginibre statistics at a coupling well below the value at which the restricted
    region does, i.e. using the whole cloud OVERESTIMATES chaos.  `keep` is ignored in this mode.

    In both cases the neighbours are found in the FULL spectrum and only the REFERENCE eigenvalues
    are restricted.  Filtering first would put an artificial edge into the point set, and the points
    beside it would report spuriously large spacings. """
function SelectReferences(quantity, values, scale, keep, window)
    if isnothing(window)
        order = sortperm(scale)                                     # smallest scale = densest
        return quantity[order[1:max(1, round(Int, keep * length(quantity)))]]
    end

    selected = quantity[real.(values) .>= window]
    isempty(selected) && error("the spectral window Re λ >= $window holds no eigenvalues")

    return selected
end


""" Complex spacing ratios of Sa, Ribeiro and Prosen (PRX 10, 021019 (2020)):

        z_i = (lambda_i^NN - lambda_i) / (lambda_i^NNN - lambda_i)

    with the nearest and next-to-nearest neighbour in the complex plane.  They live in the unit
    disc and need NO unfolding, which is why they are the more trustworthy of the two statistics
    here - the local unfolding above is a fit, this is not.  Marker values:

        2D Poisson   <|z|> = 2/3,      -<cos arg z> = 0
        Ginibre (A)  <|z|> = 0.7378,   -<cos arg z> = 0.2405

    A self-conjugate sector (q = 0, or q = pi for even L) belongs to class AI-dagger instead, whose
    values differ from the Ginibre ones - see the reference above; do not read it against the
    Ginibre numbers. """
function ComplexSpacingRatios(values; neighbours = 30, keep = 0.5, window = nothing,
        zeroTolerance = 1e-8)

    values = filter(z -> abs(z) > zeroTolerance, values)
    n = length(values)
    n > 4 * neighbours || error("too few eigenvalues for these settings")

    distances, nearest, next = NeighbourData(values, neighbours)
    ratios = [(values[nearest[i]] - values[i]) / (values[next[i]] - values[i]) for i = 1:n]

    return SelectReferences(ratios, values, distances[:, neighbours], keep, window)
end


""" Nearest-neighbour spacing density of UNCORRELATED points in the plane, normalised to <s> = 1:
    a Rayleigh law, P(s) = (pi/2) s exp(-pi s^2 / 4).  The dissipative analogue of the Poisson law,
    and what an unresolved symmetry also produces - a superposition of independent spectra looks
    uncorrelated however chaotic each piece is.  That is why §6 has to be taken seriously before
    the histogram below means anything. """
Poisson2D(s) = 0.5 * pi * s * exp(-0.25 * pi * s^2)


""" Nearest-neighbour spacing density of the Ginibre ensemble, Grobe, Haake and Sommers
    (PRL 61, 1899 (1988)).  Conditioned on an eigenvalue at the origin the moduli of the others are
    independent with Gamma(n+1, 1) distributions, so the gap probability factorises:

        E(s) = prod_{n>=1} exp(-s^2) e_n(s^2),      e_n(x) = sum_{m=0}^{n} x^m / m!
        p(s) = E(s) sum_{n>=1} 2 s^(2n+1) / (n! e_n(s^2))

    The n = 1 factor gives p(s) ~ 2 s^3: CUBIC level repulsion, the signature of a complex spectrum
    with Ginibre correlations, as against the linear repulsion of the Poisson law above and the
    linear/quadratic/quartic of Wigner-Dyson on the real line.  `p` has unit norm but not unit mean,
    so the returned density is the rescaled P(s) = c p(c s) with c = <s>_p = 1.1429. """
function GinibreSpacingUnscaled(s::Real; terms = 120)
    s <= 0 && return 0.0
    x = s * s

    logProduct = 0.0
    sum = 0.0
    e = 1.0                                      # e_0(x)
    term = 1.0                                   # x^n / n!

    for n = 1:terms
        term *= x / n                            # x^n / n!, carried along instead of formed
        e += term                                # e_n(x)
        logProduct += -x + log(e)
        sum += 2 * s * term / e                  # 2 s^(2n+1) / (n! e_n(x)), overflow-free
    end

    return exp(logProduct) * sum
end

const GINIBRE_SCALE = let grid = range(0, 6, length = 6001)
    values = [GinibreSpacingUnscaled(s) for s in grid]
    step = grid[2] - grid[1]
    sum(grid .* values) * step                   # <s> of the unscaled density
end

""" The Ginibre law normalised to <s> = 1, ready to plot against a histogram. """
GinibreSpacing(s) = GINIBRE_SCALE * GinibreSpacingUnscaled(GINIBRE_SCALE * s)


""" var(s) and P(s < 0.5) of the two reference laws, computed FROM the laws above rather than
    quoted, so the table can never drift away from the curves drawn next to the histogram.  The
    ratio markers cannot be had that way and are the published values: Ginibre (class A)
    <|z|> = 0.7378 and -<cos arg z> = 0.2405, 2D Poisson 2/3 and 0
    (Sa, Ribeiro, Prosen, PRX 10, 021019 (2020)). """
const REFERENCE = let grid = range(0, 8, length = 16001), step = 8 / 16000
    ginibre = GinibreSpacing.(grid)
    poisson = Poisson2D.(grid)
    inside = grid .<= 0.5

    (ginibreVariance = sum((grid .^ 2) .* ginibre) * step - 1,
     ginibreSmall = sum(ginibre[inside]) * step,
     ginibreRadius = 0.7378, ginibreCosine = 0.2405,
     poissonVariance = sum((grid .^ 2) .* poisson) * step - 1,
     poissonSmall = sum(poisson[inside]) * step,
     poissonRadius = 2 / 3, poissonCosine = 0.0)
end


""" One-line verdict for a set of spacings and ratios, against the two reference ensembles. """
function StatisticsSummary(spacings, ratios)
    return (meanSpacing = mean(spacings),
            varianceSpacing = var(spacings),
            smallSpacings = count(s -> s < 0.5, spacings) / length(spacings),
            meanRadius = mean(abs, ratios),
            minusCosine = -mean(cos ∘ angle, ratios),
            count = length(spacings))
end


# ---------------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------------

""" Four panels: the spectrum in the complex plane sector by sector, the nearest-neighbour spacing
    histogram against the Ginibre and Poisson laws, the same as a cumulative distribution (where
    the eye is not fooled by the binning), and the complex spacing ratios in the unit disc. """
function PlotSpectralStatistics(sectors, spacings, ratios, p::LiouvillianParameters;
        savePath = nothing, title = "", showFigure = true)

    BHL_LOAD_PLOTS || return nothing

    spectrum = plot(xlabel = "Re λ", ylabel = "Im λ", title = "Liouvillian spectrum", legend = :best)
    for s in sectors
        scatter!(spectrum, real.(s.values), imag.(s.values), ms = 1.2, msw = 0, alpha = 0.55,
                 label = @sprintf("q = 2π·%d/%d%s", s.m, p.L, s.selfConjugate ? " (self-conj.)" : ""))
    end

    grid = range(0, 3.5, length = 400)

    histogramPanel = histogram(spacings, bins = range(0, 3.5, length = 46), normalize = :pdf,
        label = "Liouvillian", color = :steelblue, lw = 0, alpha = 0.75,
        xlabel = "s", ylabel = "P(s)", title = "nearest-neighbour spacings")
    plot!(histogramPanel, grid, GinibreSpacing.(grid), lw = 2.5, color = :crimson, label = "Ginibre")
    plot!(histogramPanel, grid, Poisson2D.(grid), lw = 2.5, ls = :dash, color = :black,
          label = "2D Poisson")

    sorted = sort(spacings)
    cumulative = plot(sorted, (1:length(sorted)) ./ length(sorted), lw = 2.5, color = :steelblue,
        label = "Liouvillian", xlabel = "s", ylabel = "I(s)", title = "cumulative", legend = :bottomright)
    step = grid[2] - grid[1]
    plot!(cumulative, grid, cumsum(GinibreSpacing.(grid)) .* step, lw = 2, color = :crimson,
          label = "Ginibre")
    plot!(cumulative, grid, cumsum(Poisson2D.(grid)) .* step, lw = 2, ls = :dash, color = :black,
          label = "2D Poisson")

    angles = range(0, 2π, length = 200)
    disc = scatter(real.(ratios), imag.(ratios), ms = 1.2, msw = 0, alpha = 0.35, color = :steelblue,
        label = nothing, xlabel = "Re z", ylabel = "Im z", aspect_ratio = :equal,
        title = @sprintf("complex spacing ratios: ⟨|z|⟩ = %.3f, -⟨cos⟩ = %.3f",
                         mean(abs, ratios), -mean(cos ∘ angle, ratios)))
    plot!(disc, cos.(angles), sin.(angles), lw = 1.5, color = :black, label = nothing)

    figure = plot(spectrum, histogramPanel, cumulative, disc, layout = (2, 2), size = (1300, 1000),
                  plot_title = title)
    showFigure && display(figure)

    if !isnothing(savePath)
        savefig(figure, savePath * ".png")
        savefig(figure, savePath * ".pdf")
    end

    return figure
end


# ---------------------------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------------------------

""" The SAME unfolding pipeline applied to two point sets whose statistics are known, so that the
    table compares like with like.

    The pipeline does reproduce the two laws - at a large enough sample the Ginibre matrix lands on
    var(s) = 0.086 against the law's 0.088 - so this is not a bias correction.  It is a scatter
    measurement: a single Ginibre spectrum of a thousand eigenvalues can give anything from 0.070
    to 0.106 for var(s) and from 0.16 to 0.29 for -<cos>, which is the same order as the distance
    between the Ginibre and Poisson columns at the sample sizes one sector of a small N provides.
    Without this row it is impossible to say whether the Liouvillian sits off Ginibre or whether the
    estimator does.  `realizations` averages a few spectra to keep the yardstick itself steady. """
function PipelineReference(; dimension = 2000, realizations = 2, points = 6000, neighbours = 30,
        keep = 0.5, seed = 20240923)

    rng = Xoshiro(seed)

    ginibreSpacings = Float64[]
    ginibreRatios = ComplexF64[]
    for _ = 1:realizations
        values = eigvals(randn(rng, ComplexF64, dimension, dimension) ./ sqrt(dimension))
        append!(ginibreSpacings, NearestNeighbourSpacings(values; neighbours = neighbours, keep = keep))
        append!(ginibreRatios, ComplexSpacingRatios(values; neighbours = neighbours, keep = keep))
    end

    uniform = [2 * (rand(rng) - 0.5) + 2im * (rand(rng) - 0.5) for _ = 1:points]

    return (ginibre = StatisticsSummary(ginibreSpacings, ginibreRatios),
            poisson = StatisticsSummary(
                NearestNeighbourSpacings(uniform; neighbours = neighbours, keep = keep),
                ComplexSpacingRatios(uniform; neighbours = neighbours, keep = keep)))
end


""" Sectors that carry INDEPENDENT statistics.  q and -q are complex conjugates of each other, so
    their spacing sets are identical up to a reflection and pooling both adds nothing but a
    duplicate; one representative of each conjugate pair is kept.  Self-conjugate sectors are
    excluded here and handled separately, because they are the ones with the antiunitary
    constraint (header, §6). """
IndependentSectors(L) = [m for m = 1:(L - 1) if !SelfConjugateSector(L, m) && m <= L - m]


""" The whole calculation: spectrum, symmetry bookkeeping, statistics, figure.

    `Ns` may hold several boson numbers.  Different N are independent spectra of the same model at
    different hbar_eff = 1/N, so their (separately unfolded) spacings may be pooled - which is the
    cheap way to a decent histogram, since one sector of a single N gives only about d_N^2 / L
    eigenvalues and the bulk selection then throws away half of them.

    Dimensions grow fast: one sector is about d_N^2 / L with d_N = C(N + L - 1, L - 1), and it is
    diagonalised densely.  L = 3 gives 1452 at N = 10, 4800 at N = 14, 7803 at N = 16 and 17787 at
    N = 20 - the last needs about 5 GB and an hour. """
function Analyse(; L = 3, Ns = [12], J = 1.0, g = -20.0, η = 3.0, κ = 0.3, Γsym = 0.0, γd = 0.0,
        neighbours = 30, keep = 0.5, calibrate = true,
        showFigures = true, savePath = nothing, verbose = true)

    clean = ComplexF64[]
    conjugated = ComplexF64[]
    cleanSpacings = Float64[]
    conjugatedSpacings = Float64[]
    cleanRatios = ComplexF64[]
    conjugatedRatios = ComplexF64[]
    allSectors = []

    for N in Ns
        p = LiouvillianParameters(L, N; J = J, g = g, η = η, κ = κ, Γsym = Γsym, γd = γd)
        dimension = binomial(N + L - 1, L - 1)

        verbose && @printf("\nL = %d, N = %2d: d_N = %d, (N,N) block = %d, per Z_L sector ~ %d\n",
                           L, N, dimension, dimension^2, dimension^2 ÷ L)

        verbose && N == Ns[1] && ReflectionWarning(p)

        sectors = SectorSpectra(p; verbose = verbose)
        N == Ns[end] && (allSectors = sectors)

        zeros = sum(count(z -> abs(z) < 1e-8, s.values) for s in sectors)
        decaying = vcat([filter(z -> abs(z) > 1e-8, s.values) for s in sectors]...)
        gap = -maximum(real, decaying)

        if verbose
            slowest = sort(decaying, by = z -> -real(z))[1:min(4, length(decaying))]
            @printf("  steady states (λ = 0): %d    Liouvillian gap: %.5f\n", zeros, gap)
            println("  slowest decaying modes: ",
                    join([@sprintf("%.4f%+.4fi", real(z), imag(z)) for z in slowest], "  "))
        end

        thisSpacings = Float64[]
        thisRatios = ComplexF64[]

        for s in sectors
            spacings = NearestNeighbourSpacings(s.values; neighbours = neighbours, keep = keep)
            ratios = ComplexSpacingRatios(s.values; neighbours = neighbours, keep = keep)

            if !s.selfConjugate && s.m in IndependentSectors(L)
                append!(thisSpacings, spacings)
                append!(thisRatios, ratios)
            end

            if s.selfConjugate
                append!(conjugated, s.values)
                append!(conjugatedSpacings, spacings)
                append!(conjugatedRatios, ratios)
            elseif s.m in IndependentSectors(L)
                append!(clean, s.values)
                append!(cleanSpacings, spacings)
                append!(cleanRatios, ratios)
            end
        end

        # Per-N line: the statistics should walk towards Ginibre as hbar_eff = 1/N shrinks, which is
        # the quantum-classical correspondence the note is after.  Pooling over N hides that trend,
        # so it is printed here before the pooled table.
        if verbose && !isempty(thisSpacings)
            r = StatisticsSummary(thisSpacings, thisRatios)
            @printf("  statistics at this N: var(s) = %.4f, <|z|> = %.4f, -<cos> = %.4f (n = %d)\n",
                    r.varianceSpacing, r.meanRadius, r.minusCosine, r.count)
        end
    end

    calibration = calibrate ? PipelineReference(; neighbours = neighbours, keep = keep) : nothing

    Report(cleanSpacings, cleanRatios, conjugatedSpacings, conjugatedRatios, L, calibration;
           verbose = verbose)

    if BHL_LOAD_PLOTS && !isempty(cleanSpacings) && (showFigures || !isnothing(savePath))
        title = @sprintf("L = %d, N ∈ %s, J = %g, g = %g, η = %g, κ = %g", L, string(Ns), J, g, η, κ)
        PlotSpectralStatistics(allSectors, cleanSpacings, cleanRatios,
                               LiouvillianParameters(L, Ns[end]; J = J, g = g, η = η, κ = κ);
                               savePath = savePath, title = title, showFigure = showFigures)
    end

    return (spacings = cleanSpacings, ratios = cleanRatios, sectors = allSectors,
            selfConjugateSpacings = conjugatedSpacings, selfConjugateRatios = conjugatedRatios)
end


""" Printed when eta = 0 puts the reflection of note §6 back into the model.  It changes what the
    statistics below MEAN, and in one sector it makes them wrong, so it is not a footnote. """
function ReflectionWarning(p::LiouvillianParameters)
    ReflectionIsSymmetry(p) || return nothing

    println()
    println("  " * "!" ^ 88)
    println("  eta = 0: the ring also has the reflection j -> -j of §6, so the Z_L decomposition")
    println("  below is NOT the complete symmetry analysis.")
    println("    * The reflection superoperator maps sector q to sector -q, exactly as rho -> rho^dag")
    println("      does, so it does not split a q != 0 sector any further.")
    println("    * Their COMPOSITION maps sector q to itself.  Every sector therefore carries an")
    println("      antiunitary symmetry and none of them is in the Ginibre class A: read the table")
    println("      against class AI-dagger (Sa, Ribeiro, Prosen), not against the Ginibre row.")
    println("    * The self-conjugate sectors (q = 0, and q = pi for even L) additionally split into")
    println("      the two reflection parities, which this code does NOT resolve - their statistics")
    println("      are then a superposition of two independent spectra and will look Poisson-like")
    println("      for that reason alone, whatever the dynamics does.")
    println("  eta = 0 is also the case in which BHNumberConserving.jl finds no classical chaos at")
    println("  all (note §8.2), so it is a control rather than the case of interest.")
    println("  " * "!" ^ 88)

    return nothing
end


""" The verdict table.  The Ginibre column is the one to read for the sectors WITHOUT an
    antiunitary constraint; the self-conjugate row is printed separately and must NOT be compared
    with it (class AI-dagger, different marker values - Sa, Ribeiro and Prosen). """
function Report(spacings, ratios, conjugatedSpacings, conjugatedRatios, L, calibration;
        verbose = true)
    verbose || return nothing

    println()
    println("=" ^ 92)
    println("SPECTRAL STATISTICS")
    println("=" ^ 92)
    @printf("\nZ_L sectors used: %s of 0:%d (one per conjugate pair q <-> -q);\n",
            string(IndependentSectors(L)), L - 1)
    @printf("self-conjugate (q = -q) sectors reported separately: %s\n",
            string([m for m = 0:(L - 1) if SelfConjugateSector(L, m)]))

    println("\n  ensemble                       n      <s>     var(s)   P(s<0.5)    <|z|>   -<cos>")
    println("  " * "-" ^ 86)

    if !isempty(spacings)
        r = StatisticsSummary(spacings, ratios)
        @printf("  %-28s %5d   %6.3f   %7.4f   %7.3f   %7.4f  %7.4f\n", "Liouvillian, clean sectors",
                r.count, r.meanSpacing, r.varianceSpacing, r.smallSpacings, r.meanRadius, r.minusCosine)
    end

    if !isempty(conjugatedSpacings)
        r = StatisticsSummary(conjugatedSpacings, conjugatedRatios)
        @printf("  %-28s %5d   %6.3f   %7.4f   %7.3f   %7.4f  %7.4f\n", "  self-conjugate sectors",
                r.count, r.meanSpacing, r.varianceSpacing, r.smallSpacings, r.meanRadius, r.minusCosine)
    end

    if !isnothing(calibration)
        for (label, r) in (("Ginibre matrix, same pipeline", calibration.ginibre),
                           ("uniform points, same pipeline", calibration.poisson))
            @printf("  %-28s %5d   %6.3f   %7.4f   %7.3f   %7.4f  %7.4f\n", label,
                    r.count, r.meanSpacing, r.varianceSpacing, r.smallSpacings,
                    r.meanRadius, r.minusCosine)
        end
    end

    @printf("  %-28s %5s   %6.3f   %7.4f   %7.3f   %7.4f  %7.4f\n",
            "Ginibre law (class A)", "-", 1.0, REFERENCE.ginibreVariance, REFERENCE.ginibreSmall,
            REFERENCE.ginibreRadius, REFERENCE.ginibreCosine)
    @printf("  %-28s %5s   %6.3f   %7.4f   %7.3f   %7.4f  %7.4f\n",
            "2D Poisson law", "-", 1.0, REFERENCE.poissonVariance, REFERENCE.poissonSmall,
            REFERENCE.poissonRadius, REFERENCE.poissonCosine)

    println()
    println("  Read the Liouvillian against the two SAME PIPELINE rows: they carry the scatter of")
    println("  the estimator at a comparable sample size, which at these dimensions is the same")
    println("  order as the gap between the Ginibre and Poisson columns.  The last two columns are")
    println("  the complex spacing ratios, which need no unfolding at all and are the ones to")
    println("  trust if the two disagree.")
    println()
    println("  Ginibre statistics mean the Liouvillian is dissipatively chaotic; 2D Poisson means")
    println("  either an integrable Liouvillian or - far more often - a symmetry that has not been")
    println("  resolved, which is the failure mode this file is built to avoid.")

    return nothing
end


""" The same statistics taken from the LEADING EDGE of a Ginibre spectrum: the `count` eigenvalues
    with the largest real part, out of a matrix of dimension `dimension`.

    This is the control the window scan needs.  The region near the steady state is also the edge of
    the eigenvalue cloud, and edge statistics differ from bulk statistics in EVERY ensemble - points
    at a boundary have neighbours on one side only.  So a drift towards Poisson-looking numbers when
    the window is tightened is partly expected and says nothing on its own.  This row measures how
    much of that drift a known-chaotic spectrum shows under the same selection, at the same sample
    size. """
function EdgeReference(count_; dimension = 3000, neighbours = 30, seed = 20240923)
    rng = Xoshiro(seed)
    values = eigvals(randn(rng, ComplexF64, dimension, dimension) ./ sqrt(dimension))

    distances, nearest, next = NeighbourData(values, neighbours)
    scale = distances[:, neighbours] ./ sqrt(neighbours / pi)
    spacings = distances[:, 1] ./ scale
    ratios = [(values[nearest[i]] - values[i]) / (values[next[i]] - values[i])
              for i in eachindex(values)]

    leading = sortperm(real.(values), rev = true)[1:min(count_, length(values))]
    selected = spacings[leading]

    return StatisticsSummary(selected ./ mean(selected), ratios[leading])
end


""" Spectral statistics as a function of how close to the steady state one looks - the scan of
    Corps and Relano (arXiv:2609.18464, Fig. 2b), for this model.

    Each cut keeps the eigenvalues with Re lambda >= cut and computes the statistics on them,
    SECTOR BY SECTOR, pooling only the resulting ratios: pooling the eigenvalues themselves first
    would superimpose independent spectra and manufacture Poisson statistics out of nothing, which
    is the very failure the Z_L decomposition exists to avoid.

    Read the table from the bottom up.  The last row is the whole spectrum, dominated by the bulk
    of fast-decaying modes; the first rows are the slow region where observables actually have
    weight.  If the two disagree, the slow region is the one that carries physical meaning, and the
    whole-spectrum number is the optimistic one. """
function WindowScan(sectors, L; cuts = [-1.0, -2.0, -4.0, -8.0, -16.0], neighbours = 30,
        keep = 0.5, edge = true, verbose = true)

    usable = [s for s in sectors if !s.selfConjugate && s.m in IndependentSectors(L)]
    isempty(usable) && error("no sector without an antiunitary constraint - see the header")

    span = minimum(minimum(real, s.values) for s in usable)
    results = []

    if verbose
        @printf("  spectral window            n      var(s)    <|z|>   -<cos>\n")
        println("  " * "-" ^ 64)
    end

    for cut in vcat(cuts, nothing)
        spacings = Float64[]
        ratios = ComplexF64[]
        failed = false

        for sector in usable
            inside = count(z -> isnothing(cut) || real(z) >= cut, sector.values)
            if inside <= 4 * neighbours
                failed = true
                break
            end
            append!(spacings, NearestNeighbourSpacings(sector.values; neighbours = neighbours,
                                                       keep = keep, window = cut))
            append!(ratios, ComplexSpacingRatios(sector.values; neighbours = neighbours,
                                                 keep = keep, window = cut))
        end

        if failed
            verbose && @printf("  Re λ >= %-8.1f      too few eigenvalues in the window\n", cut)
            continue
        end

        summary = StatisticsSummary(spacings, ratios)
        push!(results, (cut = cut, summary = summary))

        if verbose
            label = isnothing(cut) ? @sprintf("all (to %.1f), bulk", span) :
                                     @sprintf("Re λ >= %.1f", cut)
            @printf("  %-22s %6d   %7.4f  %7.4f  %7.4f\n", label, summary.count,
                    summary.varianceSpacing, summary.meanRadius, summary.minusCosine)

            if edge && !isnothing(cut)
                control = EdgeReference(summary.count; neighbours = neighbours)
                @printf("    Ginibre leading edge %6d   %7.4f  %7.4f  %7.4f\n",
                        control.count, control.varianceSpacing, control.meanRadius,
                        control.minusCosine)
            end
        end
    end

    if verbose
        @printf("  %-22s %6s   %7.4f  %7.4f  %7.4f\n", "Ginibre law", "-",
                REFERENCE.ginibreVariance, REFERENCE.ginibreRadius, REFERENCE.ginibreCosine)
        @printf("  %-22s %6s   %7.4f  %7.4f  %7.4f\n", "2D Poisson law", "-",
                REFERENCE.poissonVariance, REFERENCE.poissonRadius, REFERENCE.poissonCosine)
    end

    return results
end


""" Recommended starting point, matching the classical scan: the parameters at which
    BHNumberConserving.jl finds a globally attracting strange attractor (note §8.1). """
function Demo(; Ns = [10, 11, 12], kwargs...)
    return Analyse(; L = 3, Ns = Ns, J = 1.0, g = -20.0, η = 3.0, κ = 0.3, kwargs...)
end


# ---------------------------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------------------------

""" Consistency checks of the construction, in the spirit of note §10.  The decisive ones are the
    symmetry rows: if the union of the Z_L sectors did not reproduce the whole (N, N) block, or if
    a sector leaked, the spacing statistics would come out Poisson-like for reasons that have
    nothing to do with the physics. """
function Checks(; L = 3, N = 4, verbose = true)
    results = Tuple{String, Float64}[]
    p = LiouvillianParameters(L, N; g = -20.0, η = 3.0, κ = 0.3, Γsym = 0.05, γd = 0.02)

    basis = NumberBasis(L, N)
    d = length(basis)
    H = Hamiltonian(basis, p)
    jumps = JumpOperators(basis, p)
    number = sum(Hop(basis, j, j) for j = 1:L)
    T = TranslationOperator(basis)

    push!(results, ("d_N = C(N + L - 1, L - 1)", abs(d - binomial(N + L - 1, L - 1))))
    push!(results, ("H Hermitian", norm(H - adjoint(H))))
    push!(results, ("strong U(1): [N, H] = 0 and [N, L_k] = 0 for every jump",
                    max(norm(number * H - H * number),
                        maximum(norm(number * Lk - Lk * number) for Lk in jumps))))
    push!(results, ("weak Z_L: [T, H] = 0 and [T, sum_k L_k^dag L_k] = 0",
                    max(norm(T * H - H * T),
                        norm(T * sum(adjoint(Lk) * Lk for Lk in jumps) -
                             sum(adjoint(Lk) * Lk for Lk in jumps) * T))))
    push!(results, ("T^L = 1 and T unitary",
                    max(norm(T^L - I), norm(T * adjoint(T) - I))))

    # --- momentum basis --------------------------------------------------------------------------
    data = MomentumData(basis, p)
    push!(results, ("momentum subspace dimensions sum to d_N",
                    abs(sum(length, data.ranges) - d)))

    # --- the symmetry decomposition is exact ------------------------------------------------------
    full = Matrix(LiouvillianFull(basis, p))
    fullValues = eigvals(full)
    sectorValues = vcat([eigvals(LiouvillianSector(data, m)) for m = 0:(L - 1)]...)

    Key(z) = (round(real(z), digits = 7), round(imag(z), digits = 7))
    push!(results, ("union of Z_L sectors = the whole (N, N) block, eigenvalue by eigenvalue",
                    maximum(abs, sort(fullValues, by = Key) .- sort(sectorValues, by = Key))))

    # --- Lindblad structure ------------------------------------------------------------------------
    identity = vec(Matrix{ComplexF64}(I, d, d))
    push!(results, ("trace preservation: vec(1)^dag L = 0", norm(transpose(identity) * full)))

    rng = Xoshiro(7)
    ρ = randn(rng, ComplexF64, d, d)
    image = reshape(full * vec(ρ), d, d)
    imageAdjoint = reshape(full * vec(adjoint(ρ)), d, d)
    push!(results, ("Hermiticity preservation: L(rho^dag) = L(rho)^dag",
                    maximum(abs, imageAdjoint .- adjoint(image))))

    # --- the steady state is a state ---------------------------------------------------------------
    decomposition = eigen(full)
    steady = decomposition.vectors[:, argmin(abs.(decomposition.values))]
    steady = reshape(steady, d, d)
    steady ./= tr(steady)
    push!(results, ("steady state Hermitian after normalising its trace",
                    maximum(abs, steady .- adjoint(steady))))
    push!(results, ("steady state positive (most negative eigenvalue; must be ~0)",
                    abs(min(0.0, minimum(real, eigvals(Hermitian(0.5 * (steady + adjoint(steady)))))))))
    push!(results, ("exactly one zero eigenvalue in the (N, N) block (Buca-Prosen)",
                    abs(count(z -> abs(z) < 1e-8, fullValues) - 1)))

    # --- rho -> rho^dag maps sector q to sector -q ---------------------------------------------------
    spectra = [eigvals(LiouvillianSector(data, m)) for m = 0:(L - 1)]
    worstPair = 0.0
    worstSelf = 0.0
    for m = 0:(L - 1)
        mirror = mod(L - m, L)
        difference = maximum(abs, sort(conj.(spectra[m + 1]), by = Key) .-
                                  sort(spectra[mirror + 1], by = Key))
        SelfConjugateSector(L, m) ? (worstSelf = max(worstSelf, difference)) :
                                    (worstPair = max(worstPair, difference))
    end
    push!(results, ("sector q and sector -q hold complex-conjugate spectra", worstPair))
    push!(results, ("a self-conjugate sector (q = -q) is its own conjugate", worstSelf))

    # --- the conjugation shortcut used by SectorSpectra -----------------------------------------------
    shortcut = SectorSpectra(p; verbose = false, useConjugation = true)
    explicit = SectorSpectra(p; verbose = false, useConjugation = false)
    push!(results, ("conjugation shortcut in SectorSpectra vs diagonalising every sector",
                    maximum(maximum(abs, sort(shortcut[i].values, by = Key) .-
                                         sort(explicit[i].values, by = Key)) for i = 1:L)))

    # --- the Grobe-Haake-Sommers law -----------------------------------------------------------------
    grid = range(0, 5, length = 5001)
    step = grid[2] - grid[1]
    ginibre = GinibreSpacing.(grid)
    push!(results, ("Ginibre law normalisation and mean", max(abs(sum(ginibre) * step - 1),
                                                              abs(sum(grid .* ginibre) * step - 1))))
    push!(results, ("Grobe-Haake-Sommers constant vs the published 1.1429",
                    abs(GINIBRE_SCALE - 1.1429)))
    push!(results, ("cubic repulsion: P(s)/s^3 finite and non-zero as s -> 0",
                    abs(GinibreSpacing(0.05) / 0.05^3 - GinibreSpacing(0.1) / 0.1^3)))

    # --- the statistics machinery itself, on the two ensembles it is supposed to distinguish ---------
    # A relative check: a few per cent is the sampling error of a single 1500 x 1500 matrix, so this
    # row is a per-cent row, not a 1e-15 one.
    rng = Xoshiro(20240923)
    ginibreMatrix = eigvals(randn(rng, ComplexF64, 1500, 1500) ./ sqrt(1500))
    uniform = [2 * (rand(rng) - 0.5) + 2im * (rand(rng) - 0.5) for _ = 1:4000]

    for (label, values, reference) in (("Ginibre matrix", ginibreMatrix, 0.7378),
                                       ("uniform random points", uniform, 2 / 3))
        z = ComplexSpacingRatios(values; neighbours = 30, keep = 0.5)
        push!(results, ("<|z|> recovered from $label (relative, per cent expected)",
                        abs(mean(abs, z) - reference) / reference))
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


if abspath(PROGRAM_FILE) == @__FILE__
    Checks()
    println()
    Demo(; Ns = [12, 13, 14], savePath = "liouvillian_spacings_3", showFigures = false)
end
