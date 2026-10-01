# BHNumberConservingLiouvillianSparse.jl
#
# The SLOW part of the Liouvillian spectrum at large N - item C3 of number-conserving-BH-paper-TODO.md,
# and the gap at large N of item C6 - without ever forming a dense matrix.
#
# WHY A SECOND CONSTRUCTION.  BHNumberConservingLiouvillian.jl builds each Z_L sector in the MOMENTUM
# basis of the N-boson space.  That basis is dense, so K and the jump operators become dense
# matrices and every sector is a dense (d_N^2 / L)^2 block: 5 GB at N = 20 and 107 GB at N = 30.
# Here the same sector is built in the basis of TRANSLATION ORBITS of Fock pairs,
#
#     |O, m> = N_O sum_{r=0}^{L-1} exp(-i q r) T^r |a><b| T^-r,      q = 2 pi m / L,
#
# with (a, b) the representative Fock pair of the orbit O.  The Liouvillian commutes with the
# superoperator translation (the weak Z_L of §6), and it maps a Fock pair onto O(10 L) Fock pairs,
# so the sector matrix is SPARSE with a few dozen entries per column.  The eigenvalues are those of
# LiouvillianSector(data, m) exactly (SparseChecks() below compares them).
#
# Orbit normalisation: an orbit of length l (the number of distinct pairs it holds, a divisor of L)
# contributes to sector m only if m l / L is an integer, and then N_O = sqrt(l) / L.  At L = 3 the
# only short orbit is the pair of the uniform Fock state (N/3, N/3, N/3) with itself, which exists
# for N divisible by 3 and lives in sector 0 only.
#
# SPECTRUM SLICING.  The slow region Re λ >= -w is a long strip along the right edge of the cloud
# (Im λ reaches about ±7.6 N at L = 3), and at large N it holds thousands of eigenvalues, too many
# for one Arnoldi run.  It is therefore cut into SLICES: shift-invert Arnoldi (KrylovKit) around
# shifts σ_k spread over the strip, each returning the eigenvalues nearest to its own shift.  A
# slice is TRUSTED only inside the disc |λ - σ_k| < r_k, where r_k is slightly less than the
# distance of the farthest CONVERGED eigenvalue: everything inside that disc has been found, because
# shift-invert Arnoldi converges from the largest |1 / (λ - σ)| inwards.  The union of the trusted
# discs must cover the window; Coverage() checks that and names the holes, which are closed by
# adding shifts.  The slices are independent, so they parallelise trivially - across processes on
# the laptop, or as a SLURM array on the cluster (BHNumberConservingLiouvillianSparseChimera.sh).
#
# USAGE (command line; every step is resumable and writes into RESULTS/<run>/):
#
#   julia BHNumberConservingLiouvillianSparse.jl checks                   # orbit vs momentum basis
#   julia BHNumberConservingLiouvillianSparse.jl validate  <N>            # sparse slicing vs dense
#   julia BHNumberConservingLiouvillianSparse.jl plan      <run> <N> <m> <w> [--g --eta --kappa
#                                                          --Gsym --howmany --spacing --halfwidth]
#   julia BHNumberConservingLiouvillianSparse.jl slice     <run> <k>      # one shift (SLURM task)
#   julia BHNumberConservingLiouvillianSparse.jl slices    <run> [--part i --parts n]
#                                                                         # all missing shifts, or
#                                                                         # every n-th from i
#   julia BHNumberConservingLiouvillianSparse.jl merge     <run>          # union + coverage
#   julia BHNumberConservingLiouvillianSparse.jl refine    <run>          # shifts at the holes
#   julia BHNumberConservingLiouvillianSparse.jl statistics <run>         # ⟨r⟩, -⟨cos θ⟩, var(s)
#
# `plan` writes plan.txt (model, sector, window and the list of shifts; without --spacing it first
# computes one calibration slice to choose the spacing); `slice k` computes shift k (0-based,
# = SLURM_ARRAY_TASK_ID when k is omitted) and writes slice_<k>.txt; `merge` writes spectrum.txt
# (the deduplicated eigenvalues inside the window) and reports holes in the coverage; `refine`
# appends shifts at the holes, after which `slices` and `merge` are simply run again.
#
# C3 (slow region at large N):  plan <run> <N> 1 21, then slices / merge / refine until complete;
#                               `statistics` then gives w = 8, 12, 16 (a cut needs 5 units of margin
#                               for the local unfolding - see WindowStatistics).
# C6 (gap at large N):          plan <run> <N> 0 <w> and plan <run'> <N> 1 <w> with w a little
#                               above the expected gap (1-2); merge prints the slowest mode of each
#                               sector, and the gap is the smaller of the two values.
#
# The results directory is $BH_RESULTS_DIR or ~/results/bh/number-conserving/quantum/3/sparse.

using SparseArrays
using LinearAlgebra
using Random
using Statistics
using Printf

ENV["CD_NO_PLOTS"] = get(ENV, "CD_NO_PLOTS", "true")
isdefined(Main, :LiouvillianSector) || include(joinpath(@__DIR__, "BHNumberConservingLiouvillian.jl"))

using KrylovKit


# ---------------------------------------------------------------------------------------------
# Orbit basis of Fock pairs and the sparse sector
# ---------------------------------------------------------------------------------------------

""" Translation orbits of Fock pairs (a, b), compatible with the sector m.

    `orbit[a + d (b - 1)]` is the orbit of the pair (a, b) (0 when the orbit does not contribute
    to sector m) and `shift[...]` the power s with (a, b) = T^s (representative); `size` is the
    orbit length l and `normalisation` N_O = sqrt(l) / L. """
struct OrbitBasis
    L::Int
    N::Int
    m::Int
    d::Int
    representatives::Vector{Tuple{Int, Int}}
    normalisation::Vector{Float64}
    orbit::Vector{Int32}
    shift::Vector{Int8}
end

function OrbitBasis(basis::NumberBasis, m::Integer)
    L, d = basis.L, length(basis)
    translated = [basis.index[circshift(state, 1)] for state in basis.states]     # T|a> = |t(a)>

    orbit = zeros(Int32, d * d)
    shift = zeros(Int8, d * d)
    visited = falses(d * d)
    representatives = Tuple{Int, Int}[]
    normalisation = Float64[]

    for b = 1:d, a = 1:d
        linear = a + d * (b - 1)
        visited[linear] && continue

        members = Int[]
        x, y = a, b
        for s = 0:(L - 1)
            index = x + d * (y - 1)
            visited[index] && break
            visited[index] = true
            push!(members, index)
            x, y = translated[x], translated[y]
        end

        l = length(members)
        mod(m * l, L) == 0 || continue                  # orbit does not reach this sector

        push!(representatives, (a, b))
        push!(normalisation, sqrt(l) / L)
        for (s, index) in enumerate(members)
            orbit[index] = length(representatives)
            shift[index] = s - 1
        end
    end

    return OrbitBasis(L, basis.N, m, d, representatives, normalisation, orbit, shift)
end

Base.length(orbits::OrbitBasis) = length(orbits.representatives)


""" The Liouvillian restricted to sector m, as a sparse matrix in the orbit basis.

    Column O is L|O, m> expanded back into orbit vectors: L|R_O> = sum_x c_x |x> over Fock pairs x,
    and x = T^s R_O' contributes c_x exp(i q s) N_O / N_O' to the entry (O', O).  Only the
    representative pair of each orbit is ever acted upon. """
function SparseSector(p::LiouvillianParameters, m::Integer; basis = NumberBasis(p.L, p.N),
        orbits = OrbitBasis(basis, m))

    d = length(basis)
    q = 2 * pi * m / p.L
    phases = [exp(im * q * s) for s = 0:(p.L - 1)]

    H = Hamiltonian(basis, p)
    jumps = JumpOperators(basis, p)
    K = EffectiveGenerator(H, jumps)

    rows = Int[]; cols = Int[]; values = ComplexF64[]
    sizehint!(rows, 40 * length(orbits)); sizehint!(cols, 40 * length(orbits))
    sizehint!(values, 40 * length(orbits))

    function Deposit!(column, weight, c, e, value)
        linear = c + d * (e - 1)
        target = orbits.orbit[linear]
        target == 0 && return
        push!(rows, target); push!(cols, column)
        push!(values, value * phases[orbits.shift[linear] + 1] * weight / orbits.normalisation[target])
    end

    for (column, (a, b)) in enumerate(orbits.representatives)
        weight = orbits.normalisation[column]

        # K |a><b|
        for index in nzrange(K, a)
            Deposit!(column, weight, rowvals(K)[index], b, nonzeros(K)[index])
        end
        # |a><b| K^dag : the coefficient of |a><e| is conj(K[e, b])
        for index in nzrange(K, b)
            Deposit!(column, weight, a, rowvals(K)[index], conj(nonzeros(K)[index]))
        end
        # L_k |a><b| L_k^dag
        for Lk in jumps
            for i in nzrange(Lk, a), j in nzrange(Lk, b)
                Deposit!(column, weight, rowvals(Lk)[i], rowvals(Lk)[j],
                         nonzeros(Lk)[i] * conj(nonzeros(Lk)[j]))
            end
        end
    end

    n = length(orbits)
    return sparse(rows, cols, values, n, n)          # duplicates are summed
end


# ---------------------------------------------------------------------------------------------
# Shift-invert slices
# ---------------------------------------------------------------------------------------------

""" UMFPACK settings for the sector matrices.  They are STRUCTURALLY symmetric (every pair of
    Fock pairs is connected in both directions by the Hamiltonian, even where a jump acts one
    way only), so the symmetric strategy applies, and METIS nested dissection suits their
    four-dimensional lattice structure (two points of the simplex) far better than the default:
    at N = 16 it cuts the fill of L + U by a third, and at N = 20 the default AMD/COLAMD choice
    ran past 18 GB without finishing, where METIS stays near 1 GB.  Fill grows roughly as
    n^1.5 - about 3 GB at N = 24 and 12 GB at N = 30. """
function FactorisationControl()
    control = SparseArrays.UMFPACK.get_umfpack_control(ComplexF64, Int64)
    control[6] = 3.0          # UMFPACK_STRATEGY  = UMFPACK_STRATEGY_SYMMETRIC
    control[11] = 3.0         # UMFPACK_ORDERING  = UMFPACK_ORDERING_METIS
    return control
end

""" Eigenvalues nearest to `σ`, by shift-invert Arnoldi on the sparse sector `M`.

    Returns (values, radius): the converged eigenvalues sorted by distance from σ, and the radius of
    the TRUSTED disc - `safety` times the distance of the farthest converged eigenvalue.  Inside
    that disc the list is complete.  `howmany` sets the size of the slice; the Krylov dimension is
    taken generously because the eigenvalues of (M - σ)^-1 cluster at the origin, not at the edge
    that Arnoldi converges from, and a tight subspace would waste restarts. """
function Slice(M, σ; howmany = 300, krylovdim = 2 * howmany + 50, tolerance = 1e-10,
        maxiter = 50, safety = 0.98, seed = 1, verbose = false)

    n = size(M, 1)
    time = @elapsed factor = lu(M - σ * I; control = FactorisationControl())
    verbose && @printf("    LU of %d x %d (%d nonzeros): %.1f s\n", n, n, nnz(M), time)

    x0 = randn(Xoshiro(seed), ComplexF64, n)
    howmany = min(howmany, n - 2)
    krylovdim = min(krylovdim, n - 1)

    time = @elapsed μ, _, info = eigsolve(x -> factor \ x, x0, howmany, :LM;
                                           krylovdim = krylovdim, tol = tolerance,
                                           maxiter = maxiter, eager = false, ishermitian = false)
    verbose && @printf("    Arnoldi: %d of %d converged, %d matrix-vector products, %.1f s\n",
                       info.converged, howmany, info.numops, time)

    converged = μ[1:min(info.converged, length(μ))]
    isempty(converged) && return (values = ComplexF64[], radius = 0.0)

    values = σ .+ 1 ./ converged
    order = sortperm(abs.(values .- σ))
    values = values[order]
    radius = safety * abs(values[end] - σ)

    return (values = filter(z -> abs(z - σ) < radius, values), radius = radius)
end


""" Shifts covering the strip  -w <= Re λ <= 0,  |Im λ - centre| <= halfWidth.

    A hexagonal lattice of points with nearest-neighbour distance `spacing`: discs of radius
    spacing / sqrt(3) centred on it cover the plane, so the strip is covered as soon as every
    trusted radius exceeds that value.  The first column sits at Re σ = `right`, a little to the
    right of the steady state (the sector q = 0 contains λ = 0 and σ must not coincide with it).
    The last column only has to reach -w: a column of the lattice covers a band extending
    spacing / sqrt(12) on either side of it. """
function ShiftLattice(w, halfWidth; spacing = 2.0, centre = 0.0, right = 0.25)
    shifts = ComplexF64[]
    rowHeight = spacing * sqrt(3) / 2
    columns = max(1, ceil(Int, 1 + (right + w - spacing / sqrt(12)) / rowHeight))

    for c = 0:(columns - 1)
        x = right - c * rowHeight
        offset = isodd(c) ? spacing / 2 : 0.0
        for y in (centre - halfWidth - offset):spacing:(centre + halfWidth + spacing)
            push!(shifts, complex(x, y))
        end
    end

    return shifts
end


""" Union of the slices: eigenvalues inside their own trusted discs, deduplicated (two values
    closer than `tolerance` are the same eigenvalue found by two overlapping slices). """
function MergeSlices(slices; tolerance = 1e-7)
    merged = ComplexF64[]
    for slice in slices, z in slice.values
        any(abs(z - y) < tolerance * max(1.0, abs(z)) for y in merged) || push!(merged, z)
    end
    return sort(merged, by = z -> (-real(z), imag(z)))
end


""" Points of the strip  -w <= Re λ <= 0,  |Im λ| <= halfWidth  that no trusted disc covers.  An
    empty list means the merged spectrum is COMPLETE inside the window; otherwise the holes say
    where to add shifts (Refine does that). """
function Coverage(slices, w, halfWidth; resolution = 0.1)
    centres = [s.σ for s in slices]
    radii = [s.radius for s in slices]

    holes = ComplexF64[]
    for x in range(-w, 0.0, step = resolution), y in range(-halfWidth, halfWidth, step = resolution)
        z = complex(x, y)
        any(abs(z - c) < r for (c, r) in zip(centres, radii)) || push!(holes, z)
    end

    return holes
end


""" New shifts that close the holes of a coverage: greedily, a hole becomes a shift and every hole
    within `reach` of it counts as closed.  `reach` should be well below the typical trusted radius
    (half the median of those found so far); a hole left open is caught by the next round. """
function HoleShifts(holes, reach)
    remaining = copy(holes)
    shifts = ComplexF64[]
    while !isempty(remaining)
        σ = remaining[1]
        push!(shifts, σ)
        remaining = filter(z -> abs(z - σ) > reach, remaining)
    end
    return shifts
end


# ---------------------------------------------------------------------------------------------
# Run directories: plan -> slices -> merge [-> refine -> slices -> merge] -> statistics
# ---------------------------------------------------------------------------------------------

const SPARSE_RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "quantum", "3", "sparse"))

RunPath(run) = joinpath(SPARSE_RESULTS, run)

function WritePlanFile(path, entries, shifts)
    temporary = joinpath(path, "plan.txt.part")
    open(temporary, "w") do io
        println(io, "# sparse slicing plan, BHNumberConservingLiouvillianSparse.jl")
        for (key, value) in entries
            println(io, key, "	", value)
        end
        println(io, "shifts	", join([@sprintf("%.6f%+.6fim", real(s), imag(s)) for s in shifts], ","))
    end
    mv(temporary, joinpath(path, "plan.txt"); force = true)
end

""" plan.txt: the model, the sector, the window and the shifts.

    `halfWidth` bounds the strip in Im λ; the default 8 N covers the imaginary extent of the dense
    spectra (about 7.6 N at L = 3).  `spacing = nothing` calibrates it: one slice is computed at the
    centre of the window and the lattice spacing is set to 1.5 times its trusted radius, so that
    the covering radius spacing / sqrt(3) = 0.87 r leaves a margin for the density varying along
    the strip.  At N = 30 the calibration costs one slice (half an hour); pass `spacing` from a
    smaller N scaled by the density to skip it. """
function WritePlan(run, N, m, w; L = 3, J = 1.0, g = -20.0, η = 3.0, κ = 0.3, Γsym = 0.0,
        γd = 0.0, spacing = nothing, halfWidth = 8.0 * N, howmany = 300)

    path = RunPath(run)
    mkpath(path)
    parameters = LiouvillianParameters(L, N; J = J, g = g, η = η, κ = κ, Γsym = Γsym, γd = γd)

    if isnothing(spacing)
        M = SparseSector(parameters, m)
        time = @elapsed calibration = Slice(M, complex(-w / 2, 0.0); howmany = howmany)
        spacing = 1.5 * calibration.radius
        @printf("Calibration slice at σ = %.1f: %d eigenvalues within r = %.2f (%.1f s) -> spacing %.2f
",
                -w / 2, length(calibration.values), calibration.radius, time, spacing)
    end

    shifts = ShiftLattice(w, halfWidth; spacing = spacing)
    WritePlanFile(path, (("L", L), ("N", N), ("m", m), ("w", w), ("J", J), ("g", g), ("eta", η),
                         ("kappa", κ), ("Gsym", Γsym), ("gd", γd), ("spacing", spacing),
                         ("halfWidth", halfWidth), ("howmany", howmany)), shifts)

    @printf("Plan %s: L = %d, N = %d, sector m = %d, window Re λ >= -%g, |Im λ| <= %g, %d shifts -> %s
",
            run, L, N, m, w, halfWidth, length(shifts), path)
    @printf("  SLURM: --array=0-%d
", length(shifts) - 1)

    return shifts
end

function ReadPlan(run)
    plan = Dict{String, String}()
    entries = Tuple{String, String}[]
    for line in eachline(joinpath(RunPath(run), "plan.txt"))
        (startswith(line, "#") || !occursin('	', line)) && continue
        key, value = split(line, '	'; limit = 2)
        plan[key] = value
        key == "shifts" || push!(entries, (key, value))
    end

    number(key) = parse(Float64, plan[key])
    parameters = LiouvillianParameters(parse(Int, plan["L"]), parse(Int, plan["N"]);
        J = number("J"), g = number("g"), η = number("eta"), κ = number("kappa"),
        Γsym = number("Gsym"), γd = number("gd"))
    shifts = [parse(ComplexF64, s) for s in split(plan["shifts"], ',')]

    return (parameters = parameters, m = parse(Int, plan["m"]), w = number("w"),
            halfWidth = number("halfWidth"), howmany = parse(Int, plan["howmany"]),
            shifts = shifts, entries = entries)
end

SliceFile(run, k) = joinpath(RunPath(run), "slice_$(k).txt")

""" Computes the missing shifts among `indices` (0-based).  The sector is built once per process
    and reused, so a SLURM task handling several shifts pays for the construction only once. """
function RunSlices(run, indices; verbose = true)
    plan = ReadPlan(run)
    todo = [k for k in indices if !isfile(SliceFile(run, k))]
    isempty(todo) && (println("nothing to do"); return nothing)

    time = @elapsed M = SparseSector(plan.parameters, plan.m)
    verbose && @printf("Sector m = %d at N = %d: %d x %d, %d nonzeros, built in %.1f s\n",
                       plan.m, plan.parameters.N, size(M, 1), size(M, 2), nnz(M), time)

    for k in todo
        σ = plan.shifts[k + 1]
        verbose && @printf("  shift %d: σ = %.3f%+.3fi\n", k, real(σ), imag(σ))
        time = @elapsed result = Slice(M, σ; howmany = plan.howmany, seed = k + 1, verbose = verbose)
        verbose && @printf("    %d eigenvalues in the trusted disc of radius %.3f (%.1f s)\n",
                           length(result.values), result.radius, time)

        temporary = SliceFile(run, k) * ".part"
        open(temporary, "w") do io
            @printf(io, "# sigma\t%.12e\t%.12e\n# radius\t%.12e\n", real(σ), imag(σ), result.radius)
            for z in result.values
                @printf(io, "%.15e\t%.15e\n", real(z), imag(z))
            end
        end
        mv(temporary, SliceFile(run, k); force = true)     # a killed task leaves no half file
    end

    return nothing
end

function ReadSlice(file)
    σ = 0.0im; radius = 0.0
    values = ComplexF64[]
    for line in eachline(file)
        fields = split(line, '\t')
        if startswith(line, "# sigma")
            σ = complex(parse(Float64, fields[2]), parse(Float64, fields[3]))
        elseif startswith(line, "# radius")
            radius = parse(Float64, fields[2])
        elseif !startswith(line, "#") && length(fields) == 2
            push!(values, complex(parse(Float64, fields[1]), parse(Float64, fields[2])))
        end
    end
    return (σ = σ, radius = radius, values = values)
end

""" Union of every slice file of the run, restricted to the window; writes spectrum.txt and
    reports the coverage.  Missing slices are listed, not silently skipped. """
function Merge(run; verbose = true)
    plan = ReadPlan(run)
    missing = [k for k in 0:(length(plan.shifts) - 1) if !isfile(SliceFile(run, k))]
    slices = [ReadSlice(SliceFile(run, k)) for k in 0:(length(plan.shifts) - 1)
              if isfile(SliceFile(run, k))]

    merged = filter(z -> real(z) >= -plan.w && abs(imag(z)) <= plan.halfWidth, MergeSlices(slices))
    holes = isempty(slices) ? ComplexF64[] : Coverage(slices, plan.w, plan.halfWidth)

    open(joinpath(RunPath(run), "spectrum.txt"), "w") do io
        @printf(io, "# N = %d, sector m = %d, window Re >= -%g, %d eigenvalues, %d uncovered grid points, %d missing slices\n",
                plan.parameters.N, plan.m, plan.w, length(merged), length(holes), length(missing))
        for z in merged
            @printf(io, "%.15e\t%.15e\n", real(z), imag(z))
        end
    end

    if verbose
        @printf("%s: %d of %d slices present, %d eigenvalues with Re λ >= -%g\n", run,
                length(slices), length(plan.shifts), length(merged), plan.w)
        isempty(missing) || println("  missing slices: ", join(missing, ","))
        if isempty(holes)
            println("  coverage COMPLETE: every point of the window lies in a trusted disc")
        else
            @printf("  coverage INCOMPLETE: %d grid points uncovered, e.g. %s\n", length(holes),
                    join([@sprintf("%.2f%+.2fi", real(z), imag(z)) for z in holes[1:min(5, end)]], ", "))
            println("  -> run `refine` to add shifts at the holes, then `slices` and `merge` again")
        end
        decaying = filter(z -> abs(z) > 1e-8, merged)
        if !isempty(decaying)
            @printf("  slowest decaying mode: %.6f%+.6fi  (gap of this sector %.6f)\n",
                    real(decaying[1]), imag(decaying[1]), -real(decaying[1]))
        end
    end

    return (values = merged, holes = holes, missing = missing)
end

""" Appends shifts at the holes of the current coverage to plan.txt; the new shifts get the next
    indices, so the slices already computed stay valid and only the new ones have to be run. """
function Refine(run)
    plan = ReadPlan(run)
    slices = [ReadSlice(SliceFile(run, k)) for k in 0:(length(plan.shifts) - 1)
              if isfile(SliceFile(run, k))]
    holes = Coverage(slices, plan.w, plan.halfWidth)
    if isempty(holes)
        println("coverage complete, nothing to refine")
        return 0
    end

    reach = 0.5 * median(s.radius for s in slices if s.radius > 0)
    added = HoleShifts(holes, reach)
    WritePlanFile(RunPath(run), plan.entries, vcat(plan.shifts, added))
    @printf("%d holes closed by %d new shifts (indices %d..%d); SLURM: --array=%d-%d
",
            length(holes), length(added), length(plan.shifts), length(plan.shifts) + length(added) - 1,
            length(plan.shifts), length(plan.shifts) + length(added) - 1)
    return length(added)
end

function ReadSpectrum(file)
    values = ComplexF64[]
    for line in eachline(file)
        startswith(line, "#") && continue
        fields = split(line, '\t')
        push!(values, complex(parse(Float64, fields[1]), parse(Float64, fields[2])))
    end
    return values
end


""" Ratio and spacing statistics of the merged window, for the nested windows `cuts`.

    The eigenvalues are known only inside the covered strip Re λ >= -covered, so the neighbours of a
    reference point are searched there.  The nearest and next-to-nearest neighbour (the ratios)
    lie within a spacing or two, but the local unfolding of the spacings takes the density from the
    `neighbours`-th neighbour, about sqrt(neighbours / (π ρ)) ≈ 4-5 away at L = 3.  A cut w is
    therefore used only when w + margin <= covered, with margin = 5: plan the slices to cover the
    largest cut plus 5 (-16 needs a plan with w = 21).  Validate() shows that the statistics then
    agree with the dense ones. """
function WindowStatistics(values, cuts; covered = Inf, neighbours = 30, margin = 5.0)
    values = filter(z -> abs(z) > 1e-8 && real(z) >= -covered, values)
    results = []
    for w in cuts
        w + margin <= covered || continue
        window = filter(z -> real(z) >= -w - margin, values)
        references = findall(z -> real(z) >= -w, window)
        (length(window) > 4 * neighbours && length(references) > neighbours) || continue

        distances, nearest, next = NeighbourData(window, neighbours)
        scale = distances[:, neighbours] ./ sqrt(neighbours / pi)
        spacings = (distances[:, 1] ./ scale)[references]
        ratios = [(window[nearest[i]] - window[i]) / (window[next[i]] - window[i]) for i in references]
        push!(results, (w = w, summary = StatisticsSummary(spacings ./ mean(spacings), ratios)))
    end
    return results
end


# ---------------------------------------------------------------------------------------------
# Checks and validation against the dense sectors
# ---------------------------------------------------------------------------------------------

""" The orbit construction against LiouvillianSector: same eigenvalues, sector by sector, including
    N divisible by 3 (short orbit) and η = 0. """
function SparseChecks(; verbose = true)
    Key(z) = (round(real(z), digits = 6), round(imag(z), digits = 6))
    worst = 0.0

    for (N, η, κ, Γsym, γd) in ((5, 3.0, 0.3, 0.0, 0.0), (6, 3.0, 0.3, 0.05, 0.02), (6, 0.0, 0.3, 0.0, 0.0),
                                (7, 1.0, 0.5, 0.0, 0.0))
        p = LiouvillianParameters(3, N; g = -20.0, η = η, κ = κ, Γsym = Γsym, γd = γd)
        data = MomentumData(NumberBasis(3, N), p)
        for m = 0:2
            dense = eigvals(LiouvillianSector(data, m))
            sparseValues = eigvals(Matrix(SparseSector(p, m)))
            length(dense) == length(sparseValues) ||
                error("sector sizes differ: $(length(dense)) vs $(length(sparseValues))")
            difference = maximum(abs, sort(dense, by = Key) .- sort(sparseValues, by = Key))
            worst = max(worst, difference)
            verbose && @printf("  N = %d, η = %.1f, m = %d: dimension %d, max |Δλ| = %.1e\n",
                               N, η, m, length(dense), difference)
        end
    end

    verbose && @printf("orbit basis vs momentum basis: worst eigenvalue difference %.1e\n", worst)
    return worst
end


""" The whole slicing pipeline against the dense spectrum at a size where both are available - the
    validation that C3 asks for.  Calibrated plan over Re λ >= -(w + margin), slices, refinement
    until the coverage is complete, merge; then (i) the merged strip must contain exactly the dense
    eigenvalues, nothing more, and (ii) the window statistics up to w are printed three ways: dense with the
    neighbours taken from the WHOLE spectrum (what WindowScan does at N <= 16), dense with the
    neighbours restricted to the window plus a margin (what is possible at large N), and sparse
    with the same restriction.  (i) tests the slicing, (ii) the effect of not knowing the
    eigenvalues outside the window. """
function Validate(N; m = 1, w = 12.0, g = -20.0, η = 3.0, κ = 0.3, howmany = 150, margin = 5.0)
    p = LiouvillianParameters(3, N; g = g, η = η, κ = κ)
    data = MomentumData(NumberBasis(3, N), p)
    dense = eigvals(LiouvillianSector(data, m))
    statisticsWindow = w
    w = w + margin                  # the slices cover the statistics window plus the margin
    inside = filter(z -> real(z) >= -w && abs(z) > 1e-8, dense)
    halfWidth = 8.0 * N

    M = SparseSector(p, m)
    time = @elapsed begin
        calibration = Slice(M, complex(-w / 2, 0.0); howmany = howmany)
        shifts = ShiftLattice(w, halfWidth; spacing = 1.5 * calibration.radius)
        @printf("  calibration: %d eigenvalues within r = %.2f -> %d shifts
",
                length(calibration.values), calibration.radius, length(shifts))
        slices = [merge((σ = σ,), Slice(M, σ; howmany = howmany, seed = k)) for (k, σ) in enumerate(shifts)]
        @printf("  trusted radii: min %.2f, median %.2f, max %.2f
", extrema(s.radius for s in slices)[1],
                median(s.radius for s in slices), extrema(s.radius for s in slices)[2])
        rounds = 0
        while true
            holes = Coverage(slices, w, halfWidth)
            (isempty(holes) || rounds >= 5) && break
            added = HoleShifts(holes, 0.5 * median(s.radius for s in slices if s.radius > 0))
            @printf("  refinement round %d: %d uncovered grid points, %d new shifts
", rounds + 1,
                    length(holes), length(added))
            append!(slices, [merge((σ = σ,), Slice(M, σ; howmany = howmany, seed = 1000 + k))
                             for (k, σ) in enumerate(added)])
            rounds += 1
        end
    end
    holes = Coverage(slices, w, halfWidth)
    merged = filter(z -> real(z) >= -w && abs(imag(z)) <= halfWidth && abs(z) > 1e-8, MergeSlices(slices))

    found = count(z -> minimum(abs.(merged .- z)) < 1e-7, inside)
    spurious = count(z -> minimum(abs.(dense .- z)) > 1e-7, merged)

    @printf("N = %d, sector %d, window Re λ >= -%g: dense %d, sparse %d (%d found, %d spurious); %d slices, %.1f s, %d uncovered points\n",
            N, m, w, length(inside), length(merged), found, spurious, length(slices), time, length(holes))

    cuts = [c for c in (8.0, 12.0, 16.0) if c <= statisticsWindow]
    fullScan = WindowScan(SectorSpectra(p; verbose = false, only = [m]), p.L; cuts = -cuts,
                          edge = false, verbose = false)
    for c in cuts
        # a window with too few eigenvalues for the local unfolding has no statistics at all
        full = findfirst(r -> !isnothing(r.cut) && r.cut == -c, fullScan)
        windowDense = WindowStatistics(dense, [c]; covered = w, margin = margin)
        windowSparse = WindowStatistics(merged, [c]; covered = w, margin = margin)
        if isnothing(full) || isempty(windowDense) || isempty(windowSparse)
            @printf("  w = %4.1f  too few eigenvalues for statistics
", c)
            continue
        end
        rows = (("dense, neighbours from all", fullScan[full].summary),
                ("dense, window + margin", windowDense[1].summary),
                ("sparse, window + margin", windowSparse[1].summary))
        for (label, r) in rows
            @printf("  w = %4.1f  %-28s n = %5d  <|z|> = %.4f  -<cos> = %.4f  var(s) = %.4f\n", c,
                    label, r.count, r.meanRadius, r.minusCosine, r.varianceSpacing)
        end
    end

    return (found = found, expected = length(inside), spurious = spurious, holes = length(holes))
end


# ---------------------------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------------------------

""" Positional arguments and --key value options. """
function SplitArguments(arguments)
    positional = String[]
    options = Dict{String, String}()
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

function CommandLine(arguments)
    positional, options = SplitArguments(arguments)
    isempty(positional) && (println("usage: see the header of BHNumberConservingLiouvillianSparse.jl"); return)
    command = positional[1]
    option(key, default) = haskey(options, key) ? parse(Float64, options[key]) : default

    if command == "checks"
        SparseChecks()
    elseif command == "validate"
        Validate(parse(Int, positional[2]); w = option("w", 12.0), m = round(Int, option("m", 1)),
                 howmany = round(Int, option("howmany", 150)))
    elseif command == "plan"
        N = parse(Int, positional[3])
        WritePlan(positional[2], N, parse(Int, positional[4]), parse(Float64, positional[5]);
                  g = option("g", -20.0), η = option("eta", 3.0), κ = option("kappa", 0.3),
                  Γsym = option("Gsym", 0.0), howmany = round(Int, option("howmany", 300)),
                  spacing = haskey(options, "spacing") ? parse(Float64, options["spacing"]) : nothing,
                  halfWidth = option("halfwidth", 8.0 * N))
    elseif command == "slice"
        run = positional[2]
        index = length(positional) >= 3 ? parse(Int, positional[3]) :
                parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) + parse(Int, get(ENV, "ARRAY_OFFSET", "0"))
        perTask = parse(Int, get(ENV, "SHIFTS_PER_TASK", "1"))
        total = length(ReadPlan(run).shifts)
        RunSlices(run, [k for k in (index * perTask):(index * perTask + perTask - 1) if k < total])
    elseif command == "slices"
        # --part i --parts n: this process takes the shifts k with k mod n == i, so n processes
        # started side by side on the laptop share the work without a scheduler
        parts = round(Int, option("parts", 1))
        part = round(Int, option("part", 0))
        RunSlices(positional[2], [k for k in 0:(length(ReadPlan(positional[2]).shifts) - 1)
                                  if mod(k, parts) == part])
    elseif command == "merge"
        Merge(positional[2])
    elseif command == "refine"
        Refine(positional[2])
    elseif command == "statistics"
        values = ReadSpectrum(joinpath(RunPath(positional[2]), "spectrum.txt"))
        covered = ReadPlan(positional[2]).w
        results = WindowStatistics(values, [8.0, 12.0, 16.0]; covered = covered)
        isempty(results) && println("no cut w with w + 5 <= $covered: plan a wider window")
        for r in results
            @printf("w = %4.1f: n = %5d, var(s) = %.4f, <|z|> = %.4f, -<cos> = %.4f\n", r.w,
                    r.summary.count, r.summary.varianceSpacing, r.summary.meanRadius,
                    r.summary.minusCosine)
        end
    else
        error("unknown command $command")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    CommandLine(ARGS)
end
