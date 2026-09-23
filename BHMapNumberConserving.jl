# BHMapNumberConserving.jl
#
# Parameter map of the NUMBER-CONSERVING dissipative Bose-Hubbard model of
# number-conserving-BH.md, built on the model and the Lyapunov machinery of
# BHNumberConserving.jl.  Same Distributed/pmap structure, per-trajectory seeds and resumable
# per-cell output files as BHMapDrivenDissipative.jl; what differs is that the FULL Lyapunov
# spectrum is stored, not just the largest exponent, because question 1 of note §8 is decided by
# the Kaplan-Yorke dimension and not by the sign of a single number.
#
# Section 8 of the note asks four questions.  Section8() in BHNumberConserving.jl answers all four
# on a coarse grid in a couple of minutes; this file is the production version of the same scan -
# one plane, finely resolved, many initial conditions per cell.
#
# WHICH PLANE.  J = 1 fixes the time unit and only sign(J g) matters (§6), so g < 0 throughout and
# the interesting axis is always g.  SCAN picks the second one:
#
#   :eta         (g, η) at fixed κ > 0.  The main map.  The §5 threshold g_c(η) barely moves with η
#                (it is -(9 + η²/12 + 4κ²)/2 at L = 3), but chaos does not start there: a coarse
#                scan puts the chaotic region at η >~ 1 and |g| >~ 10, so the circulation is doing
#                something the linear theory of §5 does not see.  With η = 0 the reflection j -> -j
#                of §6 is a symmetry and no strange attractor was found at all.
#   :kappa       (g, κ) at fixed η > 0.  How much phase locking the attractor needs.  With
#                MODULATION = 0 the column κ = 0 is the volume-preserving limit of §4, where no
#                attractor can exist at all - a control row, and a place where Question 2 of §8
#                expects Σλ = 0 exactly rather than a negative number.
#   :modulation  (g, modulation) at κ = 0 and fixed η - QUESTION 2.  The rates are
#                eta_j = η (1 + modulation cos(2π j / L)), so modulation > 0 makes the column sums
#                c_j = eta_j - eta_{j-1} non-zero and div F = Σ_j c_j n_j can be negative on
#                average.  modulation = 0 is the volume-preserving column, which must come out with
#                Σλ = 0 exactly - it is the control, not a result.
#                CAUTION: with modulation != 0 the uniform condensate is no longer a solution (the
#                row sums of A no longer vanish), so the threshold of §5 does not apply anywhere in
#                this plane and the whole g range has to be searched.
#
# INITIAL CONDITIONS.  Drawn uniformly on the sphere Σ_j n_j = 1, i.e. from the Fubini-Study
# measure of CP^(L-1) - no absorbing ball has to be guessed as it did in the driven model, because
# the phase space IS the sphere (§4).  Each cell is sampled with TRAJECTORIES of them, so the map
# measures a BASIN distribution and not the fate of one particular starting point; that is what
# makes question 3 (multistability) readable off the same data.
#
# OUTPUT.  One line per trajectory, 2L + 9 tab-separated columns:
#
#   1 .. 2L    the full Lyapunov spectrum, descending
#   2L + 1     λ_max of the REDUCED spectrum (the full one minus the two exact symmetry zeros)
#   2L + 2     D_KY on the reduced 2L - 2 dimensional space; fractional means strange
#   2L + 3     <div F> on the attractor
#   2L + 4     Σλ - <div F>, the residual of the trace rule of §4 (question 4; should be ~1e-9)
#   2L + 5     attractor type: 0 fixed point, 1 limit cycle, 2 torus, 3 chaotic, 4 hyperchaotic,
#              5 neutral (the whole reduced spectrum is zero - a volume-preserving invariant set,
#              not an attractor), -1 undetermined
#   2L + 6     <Σ_j Re(conj(ψ_j) ψ_{j+1})>, the bond coherence (1 on the uniform locked state)
#   2L + 7     <max_j n_j>   (1/L uniform, ~1 self-trapped)
#   2L + 8     <Σ_j n_j²>    inverse participation ratio
#   2L + 9     <Σ_j Im(conj(ψ_j) ψ_{j+1})>, the current round the ring
#
# A failed trajectory is written as a line of NaN.  In numpy:
#   data = np.loadtxt(...); spectrum = data[:, :2*L]; lambdas = data[:, 2*L]; dky = data[:, 2*L+1]
#
# PATH also receives parameters.txt, which records the grid and the §5 threshold curve, so
# analyse_map_number_conserving.py needs nothing but the output directory.

using Distributed
using Printf

workers = 16

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere using Logging
@everywhere global_logger(ConsoleLogger(stderr, Logging.Warn))

@everywhere include("BHNumberConserving.jl")

# Constants and parameters
@everywhere const L = 3                 # L >= 3: at L = 2 the reduced space is the Bloch sphere and
                                        # Poincare-Bendixson forbids chaos (note §4).  The trimer is
                                        # the minimal candidate, L = 4 the first even ring (which has
                                        # a q = pi mode, and therefore a staggered instability).
@everywhere const J = 1.0               # time unit; only sign(J g) matters, and J g < 0 is needed
@everywhere const SCAN = :eta           # :eta | :kappa | :modulation - see the header

@everywhere const κ = 0.3               # phase locking; fixed when SCAN = :eta, swept for :kappa,
                                        # and forced to 0 for :modulation (question 2)
@everywhere const η = 3.0               # circulation; fixed when SCAN = :kappa or :modulation
@everywhere const MODULATION = 0.0      # rate non-uniformity; fixed unless SCAN = :modulation

@everywhere const TRAJECTORIES = 100    # random initial conditions per cell - the basin statistics

# The exponents must be measured ON the attractor: everything before RELAXATION_TIME is the
# approach.  The reduced spectrum is then accumulated over INTEGRATION_TIME - RELAXATION_TIME.
# 4000 time units give the near-zero exponents about three digits, which is what the classification
# into fixed point / limit cycle / torus / strange needs.
@everywhere const RELAXATION_TIME = 1000.0
@everywhere const INTEGRATION_TIME = 5000.0

# A limit cycle gives exactly one zero in the reduced spectrum, a torus two, a strange attractor a
# positive exponent.  At these constants the chaotic attractors run from +0.1 to +1.8, so the gap
# above CHAOS_THRESHOLD is wide; ZERO_THRESHOLD has to sit above the finite-time scatter of a true
# zero (~1e-3 over this window) and below the smallest genuine contraction rate.
@everywhere const CHAOS_THRESHOLD = 1e-2
@everywhere const ZERO_THRESHOLD = 5e-3

# Grid.  g is the first axis in every mode.
@everywhere const G_VALUES = LinRange(-50.0, -2.0, 241)

@everywhere const Y_VALUES =
    SCAN === :eta ? LinRange(0.0, 6.0, 121) :
    SCAN === :kappa ? LinRange(0.0, 1.5, 121) :
    SCAN === :modulation ? LinRange(0.0, 1.0, 101) :
    error("SCAN must be :eta, :kappa or :modulation")

@everywhere const G_STEP = step(G_VALUES)
@everywhere const Y_STEP = step(Y_VALUES)

# Fraction of the cell over which (g, y) are randomised per trajectory: 0 pins every trajectory to
# the grid point, 1 spreads them uniformly over the cell.  Exactly as in BHMapDrivenDissipative.jl,
# this turns a point sample into a cell average and stops narrow periodic windows from aliasing
# into isolated regular holes.  Set it to 0 when the attractor COUNT of one cell is what is wanted,
# since with JITTER > 0 the trajectories of a cell no longer share their parameters and genuine
# multistability gets mixed with the variation across the cell.
@everywhere const JITTER = 1.0

@assert RELAXATION_TIME < INTEGRATION_TIME "the spectrum is accumulated on (RELAXATION_TIME, INTEGRATION_TIME)"
@assert L >= 3 "L >= 3: chaos is impossible on the 2D reduced space of the dimer (note §4)"

const PATH = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "$L", string(SCAN),
             @sprintf("J_%.3f_k_%.3f_e_%.3f_m_%.3f", J, κ, η, MODULATION)))


# The model parameters of one point of the plane.  Only the swept quantity changes; the other two
# dissipative constants stay at their fixed values.  With SCAN = :modulation the phase locking is
# switched off altogether, which is what question 2 asks about.
@everywhere function CellParameters(g, y)
    if SCAN === :eta
        return NumberConservingParameters(L; J = J, g = g, κ = κ, η = y, modulation = MODULATION)
    elseif SCAN === :kappa
        return NumberConservingParameters(L; J = J, g = g, κ = y, η = η, modulation = MODULATION)
    else
        return NumberConservingParameters(L; J = J, g = g, κ = 0.0, η = η, modulation = y)
    end
end

@everywhere const CLASS_CODES = Dict(:fixedPoint => 0, :limitCycle => 1, :torus => 2,
                                     :chaotic => 3, :hyperchaotic => 4, :neutral => 5,
                                     :undetermined => -1)

# (g, y) of one trajectory: the cell centre displaced by up to half a cell in each direction, drawn
# from the same per-trajectory rng as the initial condition, so a resumed run keeps reproducing the
# same samples.  y is clamped to its physical range (rates and the modulation are non-negative).
@everywhere function JitterParameters(g, y, rng)
    JITTER <= 0 && return (g, y)

    return (g + JITTER * G_STEP * (rand(rng) - 0.5),
            max(y + JITTER * Y_STEP * (rand(rng) - 0.5), 0.0))
end

# One trajectory.  `index` is its global index within the cell file; it is mixed into the seed, so
# a resumed run continues with fresh samples instead of repeating the ones already stored.
@everywhere function SingleTrajectory(index, g, y)
    rng = Xoshiro(hash((g, y, index)))
    gJittered, yJittered = JitterParameters(g, y, rng)
    parameters = CellParameters(gJittered, yJittered)

    result = try
        LyapunovSpectrum(RandomInitialCondition(L, rng), parameters;
                         relaxationTime = RELAXATION_TIME, integrationTime = INTEGRATION_TIME,
                         zeroThreshold = ZERO_THRESHOLD, chaosThreshold = CHAOS_THRESHOLD, rng = rng)
    catch exception
        @warn "trajectory $index at (g, y) = ($gJittered, $yJittered) failed" exception
        return fill(NaN, 2 * L + 9)
    end

    return vcat(result.spectrum,
                [result.reduced[1], result.dimension, result.divergence,
                 sum(result.spectrum) - result.divergence,
                 float(CLASS_CODES[result.classification]),
                 result.coherence, result.maximum_n, result.ipr, result.current])
end


""" All trajectories of one cell.  `indices` are the global indices still missing from the output
    file, so a resumed run neither repeats nor reuses the samples already stored. """
function LyapunovMap(g, y, indices)
    time = @elapsed result = pmap(index -> SingleTrajectory(index, g, y), indices)

    λ = [row[2 * L + 1] for row in result]
    valid = filter(isfinite, λ)
    chaotic = filter(v -> v > CHAOS_THRESHOLD, valid)
    residuals = filter(isfinite, [abs(row[2 * L + 4]) for row in result])

    println("Finished g = $g, $SCAN = $y (L = $L, J = $J, κ = $κ, η = $η, modulation = $MODULATION)")
    println("New trajectories: $(length(result)) ($(length(valid)) valid)")

    if length(valid) > 0
        @printf("chaotic fraction = %.3f\n", length(chaotic) / length(valid))
    end

    if length(chaotic) > 0
        dimensions = [row[2 * L + 2] for row in result if row[2 * L + 1] > CHAOS_THRESHOLD]
        @printf("λ_max = %.4f, D_KY(red) = %.3f (out of %d)\n",
                mean(chaotic), mean(dimensions), 2 * L - 2)
    end

    if length(residuals) > 0
        @printf("worst |Σλ - <div F>| = %.1e\n", maximum(residuals))
    end

    @printf("Elapsed time: %.1f seconds\n\n", time)

    return result
end


""" Metadata next to the results: the grid, the constants and the §5 threshold curve g_c(y), so
    that analyse_map_number_conserving.py can overlay the linear-stability boundary without
    knowing anything about the model.  The threshold is meaningful only where the uniform state is
    still a solution, i.e. for a uniform circulation - with SCAN = :modulation it is written as NaN
    for every non-zero modulation. """
function WriteMetadata()
    mkpath(PATH)

    open(joinpath(PATH, "parameters.txt"), "w") do io
        println(io, "# number-conserving dissipative Bose-Hubbard map, BHMapNumberConserving.jl")
        println(io, "L\t$L")
        println(io, "J\t$J")
        println(io, "scan\t$SCAN")
        println(io, "kappa\t$κ")
        println(io, "eta\t$η")
        println(io, "modulation\t$MODULATION")
        println(io, "trajectories\t$TRAJECTORIES")
        println(io, "relaxationTime\t$RELAXATION_TIME")
        println(io, "integrationTime\t$INTEGRATION_TIME")
        println(io, "chaosThreshold\t$CHAOS_THRESHOLD")
        println(io, "jitter\t$JITTER")
        println(io, "columns\t$(2 * L + 9)")
        println(io, "gValues\t", join(G_VALUES, ","))
        println(io, "yValues\t", join(Y_VALUES, ","))

        thresholds = map(Y_VALUES) do y
            SCAN === :modulation && y != 0 && return NaN
            SCAN === :eta ? InstabilityThreshold(L; J = J, κ = κ, η = y) :
            SCAN === :kappa ? InstabilityThreshold(L; J = J, κ = y, η = η) :
                              InstabilityThreshold(L; J = J, κ = 0.0, η = η)
        end
        println(io, "threshold\t", join(thresholds, ","))
    end
end


WriteMetadata()
println("Results in $PATH\n")

for g in G_VALUES
    for y in Y_VALUES
        file = joinpath(PATH, @sprintf("%.4f_%.4f", g, y) * ".txt")
        trajectories = isfile(file) ? countlines(file) : 0

        println("Starting g = $g, $SCAN = $y")
        println("Computed trajectories: $trajectories, to compute: $(TRAJECTORIES - trajectories)")

        if trajectories >= TRAJECTORIES
            continue
        end

        result = LyapunovMap(g, y, (trajectories + 1):TRAJECTORIES)

        isempty(result) && continue

        open(file, "a") do io
            for row in result
                println(io, join(row, "\t"))
            end
        end
    end
end
