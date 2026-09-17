using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Printf

# Cluster-aware variant of BHMapDrivenDissipative.jl for the Chimera SLURM cluster.
#
# Unlike BHMapDrivenDissipative.jl (which parallelises the trajectories of one (Δ, f) cell within
# one job via Distributed/pmap), this version runs single-threaded and relies on SLURM to
# parallelise across (Δ, f) cells instead: one array task computes a contiguous block of at most
# CELLS_PER_TASK cells sequentially, so each task needs only one CPU. Batching many cells per task
# amortises the multi-minute Julia + DifferentialEquations startup over the whole block instead of
# paying it again for every single cell.
#
# It calculates the same map: constants, grid, per-trajectory seeds and output files are those of
# BHMapDrivenDissipative.jl. Every trajectory takes its initial condition, jitter and noise from a
# seed built from (Δ, f, index) alone, so running the trajectories serially instead of through pmap
# does not change them, and a result file started by either script can be continued by the other.
#
# The full (Δ, f) grid is flattened into a single 0-based index, Δ-major / f-minor (the loop order
# of BHMapDrivenDissipative.jl). With SLURM_ARRAY_TASK_ID set, task t computes the cells with flat
# index t*CELLS_PER_TASK .. (t+1)*CELLS_PER_TASK-1 (clamped to the grid); without it the script only
# prints the grid and exits. The array therefore needs
# cld(length(DELTA_VALUES) * length(F_VALUES), CELLS_PER_TASK) tasks -- the script prints the exact
# --array range on startup. ARRAY_OFFSET shifts the task id (counted in tasks, not cells), so the
# sweep can be split across several sbatch submissions if it is larger than SLURM's MaxArraySize
# (see BHMapDrivenDissipativeChimera.sh). This script never plots, so it skips loading
# Plots/ColorSchemes/PyPlot entirely -- see the CD_NO_PLOTS note in modules/ClassicalDynamics.jl.
#
# ---------------------------------------------------------------------------------------------------
#
# (Δ, f) phase diagram of the COHERENTLY DRIVEN dissipative Bose-Hubbard model of
# BHDissipative.jl - the discrete Lugiato-Lefever / driven Kerr lattice of
# bh_dissipation_driving.md.
#
# Why (Δ, f) and not (J, E) or κ:
#
#   * The whole analytic structure of the note lives in this plane, so the map can be overlaid with
#     closed-form curves: the bistability tongue of sec. 4 (folds at
#     n_± = [2Δ ± √(Δ² - ¾κ²)]/(3g), existing for gΔ > 0 and |Δ| > (√3/2)κ) and the modulational
#     instability boundary of sec. 5 (A_k² < g²n² - κ²/4). The route to chaos of sec. 7 - uniform ->
#     patterned fixed point -> limit cycle -> chaos - is exactly an f-sweep at fixed Δ.
#   * Δ (pump frequency) and f (pump power) are the knobs an experiment actually turns; κ is set by
#     the Q factor and J by the geometry.
#   * In sec. 5's band-sliding picture the two are complementary: Δ translates the band rigidly
#     across a fixed instability window, f sets the filling n and hence the width of that window.
#   * κ makes a poor axis: it enters only as κ²/4 in the state equation and as an exact -κ/2 shift of
#     EVERY exponent (Σλ = -κL), so it mostly contributes a trivial tilt plus one threshold. Take a
#     few slices in κ instead of sweeping it finely.
#
# TWO SIGN RULES fix which quadrant is worth computing, and the defaults below satisfy both:
#
#   J g < 0   (sec. 5) or the modulational instability merely saturates into a stable stationary
#             pattern, however strong it is - no chaos. Here J = -1, g = +2, so J g = -2.
#   g Δ > 0   because the detuning that matters is Δ - g n: with g Δ > 0 filling the ring pulls the
#             mode TOWARDS resonance, which is the runaway behind the S-curve, the folds and the MI.
#             With g Δ < 0 the response is self-limiting, n never reaches the |g| n > κ/2 that sec. 5
#             needs, and every trajectory ends on a fixed point at Λ = -κ/2. Since g = +2 > 0 here,
#             THE GRID MUST RUN OVER POSITIVE Δ. This parameter set is the mirror of the note's
#             g < 0, Δ < 0 case: the exact symmetry ψ -> ψ*, (Δ, g, J, f) -> -(Δ, g, J, f) of sec. 5
#             maps it onto J = +1, g = -2, Δ < 0.
#
# Initial conditions and multistability. There is no energy shell here: the driven flow has a
# compact absorbing ball Σ_i I_i ≤ 4 f² L / κ² and settles on an attractor, forgetting where it
# started except for WHICH attractor it reaches. The note reports strong multistability, so each
# (Δ, f) point is sampled with TRAJECTORIES random initial conditions drawn uniformly in norm inside
# that ball, and the observable is the distribution over the coexisting attractors:
#
#   λ_max               the exponent on the attractor each trajectory reached
#   chaotic fraction    share of initial conditions with λ_max > CHAOS_THRESHOLD, i.e. the share of
#                       phase space draining into a strange attractor - the driven counterpart of
#                       1 - freg. A limit cycle gives λ_max = 0, a fixed point λ_max ≤ 0.
#
# Output: one line per trajectory, TWO columns
#
#   λ_max    ⟨Σ_i I_i⟩ on the attractor
#
# The second column is the filling n L, which is what makes the bistability of sec. 4 visible (the
# two branches show up as two clusters at the same (Δ, f)). Failed trajectories are written as NaN
# NaN - there is no "-1 = no initial condition" case here, since no energy shell has to be hit. In
# numpy: data = np.loadtxt(...); lambdas = data[:, 0]; chaotic = lambdas > threshold.
ENV["CD_NO_PLOTS"] = "true"

using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))

# Pulls in models/BoseHubbardFull.jl and modules/ClassicalDynamics.jl as well
include("BHDissipative.jl")

# Constants and parameters
const TRAJECTORIES = 500        # random initial conditions per (Δ, f) point; the basin statistics
const L = 3                     # allowed modes are k = 2πm/3 only, i.e. k = 0 and a doubly
                                # degenerate k = 2π/3; k = π does not exist on an odd ring. Chaos
                                # survives that (Λ up to +0.71 at these constants)
const J = -1.0
const g = 2.0                  # the note's g; BHDissipative.jl takes U = g/2. J g < 0 is required
const U = g / 2

const κ = 1.0                   # net damping; must be > 0, it is what bounds the absorbing ball
const γ = 0.0                   # dephasing; 0 keeps the attractor deterministic (and the sweep fast)
const σ = 0.0                   # finite-N additive noise

# The trajectory has to reach its attractor before the exponent is measured: everything before
# RELAXATION_TIME describes the approach, not the attractor. The exponent is then averaged over
# INTEGRATION_TIME - RELAXATION_TIME.
const RELAXATION_TIME = 1000.0
const INTEGRATION_TIME = 4000.0

# A fixed point gives λ_max ≤ 0 (exactly -κ/2 on a uniform one), a limit cycle exactly 0, while the
# chaotic attractors at these constants run from +0.26 to +0.71 - a wide, clean gap.
const CHAOS_THRESHOLD = 0.01

# Grid. Chaos tracks the UPPER FOLD of the bistability tongue, which for g = 2, κ = 1 sits at
# f = 0.354, 1.459, 2.216, 3.077, 4.031 for Δ = 1, 3, 4, 5, 6 - so the chaotic region is a diagonal
# band rather than a rectangle, and the f window has to grow with Δ. Scans at these constants put
# chaos at (Δ, f) ≈ (1, 1.0), (3, 1.0-1.5), (4, 1.5-2.0), (5, 4.0-4.5), (6, 4.5-5.0), with Λ up to
# +0.71, and it dies out by f = 5.5 everywhere in this Δ range. Δ ∈ [0, 6] × f ∈ [0, 5.5] therefore
# holds the whole tongue up to Δ = 6 (upper fold 4.031) together with the chaotic band and a margin
# above it, at a 0.1 step in both directions.
#
# Kept as LinRanges (not collected): indexing them yields exactly the values the `for Δ in
# DELTA_VALUES` loop of BHMapDrivenDissipative.jl iterates over, and with them the same file names
# and seeds.
const DELTA_VALUES = LinRange(0.0, 10.0, 501)
const F_VALUES = LinRange(0.0, 7.0, 351)
const DELTA_STEP = step(DELTA_VALUES)
const F_STEP = step(F_VALUES)

# Fraction of the cell over which (Δ, f) are randomised, per trajectory: 0 pins every trajectory to
# the exact grid point (the plain point sample), 1 spreads them uniformly over the whole cell.
#
# This matters more than it looks. Chaotic bands and periodic windows alternate along f on a scale
# of 0.02-0.1 here, i.e. at or below the 0.1 grid step, so a point sample ALIASES them: a grid point
# landing inside a periodic window reports chaotic fraction 0 while its neighbours report 1, and
# narrow chaotic bands between grid points are missed entirely. Those isolated regular "holes" are
# real periodic windows, but where they fall is an artefact of the grid. Randomising within the cell
# turns the point sample into a cell AVERAGE, which is a smooth function of (Δ, f): a window
# narrower than a cell then shows up as an intermediate fraction instead of an all-or-nothing hole.
# Same idea as the `randomize` option of SolveEnergy in modules/ClassicalDynamics.jl, applied to the
# parameter axes instead of the section coordinates. It costs nothing - the trajectories are run
# either way.
#
# Caveat: with JITTER > 0 the trajectories of one cell no longer share the same parameters, so the
# attractor count of analyse_map_driven.py mixes genuine multistability with the variation across
# the cell. Set JITTER = 0 when that particular map is what you are after.
const JITTER = 1.0

# Splitting interval of the dephasing; only used when γ > 0 or σ > 0. Accuracy needs
# 2 U max(I) noiseStep < 0.2, and under driving max(I) is bounded by the absorbing ball rather than
# by 1 - TrajectoryLyapunovDissipative checks this and warns.
const NOISE_STEP = 0.01

@assert RELAXATION_TIME < INTEGRATION_TIME "the exponent is accumulated on (RELAXATION_TIME, INTEGRATION_TIME)"
@assert κ > 0 "a driven map needs κ > 0: it is what makes the absorbing ball compact"
@assert J * g < 0 "the note finds chaos only for J g < 0; with J g > 0 the instability saturates into a stationary pattern"

const PATH = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "driven", "$L",
             @sprintf("J_%.3f_g_%.3f_k_%.3f", J, g, κ)))

""" A random point inside the absorbing ball Σ_i I_i ≤ 4 f² L / κ²: a random direction in phase space,
    scaled to a norm drawn uniformly in (0, ball]. Sampling the ball rather than starting near the
    vacuum is what makes the chaotic fraction a basin measure instead of a property of one particular
    starting point. (At f = 0 the ball collapses to the origin and the fallback scale 1 is used; that
    column of the map is the trivial undriven one, where everything decays to ψ = 0 and Λ = -κ/2.) """
function RandomInitialCondition(parameters, rng)
    L, J, U, κ, γ, σ, Δ, f = parameters

    ballRadius = 4 * f^2 * L / κ^2
    x = randn(rng, 2 * L)
    targetNorm = rand(rng) * max(ballRadius, 1.0)       # Σ_i I_i of the initial condition
    x .*= sqrt(2 * targetNorm / sum(abs2, x))

    return x
end

""" (Δ, f) of one trajectory: the cell centre displaced by up to half a cell in each direction, drawn
    from the same per-trajectory rng as the initial condition so that a resumed run keeps reproducing
    the same samples. f is clamped at 0, which only bites in the bottom row of the grid. """
function JitterParameters(parameters, rng)
    JITTER <= 0 && return parameters

    Δ = parameters.Δ + JITTER * DELTA_STEP * (rand(rng) - 0.5)
    f = max(parameters.f + JITTER * F_STEP * (rand(rng) - 0.5), 0.0)

    return merge(parameters, (Δ = Δ, f = f))
end

""" One trajectory: a random initial condition inside the absorbing ball, integrated until it settles
    on an attractor. `index` is the global index of the trajectory within its output file; it is mixed
    into both the initial condition and the noise seed, so a resumed run continues with fresh samples
    instead of repeating the ones already stored. """
function SingleTrajectory(index, cellParameters)
    rng = Xoshiro(hash((cellParameters.Δ, cellParameters.f, index)))

    # jitter first: the absorbing ball the initial condition is drawn from depends on f
    parameters = JitterParameters(cellParameters, rng)
    initialCondition = RandomInitialCondition(parameters, rng)

    _, lyapunov, _, observables = TrajectoryLyapunovDissipative(initialCondition, parameters;
        seed=hash((cellParameters.Δ, cellParameters.f, index, :noise)),
        noiseStep=NOISE_STEP,
        relaxationTime=RELAXATION_TIME,
        timeInterval=(0.0, INTEGRATION_TIME))

    if lyapunov == 0.0
        return (NaN, NaN)               # nonconvergent, flagged by TrajectoryLyapunovDissipative
    end

    # Filling of the attractor, averaged over the post-transient part of the run
    norms = last.(observables.saveval)[observables.t .>= RELAXATION_TIME]

    return (lyapunov, isempty(norms) ? NaN : mean(norms))
end

""" All trajectories of one (Δ, f) point, one after another. `indices` are the global trajectory
    indices still missing from the output file, so a resumed run neither repeats nor reuses the
    samples already stored. """
function LyapunovMap(parameters, indices)
    time = @elapsed result = [SingleTrajectory(index, parameters) for index in indices]

    _, J, U, κ, γ, σ, Δ, f = parameters

    lyapunovs = first.(result)
    valid = filter(isfinite, lyapunovs)
    chaotic = filter(v -> v > CHAOS_THRESHOLD, valid)

    println("Finished Δ = $Δ, f = $f (L = $L, J = $J, g = $g, κ = $κ, γ = $γ, σ = $σ)")
    println("Number of new trajectories: $(length(result)) ($(length(valid)) valid)")

    if length(valid) > 0
        @printf("chaotic fraction = %.3f\n", length(chaotic) / length(valid))
    end

    if length(chaotic) > 0
        print("λ_max = $(mean(chaotic))")
        if length(chaotic) > 1
            print(" ± $(std(chaotic))")
        end
        println()
    end

    norms = filter(isfinite, last.(result))
    if length(norms) > 0
        @printf("⟨Σ I⟩ = %.3f (absorbing ball 4f²L/κ² = %.3f)\n", mean(norms), 4 * f^2 * L / κ^2)
    end

    println("Elapsed time: $time seconds")
    println()

    return result
end

# One SLURM array task computes a contiguous block of at most CELLS_PER_TASK cells from the
# flattened (Δ, f) grid.
const CELLS_PER_TASK = 20

const N_F = length(F_VALUES)
const TOTAL_CELLS = length(DELTA_VALUES) * N_F
const N_TASKS = cld(TOTAL_CELLS, CELLS_PER_TASK)

# Flat 0-based index -> (Δ, f), Δ-major / f-minor.
cellAt(k) = (DELTA_VALUES[k ÷ N_F + 1], F_VALUES[k % N_F + 1])

println("Grid: $(length(DELTA_VALUES)) Δ x $N_F f = $TOTAL_CELLS cells; ",
        "$CELLS_PER_TASK cells/task -> submit with --array=0-$(N_TASKS - 1)")

if !haskey(ENV, "SLURM_ARRAY_TASK_ID")
    println("SLURM_ARRAY_TASK_ID is not set, nothing to compute (set it by hand to run one block)")
    exit(1)
end

# ARRAY_OFFSET shifts the task id (counted in tasks, not cells), so the sweep can be split across
# several sbatch submissions (see BHMapDrivenDissipativeChimera.sh).
const TASK_ID = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) + parse(Int, get(ENV, "ARRAY_OFFSET", "0"))
const START_INDEX = TASK_ID * CELLS_PER_TASK       # SLURM array indices are 0-based

if START_INDEX >= TOTAL_CELLS
    @warn "Task index past the end of the grid; nothing to compute" TASK_ID START_INDEX TOTAL_CELLS
    exit(0)
end

const STOP_INDEX = min(START_INDEX + CELLS_PER_TASK, TOTAL_CELLS) - 1
println("Task $TASK_ID: flat cell indices $START_INDEX..$STOP_INDEX ($(STOP_INDEX - START_INDEX + 1) cells)")

mkpath(PATH)

# Δ and f label the CELL here; with JITTER > 0 each trajectory samples a random point inside it.
for k in START_INDEX:STOP_INDEX
    Δ, f = cellAt(k)

    file = joinpath(PATH, @sprintf("%.3f_%.3f", Δ, f) * ".txt")
    trajectories = isfile(file) ? countlines(file) : 0

    println("Starting Δ = $Δ, f = $f")
    println("Computed trajectories: $trajectories, trajectories to compute: $(TRAJECTORIES - trajectories)")

    if trajectories >= TRAJECTORIES
        continue
    end

    parameters = DissipativeParameters((L, J, U); κ=κ, γ=γ, σ=σ, Δ=Δ, f=f)
    result = LyapunovMap(parameters, (trajectories + 1):TRAJECTORIES)

    if length(result) == 0
        continue
    end

    open(file, "a") do io
        for (lyapunov, norm) in result
            println(io, "$lyapunov\t$norm")
        end
    end
end
