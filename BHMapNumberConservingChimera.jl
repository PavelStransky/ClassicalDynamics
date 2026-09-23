using LinearAlgebra
using Random
using Statistics
using Printf

# Cluster-aware variant of BHMapNumberConserving.jl for the Chimera SLURM cluster, in the same
# relation to it as BHMapDrivenDissipativeChimera.jl is to BHMapDrivenDissipative.jl.
#
# BHMapNumberConserving.jl parallelises the TRAJECTORIES of one cell inside a single job through
# Distributed/pmap. This version does the opposite: it runs single-threaded and lets SLURM
# parallelise across CELLS instead. One array task computes a contiguous block of at most
# CELLS_PER_TASK cells sequentially, so it needs one CPU and a modest amount of memory. Batching
# many cells per task amortises the multi-minute Julia + DifferentialEquations startup over the
# whole block rather than paying it again for every cell.
#
# The map computed is identical. Constants, grid, per-trajectory seeds and output files are those of
# BHMapNumberConserving.jl, and every trajectory draws its initial condition, its jitter and its
# deviation vectors from a seed built out of (g, y, index) alone - not out of the worker id or the
# order of execution. Running the trajectories serially therefore reproduces the distributed run
# trajectory for trajectory, and a result file started by either script can be continued by the
# other. THE CONSTANTS BELOW ARE A COPY: if you change L, the grid, TRAJECTORIES or the integration
# times in one file, change them in the other, or the two will silently fill the same directory with
# results from two different models.
#
# The (g, y) grid is flattened into a single 0-based index, g-major / y-minor - the loop order of
# BHMapNumberConserving.jl, so the two scripts walk the grid the same way. With SLURM_ARRAY_TASK_ID
# set, task t computes the cells with flat index t*CELLS_PER_TASK .. (t+1)*CELLS_PER_TASK-1 (clamped
# to the grid); without it the script prints the grid and the exact --array range and exits.
# ARRAY_OFFSET shifts the task id (counted in tasks, not cells), so a sweep larger than SLURM's
# MaxArraySize can be split across several submissions - see BHMapNumberConservingChimera.sh.
#
# Restart behaviour is what makes this usable on a preemptible partition: the number of finished
# trajectories of a cell is read back from its own output file, so a requeued task picks up where it
# left off and loses at most the cell it was inside. Resubmitting the whole array after a --time
# kill skips every finished cell.
#
# BHNumberConserving.jl loads no plotting packages at all, so unlike the driven scripts there is no
# CD_NO_PLOTS to set here.
#
# ---------------------------------------------------------------------------------------------------
#
# For the physics - which plane to scan and why, what the output columns mean, how the initial
# conditions are drawn and what the §5 threshold curve in parameters.txt is for - see the header of
# BHMapNumberConserving.jl. Nothing about the model differs here.

include("BHNumberConserving.jl")


# Constants and parameters - keep in step with BHMapNumberConserving.jl (see the header)
const L = 3
const J = 1.0
const SCAN = :eta                       # :eta | :kappa | :modulation

const κ = 0.3
const η = 3.0
const MODULATION = 0.0

const TRAJECTORIES = 100

const RELAXATION_TIME = 1000.0
const INTEGRATION_TIME = 5000.0

const CHAOS_THRESHOLD = 1e-2
const ZERO_THRESHOLD = 2e-3

const G_VALUES = LinRange(-50.0, -2.0, 121)

const Y_VALUES =
    SCAN === :eta ? LinRange(0.0, 6.0, 121) :
    SCAN === :kappa ? LinRange(0.0, 1.5, 151) :
    SCAN === :modulation ? LinRange(0.0, 1.0, 101) :
    error("SCAN must be :eta, :kappa or :modulation")

const G_STEP = step(G_VALUES)
const Y_STEP = step(Y_VALUES)

const JITTER = 1.0

@assert RELAXATION_TIME < INTEGRATION_TIME "the spectrum is accumulated on (RELAXATION_TIME, INTEGRATION_TIME)"
@assert L >= 3 "L >= 3: chaos is impossible on the 2D reduced space of the dimer (note §4)"

const PATH = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "$L", string(SCAN),
             @sprintf("J_%.3f_k_%.3f_e_%.3f_m_%.3f", J, κ, η, MODULATION)))


function CellParameters(g, y)
    if SCAN === :eta
        return NumberConservingParameters(L; J = J, g = g, κ = κ, η = y, modulation = MODULATION)
    elseif SCAN === :kappa
        return NumberConservingParameters(L; J = J, g = g, κ = y, η = η, modulation = MODULATION)
    else
        return NumberConservingParameters(L; J = J, g = g, κ = 0.0, η = η, modulation = y)
    end
end

const CLASS_CODES = Dict(:fixedPoint => 0, :limitCycle => 1, :torus => 2, :chaotic => 3,
                         :hyperchaotic => 4, :neutral => 5, :undetermined => -1)

function JitterParameters(g, y, rng)
    JITTER <= 0 && return (g, y)

    return (g + JITTER * G_STEP * (rand(rng) - 0.5),
            max(y + JITTER * Y_STEP * (rand(rng) - 0.5), 0.0))
end


# One trajectory. The seed depends on (g, y, index) only, so this is the very same trajectory the
# pmap version of BHMapNumberConserving.jl would have produced for the same index.
function SingleTrajectory(index, g, y)
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


""" All remaining trajectories of one cell, computed serially: this is LyapunovMap of
    BHMapNumberConserving.jl with the pmap replaced by a comprehension. """
function LyapunovMap(g, y, indices)
    time = @elapsed result = [SingleTrajectory(index, g, y) for index in indices]

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


""" parameters.txt next to the results, as in BHMapNumberConserving.jl, but written ATOMICALLY.

    Hundreds of array tasks start within seconds of each other and every one of them would write
    this file; two overlapping writes would leave a truncated grid behind and
    analyse_map_number_conserving.py would fail on it, or worse, silently read a short grid. Writing
    to a per-task temporary and renaming it over the target makes the visible file always complete:
    rename is atomic within one filesystem, and the tasks all write the same bytes anyway. """
function WriteMetadata()
    mkpath(PATH)

    target = joinpath(PATH, "parameters.txt")
    temporary = target * ".$(get(ENV, "SLURM_ARRAY_TASK_ID", "local")).$(getpid())"

    open(temporary, "w") do io
        println(io, "# number-conserving dissipative Bose-Hubbard map, BHMapNumberConservingChimera.jl")
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

    mv(temporary, target; force = true)

    return nothing
end


# One SLURM array task computes a contiguous block of at most CELLS_PER_TASK cells from the
# flattened (g, y) grid. A cell is TRAJECTORIES full spectra; on the desktop (i9-13900HX) one
# trajectory of the default integration window costs about 1 CPU-s at L = 3, so a cell is roughly
# 2 min and a block of 20 about 40 min, before the JULIA_CPU_TARGET=generic penalty.
const CELLS_PER_TASK = 5

const N_Y = length(Y_VALUES)
const TOTAL_CELLS = length(G_VALUES) * N_Y
const N_TASKS = cld(TOTAL_CELLS, CELLS_PER_TASK)

# Flat 0-based index -> (g, y), g-major / y-minor: the loop order of BHMapNumberConserving.jl.
cellAt(k) = (G_VALUES[k ÷ N_Y + 1], Y_VALUES[k % N_Y + 1])

println("Grid: $(length(G_VALUES)) g x $N_Y $SCAN = $TOTAL_CELLS cells; ",
        "$CELLS_PER_TASK cells/task -> submit with --array=0-$(N_TASKS - 1)")
println("Results in $PATH")

if !haskey(ENV, "SLURM_ARRAY_TASK_ID")
    println("SLURM_ARRAY_TASK_ID is not set, nothing to compute (set it by hand to run one block)")
    exit(1)
end

# ARRAY_OFFSET shifts the task id (counted in tasks, not cells), so the sweep can be split across
# several sbatch submissions (see BHMapNumberConservingChimera.sh).
const TASK_ID = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) + parse(Int, get(ENV, "ARRAY_OFFSET", "0"))
const START_INDEX = TASK_ID * CELLS_PER_TASK       # SLURM array indices are 0-based

if START_INDEX >= TOTAL_CELLS
    @warn "Task index past the end of the grid; nothing to compute" TASK_ID START_INDEX TOTAL_CELLS
    exit(0)
end

const STOP_INDEX = min(START_INDEX + CELLS_PER_TASK, TOTAL_CELLS) - 1
println("Task $TASK_ID: flat cell indices $START_INDEX..$STOP_INDEX ($(STOP_INDEX - START_INDEX + 1) cells)")

WriteMetadata()

# g and y label the CELL here; with JITTER > 0 each trajectory samples a random point inside it.
for k in START_INDEX:STOP_INDEX
    g, y = cellAt(k)

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
