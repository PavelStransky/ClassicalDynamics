# BHNumberConservingLiouvillianMap.jl
#
# Driver for the dense Liouvillian calculations of number-conserving-BH-paper-TODO.md, on IPNP36
# (640 GB, where it runs them up to N = 20) or on the laptop (32 GB, up to N = 16).  It
# runs many parameter points of BHNumberConservingLiouvillian.jl in parallel (Distributed, one BLAS
# thread per worker - LAPACK's nonsymmetric eigensolver gains almost nothing from threads, so the
# parallelism belongs across points), writes one line of statistics per point into a table and,
# where the task asks for it, the eigenvalues themselves, so that the analysis items (C4 bootstrap,
# K4 unfolding, K7 edge control, C10 form factor, K6 gap fits) never have to diagonalise again.
#
#   task        item   what
#   ----------  -----  ------------------------------------------------------------------------
#   plane       C1     (g, eta) plane at N = 12 on the grid of the classical map (121 x 121), the
#                      clean sector m = 1 only - about 15 CPU-s per point, 60 CPU-h in total
#   cut-eta3    C1 C6  eta = 3, g = -3 ... -50, dense (0.1) between g_c = -5.06 and -10;
#               K8     N = 8, 10, 12, 14; sectors 0 and 1 (gap, number of steady states)
#   cut-g20     C1 C6  g = -20, eta = 0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5; N = 8, 10, 12, 14
#   eta0        C2     eta = 0, g = -4 ... -50, reflection parity resolved; N = 8 ... 14 on every
#                      integer g, N = 16 on every second one
#   reference   C4 C6  the four reference points (-20, 3), (-20, 1), (-20, 0), (-4, 3) for every
#               K3 K8  N = 6 ... 16 - the spectra behind Tables IV and V and the gap fits of K6
#   kappa       K12    kappa = 0.1 and 0.5 at the chaotic point, N = 8 ... 16
#   kappa0      C9     kappa = 0 at the chaotic point, N = 8 ... 16
#   asymmetric  C11    Gamma_- = Gamma_+ / 2 (Gsym = eta / N) at the chaotic point, N = 8 ... 16
#   dsff        C10    N = 12, 41 values of g in [-20.5, -19.5] at eta = 3: the ensemble over
#                      which the dissipative spectral form factor is averaged
#
# Usage (from the repository root; everything is resumable - points already in the table are
# skipped, so an interrupted run is simply restarted):
#
#   julia BHNumberConservingLiouvillianMap.jl <task> [--workers n] [--memory GB] [--N 8,10] [--stride k]
#
# --workers  number of worker processes (default 20, capped at the number of CPU threads - 2)
# --memory   RAM the running points may hold together, in GB (default: 70% of the machine).  Each
#            point carries an estimate (three dense copies of its largest sector, five with the
#            parity split), and a point waits until it fits, so the N = 16 points of eta0 run a
#            few at a time while the small ones fill the remaining workers.
# --N        the boson numbers, REPLACING the task's own list (the lists above are sized for the
#            laptop; IPNP36 runs e.g. `reference --N 17,18,19,20`).  A sector at N = 20 is 17787
#            wide: 5 GB per dense copy, about 2.5 CPU-h to diagonalise.
# --stride   keep every k-th value of the task's g axis (and of the η axis of the plane), e.g.
#            `plane --N 14 --stride 2` for the plane at N = 14 on a 61 x 61 grid
#
# OUTPUT in $BH_RESULTS_DIR (default ~/results/bh/number-conserving/quantum/3):
#
#   <task>.txt                  one line per point, tab separated, with a header of column names
#                               (np.genfromtxt(file, names=True) reads it):
#       N g eta kappa Gsym seconds zeros gap lead_re lead_im
#       then for every sector label S in (m1, m0, m0p, m0m) - the clean sector q = 2 pi/3, the
#       self-conjugate q = 0, and its two reflection parities where they are resolved:
#       S_dim S_real (fraction of exactly real eigenvalues, the fingerprint of an antiunitary
#       symmetry) and for every window W in (bulk, w8, w12, w16):
#       S_W_n S_W_var S_W_small S_W_r S_W_cos    (count, var(s), P(s < 1/2), <|z|>, -<cos arg z>)
#       - bulk is the densest half of the sector (keep = 0.5, K3), wX the slow region Re λ >= -X.
#       Missing sectors and windows with too few eigenvalues are NaN.  zeros counts |λ| < 1e-8 in
#       the q = 0 sector (K8: must be 1), gap = -max Re λ over the decaying modes of every sector
#       (NaN when the q = 0 sector was not computed); lead_* is the slowest decaying mode found.
#   spectra/<task>/L3_N<N>_g<g>_e<eta>_k<kappa>_s<Gsym>_<S>.txt
#                               two columns (Re, Im) per sector - for the tasks that save them

using Distributed
using Printf

const TASKS = ("plane", "cut-eta3", "cut-g20", "eta0", "reference", "kappa", "kappa0",
               "asymmetric", "dsff")

const RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "quantum", "3"))

const SECTOR_LABELS = ("m1", "m0", "m0p", "m0m")
const WINDOWS = ("bulk", "w8", "w12", "w16")
const STATISTICS = ("n", "var", "small", "r", "cos")

ColumnNames() = vcat(["N", "g", "eta", "kappa", "Gsym", "seconds", "zeros", "gap", "lead_re", "lead_im"],
                     [s * "_" * c for s in SECTOR_LABELS for c in vcat(["dim", "real"],
                         [w * "_" * x for w in WINDOWS for x in STATISTICS])])


# ---------------------------------------------------------------------------------------------
# The points of every task
# ---------------------------------------------------------------------------------------------

""" One unit of work: a parameter point at one N, which sectors to diagonalise (`only = nothing`
    means all, i.e. q = 0 and q = 2 pi / 3, the third being the conjugate of the second), whether
    to resolve the reflection parity, and whether to save the eigenvalues. """
Point(N, g, η; κ = 0.3, Γsym = 0.0, only = nothing, parity = false, save = true) =
    (N = N, g = float(g), η = float(η), κ = float(κ), Γsym = float(Γsym), only = only,
     parity = parity, save = save)

const CHAOTIC = (-20.0, 3.0)
const REFERENCE_POINTS = ((-20.0, 3.0), (-20.0, 1.0), (-20.0, 0.0), (-4.0, 3.0))

""" The points of a task.  `Ns` replaces the task's own list of boson numbers (the default lists
    are the ones sized for the laptop; IPNP36 runs the same tasks at larger N), and `stride` keeps
    every stride-th value of the task's g axis (and of the η axis of the plane) - so that, e.g., the
    plane can be done at N = 14 on every second point of the classical grid. """
function TaskPoints(task; Ns = nothing, stride = 1)
    Pick(values) = collect(values)[1:stride:end]
    Numbers(default) = isnothing(Ns) ? default : Ns

    if task == "plane"
        return [Point(N, g, η; only = [1], save = false) for N in Numbers([12])
                for g in Pick(LinRange(-50.0, -2.0, 121)) for η in Pick(LinRange(0.0, 6.0, 121))]

    elseif task == "cut-eta3"
        gs = vcat([-3.0, -4.0, -4.5, -5.0], collect(-5.1:-0.1:-10.0), collect(-11.0:-1.0:-20.0),
                  collect(-22.0:-2.0:-30.0), collect(-35.0:-5.0:-50.0))
        return [Point(N, g, 3.0) for N in Numbers([8, 10, 12, 14]) for g in Pick(gs)]

    elseif task == "cut-g20"
        return [Point(N, -20.0, η; parity = true) for N in Numbers([8, 10, 12, 14])
                for η in (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0)]

    elseif task == "eta0"
        isnothing(Ns) && stride == 1 &&
            return vcat([Point(N, g, 0.0; parity = true) for N in (8, 10, 12, 14) for g in -4.0:-1.0:-50.0],
                        [Point(16, g, 0.0; parity = true) for g in -4.0:-2.0:-50.0])
        return [Point(N, g, 0.0; parity = true) for N in Numbers([8, 10, 12, 14, 16])
                for g in Pick(-4.0:-1.0:-50.0)]

    elseif task == "reference"
        return [Point(N, g, η; parity = true) for (g, η) in REFERENCE_POINTS for N in Numbers(6:16)]

    elseif task == "kappa"
        return [Point(N, CHAOTIC...; κ = κ) for κ in (0.1, 0.5) for N in Numbers(8:16)]

    elseif task == "kappa0"
        return [Point(N, CHAOTIC...; κ = 0.0) for N in Numbers(8:16)]

    elseif task == "asymmetric"
        # Gamma_+ - Gamma_- = eta / N is kept; the added symmetric rate Gsym = eta / N makes
        # Gamma_- = eta / N = Gamma_+ / 2.  JumpOperators adds Gsym to both directions unscaled.
        return [Point(N, CHAOTIC...; Γsym = CHAOTIC[2] / N) for N in Numbers(8:16)]

    elseif task == "dsff"
        return [Point(N, g, 3.0; only = [1]) for N in Numbers([12]) for g in LinRange(-20.5, -19.5, 41)]
    end

    error("unknown task $task; one of $(join(TASKS, ", "))")
end

PointKey(N, g, η, κ, Γsym) = @sprintf("%d_%.4f_%.4f_%.4f_%.6f", N, g, η, κ, Γsym)
PointKey(point) = PointKey(point.N, point.g, point.η, point.κ, point.Γsym)

""" Memory a point holds while it runs, in GB: dense copies of its largest sector (block, LAPACK
    workspace and Hessenberg copy; the parity split adds the reflection matrix and its
    eigenvectors) plus the Julia process itself. """
function PointMemory(point)
    n = binomial(point.N + 2, 2)^2 / 3
    copies = point.parity && point.η == 0 ? 5 : 3
    return copies * 16 * n^2 / 1e9 + 0.4
end

PointCost(point) = (binomial(point.N + 2, 2)^2 / 3)^3 * (isnothing(point.only) ? 2 : length(point.only))


# ---------------------------------------------------------------------------------------------
# Worker side
# ---------------------------------------------------------------------------------------------

function Arguments(arguments)
    options = Dict("workers" => string(clamp(Sys.CPU_THREADS - 2, 1, 20)),
                   "memory" => string(round(0.7 * Sys.total_memory() / 1e9, digits = 1)),
                   "N" => "", "stride" => "1")
    task = ""
    i = 1
    while i <= length(arguments)
        if startswith(arguments[i], "--")
            options[arguments[i][3:end]] = arguments[i + 1]
            i += 2
        else
            task = arguments[i]
            i += 1
        end
    end
    return task, options
end

const PARSED = Arguments(ARGS)
const TASK = PARSED[1]
const OPTIONS = PARSED[2]
isempty(TASK) && (println("usage: julia BHNumberConservingLiouvillianMap.jl <task> [--workers n] [--memory GB] [--N 8,10] [--stride k]");
                  println("tasks: ", join(TASKS, ", ")); exit(1))

ENV["CD_NO_PLOTS"] = "true"          # inherited by the workers; no plotting package is loaded
addprocs(parse(Int, OPTIONS["workers"]))

# The packages are loaded in separate statements: a block that loaded Printf and used @printf
# would be macro-expanded before the `using` had run.
@everywhere using LinearAlgebra
@everywhere using Printf
@everywhere BLAS.set_num_threads(1)
@everywhere include(joinpath($(@__DIR__), "BHNumberConservingLiouvillian.jl"))

@everywhere begin

    const WINDOW_CUTS = Dict("bulk" => nothing, "w8" => -8.0, "w12" => -12.0, "w16" => -16.0)

    """ (n, var(s), P(s < 1/2), <|z|>, -<cos>) of one sector in one window, NaN when the window
        holds too few eigenvalues for the local unfolding (4 x neighbours). """
    function WindowRow(values, window; neighbours = 30, keep = 0.5)
        cut = WINDOW_CUTS[window]
        decaying = count(z -> abs(z) > 1e-8 && (isnothing(cut) || real(z) >= cut), values)
        (decaying <= 4 * neighbours || length(values) <= 4 * neighbours + 1) && return fill(NaN, 5)

        spacings = NearestNeighbourSpacings(values; neighbours = neighbours, keep = keep, window = cut)
        ratios = ComplexSpacingRatios(values; neighbours = neighbours, keep = keep, window = cut)
        summary = StatisticsSummary(spacings, ratios)
        return [summary.count, summary.varianceSpacing, summary.smallSpacings,
                summary.meanRadius, summary.minusCosine]
    end

    SectorLabel(s) = s.parity == 1 ? "m0p" : s.parity == -1 ? "m0m" : "m$(s.m)"

    SpectrumFile(directory, point, label) = joinpath(directory,
        @sprintf("L3_N%02d_g%.4f_e%.4f_k%.4f_s%.6f_%s.txt", point.N, point.g, point.η, point.κ,
                 point.Γsym, label))

    """ Everything about one point, as the list of values of one table line. """
    function ComputePoint(point, spectraDirectory, labels, windows)
        p = LiouvillianParameters(3, point.N; g = point.g, η = point.η, κ = point.κ,
                                  Γsym = point.Γsym)
        seconds = @elapsed sectors = SectorSpectra(p; verbose = false, only = point.only,
                                                   parity = point.parity)

        # q = 4 pi / 3 is the complex conjugate of q = 2 pi / 3 - a duplicate, not a new sample
        sectors = [s for s in sectors if s.m <= p.L - s.m || s.m == 0]

        hasZero = any(s.m == 0 for s in sectors)
        zeroCount = hasZero ? sum(count(z -> abs(z) < 1e-8, s.values) for s in sectors if s.m == 0) : -1
        decaying = vcat([filter(z -> abs(z) > 1e-8, s.values) for s in sectors]...)
        lead = decaying[argmax(real.(decaying))]
        gap = hasZero ? -real(lead) : NaN

        columns = Dict{String, Vector{Float64}}()
        for s in sectors
            label = SectorLabel(s)
            real_ = count(z -> abs(imag(z)) < 1e-9 * max(1.0, abs(z)), s.values) / length(s.values)
            columns[label] = vcat([length(s.values), real_],
                                  vcat([WindowRow(s.values, w) for w in windows]...))

            if point.save
                open(SpectrumFile(spectraDirectory, point, label), "w") do io
                    for z in s.values
                        @printf(io, "%.15e\t%.15e\n", real(z), imag(z))
                    end
                end
            end
        end

        width = 2 + 5 * length(windows)
        row = [point.N, point.g, point.η, point.κ, point.Γsym, seconds, zeroCount, gap,
               real(lead), imag(lead)]
        for label in labels
            append!(row, get(columns, label, fill(NaN, width)))
        end
        return row
    end
end


# ---------------------------------------------------------------------------------------------
# Master: queue with a memory budget, resumable table
# ---------------------------------------------------------------------------------------------

function ReadDone(table)
    done = Set{String}()
    isfile(table) || return done
    for line in eachline(table)
        startswith(line, "#") && continue
        fields = split(line, '\t')
        length(fields) < 5 && continue
        push!(done, PointKey(round(Int, parse(Float64, fields[1])), parse(Float64, fields[2]),
                             parse(Float64, fields[3]), parse(Float64, fields[4]),
                             parse(Float64, fields[5])))
    end
    return done
end

FormatValue(x) = isnan(x) ? "NaN" : isinteger(x) && abs(x) < 1e9 ? string(Int(x)) : @sprintf("%.10g", x)

function Run()
    mkpath(RESULTS)
    table = joinpath(RESULTS, "$TASK.txt")
    spectraDirectory = joinpath(RESULTS, "spectra", TASK)
    mkpath(spectraDirectory)

    points = TaskPoints(TASK; Ns = isempty(OPTIONS["N"]) ? nothing : parse.(Int, split(OPTIONS["N"], ',')),
                        stride = parse(Int, OPTIONS["stride"]))

    done = ReadDone(table)
    queue = sort([p for p in points if !(PointKey(p) in done)], by = PointCost, rev = true)
    budget = parse(Float64, OPTIONS["memory"])

    @printf("Task %s: %d points, %d already done, %d to compute on %d workers (memory budget %.1f GB)\n",
            TASK, length(points), length(points) - length(queue), length(queue), nworkers(), budget)
    @printf("Table: %s\n", table)
    isempty(queue) && return

    if !isfile(table)
        open(table, "w") do io
            println(io, "# ", join(ColumnNames(), "\t"))
        end
    end

    inflight = Ref(0.0)
    released = Condition()
    finished = Ref(0)
    total = length(queue)
    start = time()
    radiusColumn = findfirst(==("m1_bulk_r"), ColumnNames())

    @sync for worker in workers()
        @async while true
            point = nothing
            while !isempty(queue)
                index = findfirst(p -> inflight[] + PointMemory(p) <= budget || inflight[] == 0, queue)
                if isnothing(index)
                    wait(released)
                    continue
                end
                point = popat!(queue, index)
                inflight[] += PointMemory(point)
                break
            end
            isnothing(point) && break

            try
                row = remotecall_fetch(ComputePoint, worker, point, spectraDirectory,
                                       SECTOR_LABELS, WINDOWS)
                open(table, "a") do io
                    println(io, join(FormatValue.(row), "\t"))
                end
                finished[] += 1
                elapsed = time() - start
                @printf("[%5d/%d] N = %2d, g = %8.3f, η = %.3f, κ = %.2f: %6.1f s, gap %.4f, clean <|z|> %.4f  (elapsed %.2f h)\n",
                        finished[], total, point.N, point.g, point.η, point.κ, row[6], row[8],
                        row[radiusColumn], elapsed / 3600)
            catch exception
                if !(worker in workers())
                    # the worker died (typically out of memory): give the point back and retire
                    # this slot rather than letting it fail every remaining point in turn
                    @warn "worker $worker died; point returned to the queue" point
                    push!(queue, point)
                    break
                end
                @warn "point failed" point exception
            finally
                inflight[] -= PointMemory(point)
                notify(released; all = true)
            end
        end
    end

    println("done")
end

Run()
