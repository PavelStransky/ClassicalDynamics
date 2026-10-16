# BHNumberConservingLiouvillianAnalysis.jl
#
# Analysis of the Liouvillian spectra that BHNumberConservingLiouvillianMap.jl saved - nothing is
# diagonalised again.  Items of number-conserving-BH-paper-TODO.md:
#
#   tables     C4   the statistics of Tables IV and V (bulk and slow windows, clean sector) at the
#                   reference points, with two error bars: a bootstrap over the eigenvalues, and
#                   the spread of the same pipeline over an ensemble of 50 Ginibre matrices and of
#                   50 sets of uniform points of matched size (the realistic error - neighbouring
#                   ratios are correlated, which a plain bootstrap ignores and underestimates).
#   unfolding  K4   var(s) and P(s < 1/2) for neighbours = 10 ... 50 and keep = 0.3, 0.5, 0.7; the
#                   ratios need no unfolding and are printed alongside as the fixed reference.
#   edge       K7   the slow-window statistics of a Ginibre-type cloud cut by a HALF-PLANE
#                   Re λ >= -w, for elliptic clouds of aspect ratio 1 (the disk cap used so far),
#                   3, 7 (about that of the Liouvillian cloud at L = 3) and 15, at the counts of the
#                   Liouvillian windows - the edge control and its dependence on the geometry.
#   dsff       C10  the dissipative spectral form factor of Li, Prosen and Chan (PRL 127, 170602
#                   (2021)), averaged over the 41 spectra of task dsff (g ∈ [-20.5, -19.5], N = 12)
#                   and over directions in the complex plane, against Ginibre and Poisson through
#                   the same pipeline.
#   classes    A1   reference values of <|z|>, -<cos arg z> and var(s) through the same pipeline for
#                   the classes that can occur: GinUE (A), real Ginibre (AI, a K-type antiunitary -
#                   what rho -> rho^dag or reflection∘dag gives), complex symmetric (AI^dag, a
#                   transposition-type one); for AI both over the whole bulk and away from the real
#                   axis, where the eigenvalues accumulate.
#   moving     -    the statistics in a band of width W of Re λ sliding along the real axis, with an
#                   elliptic-Ginibre control cut at the same quantiles of Re λ (MovingWindow).
#                   Options --g --eta (comma lists; every combination), --kappa (default -20, 5,
#                   0.3), --width W (2), --step (W/2), --sector m1 | m0, --minimum (30),
#                   --realizations (8, 0 = none), --task (default: search every spectra directory;
#                   a point no task has saved is computed by BHNumberConservingLiouvillianMap.jl points).
#
# USAGE   julia -t 8 BHNumberConservingLiouvillianAnalysis.jl <command> [--N 14,16]
# INPUT   $BH_RESULTS_DIR/spectra/<task>/... (default ~/results/bh/number-conserving/quantum/3)
# OUTPUT  $BH_RESULTS_DIR/analysis/<command>.txt

ENV["CD_NO_PLOTS"] = "true"
isdefined(Main, :LiouvillianSector) || include(joinpath(@__DIR__, "BHNumberConservingLiouvillian.jl"))

using LinearAlgebra
using Random
using Statistics
using Printf

const ANALYSIS_RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "quantum", "3"))

const CUTS = Dict("bulk" => nothing, "w8" => -8.0, "w12" => -12.0, "w16" => -16.0)
const WINDOW_ORDER = ("bulk", "w8", "w12", "w16")

SpectrumPath(task, N, g, η; κ = 0.3, Γsym = 0.0, label = "m1") =
    joinpath(ANALYSIS_RESULTS, "spectra", task,
             @sprintf("L3_N%02d_g%.4f_e%.4f_k%.4f_s%.6f_%s.txt", N, g, η, κ, Γsym, label))

""" The saved spectrum of a point from whichever task computed it (`task = nothing`: every task
    directory is searched), or nothing. """
function FindSpectrum(task, N, g, η; κ = 0.3, label = "m1")
    tasks = isnothing(task) ? readdir(joinpath(ANALYSIS_RESULTS, "spectra")) : [task]
    for t in tasks
        file = SpectrumPath(t, N, g, η; κ = κ, label = label)
        isfile(file) && return file
    end
    return nothing
end

function ReadValues(file)
    values = ComplexF64[]
    for line in eachline(file)
        fields = split(line, '\t')
        push!(values, complex(parse(Float64, fields[1]), parse(Float64, fields[2])))
    end
    return values
end

function WriteRows(name, header, rows)
    directory = joinpath(ANALYSIS_RESULTS, "analysis")
    mkpath(directory)
    file = joinpath(directory, name)
    open(file, "w") do io
        println(io, "# ", join(header, "\t"))
        for row in rows
            println(io, join([x isa AbstractString ? x : x isa Integer ? string(x) : @sprintf("%.6g", x)
                              for x in row], "\t"))
        end
    end
    println("wrote ", file)
end


# ---------------------------------------------------------------------------------------------
# Per-eigenvalue quantities, so that resampling does not repeat the neighbour search
# ---------------------------------------------------------------------------------------------

""" Unfolded spacing, complex ratio and local scale of EVERY eigenvalue (the steady state removed),
    with the neighbours found in the full set - a selection is then only an index into these. """
function PerEigenvalue(values; neighbours = 30)
    values = filter(z -> abs(z) > 1e-8, values)
    distances, nearest, next = NeighbourData(values, neighbours)
    scale = distances[:, neighbours] ./ sqrt(neighbours / pi)
    ratios = [(values[nearest[i]] - values[i]) / (values[next[i]] - values[i]) for i in eachindex(values)]
    return (values = values, scale = scale, spacings = distances[:, 1] ./ scale, ratios = ratios)
end

""" Unfolded spacings and complex ratios of the selected eigenvalues (bulk: the `keep` densest;
    window: Re λ >= cut), exactly as NearestNeighbourSpacings / ComplexSpacingRatios select them. """
function Selected(values, cut; neighbours = 30, keep = 0.5)
    e = PerEigenvalue(values; neighbours = neighbours)
    index = isnothing(cut) ? sortperm(e.scale)[1:max(1, round(Int, keep * length(e.values)))] :
                             findall(z -> real(z) >= cut, e.values)
    return e.spacings[index], e.ratios[index]
end

Statistic(spacings, ratios) = (r = mean(abs, ratios), cos = -mean(cos ∘ angle, ratios),
                               var = var(spacings ./ mean(spacings)), n = length(ratios))

function Bootstrap(spacings, ratios; samples = 1000, rng = Xoshiro(1))
    n = length(ratios)
    draws = [Statistic(spacings[i], ratios[i]) for i in (rand(rng, 1:n, n) for _ = 1:samples)]
    return (r = std(d.r for d in draws), cos = std(d.cos for d in draws), var = std(d.var for d in draws))
end

Ginibre(n, rng) = eigvals(randn(rng, ComplexF64, n, n) ./ sqrt(n))

""" Spread of the pipeline over `realizations` Ginibre matrices of dimension n, for the same kind
    of selection (bulk keep = 0.5, or the `count` eigenvalues of largest real part - the half-plane
    cut of the disk).  Returns mean and standard deviation of each statistic. """
function Calibration(n, count_; realizations = 50, bulk = true, seed = 7)
    draws = []
    for k = 1:realizations
        rng = Xoshiro(hash((seed, n, count_, k)))
        values = Ginibre(n, rng)
        cut = bulk ? nothing : sort(real.(values), rev = true)[min(count_, n)]
        push!(draws, Statistic(Selected(values, cut)...))
    end
    return (mean = (r = mean(d.r for d in draws), cos = mean(d.cos for d in draws), var = mean(d.var for d in draws)),
            std = (r = std(d.r for d in draws), cos = std(d.cos for d in draws), var = std(d.var for d in draws)))
end

function PoissonCalibration(count_; realizations = 50, seed = 9)
    draws = []
    for k = 1:realizations
        rng = Xoshiro(hash((seed, count_, k)))
        points = [2 * (rand(rng) - 0.5) + 2im * (rand(rng) - 0.5) for _ = 1:(2 * count_)]
        push!(draws, Statistic(Selected(points, nothing)...))
    end
    return (mean = (r = mean(d.r for d in draws), cos = mean(d.cos for d in draws), var = mean(d.var for d in draws)),
            std = (r = std(d.r for d in draws), cos = std(d.cos for d in draws), var = std(d.var for d in draws)))
end


# ---------------------------------------------------------------------------------------------
# C4 tables with errors
# ---------------------------------------------------------------------------------------------

function Tables(Ns)
    rows = []
    for (g, η) in ((-20.0, 3.0), (-20.0, 1.0), (-20.0, 0.0), (-4.0, 3.0)), N in Ns
        file = SpectrumPath("reference", N, g, η)
        isfile(file) || (println("missing ", file); continue)
        values = ReadValues(file)

        for window in WINDOW_ORDER
            spacings, ratios = Selected(values, CUTS[window])
            length(ratios) > 30 || continue
            s = Statistic(spacings, ratios)
            b = Bootstrap(spacings, ratios)
            gin = window == "bulk" ? Calibration(length(values), 0) :
                                     Calibration(length(values), length(ratios); bulk = false)
            poi = PoissonCalibration(length(ratios))
            @printf("(g, η) = (%5.1f, %.1f) N = %2d %-4s n = %5d  <|z|> = %.4f ± %.4f (boot) ± %.4f (Gin)   -<cos> = %.4f ± %.4f ± %.4f   var(s) = %.4f ± %.4f ± %.4f\n",
                    g, η, N, window, s.n, s.r, b.r, gin.std.r, s.cos, b.cos, gin.std.cos, s.var, b.var, gin.std.var)
            @printf("%52s Ginibre same pipeline: %.4f            -<cos> %.4f             var %.4f;  Poisson: %.4f / %.4f / %.4f\n",
                    "", gin.mean.r, gin.mean.cos, gin.mean.var, poi.mean.r, poi.mean.cos, poi.mean.var)
            push!(rows, (g, η, N, window, s.n, s.r, b.r, gin.std.r, s.cos, b.cos, gin.std.cos, s.var, b.var,
                         gin.std.var, gin.mean.r, gin.mean.cos, gin.mean.var, poi.mean.r, poi.mean.cos, poi.mean.var,
                         poi.std.r, poi.std.cos, poi.std.var))
        end
    end
    WriteRows("tables.txt", ["g", "eta", "N", "window", "n", "r", "r_boot", "r_ginibre_spread", "cos",
                             "cos_boot", "cos_ginibre_spread", "var", "var_boot", "var_ginibre_spread",
                             "r_ginibre", "cos_ginibre", "var_ginibre", "r_poisson", "cos_poisson",
                             "var_poisson", "r_poisson_spread", "cos_poisson_spread", "var_poisson_spread"], rows)
end


# ---------------------------------------------------------------------------------------------
# K4 unfolding stability
# ---------------------------------------------------------------------------------------------

function Unfolding(Ns)
    rows = []
    for (g, η) in ((-20.0, 3.0), (-4.0, 3.0)), N in Ns
        file = SpectrumPath("reference", N, g, η)
        isfile(file) || (println("missing ", file); continue)
        values = ReadValues(file)
        for keep in (0.3, 0.5, 0.7), neighbours in (10, 20, 30, 40, 50)
            spacings, ratios = Selected(values, nothing; neighbours = neighbours, keep = keep)
            s = Statistic(spacings, ratios)
            small = count(<(0.5), spacings ./ mean(spacings)) / length(spacings)
            @printf("(g, η) = (%5.1f, %.1f) N = %2d keep %.1f neighbours %2d: var(s) = %.4f  P(s<1/2) = %.4f   <|z|> = %.4f  -<cos> = %.4f\n",
                    g, η, N, keep, neighbours, s.var, small, s.r, s.cos)
            push!(rows, (g, η, N, keep, neighbours, s.var, small, s.r, s.cos))
        end
    end
    println("the ratio columns change only through the selection (keep), never through the unfolding")
    WriteRows("unfolding.txt", ["g", "eta", "N", "keep", "neighbours", "var", "P_small", "r", "cos"], rows)
end


# ---------------------------------------------------------------------------------------------
# K7 edge control with a half-plane cut of clouds of different shape
# ---------------------------------------------------------------------------------------------

""" Eigenvalues of an elliptic Ginibre-type matrix whose cloud is `aspect` times longer along the
    imaginary axis than along the real one - the orientation of the Liouvillian cloud. """
function EllipticCloud(n, aspect, rng)
    τ = (aspect - 1) / (aspect + 1)
    H1 = randn(rng, ComplexF64, n, n); H1 = (H1 + H1') / 2
    H2 = randn(rng, ComplexF64, n, n); H2 = (H2 + H2') / 2
    J = (sqrt((1 + τ) / 2) .* H1 .+ im * sqrt((1 - τ) / 2) .* H2) ./ sqrt(n)
    return im .* eigvals(J)                 # rotate: long axis along Im
end

function Edge(Ns; realizations = 12, dimension = 1500)
    # the clouds depend on neither N nor the window - only the cut does - so they are drawn once
    aspects = (1.0, 3.0, 7.0, 15.0)
    clouds = Dict(aspect => [EllipticCloud(dimension, aspect, Xoshiro(hash((aspect, k)))) for k = 1:realizations]
                  for aspect in aspects)
    rows = []
    for N in Ns
        file = SpectrumPath("reference", N, -20.0, 3.0)
        isfile(file) || (println("missing ", file); continue)
        values = ReadValues(file)
        span = extrema(real, values)
        extent = extrema(imag, values)
        @printf("Liouvillian N = %d: Re from %.1f to %.1f, Im from %.1f to %.1f - aspect %.1f\n", N,
                span..., extent..., (extent[2] - extent[1]) / (span[2] - span[1]))

        for window in ("w8", "w12", "w16")
            spacings, ratios = Selected(values, CUTS[window])
            length(ratios) > 30 || continue
            s = Statistic(spacings, ratios)
            fraction = length(ratios) / length(values)
            @printf("  %-4s Liouvillian n = %4d (%.3f of the sector): <|z|> %.4f  -<cos> %.4f  var %.4f\n",
                    window, s.n, fraction, s.r, s.cos, s.var)
            for aspect in aspects
                draws = []
                for cloud in clouds[aspect]
                    cut = sort(real.(cloud), rev = true)[max(1, round(Int, fraction * dimension))]
                    push!(draws, Statistic(Selected(cloud, cut)...))
                end
                m = (mean(d.r for d in draws), mean(d.cos for d in draws), mean(d.var for d in draws))
                e = (std(d.r for d in draws), std(d.cos for d in draws), std(d.var for d in draws))
                @printf("       elliptic cloud, aspect %4.1f, half-plane cut: <|z|> %.4f ± %.4f  -<cos> %.4f ± %.4f  var %.4f ± %.4f\n",
                        aspect, m[1], e[1], m[2], e[2], m[3], e[3])
                push!(rows, (N, window, s.n, s.r, s.cos, s.var, aspect, m..., e...))
            end
        end
    end
    WriteRows("edge.txt", ["N", "window", "n", "r", "cos", "var", "aspect", "r_cloud", "cos_cloud",
                           "var_cloud", "r_cloud_spread", "cos_cloud_spread", "var_cloud_spread"], rows)
end


# ---------------------------------------------------------------------------------------------
# Moving window along the real axis
# ---------------------------------------------------------------------------------------------

""" Statistics of the eigenvalues in a BAND  upper - width <= Re λ <= upper  that slides along the
    real axis, upper = 0, -step, -2 step, ... to the far end of the sector.  The half-plane cuts
    Re λ >= -w of WindowScan always contain the slowest modes together with everything above the
    cut, so they cannot say WHERE along the real axis the statistics change; a band can.

    The neighbours are found in the full sector, so the band puts no artificial edge into the point
    set (SelectReferences).  What it cannot remove is the cloud's own edge: the first bands ARE the
    leading edge, where every ensemble departs from its bulk values.  The control columns measure
    that - elliptic Ginibre clouds of the sector's aspect ratio (Im extent / Re extent), cut into the
    band holding the SAME QUANTILES of Re λ, ratios pooled over `realizations` clouds.  The
    Liouvillian density is skewed along Re λ, with a long thin tail towards the steady state, so in
    the slow tail the control is a cap of a sharp-edged cloud: it is a yardstick for the edge effect,
    not a prediction of the Liouvillian's chaotic value there.  For the self-conjugate sector m0 the
    reference class is AI-dagger, not the Ginibre class of the control.

    The density varies by two orders of magnitude along the axis and so does the count per band;
    bands with fewer than `fewest` eigenvalues get NaN statistics.  The error bars are a bootstrap
    over the eigenvalues of the band, which ignores the correlation of neighbouring ratios and is
    therefore a lower bound. """
const CONTROL_CLOUDS = Dict{Tuple{Float64, Int, Int}, Any}()     # (aspect, realizations, dimension)

function MovingWindow(Ns; g = -20.0, η = 5.0, κ = 0.3, task = nothing, sector = "m1",
        width = 2.0, step = width / 2, fewest = 30, realizations = 8, dimension = 3000)

    rows = []

    for N in Ns
        file = FindSpectrum(task, N, g, η; κ = κ, label = sector)
        if isnothing(file)
            @printf("\nN = %d, g = %g, η = %g, κ = %g: no saved spectrum - compute it with\n    julia BHNumberConservingLiouvillianMap.jl points --g %g --eta %g --kappa %g --N %d\n",
                    N, g, η, κ, g, η, κ, N)
            continue
        end
        e = PerEigenvalue(ReadValues(file))
        re = real.(e.values)
        total = length(re)
        aspect = round((maximum(imag, e.values) - minimum(imag, e.values)) / (maximum(re) - minimum(re)),
                       digits = 1)

        control = realizations == 0 ? [] : get!(CONTROL_CLOUDS, (aspect, realizations, dimension)) do
            map(1:realizations) do k
                cloud = PerEigenvalue(EllipticCloud(dimension, aspect, Xoshiro(hash((aspect, k)))))
                (cloud = cloud, order = sortperm(real.(cloud.values), rev = true))
            end
        end

        @printf("\nN = %d, g = %g, η = %g, κ = %g, %s: %d eigenvalues, Re λ from %.2f to %.2f, aspect %.1f; band width %g, step %g\n    (%s)\n",
                N, g, η, κ, sector, total, minimum(re), maximum(re), aspect, width, step, file)
        println("    Re λ band             n    <|z|>          -<cos>          var(s)          P(s<1/2) | control <|z|>  -<cos>   var")

        for upper in 0.0:(-step):minimum(re)
            lower = upper - width
            index = findall(x -> lower <= x <= upper, re)
            n = length(index)
            row = Any[N, upper, lower, upper - width / 2, n]

            if n < fewest
                push!(rows, vcat(row, fill(NaN, 11)))
                continue
            end

            spacings, ratios = e.spacings[index], e.ratios[index]
            s = Statistic(spacings, ratios)
            b = Bootstrap(spacings, ratios; samples = 400)
            small = count(<(0.5), spacings ./ mean(spacings)) / n

            # the same quantiles of Re λ in every control cloud
            from, to = count(>(upper), re) / total, count(>=(lower), re) / total
            controlSpacings, controlRatios = Float64[], ComplexF64[]
            for c in control
                m = length(c.order)
                pick = c.order[(round(Int, from * m) + 1):round(Int, to * m)]
                isempty(pick) && continue
                append!(controlSpacings, c.cloud.spacings[pick] ./ mean(c.cloud.spacings[pick]))
                append!(controlRatios, c.cloud.ratios[pick])
            end
            ctrl = length(controlRatios) >= fewest ? Statistic(controlSpacings, controlRatios) :
                                                     (r = NaN, cos = NaN, var = NaN, n = length(controlRatios))

            append!(row, [s.r, b.r, s.cos, b.cos, s.var, b.var, small, ctrl.n, ctrl.r, ctrl.cos, ctrl.var])
            push!(rows, row)
            @printf("  [%6.1f, %6.1f]  %5d   %.4f ± %.4f  %.4f ± %.4f  %.4f ± %.4f  %.3f    |  %.4f  %.4f  %.4f\n",
                    lower, upper, n, s.r, b.r, s.cos, b.cos, s.var, b.var, small, ctrl.r, ctrl.cos, ctrl.var)
        end
    end

    @printf("  Ginibre (A): <|z|> %.4f  -<cos> %.4f  var(s) %.4f;  2D Poisson: %.4f  %.4f  %.4f;  AI-dagger: 0.722  0.193\n",
            REFERENCE.ginibreRadius, REFERENCE.ginibreCosine, REFERENCE.ginibreVariance,
            REFERENCE.poissonRadius, REFERENCE.poissonCosine, REFERENCE.poissonVariance)

    isempty(rows) && return nothing
    WriteRows(@sprintf("moving_g%g_e%g_k%g_%s_W%g_step%g.txt", g, η, κ, sector, width, step),
              ["N", "upper", "lower", "centre", "n", "r", "r_err", "cos", "cos_err", "var", "var_err",
               "small", "n_control", "r_control", "cos_control", "var_control"], rows)
end

""" `moving` on every combination of the comma-separated --g and --eta lists, one table each. """
function MovingCommand(Ns, options)
    Option(name, default) = haskey(options, name) ? parse(Float64, options[name]) : default
    List(name, default) = haskey(options, name) ? parse.(Float64, split(options[name], ',')) : [default]
    width = Option("width", 2.0)
    for g in List("g", -20.0), η in List("eta", 5.0)
        MovingWindow(Ns; g = g, η = η, κ = Option("kappa", 0.3),
                     task = get(options, "task", nothing), sector = get(options, "sector", "m1"),
                     width = width, step = Option("step", width / 2),
                     fewest = round(Int, Option("minimum", 30.0)),
                     realizations = round(Int, Option("realizations", 8.0)))
    end
end


# ---------------------------------------------------------------------------------------------
# C10 dissipative spectral form factor
# ---------------------------------------------------------------------------------------------

""" A spectrum prepared for the form factor: rescaled to unit mean local spacing in the bulk (the
    median local scale of the densest half), centred, and given a Gaussian WINDOW
    w_n = exp(-|z_n|^2 / (2 σ^2)) with σ a quarter of the rms extent of the bulk.  The window
    replaces a hard selection of the bulk, whose sharp artificial edge would ring in the form
    factor at τ ~ 1/size and swamp the correlation hole; the price is a smaller effective number
    of points, Σ w_n^2. """
function WindowedSpectrum(values; neighbours = 30, keep = 0.5, width = 0.25)
    values = filter(z -> abs(z) > 1e-8, values)
    distances, _, _ = NeighbourData(values, neighbours)
    scale = distances[:, neighbours] ./ sqrt(neighbours / pi)
    bulk = sortperm(scale)[1:round(Int, keep * length(values))]
    centre = mean(values[bulk])
    points = (values .- centre) ./ median(scale[bulk])
    σ = width * sqrt(mean(abs2, points[bulk]))
    return (points = points, weights = exp.(-abs2.(points) ./ (2 * σ^2)))
end

""" Connected DSFF of Li, Prosen and Chan with a window:
    K(τ) = [ < |Σ_n w_n exp(i τ ê·z_n)|^2 > - |< Σ_n w_n exp(i τ ê·z_n) >|^2 ] / < Σ_n w_n^2 >,
    averaged over the ensemble and over directions ê of the complex plane.  It tends to 1 at large
    τ (the plateau) for every ensemble; the approach to it - the correlation hole and the ramp -
    is what distinguishes Ginibre from Poisson. """
function FormFactor(spectra, τs; directions = 32)
    angles = range(0, pi, length = directions + 1)[1:end-1]
    norms = mean(sum(abs2, s.weights) for s in spectra)
    K = zeros(length(τs))
    Threads.@threads for i in eachindex(τs)
        τ = τs[i]
        total = 0.0
        for φ in angles
            e = cis(φ)
            sums = [sum(w * cis(τ * real(conj(e) * z)) for (z, w) in zip(s.points, s.weights)) for s in spectra]
            total += mean(abs2, sums) - abs2(mean(sums))
        end
        K[i] = total / directions / norms
    end
    return K
end

function DSFF()
    directory = joinpath(ANALYSIS_RESULTS, "spectra", "dsff")
    files = filter(f -> endswith(f, "_m1.txt"), readdir(directory; join = true))
    isempty(files) && error("no spectra in $directory - run BHNumberConservingLiouvillianMap.jl dsff first")
    liouvillian = [WindowedSpectrum(ReadValues(f)) for f in files]
    n = length(ReadValues(files[1]))

    ginibre = [WindowedSpectrum(Ginibre(n, Xoshiro(k))) for k in eachindex(files)]
    poisson = map(eachindex(files)) do k
        rng = Xoshiro(100 + k)
        WindowedSpectrum([complex(2 * rand(rng) - 1, 2 * rand(rng) - 1) for _ = 1:n])
    end

    τs = exp.(range(log(0.05), log(20.0), length = 120))
    rows = zip(τs, FormFactor(liouvillian, τs), FormFactor(ginibre, τs), FormFactor(poisson, τs))
    @printf("%d Liouvillian spectra of %d eigenvalues; effective points per spectrum (Σ w²) %.0f
",
            length(files), n, mean(sum(abs2, s.weights) for s in liouvillian))
    WriteRows("dsff.txt", ["tau", "K_liouvillian", "K_ginibre", "K_poisson"], collect(rows))
end


# ---------------------------------------------------------------------------------------------
# A1 reference values of the symmetry classes
# ---------------------------------------------------------------------------------------------

function Classes(; dimension = 2000, realizations = 4)
    ensembles = (("GinUE (A)", rng -> randn(rng, ComplexF64, dimension, dimension)),
                 ("real Ginibre (AI, K-type)", rng -> complex.(randn(rng, dimension, dimension))),
                 ("complex symmetric (AI†)", rng -> (X = randn(rng, ComplexF64, dimension, dimension); X + transpose(X))))
    rows = []
    for (name, draw) in ensembles
        for (selection, offAxis) in (("bulk", 0.0), ("bulk, |Im z| > 0.1 R", 0.1))
            name == "GinUE (A)" && offAxis > 0 && continue
            draws = []
            for k = 1:realizations
                values = eigvals(draw(Xoshiro(k))) ./ sqrt(dimension)
                radius = maximum(abs, values)
                spacings, ratios = Selected(values, nothing)
                values2 = values
                if offAxis > 0
                    # keep the reference points away from the real axis (neighbours from all)
                    filtered = filter(z -> abs(z) > 1e-8, values2)
                    distances, nearest, next = NeighbourData(filtered, 30)
                    scale = distances[:, 30] ./ sqrt(30 / pi)
                    allRatios = [(filtered[nearest[i]] - filtered[i]) / (filtered[next[i]] - filtered[i]) for i in eachindex(filtered)]
                    bulk = sortperm(scale)[1:round(Int, 0.5 * length(filtered))]
                    index = [i for i in bulk if abs(imag(filtered[i])) > offAxis * radius]
                    spacings = (distances[:, 1] ./ scale)[index]
                    ratios = allRatios[index]
                end
                push!(draws, Statistic(spacings, ratios))
            end
            m = (mean(d.r for d in draws), mean(d.cos for d in draws), mean(d.var for d in draws))
            e = (std(d.r for d in draws), std(d.cos for d in draws), std(d.var for d in draws))
            @printf("  %-28s %-22s <|z|> %.4f ± %.4f   -<cos> %.4f ± %.4f   var(s) %.4f ± %.4f\n",
                    name, selection, m[1], e[1], m[2], e[2], m[3], e[3])
            push!(rows, (name, selection, m..., e...))
        end
    end
    println("  published: GinUE 0.7378 / 0.2405; AI† 0.722 / 0.193 (Sá, Ribeiro, Prosen; Hamazaki et al.)")
    WriteRows("classes.txt", ["ensemble", "selection", "r", "cos", "var", "r_spread", "cos_spread", "var_spread"], rows)
end


if abspath(PROGRAM_FILE) == @__FILE__
    command = isempty(ARGS) ? "" : ARGS[1]
    options = Dict{String, String}()
    for i = 2:2:(length(ARGS) - 1)
        options[replace(ARGS[i], "--" => "")] = ARGS[i + 1]
    end
    Ns = haskey(options, "N") ? parse.(Int, split(options["N"], ',')) : [12, 14, 16]
    BLAS.set_num_threads(4)

    command == "tables" ? Tables(Ns) :
    command == "unfolding" ? Unfolding(Ns) :
    command == "edge" ? Edge(Ns) :
    command == "dsff" ? DSFF() :
    command == "classes" ? Classes() :
    command == "moving" ? MovingCommand(Ns, options) :
    println("usage: julia -t 8 BHNumberConservingLiouvillianAnalysis.jl tables | unfolding | edge | dsff | classes | moving [--N 12,14,16]")
end
