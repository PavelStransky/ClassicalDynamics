# BHNumberConservingAttractors.jl
#
# The classical calculations of number-conserving-BH-paper-TODO.md that are not Lyapunov maps (those
# run on Chimera, BHMapNumberConservingChimeraSubmit.sh) and not transients
# (BHNumberConservingTransient.jl).  All of them are cheap - minutes to an hour on the laptop.
#
#   julia -t 20 BHNumberConservingAttractors.jl <command>
#
#   bifurcation   C15  η = 3, g = -5.0 ... -10.0 in steps of 0.05: local maxima of the bond coherence
#                      C and of n_1 on the attractor, from 8 random initial conditions and from two
#                      continuation sweeps (g decreasing and increasing, each step started from the
#                      final state of the previous one, kicked by 1e-6 so that an unstable fixed
#                      point is left) - hysteresis shows up as the two sweeps disagreeing.
#                      The route to chaos (A5) is read off this diagram together with
#                      the Lyapunov cut `cut-eta3` computed on Chimera.
#   correlations  C14  autocorrelation functions on the chaotic attractor (-20, 3) of the q = 2π/3
#                      Fourier component of n_j (complex) and of C, averaged over trajectories;
#                      decay rate and frequency fitted, for comparison with the quantum gap
#                      G_∞ = 0.520 and with Im Λ of the slowest pair (sector q = 2π/3 for ñ_q, q = 0
#                      for C).
#   basins        K11  basin fractions at the multistable point (-20, 1) from 2000 initial
#                      conditions, with binomial errors, after relaxations of 2000 AND 8000 (a basin
#                      that shrinks is a long chaotic transient, not an attractor), and the attractor
#                      fingerprints (basins_g-20.00_e1.00.txt, read by `uncertainty`).
#   uncertainty   C16  uncertainty exponent of the basin boundaries at (-20, 1): the fraction f(ε)
#                      of initial conditions whose attractor changes under a perturbation of
#                      Fubini-Study size ε, f ∝ ε^α; the boundary dimension is D - α with D = 4.
#   simplex       C17  data for the attractor figure: the limit cycle(s) at g = -6, the strange
#                      attractor at g = -20 with its visit density on the population simplex, the
#                      three attractors at η = 1, and a Poincaré section at g = -20
#                      (plot with analyse_attractors_number_conserving.py).
#   robustness    K9   the g = -6 limit cycle with longer windows (1000, 20000) and its exact Floquet
#                      multipliers; the reduced contraction rate of the regular reference points
#                      (the FP entry of Table III);
#                 K10  the L = 4 hyperchaos with three windows and three integrator tolerances
#                      (`robustness <g> <η>` sets the L = 4 point; default (-20, 3)).
#
# OUTPUT in $BH_RESULTS_DIR (default ~/results/bh/number-conserving/classical/3/attractors).

include(joinpath(@__DIR__, "BHNumberConservingTransient.jl"))     # also brings BHNumberConserving.jl

using ProgressMeter

const ATTRACTOR_RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "classical", "3", "attractors"))

Ring(g, η; L = 3, κ = 0.3) = NumberConservingParameters(L; J = 1.0, g = g, κ = κ, η = η)

BondCoherence(ψ) = sum(real(conj(ψ[j]) * ψ[mod1(j + 1, length(ψ))]) for j in eachindex(ψ))
Current(ψ) = sum(imag(conj(ψ[j]) * ψ[mod1(j + 1, length(ψ))]) for j in eachindex(ψ))
FourierDensity(ψ) = sum(abs2(ψ[j]) * exp(-2im * pi * (j - 1) / length(ψ)) for j in eachindex(ψ))

""" Samples of the state every `step` over (start, stop], starting from ψ0 at t = 0. """
function Sample(ψ0, parameters, start, stop, step)
    solution = solve(PlainProblem(ψ0, parameters, stop), DP8(); reltol = 1e-10, abstol = 1e-10,
                     saveat = start:step:stop, maxiters = 10^8)
    return solution.u
end

function WriteTable(file, header, rows)
    mkpath(dirname(file))
    open(file, "w") do io
        println(io, "# ", join(header, "\t"))
        for row in rows
            println(io, join([x isa Integer ? string(x) : @sprintf("%.10g", x) for x in row], "\t"))
        end
    end
    println("wrote ", file)
end


# ---------------------------------------------------------------------------------------------
# C15  bifurcation diagram
# ---------------------------------------------------------------------------------------------

""" Local maxima of C and of n_1 along the attractor reached from ψ0: after `relaxation`, over a
    window of length `window`.  Maxima are located by the downward zero crossings of the time
    derivatives, dC/dt and dn_1/dt, computed from the vector field - exact, no sampling. """
function LocalMaxima(ψ0, parameters; relaxation = 2000.0, window = 600.0, maximum_ = 300)
    L = parameters.L
    work = NumberConservingWorkspace(L)
    field = zeros(ComplexF64, L)
    function Derivatives(u)
        UpdateWorkspace!(work, u, parameters)
        VectorField!(field, u, parameters, work)
        dC = sum(real(conj(field[j]) * u[mod1(j + 1, L)] + conj(u[j]) * field[mod1(j + 1, L)]) for j = 1:L)
        dn = 2 * real(conj(u[1]) * field[1])
        return dC, dn
    end

    maximaC = Float64[]
    maximaN = Float64[]
    active(t) = t > relaxation
    # a maximum is a DOWNWARD zero crossing of the derivative: affect! (upward) does nothing
    function Record!(integrator, index)
        active(integrator.t) || return
        u = integrator.u
        index == 1 && length(maximaC) < maximum_ && push!(maximaC, BondCoherence(u) / sum(abs2, u))
        index == 2 && length(maximaN) < maximum_ && push!(maximaN, abs2(u[1]) / sum(abs2, u))
    end
    # On a fixed point both derivatives are rounding noise (~1e-17) of either sign, so that every
    # step "crosses zero" right at its start; the two callbacks then retrigger each other and the
    # time stops advancing until maxiters.  The offset keeps the sign definite there.  A true
    # maximum is displaced by it by ~1e-12 in time, far below the integration tolerance.
    offset = 1e-12
    callback = CallbackSet(
        ContinuousCallback((u, t, integrator) -> Derivatives(u)[1] + offset, nothing,
                           integrator -> Record!(integrator, 1); save_positions = (false, false)),
        ContinuousCallback((u, t, integrator) -> Derivatives(u)[2] + offset, nothing,
                           integrator -> Record!(integrator, 2); save_positions = (false, false)))

    solution = solve(PlainProblem(ψ0, parameters, relaxation + window), DP8(); reltol = 1e-10,
                     abstol = 1e-10, callback = callback, save_everystep = false, maxiters = 10^8)
    final = solution.u[end]

    # a fixed point has no maxima: its constant values stand for them
    isempty(maximaC) && push!(maximaC, BondCoherence(final))
    isempty(maximaN) && push!(maximaN, abs2(final[1]))
    return maximaC, maximaN, final
end

function Bifurcation(; η = 3.0, gs = collect(-2.0:-0.01:-10.0), randoms = 100, kick = 1e-6)
    rows = []
    lock_ = ReentrantLock()

    randomProgress = Progress(length(gs) * randoms; desc = "random initial conditions ")
    Threads.@threads for g in gs
        parameters = Ring(g, η)
        for r = 1:randoms
            rng = Xoshiro(hash((g, η, r, "bifurcation")))
            maximaC, maximaN, _ = LocalMaxima(RandomInitialCondition(3, rng), parameters)
            lock(lock_) do
                append!(rows, [(g, r, 1, m) for m in maximaC])
                append!(rows, [(g, r, 2, m) for m in maximaN])
            end
            next!(randomProgress)
        end
    end

    # continuation, both directions; sources -1 (g decreasing) and -2 (g increasing)
    sweepProgress = Progress(2 * length(gs); desc = "continuation sweeps ")
    for (source, sweep) in ((-1, gs), (-2, reverse(gs)))
        rng = Xoshiro(source)
        ψ = UniformInitialCondition(3; amplitude = 1e-3, rng = rng)
        for g in sweep
            # The carried state is kicked at every step.  On the stable side it converges onto the
            # uniform state EXACTLY (identical amplitudes), and the vector field keeps identical
            # amplitudes identical: with nothing to seed the instability, the sweep would follow
            # the unstable fixed point beyond the threshold for ever.
            ψ = ψ .+ kick .* randn(rng, ComplexF64, 3)
            maximaC, maximaN, ψ = LocalMaxima(ψ, Ring(g, η); relaxation = 1000.0)
            append!(rows, [(g, source, 1, m) for m in maximaC])
            append!(rows, [(g, source, 2, m) for m in maximaN])
            next!(sweepProgress)
        end
    end

    sort!(rows)
    WriteTable(joinpath(ATTRACTOR_RESULTS, @sprintf("bifurcation_e%.2f.txt", η)),
               ["g", "source", "observable", "maximum"], rows)
    println("  source r > 0: random initial condition r; -1: continuation towards more negative g; ",
            "-2: towards less negative g.  observable 1 = bond coherence C, 2 = n_1")
end


# ---------------------------------------------------------------------------------------------
# C14  correlation decay
# ---------------------------------------------------------------------------------------------

""" Normalised autocorrelation <δA(t + s) conj(δA(t))> / <|δA|^2> of a (complex) time series for
    lags 0 ... maxLag, by direct summation (no FFT package needed at these lengths). """
function Autocorrelation(series, maxLag)
    x = series .- mean(series)
    n = length(x)
    c = zeros(ComplexF64, maxLag + 1)
    Threads.@threads for lag = 0:maxLag
        total = zero(ComplexF64)
        @inbounds for t = 1:(n - lag)
            total += x[t + lag] * conj(x[t])
        end
        c[lag + 1] = total / (n - lag)
    end
    return c ./ real(c[1])
end

""" Decay rate from a straight-line fit of log|C(s)| over the lags where 0.02 < |C| < 0.7, and
    frequency from the slope of the unwrapped phase (complex series) or, for a real series, from
    the spacing of the zero crossings. """
function FitDecay(lags, c)
    selection = findall(i -> 0.02 < abs(c[i]) < 0.7 && lags[i] > 0, eachindex(c))
    # stop at the first time the envelope reaches the noise floor
    last = findfirst(i -> abs(c[i]) < 0.02, eachindex(c))
    isnothing(last) || (selection = filter(i -> i < last, selection))
    length(selection) < 3 && return (rate = NaN, frequency = NaN)

    X = hcat(ones(length(selection)), lags[selection])
    slope = (X \ log.(abs.(c[selection])))[2]

    frequency = NaN
    if any(z -> abs(imag(z)) > 1e-3 * abs(z), c)
        phase = unwrap_(angle.(c[1:(isnothing(last) ? length(c) : last)]))
        k = 1:length(phase)
        frequency = -(hcat(ones(length(k)), lags[k]) \ phase)[2]
    else
        crossings = [lags[i] for i = 2:length(c) if sign(real(c[i])) != sign(real(c[i - 1]))]
        length(crossings) > 2 && (frequency = pi / mean(diff(crossings)))
    end

    return (rate = -slope, frequency = frequency)
end

function unwrap_(φ)
    out = copy(φ)
    for i = 2:length(out)
        δ = out[i] - out[i - 1]
        out[i] -= 2pi * round(δ / (2pi))
    end
    return out
end

function Correlations(; g = -20.0, η = 3.0, trajectories = 8, relaxation = 1000.0, length_ = 20000.0,
        step = 0.05, maxLag = 50.0)
    parameters = Ring(g, η)
    lags = collect(0:step:maxLag)
    density = zeros(ComplexF64, length(lags))
    coherence = zeros(ComplexF64, length(lags))

    for r = 1:trajectories
        rng = Xoshiro(hash((g, η, r, "correlations")))
        states = Sample(RandomInitialCondition(3, rng), parameters, relaxation, relaxation + length_, step)
        density .+= Autocorrelation(FourierDensity.(states), length(lags) - 1) ./ trajectories
        coherence .+= Autocorrelation(complex.(BondCoherence.(states)), length(lags) - 1) ./ trajectories
        @printf("  trajectory %d of %d\n", r, trajectories)
    end

    WriteTable(joinpath(ATTRACTOR_RESULTS, @sprintf("correlations_g%.2f_e%.2f.txt", g, η)),
               ["lag", "Re_C_nq", "Im_C_nq", "C_coherence"],
               [(lags[i], real(density[i]), imag(density[i]), real(coherence[i])) for i in eachindex(lags)])

    a = FitDecay(lags, density)
    b = FitDecay(lags, coherence)
    @printf("\n  ñ_q (q = 2π/3):  decay rate %.4f, frequency %.4f   -> compare Λ of the slowest pair in sector q = 2π/3\n",
            a.rate, a.frequency)
    @printf("  C (q = 0):       decay rate %.4f, frequency %.4f   -> compare the slowest pair in sector q = 0\n",
            b.rate, b.frequency)
    println("  quantum gap extrapolated in the draft: G_∞ = 0.520 at the chaotic point")
end


# ---------------------------------------------------------------------------------------------
# K11 basins and C16 uncertainty exponent
# ---------------------------------------------------------------------------------------------

""" Gauge- and translation-invariant fingerprint of the attractor reached from ψ0: time averages and
    standard deviations of (C, max_j n_j, current) over a window after the relaxation, the reduced
    finite-time Lyapunov exponent over the same window, and from it the TYPE - 0 fixed point
    (λ < -0.005: only a fixed point contracts in every reduced direction), 1 regular (limit cycle or
    torus, λ ≈ 0), 2 chaotic (λ > 0.01).  The type is a
    hard criterion when grouping, as the classification is in CountAttractors of
    BHNumberConserving.jl: without it the wide time fluctuation of a chaotic attractor lets it
    swallow a limit cycle lying inside its range of averages.  `final` is the last state, so a
    later, longer relaxation can continue from it. """
function Fingerprint(ψ0, parameters; relaxation = 2000.0, window = 500.0, step = 0.5)
    L = parameters.L
    ψ = relaxation > 0 ? Evolve(ψ0, parameters, relaxation) : ψ0 ./ norm(ψ0)

    tracker = TransientTracker(parameters, NumberConservingWorkspace(L), Any[], 0.0, Inf,
                               fill(0.0, 3), fill(Inf, 3), NaN, 0, 0.0)
    δ = randn(Xoshiro(hash(ψ)), ComplexF64, L)
    ReducedDeviation!(δ, ψ)
    u0 = vcat(ψ, δ ./ norm(δ))
    solution = solve(ODEProblem(TransientField!, u0, (0.0, window), tracker), DP8();
                     reltol = 1e-10, abstol = 1e-10, saveat = 0.0:step:window, maxiters = 10^8,
                     callback = PeriodicCallback(Renormalise!, 1.0; save_positions = (false, false)))
    final = solution.u[end]
    λ = (tracker.logGrowth + log(ReducedDeviation!(copy(final[(L + 1):end]), final[1:L]))) / window

    values = [(BondCoherence(u[1:L]), maximum(abs2, u[1:L]), Current(u[1:L])) for u in solution.u]
    means = [mean(v[k] for v in values) for k = 1:3]
    spreads = [std(v[k] for v in values) for k = 1:3]
    type = λ > 0.01 ? 2 : λ < -0.005 ? 0 : 1

    return (means = means, spreads = spreads, λ = λ, type = type, final = final[1:L])
end

const TYPE_NAMES = Dict(0 => "fixed point", 1 => "regular", 2 => "chaotic")

SameAttractor(a, b; tolerance = 0.02) =
    a.type == b.type && all(abs.(a.means .- b.means) .<= tolerance .+ 0.5 .* (a.spreads .+ b.spreads))

""" Greedy grouping of fingerprints: same type, and every mean within `tolerance` plus half the
    time fluctuation.  The count is a lower bound, as in CountAttractors. """
function GroupFingerprints(fingerprints; tolerance = 0.02)
    centroids = []
    labels = zeros(Int, length(fingerprints))
    for (i, f) in enumerate(fingerprints)
        k = findfirst(c -> SameAttractor(f, c; tolerance = tolerance), centroids)
        if isnothing(k)
            push!(centroids, f)
            k = length(centroids)
        end
        labels[i] = k
    end
    return labels, centroids
end

""" Nearest attractor of the same type, in units of the allowed deviation; 0 when none is within
    reach (a new or unresolved attractor - counted separately, never forced into a class). """
function Classify(f, centroids; tolerance = 0.02, reach = 3.0)
    best, index = Inf, 0
    for (k, c) in enumerate(centroids)
        c.type == f.type || continue
        d = maximum(abs.(f.means .- c.means) ./ (tolerance .+ 0.5 .* (f.spreads .+ c.spreads)))
        d < best && ((best, index) = (d, k))
    end
    return best <= reach ? index : 0
end

BasinFile(g, η) = joinpath(ATTRACTOR_RESULTS, @sprintf("basins_g%.2f_e%.2f.txt", g, η))

""" Basin fractions with binomial errors (K11), measured TWICE: after a relaxation of 2000 (the
    one of the maps and of the draft's Table) and after 8000.  If a basin shrinks between the two,
    part of what looked like an attractor is a long chaotic transient - which also decides how
    the uncertainty exponent below has to be read. """
function Basins(; g = -20.0, η = 1.0, samples = 2000, relaxations = (2000.0, 8000.0), window = 500.0)
    parameters = Ring(g, η)
    early = Vector{Any}(undef, samples)
    late = Vector{Any}(undef, samples)
    Threads.@threads for i = 1:samples
        early[i] = Fingerprint(RandomInitialCondition(3, Xoshiro(hash((g, η, i, "basins")))), parameters;
                               relaxation = relaxations[1], window = window)
        late[i] = Fingerprint(early[i].final, parameters;
                              relaxation = relaxations[2] - relaxations[1] - window, window = window)
    end

    rows = []
    centroids = nothing
    for (relaxation, fingerprints) in zip(relaxations, (early, late))
        labels, groups = GroupFingerprints(fingerprints)
        relaxation == relaxations[1] && (centroids = groups)

        @printf("(g, η) = (%.1f, %.1f), relaxation %g: %d initial conditions, %d attractors\n",
                g, η, relaxation, samples, length(groups))
        println("   k   type          basin ± binomial   <C>      max n    current    λ       (spreads)")
        for (k, c) in enumerate(groups)
            p = count(==(k), labels) / samples
            error_ = sqrt(p * (1 - p) / samples)
            @printf("  %2d   %-12s  %.4f ± %.4f   %+.4f  %.4f  %+.4f  %+.4f  (%.3f %.3f %.3f)\n", k,
                    TYPE_NAMES[c.type], p, error_, c.means..., c.λ, c.spreads...)
            push!(rows, (relaxation, k, c.type, p, error_, c.means..., c.spreads..., c.λ))
        end
    end
    @printf("  for comparison: 96 initial conditions give a binomial error of up to ±%.3f\n", sqrt(0.25 / 96))

    WriteTable(BasinFile(g, η), ["relaxation", "attractor", "type", "basin", "error", "C", "max_n",
                                 "current", "sC", "smax_n", "scurrent", "lambda"], rows)
    return centroids
end

""" The attractors after the SHORT relaxation, as written by Basins - the reference for Classify. """
function ReadCentroids(g, η)
    file = BasinFile(g, η)
    isfile(file) || return Basins(; g = g, η = η)
    centroids = []
    first = nothing
    for line in eachline(file)
        startswith(line, "#") && continue
        x = parse.(Float64, split(line, '\t'))
        isnothing(first) && (first = x[1])
        x[1] == first || continue
        push!(centroids, (means = x[6:8], spreads = x[9:11], λ = x[12], type = round(Int, x[3])))
    end
    return centroids
end

""" Uncertainty exponent (C16): f(ε) = share of initial conditions whose attractor changes under a
    Fubini-Study perturbation of size ε, fitted as f ∝ ε^α.  Every initial condition is classified
    after the relaxation of the basin table (2000); if Basins shows a basin that shrinks with a
    longer relaxation, rerun with a longer one (keyword `relaxation`) - a transient mistaken for an
    attractor gives α ≈ 0 for a reason that has nothing to do with the boundary geometry. """
function Uncertainty(; g = -20.0, η = 1.0, samples = 2000, exponents = -9.0:0.5:-2.0,
        relaxation = 2000.0)
    parameters = Ring(g, η)
    centroids = ReadCentroids(g, η)

    bases = [RandomInitialCondition(3, Xoshiro(hash((g, η, i, "uncertainty")))) for i = 1:samples]
    baseLabels = zeros(Int, samples)
    Threads.@threads for i = 1:samples
        baseLabels[i] = Classify(Fingerprint(bases[i], parameters; relaxation = relaxation), centroids)
    end
    @printf("  %d of %d reference initial conditions classified\n", count(>(0), baseLabels), samples)

    rows = []
    for e in exponents
        ε = 10.0^e
        uncertain = zeros(Bool, samples)
        valid = zeros(Bool, samples)
        Threads.@threads for i = 1:samples
            rng = Xoshiro(hash((g, η, i, e, "perturbation")))
            label = Classify(Fingerprint(AtDistance(bases[i], ε, rng), parameters; relaxation = relaxation),
                             centroids)
            valid[i] = label > 0 && baseLabels[i] > 0
            uncertain[i] = valid[i] && label != baseLabels[i]
        end
        n = count(valid)
        f = count(uncertain) / max(n, 1)
        push!(rows, (ε, f, sqrt(f * (1 - f) / max(n, 1)), n))
        @printf("  ε = %.1e: f = %.4f ± %.4f (%d valid pairs)\n", ε, rows[end][2], rows[end][3], n)
    end

    fit = [r for r in rows if r[2] > 0]
    if length(fit) >= 3
        x = log10.([r[1] for r in fit])
        y = log10.([r[2] for r in fit])
        X = hcat(ones(length(x)), x)
        coefficients = X \ y
        residual = y .- X * coefficients
        σ = sqrt(sum(abs2, residual) / (length(x) - 2) * inv(X' * X)[2, 2])
        @printf("\n  uncertainty exponent α = %.3f ± %.3f; boundary dimension D - α = %.3f (D = 4, the reduced space)\n",
                coefficients[2], σ, 4 - coefficients[2])
        println("  α close to 1 means smooth basin boundaries, α well below 1 fractal ones")
    end

    WriteTable(joinpath(ATTRACTOR_RESULTS, @sprintf("uncertainty_g%.2f_e%.2f.txt", g, η)),
               ["epsilon", "f", "error", "pairs"], rows)
end


# ---------------------------------------------------------------------------------------------
# C17  attractors on the simplex
# ---------------------------------------------------------------------------------------------

function SaveTrajectory(file, states; header = ["t", "n1", "n2", "n3", "C", "current"], step = 1.0)
    rows = [((i - 1) * step, abs2.(ψ)..., BondCoherence(ψ), Current(ψ)) for (i, ψ) in enumerate(states)]
    WriteTable(file, header, rows)
end

function Simplex()
    directory = joinpath(ATTRACTOR_RESULTS, "simplex")

    # every distinct attractor at the limit-cycle point and at the multistable point, one trajectory each
    for (g, η, span) in ((-6.0, 3.0, 30.0), (-20.0, 1.0, 100.0))
        parameters = Ring(g, η)
        fingerprints = Vector{Any}(undef, 48)
        finals = Vector{Any}(undef, 48)
        Threads.@threads for i = 1:48
            ψ = RandomInitialCondition(3, Xoshiro(hash((g, η, i, "simplex"))))
            finals[i] = Evolve(ψ, parameters, 2000.0)
            fingerprints[i] = Fingerprint(finals[i], parameters; relaxation = 0.0, window = 200.0)
        end
        labels, centroids = GroupFingerprints(fingerprints)
        for k in eachindex(centroids)
            i = findfirst(==(k), labels)
            states = Sample(finals[i], parameters, 0.0, span, 0.01)
            SaveTrajectory(joinpath(directory, @sprintf("attractor_g%.2f_e%.2f_%d.txt", g, η, k)), states;
                           step = 0.01)
        end
        @printf("(g, η) = (%.1f, %.1f): %d attractors saved\n", g, η, length(centroids))
    end

    # strange attractor at the chaotic point: visit density of the populations
    parameters = Ring(-20.0, 3.0)
    bins = 200
    histogram = zeros(Int, bins, bins)
    states = Sample(RandomInitialCondition(3, Xoshiro(3)), parameters, 1000.0, 21000.0, 0.02)
    for ψ in states
        n = abs2.(ψ) ./ sum(abs2, ψ)
        x = n[2] + 0.5 * n[3]                     # ternary coordinates of (n1, n2, n3)
        y = sqrt(3) / 2 * n[3]
        i = clamp(floor(Int, x * bins) + 1, 1, bins)
        j = clamp(floor(Int, y / (sqrt(3) / 2) * bins) + 1, 1, bins)
        histogram[i, j] += 1
    end
    WriteTable(joinpath(directory, "density_g-20.00_e3.00.txt"), ["i", "j", "count"],
               [(i, j, histogram[i, j]) for i = 1:bins, j = 1:bins if histogram[i, j] > 0])
    SaveTrajectory(joinpath(directory, "trajectory_g-20.00_e3.00.txt"), states[1:5:min(end, 50000)];
                   step = 0.1)

    # Poincaré section at the chaotic point: upward crossings of n_1 = n_2
    crossings = SectionCrossings(states[end], parameters, 1, 50000.0)
    WriteTable(joinpath(directory, "section_g-20.00_e3.00.txt"), ["n3", "C", "current"],
               [(abs2(ψ[3]), BondCoherence(ψ), Current(ψ)) for ψ in crossings])
end


# ---------------------------------------------------------------------------------------------
# K9, K10  robustness of the classification
# ---------------------------------------------------------------------------------------------

function PrintSpectrum(label, r)
    @printf("  %-34s reduced %s  ± %s   class %s\n", label,
            join([@sprintf("%+.5f", x) for x in r.reduced], " "),
            @sprintf("%.1e", maximum(r.uncertainty)), r.classification)
end

""" K9 and K10.  `hyperchaos` is the L = 4 point of K10; (-20, 3) is the default, and the draft's
    hyperchaotic point should be passed if it is another one (command line: robustness <g> <η>).
    At (-20, 3) the test runs found ONE positive exponent (+0.233) and the flow zero, stable
    across windows and tolerances - chaotic, not hyperchaotic. """
function Robustness(; hyperchaos = (-20.0, 3.0))
    println("K9: the limit cycle at (g, η) = (-6, 3)")
    parameters = Ring(-6.0, 3.0)
    for (relaxation, stop) in ((1000.0, 5000.0), (1000.0, 20000.0), (2000.0, 40000.0))
        results = Vector{Any}(undef, 4)
        Threads.@threads for i = 1:4
            rng = Xoshiro(hash((i, "K9")))
            results[i] = LyapunovSpectrum(RandomInitialCondition(3, rng), parameters;
                                          relaxationTime = relaxation, integrationTime = stop,
                                          rng = rng, zeroThreshold = 2e-3)
        end
        for (i, r) in enumerate(results)
            PrintSpectrum(@sprintf("window (%g, %g), IC %d", relaxation, stop, i), r)
        end
    end
    println("  exact relaxation times from the Floquet multipliers of every cycle found:")
    for a in FindAttractors(parameters)
        a.kind === :cycle && a.basin[] > 0 && @printf("    cycle with %d section points: τ = %.4f, i.e. slowest Floquet exponent %.5f\n",
                                                      length(a.points), a.τ, -1 / a.τ)
    end

    println("\nK9: reduced contraction rate at the regular reference points (Table III, FP entry)")
    for (g, η) in ((-20.0, 0.0), (-4.0, 3.0))
        parameters = Ring(g, η)
        for a in FindAttractors(parameters; reflect = η == 0)
            (a.kind === :fixed && a.basin[] > 0) || continue
            r = LyapunovSpectrum(a.state, parameters; relaxationTime = 100.0, integrationTime = 5100.0,
                                 rng = Xoshiro(1))
            @printf("  (g, η) = (%.0f, %.0f): max Re μ of the reduced Jacobian %.5f; reduced λ_max %.5f\n",
                    g, η, -1 / a.τ, r.reduced[1])
        end
    end

    @printf("\nK10: L = 4 at (g, η) = (%.1f, %.1f) - two positive reduced exponents would be hyperchaos\n",
            hyperchaos...)
    parameters = Ring(hyperchaos...; L = 4)
    combinations = [(w, tol, i) for w in ((1000.0, 5000.0), (1000.0, 20000.0), (2000.0, 40000.0))
                    for tol in (1e-8, 1e-10, 1e-12) for i = 1:3]
    results = Vector{Any}(undef, length(combinations))
    Threads.@threads for k in eachindex(combinations)
        (window, tol, i) = combinations[k]
        rng = Xoshiro(hash((i, "K10")))
        results[k] = LyapunovSpectrum(RandomInitialCondition(4, rng), parameters;
                                      relaxationTime = window[1], integrationTime = window[2],
                                      tolerance = tol, rng = rng, zeroThreshold = 2e-3)
    end
    for (k, (window, tol, i)) in enumerate(combinations)
        PrintSpectrum(@sprintf("window (%g, %g), tol %.0e, IC %d", window..., tol, i), results[k])
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    command = isempty(ARGS) ? "" : ARGS[1]
    @printf("Threads: %d\n", Threads.nthreads())
    if command == "bifurcation"
        Bifurcation()
    elseif command == "correlations"
        Correlations()
    elseif command == "basins"
        Basins()
    elseif command == "uncertainty"
        Uncertainty()
    elseif command == "simplex"
        Simplex()
    elseif command == "robustness"
        Robustness(; hyperchaos = length(ARGS) >= 3 ? (parse(Float64, ARGS[2]), parse(Float64, ARGS[3])) :
                                                      (-20.0, 3.0))
    else
        println("usage: julia -t <threads> BHNumberConservingAttractors.jl ",
                "bifurcation | correlations | basins | uncertainty | simplex | robustness")
    end
end
