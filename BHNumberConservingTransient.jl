# BHNumberConservingTransient.jl
#
# Transient chaos in the classical number-conserving dissipative Bose-Hubbard ring - item C13 of
# number-conserving-BH-paper-TODO.md, following Corps and Relano: how long does a trajectory take to
# reach its attractor, and how sensitively does that time depend on where it started?
#
#   t_crit     the first time the trajectory is within TOLERANCE = 1e-4 of the attractor, measured
#              with the FUBINI-STUDY distance  d(ψ, φ) = arccos |<ψ|φ>|  - the distance on CP^(L-1),
#              blind to the global phase that the flow rotates for ever, so a relative fixed point
#              really is a point.  A Euclidean distance would never converge.
#   τ          the linear relaxation time of the attractor from the REDUCED linearisation (the two
#              exact symmetry zeros removed): 1 / |max Re μ| of the Jacobian in the co-rotating
#              frame for a relative fixed point, T / |log |μ|| with μ the largest non-trivial
#              Floquet multiplier (monodromy over one period T) for a limit cycle.  t_crit / τ is
#              what is compared across parameters: close to a bifurcation τ diverges and t_crit
#              with it, for no chaotic reason at all.
#   σ_d        spread of t_crit / τ over initial conditions at a FIXED distance d from the attractor
#              (a sphere around it in the Fubini-Study metric), for several d.
#   σ_ε        spread of t_crit / τ over perturbations of size ε of ONE initial condition: for a
#              regular transient it shrinks in proportion to ε, for a chaotic transient (a chaotic
#              saddle on the way) it stays finite however small ε is - the signature sought.
#   FTLE       finite-time Lyapunov exponent over [0, t_crit] of the REDUCED tangent dynamics (the
#              deviation vector is projected off the U(1) and norm directions ψ and iψ at every
#              renormalisation, otherwise the neutral symmetry directions win in the end and every
#              FTLE tends to 0).  A positive FTLE over the transient is transient chaos.
#
# ATTRACTORS.  Found at every point from REFERENCES long trajectories.  A relative fixed point is a
# single state; a limit cycle is represented by its crossings of a Poincare section (one of the
# three surfaces n_i = n_j, upward), where the trajectory's own crossings are compared with them.
# Symmetry images (the L translations, and at η = 0 the reflection) are added, because a trajectory
# can land on any copy.  An attractor with more than MAX_SECTION_POINTS distinct crossings is a
# torus or strange: no 1e-4 neighbourhood can be resolved on it, and the trajectories reaching it
# are reported with t_crit = NaN (and counted); so is a point whose MAIN attractor is of that kind,
# without running its trajectories.  In the two regions of the TODO item the attractors are fixed
# points (η = 0) and mostly limit cycles (η = 3); the tests met a torus or strange attractor from
# g ≈ -7 on at η = 3, where the cut ends.
#
# USAGE (multithreaded over the trajectories of one point):
#
#   julia -t 20 BHNumberConservingTransient.jl eta0             # η = 0, g = -4 ... -50
#   julia -t 20 BHNumberConservingTransient.jl eta3             # η = 3, g = g_c ... -8
#   julia -t 20 BHNumberConservingTransient.jl point <g> <η>    # one point, printed in detail
#
# OUTPUT in $BH_RESULTS_DIR (default ~/results/bh/number-conserving/classical/3/transient):
#   summary_<task>.txt       one line per g (header of column names; np.genfromtxt names=True);
#                            fixed_points and cycles count distinct attractors (with their symmetry
#                            images), complex the reference trajectories that ended on a torus or a
#                            strange set (those are not grouped); all-NaN columns mean that the main
#                            attractor has no resolvable 1e-4 neighbourhood
#   <task>/g<g>_e<η>.txt     every trajectory: kind (d / eps), value, base, t_crit, t_crit/τ,
#                            FTLE, attractor index
# Both are resumable: a point with a raw file is skipped.

include(joinpath(@__DIR__, "BHNumberConserving.jl"))

using Printf
using Random
using Statistics
using LinearAlgebra

const L_RING = 3
const κ_RING = 0.3

const TOLERANCE = 1e-4
const MAXIMUM_TIME = 20000.0
const REFERENCES = 24
const RELAXATION = 3000.0
const MAX_SECTION_POINTS = 12

const DISTANCES = (0.05, 0.2, 0.5, 1.0)          # σ_d; π/2 is the diameter of CP^(L-1)
const SAMPLES_D = 100
const EPSILONS = (1e-8, 1e-6, 1e-4, 1e-2)        # σ_ε
const BASES = 5
const SAMPLES_ε = 40

const RESULTS = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "number-conserving", "classical", "3", "transient"))


# ---------------------------------------------------------------------------------------------
# Geometry of CP^(L-1)
# ---------------------------------------------------------------------------------------------

FubiniStudy(ψ, φ) = acos(clamp(abs(dot(ψ, φ)) / sqrt(real(dot(ψ, ψ)) * real(dot(φ, φ))), 0.0, 1.0))

""" A state at Fubini-Study distance exactly d from ψ, in a uniformly random direction. """
function AtDistance(ψ, d, rng)
    ψ = ψ ./ norm(ψ)
    χ = randn(rng, ComplexF64, length(ψ))
    χ .-= dot(ψ, χ) .* ψ
    χ ./= norm(χ)
    return cos(d) .* ψ .+ sin(d) .* χ
end

""" The images of a state under the symmetries of the ring: the L translations, and with them the
    reflection j -> 2 - j when it is a symmetry (η = 0, note §6). """
function Images(ψ, reflect)
    L = length(ψ)
    images = [circshift(ψ, r) for r = 0:(L - 1)]
    if reflect
        mirrored = [ψ[mod1(2 - j, L)] for j = 1:L]
        append!(images, [circshift(mirrored, r) for r = 0:(L - 1)])
    end
    return images
end

# The three Poincare sections n_i - n_j = 0, crossed upwards.
const SECTIONS = ((1, 2), (2, 3), (1, 3))
SectionValue(ψ, k) = abs2(ψ[SECTIONS[k][1]]) - abs2(ψ[SECTIONS[k][2]])


# ---------------------------------------------------------------------------------------------
# Plain trajectories and the attractors
# ---------------------------------------------------------------------------------------------

function PlainProblem(ψ0, parameters, time)
    work = NumberConservingWorkspace(parameters.L)
    function Field!(du, u, p, t)
        UpdateWorkspace!(work, u, parameters)
        VectorField!(du, u, parameters, work)
    end
    return ODEProblem(Field!, ComplexF64.(ψ0 ./ norm(ψ0)), (0.0, float(time)))
end

Evolve(ψ0, parameters, time) =
    solve(PlainProblem(ψ0, parameters, time), DP8(); reltol = 1e-10, abstol = 1e-10,
          save_everystep = false, maxiters = 10^8).u[end]

""" ||F(ψ) + i ω ψ|| minimised over the phase velocity ω: zero on a relative fixed point. """
function Stationarity(ψ, parameters)
    field = VectorField(ψ, parameters)
    ω = -imag(dot(ψ, field)) / real(dot(ψ, ψ))
    return norm(field .+ im * ω .* ψ), ω
end

""" Upward crossings of the section k over `time`, starting on the attractor at ψ. """
function SectionCrossings(ψ, parameters, k, time)
    crossings = Vector{ComplexF64}[]
    callback = ContinuousCallback((u, t, integrator) -> SectionValue(u, k),
                                  integrator -> push!(crossings, integrator.u ./ norm(integrator.u)),
                                  nothing; save_positions = (false, false))
    solve(PlainProblem(ψ, parameters, time), DP8(); reltol = 1e-11, abstol = 1e-11,
          callback = callback, save_everystep = false, maxiters = 10^8)
    return crossings
end

""" Crossings merged when closer than `tolerance` (a limit cycle crosses the section at a few
    points, every period again). """
function DistinctPoints(points; tolerance = 1e-6)
    distinct = Vector{ComplexF64}[]
    for p in points
        any(FubiniStudy(p, q) < tolerance for q in distinct) || push!(distinct, p)
    end
    return distinct
end

""" Every attractor of the point, with its symmetry images.  Each is a NamedTuple
    (kind = :fixed | :cycle | :complex, state, points, section, basin, τ). """
function FindAttractors(parameters; references = REFERENCES, seed = 1, reflect = false)
    finals = Vector{Vector{ComplexF64}}(undef, references)
    Threads.@threads for i = 1:references
        rng = Xoshiro(hash((seed, i, parameters.g)))
        ψ0 = i == 1 ? UniformInitialCondition(parameters.L; amplitude = 1e-3, rng = rng) :
                      RandomInitialCondition(parameters.L, rng)
        finals[i] = Evolve(ψ0, parameters, RELAXATION)
    end

    attractors = []
    for ψ in finals
        residual, _ = Stationarity(ψ, parameters)

        if residual < 1e-7
            # already known (as itself or as an image)?
            known = findfirst(a -> a.kind === :fixed && FubiniStudy(ψ, a.state) < 1e-5, attractors)
            if isnothing(known)
                τ = FixedPointTime(ψ, parameters)
                for image in DistinctPoints(Images(ψ, reflect); tolerance = 1e-5)
                    push!(attractors, (kind = :fixed, state = image, points = [image], section = 0,
                                       basin = Ref(0), τ = τ))
                end
                known = findfirst(a -> a.kind === :fixed && FubiniStudy(ψ, a.state) < 1e-5, attractors)
            end
            attractors[known].basin[] += 1
            continue
        end

        # not stationary: is it on a limit cycle already known?
        known = nothing
        for (index, a) in enumerate(attractors)
            a.kind === :cycle || continue
            crossings = SectionCrossings(ψ, parameters, a.section, 200.0)
            if !isempty(crossings) && any(FubiniStudy(c, p) < 1e-4 for c in crossings for p in a.points)
                known = index
                break
            end
        end
        if !isnothing(known)
            attractors[known].basin[] += 1
            continue
        end

        section = 0
        points = Vector{ComplexF64}[]
        for k in eachindex(SECTIONS)
            crossings = SectionCrossings(ψ, parameters, k, 400.0)
            if length(crossings) >= 2
                section = k
                points = DistinctPoints(crossings)
                break
            end
        end

        if section == 0 || length(points) > MAX_SECTION_POINTS
            push!(attractors, (kind = :complex, state = ψ, points = points, section = section,
                               basin = Ref(1), τ = NaN))
            continue
        end

        τ = CycleTime(points[1], section, parameters)
        if isnan(τ)
            # The section points never recur exactly: a torus, or a strange set whose crossings
            # happened to cluster.  The Lyapunov spectrum decides; only a genuine limit cycle keeps
            # a 1e-4 neighbourhood that trajectories can be timed into.
            spectrum = LyapunovSpectrum(ψ, parameters; relaxationTime = 10.0, integrationTime = 3010.0,
                                        rng = Xoshiro(7), zeroThreshold = 2e-3)
            contracting = filter(<(-2e-3), spectrum.reduced)
            if spectrum.classification !== :limitCycle || isempty(contracting)
                push!(attractors, (kind = :complex, state = ψ, points = points, section = section,
                                   basin = Ref(1), τ = NaN))
                continue
            end
            τ = 1 / abs(maximum(contracting))
        end
        first = length(attractors) + 1
        for image in Images(ψ, reflect)
            # a translated or reflected cycle need not cross the same section: take the first of
            # the three that it does cross
            imagePoints, imageSection = Vector{ComplexF64}[], 0
            for k in circshift(collect(eachindex(SECTIONS)), 1 - section)
                crossings = SectionCrossings(image, parameters, k, 400.0)
                if length(crossings) >= 2
                    imagePoints, imageSection = DistinctPoints(crossings), k
                    break
                end
            end
            imageSection == 0 && continue
            # an image that is the same cycle adds nothing
            any(a.kind === :cycle && a.section == imageSection &&
                any(FubiniStudy(p, q) < 1e-5 for p in imagePoints for q in a.points)
                for a in attractors) && continue
            push!(attractors, (kind = :cycle, state = image, points = imagePoints,
                               section = imageSection, basin = Ref(0), τ = τ))
        end
        attractors[first].basin[] += 1
    end

    return attractors
end

""" 1 / |max Re μ| over the reduced linearisation of a relative fixed point, in the co-rotating
    frame.  The two eigenvalues closest to zero are the exact symmetry zeros (phase and norm) and
    are dropped.  NaN if the point is linearly unstable (then it is not the attractor). """
function FixedPointTime(ψ, parameters)
    _, ω = Stationarity(ψ, parameters)
    L = parameters.L
    rotation = [zeros(L, L) -Matrix(I, L, L); Matrix(I, L, L) zeros(L, L)]
    μ = eigvals(JacobianMatrix(ψ, parameters) .+ ω .* rotation)
    reduced = μ[sortperm(abs.(μ))[3:end]]
    rate = maximum(real, reduced)
    return rate < 0 ? 1 / abs(rate) : NaN
end

""" Relaxation time onto a limit cycle from its Floquet multipliers, and the period.

    The period T is the time of the first return of the section point p0 to itself.  The monodromy
    matrix is the real 2L x 2L tangent propagator over one period, with the global phase the orbit
    has picked up (ψ(T) = exp(iθ) p0 - a RELATIVE periodic orbit) rotated back.  Three multipliers
    are then trivially 1 - the flow direction, the U(1) phase and the conserved norm - and are
    dropped; the largest of the rest gives the slowest contraction, τ = T / |log |μ||.  Exact up to
    the integration tolerance, unlike a finite-time Lyapunov estimate, which near a marginal cycle
    (g = -6, see K9) cannot separate a slow contraction from zero. """
function CycleTime(p0, section, parameters; maximumCrossings = 50)
    L = parameters.L
    p0 = p0 ./ norm(p0)

    period = NaN
    count_ = Ref(0)
    function Returned!(integrator)
        count_[] += 1
        if FubiniStudy(view(integrator.u, 1:L), p0) < 1e-7 && integrator.t > 1e-6
            period = integrator.t
            terminate!(integrator)
        elseif count_[] >= maximumCrossings
            terminate!(integrator)
        end
    end
    callback = ContinuousCallback((u, t, integrator) -> SectionValue(u, section), Returned!, nothing;
                                  save_positions = (false, false))
    solve(PlainProblem(p0, parameters, 1e4), DP8(); reltol = 1e-12, abstol = 1e-12,
          callback = callback, save_everystep = false, maxiters = 10^8)
    isfinite(period) || return NaN

    work = NumberConservingWorkspace(L)
    function Monodromy!(du, u, p, t)
        ψ = view(u, 1:L)
        UpdateWorkspace!(work, ψ, parameters)
        VectorField!(view(du, 1:L), ψ, parameters, work)
        for k = 1:(2 * L)
            TangentField!(view(du, (k * L + 1):((k + 1) * L)), view(u, (k * L + 1):((k + 1) * L)), ψ,
                          parameters, work)
        end
    end

    u0 = zeros(ComplexF64, L * (2 * L + 1))
    u0[1:L] .= p0
    for k = 1:(2 * L)
        u0[k * L + (k <= L ? k : k - L)] = k <= L ? 1.0 : im
    end
    uT = solve(ODEProblem(Monodromy!, u0, (0.0, period)), DP8(); reltol = 1e-12, abstol = 1e-12,
               save_everystep = false, maxiters = 10^8).u[end]

    rotation = exp(-im * angle(dot(p0, uT[1:L])))
    M = zeros(2 * L, 2 * L)
    for k = 1:(2 * L)
        δ = uT[(k * L + 1):((k + 1) * L)] .* rotation
        M[1:L, k] .= real.(δ)
        M[(L + 1):end, k] .= imag.(δ)
    end

    μ = eigvals(M)
    nontrivial = μ[sortperm(abs.(μ .- 1))[4:end]]
    rate = log(maximum(abs, nontrivial)) / period
    return rate < 0 ? 1 / abs(rate) : NaN
end


# ---------------------------------------------------------------------------------------------
# One transient: t_crit and the FTLE
# ---------------------------------------------------------------------------------------------

mutable struct TransientTracker
    parameters::NumberConservingParameters
    work::NumberConservingWorkspace
    attractors::Vector{Any}
    lastTime::Float64
    lastDistance::Float64
    lastCrossing::Vector{Float64}          # per section: time and distance of the previous crossing
    lastCrossingDistance::Vector{Float64}
    crossTime::Float64
    reached::Int
    logGrowth::Float64
end

function TransientField!(du, u, tracker, t)
    L = tracker.parameters.L
    ψ = view(u, 1:L)
    UpdateWorkspace!(tracker.work, ψ, tracker.parameters)
    VectorField!(view(du, 1:L), ψ, tracker.parameters, tracker.work)
    TangentField!(view(du, (L + 1):(2 * L)), view(u, (L + 1):(2 * L)), ψ, tracker.parameters,
                  tracker.work)
    return nothing
end

""" The deviation vector off the symmetry directions ψ and iψ (real inner product), and its norm
    - the reduced tangent dynamics. """
function ReducedDeviation!(δ, ψ)
    ψn = ψ ./ norm(ψ)
    for direction in (ψn, im .* ψn)
        δ .-= RealDot(direction, δ) .* direction
    end
    return norm(δ)
end

function Renormalise!(integrator)
    tracker = integrator.p
    L = tracker.parameters.L
    δ = view(integrator.u, (L + 1):(2 * L))
    growth = ReducedDeviation!(δ, view(integrator.u, 1:L))
    tracker.logGrowth += log(growth)
    δ ./= growth
    u_modified!(integrator, true)
end

function CheckFixedPoints!(integrator)
    tracker = integrator.p
    ψ = view(integrator.u, 1:tracker.parameters.L)
    distance, index = Inf, 0
    for (k, a) in enumerate(tracker.attractors)
        a.kind === :fixed || continue
        d = FubiniStudy(ψ, a.state)
        d < distance && ((distance, index) = (d, k))
    end

    if distance < TOLERANCE
        t = integrator.t
        # the approach is exponential, so log(distance) is interpolated linearly in time
        if isfinite(tracker.lastDistance) && tracker.lastDistance > TOLERANCE
            t = tracker.lastTime + (t - tracker.lastTime) *
                log(tracker.lastDistance / TOLERANCE) / log(tracker.lastDistance / distance)
        end
        tracker.crossTime = t
        tracker.reached = index
        terminate!(integrator)
    end
    tracker.lastTime = integrator.t
    tracker.lastDistance = distance
end

function CheckCycles!(integrator, section)
    tracker = integrator.p
    ψ = view(integrator.u, 1:tracker.parameters.L)
    distance, index = Inf, 0
    for (k, a) in enumerate(tracker.attractors)
        (a.kind === :cycle && a.section == section && !isempty(a.points)) || continue
        d = minimum(FubiniStudy(ψ, p) for p in a.points)
        d < distance && ((distance, index) = (d, k))
    end
    index == 0 && return

    if distance < TOLERANCE
        t = integrator.t
        previous, previousDistance = tracker.lastCrossing[section], tracker.lastCrossingDistance[section]
        if isfinite(previousDistance) && previousDistance > TOLERANCE
            t = previous + (t - previous) * log(previousDistance / TOLERANCE) / log(previousDistance / distance)
        end
        tracker.crossTime = t
        tracker.reached = index
        terminate!(integrator)
    end
    tracker.lastCrossing[section] = integrator.t
    tracker.lastCrossingDistance[section] = distance
end

""" t_crit, the attractor reached and the reduced FTLE over [0, t_crit] of one initial state. """
function Transient(ψ0, parameters, attractors, rng)
    L = parameters.L
    u0 = zeros(ComplexF64, 2 * L)
    u0[1:L] .= ψ0 ./ norm(ψ0)
    δ = randn(rng, ComplexF64, L)
    ReducedDeviation!(δ, u0[1:L])
    u0[(L + 1):end] .= δ ./ norm(δ)

    tracker = TransientTracker(parameters, NumberConservingWorkspace(L), attractors, 0.0, Inf,
                               fill(0.0, 3), fill(Inf, 3), NaN, 0, 0.0)

    callbacks = Any[PeriodicCallback(Renormalise!, 1.0; save_positions = (false, false))]
    any(a.kind === :fixed for a in attractors) &&
        push!(callbacks, PeriodicCallback(CheckFixedPoints!, 0.05; save_positions = (false, false)))
    for section in unique(a.section for a in attractors if a.kind === :cycle)
        push!(callbacks, ContinuousCallback((u, t, integrator) -> SectionValue(view(u, 1:L), section),
                                            integrator -> CheckCycles!(integrator, section), nothing;
                                            save_positions = (false, false)))
    end

    problem = ODEProblem(TransientField!, u0, (0.0, MAXIMUM_TIME), tracker)
    solution = solve(problem, DP8(); reltol = 1e-10, abstol = 1e-10, callback = CallbackSet(callbacks...),
                     save_everystep = false, maxiters = 10^8)

    # growth since the last renormalisation
    final = solution.u[end]
    δ = copy(final[(L + 1):end])
    elapsed = solution.t[end]
    ftle = (tracker.logGrowth + log(ReducedDeviation!(δ, final[1:L]))) / max(elapsed, 1e-12)

    return (t = tracker.crossTime, attractor = tracker.reached, ftle = ftle)
end


# ---------------------------------------------------------------------------------------------
# One parameter point
# ---------------------------------------------------------------------------------------------

""" Every trajectory of one point: SAMPLES_D initial conditions on each Fubini-Study sphere of
    radius d around the main attractor (σ_d), and SAMPLES_ε perturbations of size ε of BASES
    random initial conditions (σ_ε). """
function TransientPoint(g, η; verbose = true, seed = 20260923)
    parameters = NumberConservingParameters(L_RING; J = 1.0, g = g, κ = κ_RING, η = η)
    reflect = η == 0
    attractors = FindAttractors(parameters; reflect = reflect)

    main = argmax([a.basin[] for a in attractors])
    centre = attractors[main]

    if verbose
        @printf("g = %.3f, η = %.3f: %d attractors (with images)\n", g, η, length(attractors))
        for (k, a) in enumerate(attractors)
            a.basin[] > 0 || continue
            @printf("  %2d  %-8s basin %2d/%d  τ = %8.3f  section points %d\n", k, a.kind,
                    a.basin[], REFERENCES, a.τ, length(a.points))
        end
    end

    # A torus or a strange attractor has no resolvable 1e-4 neighbourhood: there is no t_crit to
    # measure, and every trajectory would run to MAXIMUM_TIME.  The point is recorded as such.
    if centre.kind === :complex
        verbose && println("  the main attractor is a torus or strange - no t_crit, point recorded as NaN")
        return (parameters = parameters, attractors = attractors, rows = [], main = main)
    end

    # the centre of the distance spheres: the fixed point itself, or a point of the cycle
    function CentreState(rng)
        centre.kind === :fixed && return centre.state
        return Evolve(centre.state, parameters, 50.0 * rand(rng))
    end

    jobs = []
    for d in DISTANCES, i = 1:SAMPLES_D
        push!(jobs, (kind = "d", value = d, base = 0, index = i))
    end
    for b = 1:BASES, ε in EPSILONS, i = 1:SAMPLES_ε
        push!(jobs, (kind = "eps", value = ε, base = b, index = i))
    end

    rows = Vector{Any}(undef, length(jobs))
    Threads.@threads for j in eachindex(jobs)
        job = jobs[j]
        rng = Xoshiro(hash((seed, g, η, job.kind, job.value, job.base, job.index)))
        if job.kind == "d"
            ψ0 = AtDistance(CentreState(rng), job.value, rng)
        else
            base = RandomInitialCondition(L_RING, Xoshiro(hash((seed, g, η, "base", job.base))))
            ψ0 = AtDistance(base, job.value, rng)
        end
        result = Transient(ψ0, parameters, attractors, rng)
        τ = result.attractor > 0 ? attractors[result.attractor].τ : NaN
        rows[j] = (job..., t = result.t, tnorm = result.t / τ, ftle = result.ftle,
                   attractor = result.attractor)
    end

    return (parameters = parameters, attractors = attractors, rows = rows, main = main)
end


""" The summary line of one point: per d, mean / std of t_crit/τ, converged share, mean FTLE and
    the share of positive FTLEs; per ε, the median over the bases of σ_ε and the mean t_crit/τ. """
function Summary(g, η, point)
    rows = point.rows
    a = point.attractors[point.main]
    fixed = count(x -> x.kind === :fixed && x.basin[] > 0, point.attractors)
    cycles = count(x -> x.kind === :cycle && x.basin[] > 0, point.attractors)
    complex_ = count(x -> x.kind === :complex, point.attractors)

    values = Any[g, η, fixed, cycles, complex_, a.τ]
    for d in DISTANCES
        selection = [r for r in rows if r.kind == "d" && r.value == d]
        converged = [r for r in selection if isfinite(r.tnorm)]
        tnorm = [r.tnorm for r in converged]
        ftle = [r.ftle for r in converged]
        append!(values, [isempty(tnorm) ? NaN : mean(tnorm), length(tnorm) > 1 ? std(tnorm) : NaN,
                         length(converged) / length(selection),
                         isempty(ftle) ? NaN : mean(ftle),
                         isempty(ftle) ? NaN : count(>(1e-3), ftle) / length(ftle)])
    end
    for ε in EPSILONS
        spreads = Float64[]
        means = Float64[]
        for b = 1:BASES
            tnorm = [r.tnorm for r in rows if r.kind == "eps" && r.value == ε && r.base == b && isfinite(r.tnorm)]
            length(tnorm) > 1 && push!(spreads, std(tnorm))
            isempty(tnorm) || push!(means, mean(tnorm))
        end
        append!(values, [isempty(spreads) ? NaN : median(spreads), isempty(means) ? NaN : mean(means)])
    end
    return values
end

function SummaryNames()
    names = ["g", "eta", "fixed_points", "cycles", "complex", "tau"]
    for d in DISTANCES
        tag = @sprintf("d%g", d)
        append!(names, [tag * "_mean", tag * "_sigma", tag * "_converged", tag * "_ftle", tag * "_ftle_positive"])
    end
    for ε in EPSILONS
        tag = @sprintf("eps%.0e", ε)
        append!(names, [tag * "_sigma", tag * "_mean"])
    end
    return names
end


function RunTask(task, gs, η)
    directory = joinpath(RESULTS, task)
    mkpath(directory)
    summaryFile = joinpath(RESULTS, "summary_$task.txt")
    isfile(summaryFile) || open(io -> println(io, "# ", join(SummaryNames(), "\t")), summaryFile, "w")

    for g in gs
        raw = joinpath(directory, @sprintf("g%.4f_e%.4f.txt", g, η))
        if isfile(raw)
            @printf("g = %.3f done, skipped\n", g)
            continue
        end

        time = @elapsed point = TransientPoint(g, η)
        values = Summary(g, η, point)

        temporary = raw * ".part"
        open(temporary, "w") do io
            println(io, "# kind\tvalue\tbase\tt_crit\tt_over_tau\tftle\tattractor")
            for r in point.rows
                @printf(io, "%s\t%.3e\t%d\t%.8g\t%.8g\t%.8g\t%d\n", r.kind, r.value, r.base, r.t,
                        r.tnorm, r.ftle, r.attractor)
            end
        end
        mv(temporary, raw; force = true)
        open(io -> println(io, join([v isa Integer ? string(v) : @sprintf("%.6g", v) for v in values], "\t")),
             summaryFile, "a")

        @printf("  t_crit/τ at d = %.2f: %.2f ± %.2f;  σ_ε at ε = 1e-8: %.3g;  mean FTLE (d = %.2f): %+.4f  (%.1f s)\n",
                DISTANCES[end], values[6 + 5 * (length(DISTANCES) - 1) + 1],
                values[6 + 5 * (length(DISTANCES) - 1) + 2], values[6 + 5 * length(DISTANCES) + 1],
                DISTANCES[end], values[6 + 5 * (length(DISTANCES) - 1) + 4], time)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    task = isempty(ARGS) ? "" : ARGS[1]
    @printf("Threads: %d\n", Threads.nthreads())

    if task == "eta0"
        RunTask("eta0", -4.0:-1.0:-50.0, 0.0)
    elseif task == "eta3"
        threshold = InstabilityThreshold(L_RING; J = 1.0, κ = κ_RING, η = 3.0)
        @printf("η = 3: g_c = %.4f\n", threshold)
        RunTask("eta3", collect((floor(threshold * 10) / 10):-0.1:-8.0), 3.0)
    elseif task == "point"
        g, η = parse(Float64, ARGS[2]), parse(Float64, ARGS[3])
        point = TransientPoint(g, η)
        columnNames = SummaryNames()
        for (name, value) in zip(columnNames, Summary(g, η, point))
            @printf("  %-22s %s\n", name, value isa Integer ? string(value) : @sprintf("%.6g", value))
        end
    else
        println("usage: julia -t <threads> BHNumberConservingTransient.jl eta0 | eta3 | point <g> <η>")
    end
end
