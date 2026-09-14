using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Distributed
using Printf

# (g, d) phase diagram of the INCOHERENTLY PUMPED dissipative Bose-Hubbard model of BHPumped.jl -
# the discrete complex Ginzburg-Landau lattice of pump_bh_cgle.py. Same Distributed/pmap structure,
# per-trajectory seeds and resumable output files as BHMapDrivenDissipative.jl; the swept plane, the
# initial conditions and the meaning of the exponent are what differ.
#
# Why (g, d) - and why that is the WHOLE phase diagram:
#
#   * The flow dψ_j/dt = (d + iJ) Lap ψ_j + (P/2) ψ_j - (i g + Γ/2)|ψ_j|² ψ_j has five constants but
#     only THREE dimensionless ones. Rescaling time by P and the amplitude by √(P/Γ),
#     τ = P t and φ = ψ √(Γ/P), turns it into
#         dφ/dτ = ((d + iJ)/P) Lap φ + φ/2 - (i g/Γ + 1/2) |φ|² φ,
#     so only J/P, d/P and g/Γ matter. Fixing J = P = Γ = 1, as pump_bh_cgle.py does, is therefore
#     no loss of generality at all: the (g, d) plane below IS the phase diagram of the model, and
#     every other parameter set is one of its points in disguise.
#   * Both analytic boundaries live in exactly this plane and can be overlaid on the map in closed
#     form: the Benjamin-Feir line g = -d Γ/(2J) of the continuum, and the edge of the region where
#     at least one of the L plane waves is stable - its lower branch d*(g) = CriticalDiffusion(J, g,
#     P, Γ; L = L), and for g ≳ -1 an upper, re-entrant branch beyond which a wave is stable again.
#   * The route the note describes - uniform state -> plane waves of every Q -> no stable wave at all
#     -> strange attractor - is exactly a d-sweep at fixed g.
#
# TWO RULES fix which quadrant is worth computing, and the grid below satisfies both:
#
#   J g < 0   or the Benjamin-Feir criterion J g < -d Γ/2 cannot hold for any d ≥ 0 and the uniform
#             state is stable however strong the nonlinearity is. Measured: at J = +1 and g = +0.5
#             or +1 every one of 8 random initial conditions ends on a regular attractor, at every d
#             from 0 to 0.8. Hence J = +1 here and THE GRID RUNS OVER NEGATIVE g.
#   d > 0     because J g < 0 is necessary but NOT sufficient on a lattice: at d = 0 the short waves
#             |Q| > π/2 have an effective mass J cos Q of the opposite sign and stay stable, and
#             every trajectory ends on one of them. Measured: Λ = 0 to 1e-5 for all of
#             g ∈ [-6, -0.5] at d = 0. The whole d = 0 edge of the map is therefore regular. Chaos
#             takes over where the diffusion has damped the last stable wave, above d* ≈ 0.071-0.073
#             on this ring - but for -1.7 ≤ g ≤ -0.4 it already coexists with the stable waves below
#             d*, in up to 93 % of the initial conditions of a cell.
#
# Ring size. Unlike the driven model, where L = 3 already sustains chaos, the CGLE needs enough
# modes: the maximum of Λ over 8 random initial conditions at g = -1, J = P = Γ = 1 is
#
#   d           0.00    0.05    0.10    0.15    0.20    0.30    0.40    0.60
#   L = 3, 4      0       0       0       0       0       0       0       0
#   L = 5         0    +0.14   +0.10       0       0       0       0       0
#   L = 6         0    +0.17   +0.15   +0.09   +0.09       0       0       0
#   L = 8         0    +0.07   +0.18   +0.16   +0.13   +0.12       0       0
#   L = 12        0    +0.18   +0.21   +0.18   +0.17   +0.15       0    +0.05
#   L = 16        0    +0.19   +0.20   +0.18   +0.17   +0.14   +0.13   +0.06
#
# L = 8 is the smallest ring whose picture is already the large-L one (L = 3 and 4 have no chaos
# anywhere, L = 5 and 6 only in a sliver), and it costs about 0.1 s per trajectory - the whole grid
# below is then some 17 hours on 10 workers, and it is resumable.
#
# Initial conditions and multistability. There is no energy shell and no absorbing ball here; what
# bounds the flow is the saturation of the pump itself, through the trapping region
# P/Γ ≤ Σ_i I_i ≤ L P/Γ of BHPumped.jl, whose upper end is the uniform state. Each (g, d) point is
# therefore sampled with TRAJECTORIES random initial conditions of density drawn uniformly in
# (0, L P/Γ], and the observable is the distribution over the coexisting attractors:
#
#   λ_max               the exponent on the attractor each trajectory reached
#   chaotic fraction    share of initial conditions with λ_max > CHAOS_THRESHOLD, i.e. the share of
#                       phase space draining into a strange attractor.
#
# Output: one line per trajectory, TWO columns
#
#   λ_max    ⟨Σ_i I_i⟩ on the attractor
#
# The second column is the filling n L, and it is far more informative here than the driven map's
# bistability branches: a trajectory that locks onto the plane wave of wavenumber Q sits at exactly
# Σ_i I_i = L (P - 2 d D_Q)/Γ with D_Q = 2(1 - cos Q), so wherever the regular attractors are plane
# waves the second column says WHICH wave the trajectory fell onto (up to the mirror image -Q), and
# its spread across the TRAJECTORIES samples is a direct picture of how the basins are shared out -
# analyse_map_pumped.py reads it exactly that way. Three limits, all visible in the complete map: the
# regular attractors are plane waves only where one is stable (below d*, and beyond the re-entrant
# edge) and match no filling elsewhere, e.g. in the regular fingers above the chaotic wedge; at d = 0
# every wave has the same filling L P/Γ; and on L = 8 the filling of Q = π/2 is shared by a family of
# uniform-modulus superpositions of ±π/2, some of them stable where the pure wave is not. Failed
# trajectories are
# written as NaN NaN; there is no "-1 = no initial condition" case, since no energy shell has to
# be hit. In numpy: data = np.loadtxt(...); lambdas = data[:, 0]; chaotic = lambdas > threshold.
#
# NOTE ON THE THRESHOLD. U(1) is unbroken without a drive, so every attractor carries the exact zero
# of the global phase mode. A fixed point and a limit cycle therefore BOTH give λ_max = 0 here, not
# the -κ/2 of the driven model - a regular point of this map is a zero, not a negative number, and
# the classification rests on telling that zero from the chaotic exponents, which reach +0.86 deep
# inside the wedge. Measured scatter of the zero over 10 random initial conditions at
# INTEGRATION_TIME = 2500: below 2e-3, so the threshold 1e-2 below clears the numerical zero by a
# factor of five. What no threshold can do is sharpen the EDGE of the wedge, where λ_max goes to
# zero continuously as the attractor stops being strange; the cells there are classified by where
# λ_max happens to cross 1e-2, and that is a property of the transition, not of the integration.

workers = 10

if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

# Set on the workers themselves: addprocs copies the environment at spawn time, so assigning to
# ENV on the master afterwards would no longer reach them and they would each load Plots/PyPlot.
@everywhere ENV["CD_NO_PLOTS"] = "true"

@everywhere using Logging
@everywhere global_logger(ConsoleLogger(stderr, Logging.Warn))

# Pulls in models/BoseHubbardFull.jl and modules/ClassicalDynamics.jl as well
@everywhere include("BHPumped.jl")

# Constants and parameters
const TRAJECTORIES = 500        # random initial conditions per (g, d) point; the basin statistics
const L = 8                     # see the table in the header: the smallest ring with the large-L
                                # picture. Allowed modes are k = 2πm/8, k = π among them.

# J = P = Γ = 1 is not a choice of units among many - by the rescaling in the header it is the
# GENERAL case, with g and d carrying everything that is left.
const J = 1.0                   # J g < 0 is required, hence the negative g grid below
const P = 1.0                   # net linear gain γ_pump - γ_loss, the -κ of BHDissipative.jl
const Γ = 1.0                   # two-body loss; it alone saturates the pump

# The constants below are read by SingleTrajectory, which runs on the workers, so they have to be
# defined there as well - a plain `const` would only ever exist on the master.

# The trajectory has to reach its attractor before the exponent is measured: everything before
# RELAXATION_TIME describes the approach, not the attractor. The exponent is then averaged over
# INTEGRATION_TIME - RELAXATION_TIME. (pump_bh_cgle.py uses the same 400 + 2500.)
@everywhere const RELAXATION_TIME = 500.0
@everywhere const INTEGRATION_TIME = 2500.0

# A fixed point and a limit cycle both give exactly 0 here (the U(1) phase mode), and the chaotic
# attractors reach +0.86 deep inside the chaotic wedge - see the note in the header.
@everywhere const CHAOS_THRESHOLD = 0.01

# Grid. Chaos fills a wedge that opens towards strong interaction. Measured on the complete map: its
# lower edge is the lattice line d*(g) ≈ 0.071...0.073 of this ring (the 0.08...0.09 of
# pump_bh_cgle.py is the L = 256 value), with chaos leaking below it for -1.7 ≤ g ≤ -0.4; the
# chaotic fraction stays above 1/2 without a break up to d ≈ 0.36 at g = -1, 0.38 at g = -2, 0.54
# at g = -3 and all the way to d = 1 at g = -6, with regular fingers cutting in above that; and it
# is regular at d = 0, nowhere above 1/2 to the right of g = -0.5, and regular again beyond the
# re-entrant edge where a plane wave turns stable (d = 0.92 at g = -1 down to 0.19 at g = -0.4).
# g ∈ [-6, 0] × d ∈ [0, 1] holds all of it together with a regular margin.
@everywhere const G_VALUES = LinRange(-6.0, 0.0, 61)
@everywhere const D_VALUES = LinRange(0.0, 1.0, 51)
@everywhere const G_STEP = step(G_VALUES)
@everywhere const D_STEP = step(D_VALUES)

# Fraction of the cell over which (g, d) are randomised, per trajectory: 0 pins every trajectory to
# the exact grid point (the plain point sample), 1 spreads them uniformly over the whole cell.
#
# The same argument as in BHMapDrivenDissipative.jl: periodic windows alternate with chaotic bands
# on a scale comparable to the grid step, so a point sample ALIASES them - a grid point inside a
# window reports chaotic fraction 0 while its neighbours report 1, and a window narrower than a cell
# is either missed or reported as a hole, depending on where it falls. Randomising within the cell
# turns the point sample into a cell AVERAGE, a smooth function of (g, d) in which a narrow window
# shows up as an intermediate fraction. It costs nothing - the trajectories are run either way.
#
# Caveat: with JITTER > 0 the trajectories of one cell no longer share the same parameters, so an
# attractor count mixes genuine multistability with the variation across the cell. Set JITTER = 0
# when that particular map is what you are after.
@everywhere const JITTER = 0.0

@assert RELAXATION_TIME < INTEGRATION_TIME "the exponent is accumulated on (RELAXATION_TIME, INTEGRATION_TIME)"
@assert P > 0 "a pumped map needs P > 0: without a net gain the vacuum is stable and everything decays to ψ = 0"
@assert Γ > 0 "the two-body loss is the only saturation; with Γ = 0 nothing bounds the amplitude"
@assert all(g -> J * g <= 0, G_VALUES) "chaos needs J g < 0; with J g > 0 every trajectory ends on a stable uniform state"

const PATH = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "pumped", "$L",
             @sprintf("J_%.3f_P_%.3f_G_%.3f", J, P, Γ)))

# (g, d) of one trajectory: the cell centre displaced by up to half a cell in each direction, drawn
# from the same per-trajectory rng as the initial condition so that a resumed run keeps reproducing
# the same samples. g is clamped at 0 and d at 0, which only bites on the two edges of the grid.
# The parameter tuple carries U = g/2, so that is what is written back.
@everywhere function JitterParameters(parameters, rng)
    JITTER <= 0 && return parameters

    g = min(2 * parameters.U + JITTER * G_STEP * (rand(rng) - 0.5), 0.0)
    d = max(parameters.d + JITTER * D_STEP * (rand(rng) - 0.5), 0.0)

    return merge(parameters, (U = 0.5 * g, d = d))
end

# One trajectory: a random initial condition inside the trapping region, integrated until it settles
# on an attractor. `index` is the global index of the trajectory within its output file; it is mixed
# into the initial condition and into the seed of the deviation vector, so a resumed run continues
# with fresh samples instead of repeating the ones already stored.
# (A docstring cannot be attached to an @everywhere expression, hence the plain comment.)
@everywhere function SingleTrajectory(index, cellParameters)
    g = 2 * cellParameters.U
    rng = Xoshiro(hash((g, cellParameters.d, index)))

    parameters = JitterParameters(cellParameters, rng)
    initialCondition = PumpedInitialCondition(parameters, rng)

    _, lyapunov, _, _, observables = TrajectoryLyapunovPumped(initialCondition, parameters;
        seed=hash((g, cellParameters.d, index, :deviation)),
        relaxationTime=RELAXATION_TIME,
        timeInterval=(0.0, INTEGRATION_TIME))

    if isnan(lyapunov)
        return (NaN, NaN)               # nonconvergent, flagged by TrajectoryLyapunovPumped
    end

    # Density of the attractor, averaged over the post-transient part of the run
    norms = last.(observables.saveval)[observables.t .>= RELAXATION_TIME]

    return (lyapunov, isempty(norms) ? NaN : mean(norms))
end

""" All trajectories of one (g, d) point. `indices` are the global trajectory indices still missing
    from the output file, so a resumed run neither repeats nor reuses the samples already stored. """
function LyapunovMap(parameters, indices)
    time = @elapsed result = pmap(index -> SingleTrajectory(index, parameters), indices)

    g = 2 * parameters.U
    d = parameters.d

    lyapunovs = first.(result)
    valid = filter(isfinite, lyapunovs)
    chaotic = filter(v -> v > CHAOS_THRESHOLD, valid)

    println("Finished g = $g, d = $d (L = $L, J = $J, P = $P, Γ = $Γ)")
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
        @printf("⟨Σ I⟩ = %.3f (uniform state L P/Γ = %.3f)\n", mean(norms), L * P / Γ)
    end

    println("Elapsed time: $time seconds")
    println()

    return result
end

# g and d label the CELL here; with JITTER > 0 each trajectory samples a random point inside it.
for g in G_VALUES
    for d in D_VALUES
        mkpath(PATH)
        file = PATH * "/" * @sprintf("%.3f_%.3f", g, d) * ".txt"
        if isfile(file)
            trajectories = countlines(file)
        else
            trajectories = 0
        end

        println("Starting g = $g, d = $d")
        println("Computed trajectories: $trajectories, trajectories to compute: $(TRAJECTORIES - trajectories)")

        if trajectories >= TRAJECTORIES
            continue
        end

        parameters = PumpedParametersG(L, J, g, P, Γ; d=d)
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
end
