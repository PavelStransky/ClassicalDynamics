using DifferentialEquations
using LinearAlgebra
using Random
using Statistics
using Distributed
using Printf

# Dissipative counterpart of BHMap.jl: the same (J, E) sweep of the Lyapunov exponent, but for the
# OPEN Bose-Hubbard model of BHDissipative.jl - net damping κ, dephasing γ and the finite-N
# additive noise σ. Trajectories are parallelised within one job via Distributed/pmap, exactly as
# in BHMap.jl, and the output files keep the format expected by AnalyseLyapunov.py.
#
# Differences from BHMap.jl, all forced by the physics of the open system:
#
#   * ENERGY is the energy of the INITIAL CONDITION only. With κ ≠ 0 or γ ≠ 0 it is not conserved,
#     so the map reads "Lyapunov exponent as a function of the initial energy shell".
#   * There is no early stopping on convergence. The noise keeps the running exponent fluctuating,
#     so the plateau test of ClassicalDynamics.jl is unreliable; every trajectory is integrated for
#     the fixed INTEGRATION_TIME instead. This also makes the cost per point predictable.
#   * The number written to the file is the INTRINSIC exponent Λ + κ/2. The damping contracts every
#     phase-space direction equally and shifts every exponent by exactly -κ/2; removing it keeps the
#     stored values positive, which is what the `data[data > 0]` filter of AnalyseLyapunov.py and
#     the sentinels below assume. Set SAVE_INTRINSIC_EXPONENT = false to store the raw Λ.
#   * Every trajectory gets its own noise seed, derived deterministically from (J, E, index), so a
#     resumed run continues with fresh realisations instead of repeating the ones already stored.
#
# Sentinels in the output files, unchanged from BHMap.jl:  -1 = no initial condition found on the
# shell,  0 = nonconvergent / failed trajectory,  > 0 = a valid exponent.

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
@everywhere include("BHDissipative.jl")

# Random.seed!(1234)

# Constants and parameters
const TRAJECTORIES = 1000
const U = 1.0            # Float64 so modelParameters is a concrete NamedTuple -> type-stable equations of motion
const L = 4

# Dissipative rates, the counterpart of U above - see DissipativeParameters in BHDissipative.jl.
#   κ = γ_loss - γ_pump   net damping;  κ = 0 keeps Σ I_i = 1 and is the cleanest comparison with BHMap.jl
#   γ                     dephasing rate
#   σ                     finite-N (truncated Wigner) additive noise; 0 is the N -> infinity limit
const κ = 0.0
const γ = 0.1
const σ = 0.0

# The four constants below are read by SingleTrajectory, which runs on the workers, so they have to
# be defined there as well - a plain `const` would only ever exist on the master.

# Integration time of every trajectory. With κ > 0 the norm decays as e^{-κt}, the nonlinearity dies
# with it and the exponent drifts to the trivial -κ/2, so do not integrate much beyond a few 1/κ -
# use WINDOW below to read the exponent while the ring is still populated.
@everywhere const INTEGRATION_TIME = 2000.0

# nothing        -> the exponent is the trailing average of the whole run (right for κ = 0)
# (t1, t2)       -> the finite-time exponent on [t1, t2] via WindowLyapunov (right for κ > 0,
#                   e.g. (0.2/κ, 2/κ), while Σ I_i is still of order one)
@everywhere const WINDOW = nothing

# Splitting interval of the dephasing. Accuracy needs 2 U max(I) noiseStep < 0.2 and max(I) <= Σ I_i
# = 1, so anything below 0.1/U is safe here; the cost is proportional to 1/noiseStep.
@everywhere const NOISE_STEP = 0.05

@everywhere const SAVE_INTRINSIC_EXPONENT = true

const PATH = get(ENV, "BH_RESULTS_DIR",
    joinpath(homedir(), "results", "bh", "dissipative", "$L", @sprintf("k%.3f_g%.3f_s%.3f", κ, γ, σ)))

# One trajectory: a random initial condition on the given energy shell, integrated with the
# dissipative dynamics. `index` is the global index of the trajectory within its output file and is
# mixed into the noise seed, so every stored exponent comes from an independent realisation.
# (A docstring cannot be attached to an @everywhere expression, hence the plain comment.)
@everywhere function SingleTrajectory(index, energy, parameters, initialConditionEnergyTolerance)
    initialCondition = InitialCondition(energy, parameters, initialConditionEnergyTolerance)

    if initialCondition === nothing
        return -1.0
    end

    seed = hash((parameters.J, parameters.U, parameters.κ, parameters.γ, parameters.σ, energy, index))

    _, lyapunov, lyapunovs, _ = TrajectoryLyapunovDissipative(initialCondition, parameters;
        seed=seed,
        noiseStep=NOISE_STEP,
        timeInterval=(0.0, INTEGRATION_TIME))

    if lyapunov == 0.0
        return 0.0              # nonconvergent, flagged by TrajectoryLyapunovDissipative
    end

    if WINDOW !== nothing
        lyapunov = WindowLyapunov(lyapunovs, WINDOW[1], WINDOW[2])
    end

    return SAVE_INTRINSIC_EXPONENT ? lyapunov + 0.5 * parameters.κ : lyapunov
end

""" All trajectories of one (J, E) point. `indices` are the global trajectory indices still missing
    from the output file, so a resumed run neither repeats nor reuses the seeds already stored. """
function LyapunovMap(parameters, energy, indices; initialConditionEnergyTolerance=0.0001)
    time = @elapsed result = pmap(index -> SingleTrajectory(index, energy, parameters, initialConditionEnergyTolerance), indices)

    _, J, U, κ, γ, σ = parameters

    nonzero = filter(x -> x > 0, result)
    positive = filter(x -> x > 0.001, result)

    println("Finished J = $J, U = $U, κ = $κ, γ = $γ, σ = $σ, E = $energy");
    println("Number of new trajectories: $(length(result)) ($(length(positive)) unstable)");

    if length(positive) > 0
        print("Λ = $(mean(positive))")
        if length(positive) > 1
            print(" ± $(var(positive))")
        end
        println()
    end

    if length(nonzero) > 0
        println("freg = ", 1 - length(positive) / length(nonzero))
    end

    println("Elapsed time: $time seconds")
    println()

    return result, positive
end

for j in LinRange(0.0, 1.0, 51)
    for energy in LinRange(-0.5, 1.5, 101)
        mkpath(PATH)
        file = PATH * "/" * @sprintf("%.3f_%.3f_%.3f", j, U, energy) * ".txt"
        if isfile(file)
            trajectories = countlines(file)
        else
            trajectories = 0
        end

        println("Starting J = $j, U = $U, κ = $κ, γ = $γ, σ = $σ, E = $energy");
        println("Computed trajectories: $trajectories, trajectories to compute: $(TRAJECTORIES - trajectories)");

        if trajectories >= TRAJECTORIES
            continue
        end

        parameters = DissipativeParameters((L, j, U); κ=κ, γ=γ, σ=σ)
        lyapunovs, positive = LyapunovMap(parameters, energy, (trajectories + 1):TRAJECTORIES)

        if length(lyapunovs) == 0
            continue
        end

        open(file, "a") do io
            for lyapunov in lyapunovs
                println(io, lyapunov)
            end
        end
    end
end
