"""Heatmaps of the (g, d) sweep produced by BHMapPumped.jl.

The pumped counterpart of analyse_map_driven.py, for the incoherently pumped Bose-Hubbard ring with
two-body loss of BHPumped.jl and pump_bh_cgle.py. The file layout is the same:

  * the plane is (g, d) and each file is named g_d.txt;
  * each file has TWO columns, one line per trajectory,

        lambda_max    <sum_i I_i> on the attractor

    written from a random initial condition inside the trapping region. A failed trajectory is
    NaN NaN, and a cell whose file is missing or all-NaN is drawn white - which is also how a map
    that is still running looks.

Two properties of this model change how the two columns are read:

  * U(1) is unbroken, so every attractor carries the exact zero of the global phase mode. A fixed
    point and a limit cycle BOTH give lambda_max = 0 here, not a negative number, and chaos is
    lambda_max > THRESHOLD_CHAOS.
  * A trajectory on the plane wave psi_j = A e^{i(Qj - wt)}, Q = 2 pi m / L, has the exactly known
    filling

        sum_i I_i = L (P - 2 d D_Q) / Gamma,    D_Q = 2 (1 - cos Q),

    so the second column IDENTIFIES the wave: |m| is read off it (a wave and its mirror image -m
    share the filling and are counted once). Regular attractors that are not plane waves match no
    filling and stay unidentified. At d = 0 every wave has the same filling L P / Gamma, so there a
    trajectory is known to be on a plane wave but its |m| is not.

    One family shares a filling without being a single plane wave. When L is divisible by 4 the pair
    Q = +-pi/2 (|m| = L/4) spans a subspace the cubic term cannot leave, since 3 pi/2 = -pi/2 mod
    2 pi, so every a e^{i pi j/2} + b e^{-i pi j/2} with a uniform modulus (Re a b* = 0) is an exact
    stationary state with exactly the filling of the pure wave. |m| = L/4 stands for the whole
    family, and its mixed members have a stability of their own: re-running trajectories of the
    L = 8 map that the check below flags gave final states with 27-39 % of the power in one of the
    two waves, two zero exponents (the phase and the family) and all the others negative - stable,
    at points where the pure wave grows at +0.06 to +0.23.

The trajectories at one (g, d) are samples of the coexisting attractors, so the maps are basin
statistics:

    chaotic fraction      share of initial conditions reaching a strange attractor
    lambda (chaotic)      mean exponent over those, i.e. how strong the chaos is where it exists
    lambda (max)          strongest attractor found at the point
    filling n             <sum_i I_i> / L; the uniform state has n = P / Gamma
    plane-wave fraction   share of initial conditions ending on a plane wave
    plane waves by |m|    the same, split by wavenumber - one panel per |m|
    attractors            distinct coexisting attractors: one per |m| found, plus heuristic clusters
                          among the trajectories on no plane wave

The overlay plays the role of the fold tongue of analyse_map_driven.py, again in closed form:

  * the edge of the region where at least one of the L plane waves of the ring is linearly stable,
    from the Bogoliubov rates of pw_growth in pump_bh_cgle.py (PlaneWaveGrowth in BHPumped.jl);
  * the continuum Benjamin-Feir line J g = -d Gamma / 2, below which the uniform state is unstable.

The first is evaluated on a fine grid and contoured rather than taken from d_star of
pump_bh_cgle.py. d_star bisects, which assumes that stability, once lost, never returns as d grows -
and on a ring of 8 sites that is false: for -1.5 <= g <= -0.5 a wave becomes stable again at larger
d, the bisection lands in that upper window and reports that a stable wave always exists, and the
window without one disappears from the plot. The plane-wave fraction gives an independent check of
that edge: no trajectory can end on a plane wave where none is stable, and the script counts the
cells where the data says otherwise. Two kinds of cell are expected among them: the |m| = L/4 family
above, which the stability of the pure wave does not cover and which is therefore reported apart,
and cells a hair's breadth from an edge, where a wave is unstable but grows so slowly (1e-4 at 2e-4
from the edge) that the trajectory has not left it by the end of the run. Everything but the family
is listed cell by cell.
"""

import os
from fractions import Fraction

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

FIGSIZE = (9.2, 6.4)

CMAP = plt.get_cmap("viridis").copy()
CMAP.set_bad(color="white")

# Must match the constants of BHMapPumped.jl
L = 8
J = 1.0
P = 1.0
GAMMA = 1.0

# Same value as CHAOS_THRESHOLD in the Julia script: a fixed point and a limit cycle both give
# lambda_max = 0 (the U(1) phase mode), with a numerical scatter below 2e-3.
THRESHOLD_CHAOS = 0.01

# A regular trajectory is on the plane wave |m| when its filling agrees with L (P - 2 d D_Q) / Gamma
# to within this absolute tolerance. On a plane wave the stored filling is the time average of a
# constant and is exact to ~1e-7, while two different |m| are at least 2 L d min|D_Q - D_Q'| apart -
# 0.19 already at the first nonzero grid value d = 0.02 - so the window is wide on both sides.
TOLERANCE_PLANE_WAVE = 1e-4

# Two attractors that are not plane waves count as the same one when both coordinates agree to
# within these relative tolerances - the heuristic of analyse_map_driven.py.
TOLERANCE_LAMBDA = 0.05
TOLERANCE_NORM = 0.05

# Grid of BHMapPumped.jl
minG, maxG, numG = -6.0, 0.0, 61
minD, maxD, numD = 0.0, 1.0, 51

# Resolution of the analytic overlay. The onset of the lattice escape sits at d ~ 0.07, so d needs a
# much finer step than the data grid, and the upper edge runs so steeply through the plane that a
# coarse g step draws it as a staircase. About 3 s.
numGAnalytic, numDAnalytic = 961, 1001

PATH = f"c:/Users/micro/results/bh/pumped/{L}/J_{J:.3f}_P_{P:.3f}_G_{GAMMA:.3f}/"

STEM = f"pumped_{L}"

gs = np.linspace(minG, maxG, numG)
ds = np.linspace(minD, maxD, numD)

MODES = L // 2 + 1                                              # |m| = 0 ... L/2
FAMILY = L // 4 if L % 4 == 0 else None                         # the +-pi/2 family, see the docstring
DQ = 2 * (1 - np.cos(2 * np.pi * np.arange(MODES) / L))         # D_Q of each |m|


def data_file(g, d):
    """Name of the file written by BHMapPumped.jl for this grid point."""
    return f"{PATH}{g:.3f}_{d:.3f}.txt"


def plane_wave_modes(d, lambdas, norms):
    """|m| of the plane wave each trajectory ended on.

    Returns an int array: |m| for an identified wave, -1 for a trajectory on no plane wave (chaotic,
    or regular but not a plane wave), and -2 for one on a plane wave whose |m| cannot be resolved
    because several waves share its filling (at d = 0). |m| = FAMILY collects the uniform-modulus
    superpositions of +-pi/2 as well as the pure wave."""
    fillings = L * (P - 2 * d * DQ) / GAMMA                 # nonpositive: the wave does not exist
    match = (np.abs(norms[:, None] - fillings[None, :]) < TOLERANCE_PLANE_WAVE) & (fillings > 0)
    match &= (lambdas <= THRESHOLD_CHAOS)[:, None]

    matches = match.sum(axis=1)
    modes = np.where(matches == 1, np.argmax(match, axis=1), -1)
    modes[matches > 1] = -2

    return modes


def count_attractors(lambdas, norms, modes):
    """Number of distinct coexisting attractors among the sampled initial conditions, counted modulo
    the symmetries (the U(1) circle of every attractor, and the mirror image of a plane wave).

    Every |m| found is one attractor, all unresolved plane waves together count as one (so the count
    is a lower bound at d = 0), and the trajectories on no plane wave are grouped by the greedy
    single-pass clustering of analyse_map_driven.py: a trajectory joins an existing cluster when
    both coordinates match its representative to within the relative tolerances above, otherwise
    it opens a new one."""
    waves = len(np.unique(modes[modes >= 0])) + (1 if np.any(modes == -2) else 0)

    clusters = []
    for lam, norm in zip(lambdas[modes == -1], norms[modes == -1]):
        for clusterLambda, clusterNorm in clusters:
            closeLambda = abs(lam - clusterLambda) <= TOLERANCE_LAMBDA * max(1.0, abs(clusterLambda))
            closeNorm = abs(norm - clusterNorm) <= TOLERANCE_NORM * max(1.0, abs(clusterNorm))
            if closeLambda and closeNorm:
                break
        else:
            clusters.append((lam, norm))

    return waves + len(clusters)


shape = (numG, numD)
chaotic_fraction = np.full(shape, np.nan)
lyapunov_chaotic = np.full(shape, np.nan)
lyapunov_max = np.full(shape, np.nan)
filling = np.full(shape, np.nan)
plane_wave_fraction = np.full(shape, np.nan)
mode_fraction = np.full(shape + (MODES,), np.nan)
attractors = np.full(shape, np.nan)

missing = 0
for gi, g in enumerate(gs):
    for di, d in enumerate(ds):
        file = data_file(g, d)
        if not os.path.isfile(file):
            missing += 1                    # not computed yet - a running map is mostly this
            continue

        try:
            data = np.atleast_2d(np.loadtxt(file))

        except Exception as e:
            print(f"Error loading g = {g:.3f}, d = {d:.3f}: {e}")
            continue

        if data.size == 0:
            continue

        lambdas = data[:, 0]
        norms = data[:, 1]

        # NaN marks a trajectory the integrator could not finish; a cell with nothing left is
        # missing data and stays white.
        valid = np.isfinite(lambdas) & np.isfinite(norms)
        if not np.any(valid):
            continue

        lambdas = lambdas[valid]
        norms = norms[valid]

        chaotic = lambdas[lambdas > THRESHOLD_CHAOS]
        modes = plane_wave_modes(d, lambdas, norms)

        chaotic_fraction[gi, di] = len(chaotic) / len(lambdas)
        lyapunov_max[gi, di] = np.max(lambdas)
        filling[gi, di] = np.mean(norms) / L
        plane_wave_fraction[gi, di] = np.mean(modes != -1)
        attractors[gi, di] = count_attractors(lambdas, norms, modes)

        # |m| is unresolvable at d = 0: leave the cell white in every panel rather than guess
        if not np.any(modes == -2):
            mode_fraction[gi, di] = [np.mean(modes == m) for m in range(MODES)]

        # left as NaN (white) where no initial condition reached a strange attractor
        if len(chaotic) > 0:
            lyapunov_chaotic[gi, di] = np.mean(chaotic)

print(f"{numG * numD - missing} of {numG * numD} cells computed, {missing} missing (drawn white)")


def stable_plane_waves(g, d):
    """Number of the L plane waves of the ring that exist and are linearly stable, at one g and an
    array of d.

    The analytic 2x2 Bogoliubov blocks of pw_growth in pump_bh_cgle.py, vectorised over d, over the
    L waves Q and over the L perturbation wavenumbers k: the wave Q mixes only Q + k with Q - k, and
    its largest rate is the largest real part of [tr +- sqrt(tr^2 - 4 det)] / 2 over all k. The
    phase mode makes that rate exactly 0 on a stable wave, hence the tolerance."""
    d = np.asarray(d, dtype=float)[:, None, None]
    Q = (2 * np.pi * np.arange(L) / L)[:, None]
    k = (2 * np.pi * np.arange(L) / L)[None, :]

    DQwave = 2 * (1 - np.cos(Q))
    amplitude2 = (P - 2 * d * DQwave) / GAMMA

    W = d + 1j * J
    c = (1j * g + 0.5 * GAMMA) * amplitude2
    A11 = W * (DQwave - 2 * (1 - np.cos(Q + k))) - c
    A22 = np.conj(W) * (DQwave - 2 * (1 - np.cos(Q - k))) - np.conj(c)

    trace = A11 + A22
    determinant = A11 * A22 - c * np.conj(c)
    root = np.sqrt(trace * trace - 4 * determinant + 0j)
    growth = np.maximum((0.5 * (trace + root)).real, (0.5 * (trace - root)).real).max(axis=-1)

    return np.sum((amplitude2[..., 0] > 0) & (growth < 1e-10), axis=-1)


gsAnalytic = np.linspace(minG, maxG, numGAnalytic)
dsAnalytic = np.linspace(minD, maxD, numDAnalytic)
stable_count = np.array([stable_plane_waves(g, dsAnalytic) for g in gsAnalytic])

# Consistency check of the simulation, of the identification and of the overlay at once: a trajectory
# cannot end on a plane wave where none of them is stable. The |m| = FAMILY family is exempt, since the
# stability of the pure wave does not cover its mixed members, and is only counted; every other
# flagged cell is listed.
stable_at_cells = np.array([stable_plane_waves(g, ds) for g in gs])
flagged = [(gi, di, m) for gi in range(numG) for di in range(numD) for m in range(MODES)
           if stable_at_cells[gi, di] == 0 and mode_fraction[gi, di, m] > 0]

if FAMILY is not None:
    family = sum(1 for _, _, m in flagged if m == FAMILY)
    print(f"cells reaching the |m| = {FAMILY} family where the pure plane waves are all unstable: {family} "
          f"(expected - its mixed members are stable there)")

others = [(gi, di, m) for gi, di, m in flagged if m != FAMILY]
print(f"cells reaching any other plane wave where none is stable: {len(others)} (should be 0 away from an edge)")
for gi, di, m in others:
    print(f"   g = {gs[gi]:.3f}, d = {ds[di]:.3f}, |m| = {m}: {mode_fraction[gi, di, m]:.3f} of the trajectories")


# The white line gets a dark casing so that it stays visible over white (missing) cells as well
CASING = [patheffects.Stroke(linewidth=3.0, foreground="black"), patheffects.Normal()]


def PlotAnalytics(ax=None, add_legend=True):
    """Overlay the closed-form structure: the edge of the stable plane waves of the ring and the
    continuum Benjamin-Feir line."""
    if ax is None:
        ax = plt.gca()

    if np.any(stable_count > 0) and np.any(stable_count == 0):
        ax.contour(gsAnalytic, dsAnalytic, np.transpose(stable_count), levels=[0.5],
                   colors="black", linewidths=3.0)
        ax.contour(gsAnalytic, dsAnalytic, np.transpose(stable_count), levels=[0.5],
                   colors="white", linewidths=1.6)

    # J g = -d Gamma / 2; below it (smaller d) the uniform state is Benjamin-Feir unstable
    ax.plot(gsAnalytic, -2 * J * gsAnalytic / GAMMA, color="red", lw=1.4, ls="dashed")

    if add_legend:
        handles = [plt.Line2D([0], [0], color="white", lw=1.6, path_effects=CASING,
                              label=f"edge of stable plane waves (L = {L})"),
                   plt.Line2D([0], [0], color="red", lw=1.4, ls="dashed",
                              label="Benjamin-Feir (continuum)")]
        ax.legend(handles=handles, loc="upper left", fontsize=8)


def SaveAndShow(name):
    plt.savefig(f"{name}_{STEM}.png", dpi=150)
    plt.savefig(f"{name}_{STEM}.pdf")
    plt.show()


def PlotMap(values, label, title, name, cmap=CMAP, vmin=None, vmax=None):
    """Colormesh of `values` over the (g, d) plane with the analytic curves on top, saved as
    name_STEM.png / .pdf."""
    plt.figure(figsize=FIGSIZE)
    plt.pcolormesh(gs, ds, np.transpose(values), cmap=cmap, shading="auto", vmin=vmin, vmax=vmax)
    plt.colorbar(label=label)
    PlotAnalytics()
    plt.title(f"{title}, L = {L}, J = {J:g}, P = {P:g}, $\\Gamma$ = {GAMMA:g}")
    plt.xlabel("g")
    plt.ylabel("d")
    plt.xlim(gs[0], gs[-1])
    plt.ylim(ds[0], ds[-1])
    plt.tight_layout()
    SaveAndShow(name)


def WavenumberLabel(m):
    """Q = 2 pi m / L written as a multiple of pi, e.g. 3pi/4."""
    fraction = Fraction(2 * m, L)
    if fraction == 0:
        return "$Q = 0$"

    numerator = "" if fraction.numerator == 1 else str(fraction.numerator)
    denominator = "" if fraction.denominator == 1 else f"/{fraction.denominator}"
    return f"$Q = {numerator}\\pi{denominator}$"


def PlotModes(name):
    """Small multiples of the share of initial conditions ending on the plane wave |m|, one panel per
    |m| on a common 0-1 scale, saved as name_STEM.png / .pdf."""
    figure, axes = plt.subplots(1, MODES, figsize=(3.0 * MODES + 1.4, 4.8), sharex=True, sharey=True,
                                layout="constrained")

    for m, ax in enumerate(axes):
        mesh = ax.pcolormesh(gs, ds, np.transpose(mode_fraction[:, :, m]), cmap=CMAP, shading="auto",
                             vmin=0.0, vmax=1.0)
        PlotAnalytics(ax, add_legend=(m == 0))
        family = "\nand the uniform $\\pm Q$ mixtures" if m == FAMILY else ""
        ax.set_title(f"|m| = {m},  {WavenumberLabel(m)}{family}")
        ax.set_xlabel("g")
        ax.set_xlim(gs[0], gs[-1])
        ax.set_ylim(ds[0], ds[-1])

    axes[0].set_ylabel("d")
    figure.colorbar(mesh, ax=axes, label="fraction of initial conditions")
    figure.suptitle(f"Plane wave reached, by wavenumber, L = {L}, J = {J:g}, P = {P:g}, "
                    f"$\\Gamma$ = {GAMMA:g}  (d = 0 is white: |m| unresolvable there)")
    SaveAndShow(name)


PlotMap(chaotic_fraction, "chaotic fraction", "Fraction of initial conditions reaching chaos",
        "chaotic_fraction", vmin=0.0, vmax=1.0)
PlotMap(lyapunov_chaotic, "$\\Lambda$", "Lyapunov exponent on the chaotic attractors",
        "lyapunov_chaotic")
PlotMap(lyapunov_max, "$\\Lambda_{max}$", "Largest Lyapunov exponent found", "lyapunov_max")
PlotMap(filling, "n", "Filling per site on the attractor", "filling")
PlotMap(plane_wave_fraction, "plane-wave fraction", "Fraction of initial conditions ending on a plane wave",
        "plane_wave_fraction", vmin=0.0, vmax=1.0)
PlotModes("plane_wave_modes")
PlotMap(attractors, "attractors", "Number of coexisting attractors", "attractors")
