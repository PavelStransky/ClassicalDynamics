"""Heatmaps of the (Delta, f) sweep produced by BHMapDrivenDissipative.jl.

The driven counterpart of analyse_map.py. Two things differ from the undriven case:

  * the plane is (Delta, f) instead of (J, E), and each file is named Delta_f.txt;
  * each file has TWO columns, one line per trajectory,

        lambda_max    <sum_i I_i> on the attractor

    written from a random initial condition inside the absorbing ball. There is no "-1 = no
    initial condition" sentinel here (no energy shell has to be hit); a failed trajectory is
    NaN NaN, and a cell whose file is missing or all-NaN is drawn white.

Because the driven flow settles on an ATTRACTOR and the note reports strong multistability, the
trajectories at one (Delta, f) are samples of the coexisting attractors rather than repeats of one
number. The maps below are therefore basin statistics:

    chaotic fraction   share of initial conditions reaching a strange attractor - the driven
                       counterpart of 1 - freg in analyse_map.py
    lambda (chaotic)   mean exponent over those, i.e. how strong the chaos is where it exists
    lambda (max)       strongest attractor found at the point
    filling n          <sum_i I_i> / L, which is what makes the bistability visible
    attractors         heuristic count of distinct coexisting attractors

Chaotic bands and periodic windows alternate along f on a scale of 0.02-0.1, i.e. at or below the
0.1 grid step, so a point sample aliases them: isolated cells report a chaotic fraction of 0 (a
real periodic window) while their neighbours report 1, and chaotic bands falling between grid
points are missed altogether. BHMapDrivenDissipative.jl therefore has a JITTER constant that
randomises (Delta, f) within each cell, turning every cell into an average rather than a point
sample; with JITTER = 1 the maps below are smooth and a narrow window reads as an intermediate
fraction. Note that the attractor count is then a mixture of genuine multistability and the
variation across the cell - set JITTER = 0 if that map is what you want.

The overlay plays the role that the stationary-point CSV plays in analyse_map.py: instead of
ESQPT lines it draws the analytic structure of bh_dissipation_driving.md - the fold (saddle-node)
tongue of sec. 4 and the modulational-instability boundary of sec. 5 - both in closed form, so no CSV
is needed.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

FIGSIZE = (9.2, 6.4)

CMAP = plt.get_cmap("viridis").copy()
# CMAP = plt.get_cmap("gist_rainbow").copy()
CMAP.set_bad(color="white")

# Must match the constants of BHMapDrivenDissipative.jl
L = 2
J = -1.0
G = 2.0             # the note's g; the Julia side takes U = g/2
KAPPA = 1.0

# Same value as CHAOS_THRESHOLD in the Julia script: a fixed point gives lambda <= 0, a limit cycle
# exactly 0, chaos here runs from +0.26 to +0.71.
THRESHOLD_CHAOS = 0.01

# Two attractors count as the same one when both coordinates agree to within these relative
# tolerances. Loose enough to keep a single chaotic attractor together (its filling and exponent
# fluctuate from run to run), tight enough to separate the distinct branches. A heuristic.
TOLERANCE_LAMBDA = 0.05
TOLERANCE_NORM = 0.05

# Grid of BHMapDrivenDissipative.jl
minDelta, maxDelta, numDelta = 0.5, 3.0, 501
minF, maxF, numF = 0.5, 2.5, 401

PATH = f"c:/Users/micro/results/bh/driven/{L}/J_{J:.3f}_g_{G:.3f}_k_{KAPPA:.3f}/"

# STEM = f"driven_{L}"
STEM = f"driven_detail_{L}"

deltas = np.linspace(minDelta, maxDelta, numDelta)
fs = np.linspace(minF, maxF, numF)


def data_file(delta, f):
    """Name of the file written by BHMapDrivenDissipative.jl for this grid point."""
    return f"{PATH}{delta:.3f}_{f:.3f}.txt"


def count_attractors(lambdas, norms):
    """Number of distinct coexisting attractors among the sampled initial conditions.

    Greedy single-pass clustering in the (lambda, sum I) plane: a trajectory joins an existing
    cluster when both coordinates match its representative to within the relative tolerances
    above, otherwise it opens a new one."""
    clusters = []

    for lam, norm in zip(lambdas, norms):
        for clusterLambda, clusterNorm in clusters:
            closeLambda = abs(lam - clusterLambda) <= TOLERANCE_LAMBDA * max(1.0, abs(clusterLambda))
            closeNorm = abs(norm - clusterNorm) <= TOLERANCE_NORM * max(1.0, abs(clusterNorm))
            if closeLambda and closeNorm:
                break
        else:
            clusters.append((lam, norm))

    return len(clusters)


chaotic_fraction = np.full((numDelta, numF), np.nan)
lyapunov_chaotic = np.full((numDelta, numF), np.nan)
lyapunov_max = np.full((numDelta, numF), np.nan)
filling = np.full((numDelta, numF), np.nan)
attractors = np.full((numDelta, numF), np.nan)

for di, delta in enumerate(deltas):
    for fi, f in enumerate(fs):
        try:
            data = np.atleast_2d(np.loadtxt(data_file(delta, f)))

        except Exception as e:
            print(f"Error loading Delta = {delta:.3f}, f = {f:.3f}: {e}")
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

        chaotic_fraction[di, fi] = len(chaotic) / len(lambdas)
        lyapunov_max[di, fi] = np.max(lambdas)
        filling[di, fi] = np.mean(norms) / L
        attractors[di, fi] = count_attractors(lambdas, norms)

        # left as NaN (white) where no initial condition reached a strange attractor
        if len(chaotic) > 0:
            lyapunov_chaotic[di, fi] = np.mean(chaotic)


def fold_curves(delta):
    """The two saddle-node (fold) lines of the uniform branch, note sec. 4:

        n_pm = [2 Delta +- sqrt(Delta^2 - 3/4 kappa^2)] / (3 g),
        f^2  = n [(Delta - g n)^2 + kappa^2 / 4]

    They exist only for g Delta > 0 and |Delta| > (sqrt(3)/2) kappa, and meet at that cusp; the
    region between them is bistable. NaN outside."""
    lower = np.full(len(delta), np.nan)
    upper = np.full(len(delta), np.nan)

    for i, d in enumerate(delta):
        radicand = d * d - 0.75 * KAPPA ** 2
        if radicand <= 0 or G * d <= 0:
            continue

        root = np.sqrt(radicand)
        ns = [(2 * d + s * root) / (3 * G) for s in (1, -1)]
        values = sorted(np.sqrt(n * ((d - G * n) ** 2 + KAPPA ** 2 / 4)) for n in ns if n > 0)

        if len(values) == 2:
            lower[i], upper[i] = values

    return lower, upper


def modulational_growth_rate(delta, f):
    """Largest modulational-instability exponent at (Delta, f), note sec. 5:

        lambda_k = -kappa/2 + sqrt(g^2 n^2 - A_k^2),  A_k = eps_k - Delta + 2 g n,
        eps_k = 2 J (1 - cos k),  k = 2 pi m / L,

    with n running over the real positive roots of the state equation. MI means a k != 0 mode
    grows while k = 0 is stable, so only states with lambda_0 <= 0 enter, and for them only
    m >= 1; the middle branch, unstable at k = 0, is the saddle of the fold, not MI. Positive
    means at least one uniform state is modulationally unstable (a stable one may coexist), so the
    zero contour is the onset of pattern formation."""
    roots = np.roots([G ** 2, -2 * G * delta, delta ** 2 + KAPPA ** 2 / 4, -f ** 2])

    largest = np.nan
    for root in roots:
        if abs(root.imag) > 1e-9 or root.real <= 0:
            continue

        n = root.real
        rates = []
        for m in range(L):
            A = 2 * J * (1 - np.cos(2 * np.pi * m / L)) - delta + 2 * G * n
            radicand = G ** 2 * n ** 2 - A ** 2
            rates.append(-KAPPA / 2 + (np.sqrt(radicand) if radicand > 0 else 0.0))

        # k = 0 unstable: the saddle branch of the S-curve
        if rates[0] > 0:
            continue

        rate = max(rates[1:])
        if np.isnan(largest) or rate > largest:
            largest = rate

    return largest


growth = np.array([[modulational_growth_rate(d, f) for f in fs] for d in deltas])


def PlotAnalytics(ax=None, add_legend=True):
    """Overlay the closed-form structure of the note: the fold tongue of sec. 4 and the
    modulational-instability boundary of sec. 5."""
    if ax is None:
        ax = plt.gca()

    lower, upper = fold_curves(deltas)
    ax.plot(deltas, lower, color="white", lw=1.6, label="fold (saddle-node)")
    ax.plot(deltas, upper, color="white", lw=1.6)

    # zero contour of the analytic Bogoliubov rate = onset of modulational instability
    if np.any(np.isfinite(growth)):
        ax.contour(deltas, fs, np.transpose(growth), levels=[0.0],
                   colors="red", linewidths=1.4, linestyles="dashed")

    if add_legend:
        handles = [plt.Line2D([0], [0], color="white", lw=1.6, label="fold (bistability)"),
                   plt.Line2D([0], [0], color="red", lw=1.4, ls="dashed", label="MI onset")]
        ax.legend(handles=handles, loc="upper left", fontsize=8)


def PlotMap(values, label, title, name, cmap=CMAP, vmin=None, vmax=None):
    """Colormesh of `values` over the (Delta, f) plane with the analytic curves on top, saved as
    name_STEM.png / .pdf."""
    plt.figure(figsize=FIGSIZE)
    plt.pcolormesh(deltas, fs, np.transpose(values), cmap=cmap, shading="auto",
                   vmin=vmin, vmax=vmax)
    plt.colorbar(label=label)
    PlotAnalytics()
    plt.title(f"{title}, L = {L}, J = {J:g}, g = {G:g}, $\\kappa$ = {KAPPA:g}")
    plt.xlabel("$\\Delta$")
    plt.ylabel("f")
    plt.xlim(deltas[0], deltas[-1])
    plt.ylim(fs[0], fs[-1])
    plt.tight_layout()
    plt.savefig(f"{name}_{STEM}.png", dpi=150)
    plt.savefig(f"{name}_{STEM}.pdf")
    plt.show()

def PlotMapLyapunov(values, label, title, name):
    """Colormesh of `values` over the (Delta, f) plane with the analytic curves on top, saved as
    name_STEM.png / .pdf."""
    plt.figure(figsize=FIGSIZE)
    norm = TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=3.0)
    plt.pcolormesh(deltas, fs, np.transpose(values), cmap='seismic', shading="auto",
                   norm=norm)
    plt.colorbar(label=label)
    PlotAnalytics()
    plt.title(f"{title}, L = {L}, J = {J:g}, g = {G:g}, $\\kappa$ = {KAPPA:g}")
    plt.xlabel("$\\Delta$")
    plt.ylabel("f")
    plt.xlim(deltas[0], deltas[-1])
    plt.ylim(fs[0], fs[-1])
    plt.tight_layout()
    plt.savefig(f"{name}_{STEM}.png", dpi=150)
    plt.savefig(f"{name}_{STEM}.pdf")
    plt.show()

PlotMap(chaotic_fraction, "chaotic fraction", "Fraction of initial conditions reaching chaos",
        "chaotic_fraction", vmin=0.0, vmax=1.0)
PlotMap(lyapunov_chaotic, "$\\Lambda$", "Lyapunov exponent on the chaotic attractors",
        "lyapunov_chaotic")
PlotMapLyapunov(lyapunov_max, "$\\Lambda_{max}$", "Largest Lyapunov exponent found", "lyapunov_max")
PlotMap(filling, "n", "Filling per site on the attractor", "filling")
PlotMap(attractors, "attractors", "Number of coexisting attractors", "attractors", vmax=10)
