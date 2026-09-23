"""Heatmaps of the sweep produced by BHMapNumberConserving.jl.

The number-conserving counterpart of analyse_map_driven.py, for the model of
number-conserving-BH.md. Three things differ from the driven case:

  * the plane is (g, y) where y is whichever of eta, kappa or modulation the Julia script swept,
    and the grid, the constants and the section-5 threshold curve all come from parameters.txt in
    the results directory, so nothing has to be kept in sync by hand;
  * each file has 2L + 9 columns, one line per trajectory - the FULL Lyapunov spectrum followed by
    the reduced largest exponent, the Kaplan-Yorke dimension of the reduced 2L - 2 dimensional
    space, the mean divergence, the trace-rule residual, the attractor type and four invariant
    averages. The layout is documented at the top of BHMapNumberConserving.jl;
  * a strange attractor is identified by a positive exponent AND a fractional reduced D_KY, so the
    dimension map is the one that answers question 1 of section 8 - a positive exponent alone
    cannot tell a strange attractor from a torus that has not converged.

The maps are basin statistics: every cell is sampled with many initial conditions drawn uniformly
on the sphere sum_j n_j = 1 (the Fubini-Study measure of CP^(L-1)), so

    chaotic fraction   share of initial conditions reaching a strange attractor
    lambda (chaotic)   mean exponent over those - how strong the chaos is where it exists
    D_KY (chaotic)     mean reduced Kaplan-Yorke dimension over those, out of 2L - 2
    attractors         number of distinct coexisting attractors (question 3)
    locked fraction    share ending on the uniform phase-locked condensate, the dark state of the
                       kappa channel, which is what a Liouvillian metastable state would look like
    trace residual     log10 |sum lambda - <div F>|, the check of question 4 - this map should be
                       uniformly near -9 and is a diagnostic, not physics

The overlay is the modulational-instability threshold of section 5, g_c(y), read from
parameters.txt. Past it the uniform dark state is unstable; the point of the maps is that chaos
starts well beyond that line, not at it.
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, BoundaryNorm, ListedColormap

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

FIGSIZE = (9.2, 6.4)

CMAP = plt.get_cmap("viridis").copy()
CMAP.set_bad(color="white")

# Same value as CHAOS_THRESHOLD in the Julia script.
THRESHOLD_CHAOS = 1e-2

# Two trajectories count as the same attractor when all three invariant averages agree to within
# this absolute tolerance, widened by the time fluctuation of the observable on the attractor (a
# single strange attractor would otherwise split into as many clusters as there are trajectories).
# The same heuristic as CountAttractors in BHNumberConserving.jl; only the spreads are missing
# here, so the tolerance is a flat number and the count is a lower bound.
TOLERANCE = 0.08

PATH = sys.argv[1] if len(sys.argv) > 1 else \
    os.path.join(os.path.expanduser("~"), "results", "bh", "number-conserving", "3", "eta",
                 "J_1.000_k_0.300_e_3.000_m_0.000")

AXIS_LABELS = {"eta": "$\\eta$", "kappa": "$\\kappa$", "modulation": "rate modulation"}
CLASS_NAMES = {-1: "undetermined", 0: "fixed point", 1: "limit cycle", 2: "torus",
               3: "chaotic", 4: "hyperchaotic", 5: "neutral (no contraction)"}


def read_metadata(path):
    """Grid, constants and the section-5 threshold curve, as written by WriteMetadata()."""
    metadata = {}

    with open(os.path.join(path, "parameters.txt")) as handle:
        for line in handle:
            if line.startswith("#") or "\t" not in line:
                continue
            key, value = line.rstrip("\n").split("\t", 1)
            metadata[key] = value

    for key in ("L", "trajectories", "columns"):
        metadata[key] = int(metadata[key])
    for key in ("J", "kappa", "eta", "modulation", "chaosThreshold", "jitter"):
        metadata[key] = float(metadata[key])
    for key in ("gValues", "yValues", "threshold"):
        metadata[key] = np.array([float(v) for v in metadata[key].split(",")])

    return metadata


META = read_metadata(PATH)
L = META["L"]
SCAN = META["scan"]
gs = META["gValues"]
ys = META["yValues"]
STEM = f"nc_{SCAN}_{L}"

# Column layout of BHMapNumberConserving.jl
COLUMN_LAMBDA = 2 * L            # largest exponent of the reduced spectrum
COLUMN_DIMENSION = 2 * L + 1     # D_KY on the reduced 2L - 2 space
COLUMN_DIVERGENCE = 2 * L + 2
COLUMN_RESIDUAL = 2 * L + 3
COLUMN_CLASS = 2 * L + 4
COLUMN_COHERENCE = 2 * L + 5
COLUMN_MAXN = 2 * L + 6
COLUMN_IPR = 2 * L + 7
COLUMN_CURRENT = 2 * L + 8


def data_file(g, y):
    """Name of the file written by BHMapNumberConserving.jl for this grid point."""
    return os.path.join(PATH, f"{g:.4f}_{y:.4f}.txt")


def count_attractors(data):
    """Distinct coexisting attractors among the sampled initial conditions.

    Greedy single-pass clustering on (attractor type, bond coherence, max n, current) - all of
    them gauge- and translation-invariant, so the L states related by the Z_L symmetry of
    section 6 count as one attractor. The largest exponent is deliberately left out: its
    finite-time scatter on a strange attractor is several per cent and would split one attractor
    into many."""
    clusters = []

    for row in data:
        fingerprint = (row[COLUMN_CLASS], row[COLUMN_COHERENCE], row[COLUMN_MAXN],
                       row[COLUMN_CURRENT])
        for cluster in clusters:
            if cluster[0] == fingerprint[0] and \
                    max(abs(a - b) for a, b in zip(fingerprint[1:], cluster[1:])) <= TOLERANCE:
                break
        else:
            clusters.append(fingerprint)

    return len(clusters)


shape = (len(gs), len(ys))
chaotic_fraction = np.full(shape, np.nan)
lyapunov_chaotic = np.full(shape, np.nan)
lyapunov_max = np.full(shape, np.nan)
dimension_chaotic = np.full(shape, np.nan)
attractors = np.full(shape, np.nan)
locked_fraction = np.full(shape, np.nan)
divergence = np.full(shape, np.nan)
residual = np.full(shape, np.nan)
dominant_class = np.full(shape, np.nan)

for gi, g in enumerate(gs):
    for yi, y in enumerate(ys):
        try:
            data = np.atleast_2d(np.loadtxt(data_file(g, y)))
        except Exception as error:
            print(f"Error loading g = {g:.4f}, {SCAN} = {y:.4f}: {error}")
            continue

        if data.size == 0:
            continue

        # NaN marks a trajectory the integrator could not finish; a cell with nothing left is
        # missing data and stays white.
        data = data[np.isfinite(data[:, COLUMN_LAMBDA])]
        if len(data) == 0:
            continue

        lambdas = data[:, COLUMN_LAMBDA]
        chaotic = lambdas > THRESHOLD_CHAOS

        chaotic_fraction[gi, yi] = np.mean(chaotic)
        lyapunov_max[gi, yi] = np.max(lambdas)
        attractors[gi, yi] = count_attractors(data)
        divergence[gi, yi] = np.mean(data[:, COLUMN_DIVERGENCE])
        residual[gi, yi] = np.log10(np.max(np.abs(data[:, COLUMN_RESIDUAL])) + 1e-16)

        # the uniform phase-locked condensate: every bond fully coherent, every site at n = 1/L
        locked = (np.abs(data[:, COLUMN_COHERENCE] - 1) < 0.05) & \
                 (np.abs(data[:, COLUMN_MAXN] - 1 / L) < 0.05)
        locked_fraction[gi, yi] = np.mean(locked)

        classes, counts = np.unique(data[:, COLUMN_CLASS], return_counts=True)
        dominant_class[gi, yi] = classes[np.argmax(counts)]

        # left as NaN (white) where no initial condition reached a strange attractor
        if np.any(chaotic):
            lyapunov_chaotic[gi, yi] = np.mean(lambdas[chaotic])
            dimension_chaotic[gi, yi] = np.mean(data[chaotic, COLUMN_DIMENSION])


def PlotThreshold(ax=None, add_legend=True):
    """The modulational-instability threshold g_c(y) of section 5, read from parameters.txt: the
    uniform dark state is linearly unstable to the left of it (g < g_c). It is the only analytic
    structure the note has, and the maps show that chaos begins well past it."""
    if ax is None:
        ax = plt.gca()

    threshold = META["threshold"]
    if not np.any(np.isfinite(threshold)):
        return

    ax.plot(threshold, ys, color="white", lw=1.8, label="$\\S$5 threshold $g_c$")

    if add_legend:
        handles = [plt.Line2D([0], [0], color="white", lw=1.8,
                              label="MI threshold $g_c$ (uniform state)")]
        ax.legend(handles=handles, loc="upper left", fontsize=8)


def Title(text):
    return (f"{text}, L = {L}, J = {META['J']:g}, "
            + (f"$\\kappa$ = {META['kappa']:g}" if SCAN != "kappa" else "")
            + (f", $\\eta$ = {META['eta']:g}" if SCAN != "eta" else "")
            + (f", mod = {META['modulation']:g}" if SCAN != "modulation" else ""))


def PlotMap(values, label, title, name, cmap=CMAP, vmin=None, vmax=None, norm=None):
    """Colormesh of `values` over the (g, y) plane with the section-5 threshold on top, saved as
    name_STEM.png / .pdf."""
    plt.figure(figsize=FIGSIZE)
    plt.pcolormesh(gs, ys, np.transpose(values), cmap=cmap, shading="auto",
                   vmin=vmin, vmax=vmax, norm=norm)
    plt.colorbar(label=label)
    PlotThreshold()
    plt.title(Title(title))
    plt.xlabel("$g$")
    plt.ylabel(AXIS_LABELS.get(SCAN, SCAN))
    plt.xlim(gs[0], gs[-1])
    plt.ylim(ys[0], ys[-1])
    plt.tight_layout()
    plt.savefig(f"{name}_{STEM}.png", dpi=150)
    plt.savefig(f"{name}_{STEM}.pdf")
    plt.show()


def PlotClasses():
    """Which attractor type the majority of initial conditions reaches - the qualitative answer to
    question 1 of section 8, read directly rather than through a threshold on one exponent."""
    codes = [-1, 0, 1, 2, 3, 4, 5]
    colors = ["0.6", "#4477aa", "#66ccee", "#228833", "#ee6677", "#aa3377", "#ccbb44"]

    plt.figure(figsize=FIGSIZE)
    plt.pcolormesh(gs, ys, np.transpose(dominant_class), shading="auto",
                   cmap=ListedColormap(colors),
                   norm=BoundaryNorm([c - 0.5 for c in codes] + [codes[-1] + 0.5], len(codes)))
    colorbar = plt.colorbar(ticks=codes)
    colorbar.ax.set_yticklabels([CLASS_NAMES[c] for c in codes])
    PlotThreshold()
    plt.title(Title("dominant attractor type"))
    plt.xlabel("$g$")
    plt.ylabel(AXIS_LABELS.get(SCAN, SCAN))
    plt.xlim(gs[0], gs[-1])
    plt.ylim(ys[0], ys[-1])
    plt.tight_layout()
    plt.savefig(f"attractor_type_{STEM}.png", dpi=150)
    plt.savefig(f"attractor_type_{STEM}.pdf")
    plt.show()


PlotMap(chaotic_fraction, "chaotic fraction", "fraction of initial conditions reaching chaos",
        "chaotic_fraction", vmin=0, vmax=1)
PlotMap(lyapunov_max, "$\\lambda_\\max$", "largest exponent found", "lyapunov_max",
        cmap=plt.get_cmap("coolwarm"),
        norm=TwoSlopeNorm(vcenter=0.0, vmin=min(-1e-6, np.nanmin(lyapunov_max)),
                          vmax=max(1e-6, np.nanmax(lyapunov_max))))
PlotMap(lyapunov_chaotic, "$\\langle\\lambda_\\max\\rangle$ over chaotic trajectories",
        "strength of the chaos", "lyapunov_chaotic")
PlotMap(dimension_chaotic, f"$D_{{KY}}$ (reduced, out of {2 * L - 2})",
        "Kaplan-Yorke dimension of the strange attractor", "dimension",
        vmin=0, vmax=2 * L - 2)
PlotMap(attractors, "coexisting attractors", "multistability", "attractors", vmin=1)
PlotMap(locked_fraction, "basin fraction", "share ending on the uniform locked condensate",
        "locked_fraction", vmin=0, vmax=1)
PlotMap(divergence, "$\\langle$div $F\\rangle$", "mean phase-space contraction", "divergence",
        cmap=plt.get_cmap("coolwarm"),
        norm=TwoSlopeNorm(vcenter=0.0, vmin=min(-1e-6, np.nanmin(divergence)),
                          vmax=max(1e-6, np.nanmax(divergence))))
PlotMap(residual, "$\\log_{10}|\\Sigma\\lambda - \\langle$div$F\\rangle|$",
        "trace-rule residual (diagnostic)", "trace_residual")
PlotClasses()
