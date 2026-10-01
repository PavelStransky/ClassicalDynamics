"""One-dimensional cuts produced by BHMapNumberConservingChimera.jl (presets cut-eta3, cut-eta0,
cut-g20 and L4-eta3 of BHMapNumberConservingChimeraSubmit.sh).

The companion of analyse_map_number_conserving.py for a grid in which one of the two axes has a
single point.  The file format is the same (2L + 9 columns, one line per trajectory, layout in the
header of BHMapNumberConserving.jl); only the presentation differs: every trajectory is drawn as a
point against the swept parameter, so coexisting attractors, periodic windows and the onset of
chaos are visible directly rather than through a cell average.

Figures (saved as cut_<quantity>_<stem>.png/.pdf):

    lyapunov     reduced lambda_max of every trajectory, with the chaotic threshold
    dimension    reduced Kaplan-Yorke dimension of every trajectory
    classes      share of the initial conditions reaching each attractor type
    invariants   bond coherence, max n and current of every trajectory - a bifurcation diagram
                 built from time averages, complementary to the local maxima of
                 BHNumberConservingAttractors.jl bifurcation

and a table cut_<stem>.txt with one row per parameter value:

    x  chaotic_fraction  lambda_max  mean_lambda_chaotic  D_KY_chaotic  attractors  dominant_class

which is what A3 (where quantum chaos sets in relative to classical chaos) and A5 (route to chaos)
read.

The preset ray-J (scan J) is the experimental path: g, eta and kappa fixed, the hopping J swept on
a logarithmic axis.  Its exponents are stored in units of J; the extra figure lyapunov_physical
shows J lambda, the exponent in the time unit of the fixed rates, and the dashed lines mark the J
at which the fixed g crosses the modulational-instability threshold g_c(J).

Usage:  python analyse_cut_number_conserving.py <results directory> [--no-show]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

THRESHOLD_CHAOS = 1e-2
TOLERANCE = 0.08

CLASS_NAMES = {-1: "undetermined", 0: "fixed point", 1: "limit cycle", 2: "torus",
               3: "chaotic", 4: "hyperchaotic", 5: "neutral"}
CLASS_COLORS = {-1: "0.6", 0: "#4477aa", 1: "#66ccee", 2: "#228833", 3: "#ee6677",
                4: "#aa3377", 5: "#ccbb44"}


def read_metadata(path):
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


def count_attractors(data, L):
    """Same greedy clustering as analyse_map_number_conserving.py."""
    clusters = []
    for row in data:
        fingerprint = (row[2 * L + 4], row[2 * L + 5], row[2 * L + 6], row[2 * L + 8])
        for cluster in clusters:
            if cluster[0] == fingerprint[0] and \
                    max(abs(a - b) for a, b in zip(fingerprint[1:], cluster[1:])) <= TOLERANCE:
                break
        else:
            clusters.append(fingerprint)
    return len(clusters)


def main():
    arguments = [a for a in sys.argv[1:] if not a.startswith("--")]
    show = "--no-show" not in sys.argv
    if not arguments:
        print(__doc__)
        sys.exit(1)

    path = arguments[0]
    meta = read_metadata(path)
    L = meta["L"]
    gs, ys = meta["gValues"], meta["yValues"]

    ray = meta["scan"] == "J"
    if len(ys) == 1:
        xs, label, fixed = gs, "$g$", f"$\\eta$ = {ys[0]:g}"
        cells = [(g, ys[0]) for g in gs]
        marks = [meta["threshold"][0]]
    elif len(gs) == 1:
        label = {"eta": "$\\eta$", "kappa": "$\\kappa$",
                 "J": "$J$  ($g$, $\\eta$, $\\kappa$ fixed)"}.get(meta["scan"], meta["scan"])
        xs, fixed = ys, f"$g$ = {gs[0]:g}" + (f", $\\eta$ = {meta['eta']:g}" if ray else "")
        cells = [(gs[0], y) for y in ys]
        # where the fixed g crosses the MI threshold g_c(y) along the cut (twice along the J ray)
        excess = gs[0] - meta["threshold"]
        marks = [ys[i] - excess[i] * (ys[i + 1] - ys[i]) / (excess[i + 1] - excess[i])
                 for i in range(len(ys) - 1)
                 if np.isfinite(excess[i]) and np.isfinite(excess[i + 1]) and excess[i] * excess[i + 1] < 0]
    else:
        sys.exit("not a one-dimensional cut - use analyse_map_number_conserving.py")

    stem = os.path.basename(os.path.normpath(path)) + f"_L{L}_" + \
        ("g" if len(ys) == 1 else "y") + "cut"
    title = f"L = {L}, $\\kappa$ = {meta['kappa']:g}, {fixed}, jitter = {meta['jitter']:g}"

    points = {k: [] for k in ("x", "lambda", "dimension", "class", "coherence", "maxn", "current")}
    table = []

    for x, (g, y) in zip(xs, cells):
        name = os.path.join(path, f"{g:.4f}_{y:.4f}.txt")
        if not os.path.isfile(name):
            continue
        data = np.atleast_2d(np.loadtxt(name))
        data = data[np.isfinite(data[:, 2 * L])] if data.size else data
        if data.size == 0:
            continue

        lambdas = data[:, 2 * L]
        chaotic = lambdas > THRESHOLD_CHAOS
        classes = data[:, 2 * L + 4].astype(int)

        points["x"].extend([x] * len(data))
        points["lambda"].extend(lambdas)
        points["dimension"].extend(data[:, 2 * L + 1])
        points["class"].extend(classes)
        points["coherence"].extend(data[:, 2 * L + 5])
        points["maxn"].extend(data[:, 2 * L + 6])
        points["current"].extend(data[:, 2 * L + 8])

        values, counts = np.unique(classes, return_counts=True)
        table.append((x, np.mean(chaotic), np.max(lambdas),
                      np.mean(lambdas[chaotic]) if np.any(chaotic) else np.nan,
                      np.mean(data[chaotic, 2 * L + 1]) if np.any(chaotic) else np.nan,
                      count_attractors(data, L), values[np.argmax(counts)]))

    if not table:
        sys.exit(f"no result files in {path}")

    for key in points:
        points[key] = np.array(points[key])
    table = np.array(table)

    np.savetxt(f"cut_{stem}.txt", table, fmt="%.6g", delimiter="\t",
               header="x\tchaotic_fraction\tlambda_max\tmean_lambda_chaotic\tD_KY_chaotic\t"
                      "attractors\tdominant_class")
    print(f"{len(table)} of {len(xs)} parameter values present; table in cut_{stem}.txt")

    onset = table[table[:, 1] > 0.5, 0]
    if len(onset):
        print(f"chaotic fraction > 1/2 between x = {onset.min():.4g} and {onset.max():.4g}; "
              + (f"the uniform state is unstable for g < g_c = {marks[0]:.3f}" if len(ys) == 1 else
                 ("the fixed g crosses the MI threshold at x = " + ", ".join(f"{m:.4g}" for m in marks)
                  if marks else "the fixed g does not cross the MI threshold in this range")))

    def finish(name, ylabel):
        for k, m in enumerate(marks):
            if np.isfinite(m):
                plt.axvline(m, color="k", ls="--", lw=1, label="MI threshold" if k == 0 else None)
        if ray:
            plt.xscale("log")
        plt.xlabel(label)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend(loc="best", fontsize=8)
        plt.tight_layout()
        plt.savefig(f"cut_{name}_{stem}.png", dpi=150)
        plt.savefig(f"cut_{name}_{stem}.pdf")
        if show:
            plt.show()
        plt.close()

    plt.figure(figsize=(9, 5))
    plt.scatter(points["x"], points["lambda"], s=2, c="#4477aa", alpha=0.4, label="trajectories")
    plt.plot(table[:, 0], table[:, 2], color="#ee6677", lw=1, label="max over trajectories")
    plt.axhline(THRESHOLD_CHAOS, color="0.5", lw=0.8)
    finish("lyapunov", "$\\lambda_\\max$ (reduced)" + (", units of $J$" if ray else ""))

    if ray:
        # the exponents are stored in units of 1/J; in the time unit of the fixed rates they are λ J
        plt.figure(figsize=(9, 5))
        plt.scatter(points["x"], points["lambda"] * points["x"], s=2, c="#4477aa", alpha=0.4,
                    label="trajectories")
        plt.plot(table[:, 0], table[:, 2] * table[:, 0], color="#ee6677", lw=1,
                 label="max over trajectories")
        finish("lyapunov_physical", "$J\\,\\lambda_\\max$ (time unit of the fixed rates)")

    plt.figure(figsize=(9, 5))
    plt.scatter(points["x"], points["dimension"], s=2, c="#228833", alpha=0.4, label="trajectories")
    finish("dimension", f"$D_{{KY}}$ (reduced, of {2 * L - 2})")

    plt.figure(figsize=(9, 5))
    order = np.argsort(table[:, 0])
    xsorted = table[order, 0]
    bottom = np.zeros(len(xsorted))
    for code, name in CLASS_NAMES.items():
        share = []
        for x in xsorted:
            selection = points["class"][points["x"] == x]
            share.append(np.mean(selection == code) if len(selection) else 0.0)
        share = np.array(share)
        if np.any(share > 0):
            plt.fill_between(xsorted, bottom, bottom + share, step="mid",
                             color=CLASS_COLORS[code], label=name, lw=0)
            bottom += share
    plt.ylim(0, 1)
    finish("classes", "share of initial conditions")

    figure, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for ax, key, name in zip(axes, ("coherence", "maxn", "current"),
                             ("bond coherence", "max$_j$ $n_j$", "current")):
        ax.scatter(points["x"], points[key], s=2, c=[CLASS_COLORS.get(c, "k") for c in points["class"]],
                   alpha=0.5)
        ax.set_ylabel(name)
        for m in marks:
            if np.isfinite(m):
                ax.axvline(m, color="k", ls="--", lw=1)
        if ray:
            ax.set_xscale("log")
    axes[-1].set_xlabel(label)
    axes[0].set_title(title + " - time averages, coloured by attractor type")
    plt.tight_layout()
    plt.savefig(f"cut_invariants_{stem}.png", dpi=150)
    plt.savefig(f"cut_invariants_{stem}.pdf")
    if show:
        plt.show()
    plt.close()


if __name__ == "__main__":
    main()
