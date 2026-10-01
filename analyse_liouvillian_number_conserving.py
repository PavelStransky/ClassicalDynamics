"""Figures and fits for the Liouvillian tables written by BHNumberConservingLiouvillianMap.jl.

Items of number-conserving-BH-paper-TODO.md:

    plane   C1 (Fig. 6)  <|z|>, -<cos arg z> and var(s) of the clean sector over the (g, eta) plane
                         for every N in the table (N = 12 on the full grid, N = 14 on the thinned
                         grid of IPNP36), for the bulk and the slow windows w = 12 and 16, with the
                         modulational-instability threshold g_c(eta) and - if a classical map
                         directory is given - the contour where half of the classical initial
                         conditions are chaotic.
    cuts    C1, A3       the same statistics along eta = 3 (against g) and g = -20 (against eta)
                         for every N, and the CROSSOVER: for each N the g at which -<cos> passes the
                         midpoint between 2D Poisson (0) and Ginibre (0.2405) - the number A3 needs,
                         tabulated against g_c and the onset of classical chaos.
    gap     C6, K6       the gap against N at every point of the tables, and its extrapolation
                         G_inf with three fit forms (a + b/N + c/N^2, a + b N^-alpha,
                         a + b exp(-c N)) and three fit ranges (N >= 8, 10, 12); the spread of the
                         nine values is the uncertainty K6 asks for.
    eta0    C2           statistics of the parity-resolved sectors along eta = 0, and the fraction
                         of exactly real eigenvalues in each (the fingerprint of the antiunitary
                         symmetry, which A1 has to interpret).

Usage:
    python analyse_liouvillian_number_conserving.py plane [--classical <classical map dir>]
    python analyse_liouvillian_number_conserving.py cuts | gap | eta0
Options: --results <dir> (default ~/results/bh/number-conserving/quantum/3), --no-show
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

GINIBRE = {"r": 0.7378, "cos": 0.2405, "var": 0.0875}
POISSON = {"r": 2.0 / 3.0, "cos": 0.0, "var": 4.0 / np.pi - 1.0}
KAPPA = 0.3
LABELS = {"r": r"$\langle|z|\rangle$", "cos": r"$-\langle\cos\theta\rangle$", "var": "var(s)"}


def option(name, default=None):
    if name in sys.argv:
        return sys.argv[sys.argv.index(name) + 1]
    return default


RESULTS = option("--results", os.path.join(os.path.expanduser("~"), "results", "bh",
                                             "number-conserving", "quantum", "3"))
SHOW = "--no-show" not in sys.argv


def threshold(eta, kappa=KAPPA):
    """Section-5 threshold at L = 3, J = 1: the uniform state is unstable for g < g_c."""
    return -(9 + eta ** 2 / 12 + 4 * kappa ** 2) / 2


def load(task):
    data = np.genfromtxt(os.path.join(RESULTS, f"{task}.txt"), names=True)
    return np.atleast_1d(data)


def finish(name):
    plt.tight_layout()
    plt.savefig(f"{name}.png", dpi=150)
    plt.savefig(f"{name}.pdf")
    if SHOW:
        plt.show()
    plt.close()


def plane():
    everything = load("plane")
    classical = option("--classical")
    contour = None
    if classical:
        contour = classical_chaotic_fraction(classical)

    # one set of figures per N (N = 12 on the full grid; IPNP36 adds N = 14 on a thinned one)
    for N in np.unique(everything["N"]).astype(int):
        plane_at(everything[everything["N"] == N], N, contour)


def plane_at(data, N, contour):
    gs, etas = np.unique(data["g"]), np.unique(data["eta"])
    for window in ("bulk", "w12", "w16"):
        for statistic in ("r", "cos", "var"):
            grid = np.full((len(gs), len(etas)), np.nan)
            for row in data:
                grid[np.searchsorted(gs, row["g"]), np.searchsorted(etas, row["eta"])] = \
                    row[f"m1_{window}_{statistic}"]
            plt.figure(figsize=(9, 6))
            vmin, vmax = sorted((POISSON[statistic], GINIBRE[statistic]))
            plt.pcolormesh(gs, etas, grid.T, shading="auto", cmap="viridis",
                           vmin=vmin - 0.3 * (vmax - vmin), vmax=vmax + 0.3 * (vmax - vmin))
            bar = plt.colorbar(label=LABELS[statistic])
            bar.ax.axhline(GINIBRE[statistic], color="r", lw=2)
            bar.ax.axhline(POISSON[statistic], color="w", lw=2)
            plt.plot(threshold(etas), etas, "w-", lw=1.5, label="$g_c(\\eta)$")
            if contour is not None:
                plt.contour(contour[0], contour[1], contour[2].T, levels=[0.5], colors="r",
                            linewidths=1.2)
                plt.plot([], [], "r-", label="classical chaotic fraction 1/2")
            plt.xlabel("$g$")
            plt.ylabel("$\\eta$")
            plt.legend(loc="upper left", fontsize=8)
            plt.title(f"N = {N}, clean sector, {window}: {LABELS[statistic]} "
                      f"(Ginibre {GINIBRE[statistic]:.3f}, Poisson {POISSON[statistic]:.3f})")
            finish(f"quantum_plane_N{N}_{window}_{statistic}")


def classical_chaotic_fraction(path):
    """Chaotic fraction of a classical map directory (format of BHMapNumberConserving.jl)."""
    meta = {}
    with open(os.path.join(path, "parameters.txt")) as handle:
        for line in handle:
            if "\t" in line and not line.startswith("#"):
                key, value = line.rstrip("\n").split("\t", 1)
                meta[key] = value
    L = int(meta["L"])
    gs = np.array([float(v) for v in meta["gValues"].split(",")])
    ys = np.array([float(v) for v in meta["yValues"].split(",")])
    fraction = np.full((len(gs), len(ys)), np.nan)
    for i, g in enumerate(gs):
        for j, y in enumerate(ys):
            name = os.path.join(path, f"{g:.4f}_{y:.4f}.txt")
            if os.path.isfile(name):
                d = np.atleast_2d(np.loadtxt(name))
                d = d[np.isfinite(d[:, 2 * L])]
                if len(d):
                    fraction[i, j] = np.mean(d[:, 2 * L] > 1e-2)
    return gs, ys, fraction


def crossover(x, y, level):
    """First x (scanning from the regular end) where y crosses `level`, by linear interpolation."""
    order = np.argsort(-x)            # from g close to 0 towards more negative g
    x, y = x[order], y[order]
    for i in range(1, len(x)):
        if np.isfinite(y[i - 1]) and np.isfinite(y[i]) and (y[i - 1] - level) * (y[i] - level) <= 0:
            return x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])
    return np.nan


def cuts():
    for task, axis, xlabel in (("cut-eta3", "g", "$g$"), ("cut-g20", "eta", "$\\eta$")):
        data = load(task)
        Ns = np.unique(data["N"]).astype(int)
        figure, axes = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
        for column, window in enumerate(("bulk", "w12", "w16")):
            for row, statistic in enumerate(("r", "cos", "var")):
                ax = axes[row, column]
                for N in Ns:
                    select = data["N"] == N
                    order = np.argsort(data[axis][select])
                    ax.plot(data[axis][select][order], data[f"m1_{window}_{statistic}"][select][order],
                            "o-", ms=2.5, lw=1, label=f"N = {N}")
                ax.axhline(GINIBRE[statistic], color="r", ls="--", lw=1)
                ax.axhline(POISSON[statistic], color="k", ls=":", lw=1)
                if axis == "g":
                    ax.axvline(threshold(3.0), color="0.5", lw=1)
                ax.set_ylabel(LABELS[statistic])
                if row == 0:
                    ax.set_title(window)
                if row == 2:
                    ax.set_xlabel(xlabel)
        axes[0, 0].legend(fontsize=7)
        finish(f"quantum_{task}")

        if axis == "g":
            level = 0.5 * GINIBRE["cos"]
            print(f"\nA3: g at which -<cos> passes {level:.4f} (half way Poisson -> Ginibre); "
                  f"g_c = {threshold(3.0):.3f}")
            print("     N    bulk      w12      w16")
            for N in Ns:
                select = data["N"] == N
                values = [crossover(data["g"][select], data[f"m1_{w}_cos"][select], level)
                          for w in ("bulk", "w12", "w16")]
                print(f"  {N:4d}  " + "  ".join(f"{v:7.2f}" for v in values))


def fit_forms():
    return {
        "a + b/N + c/N^2": (lambda N, a, b, c: a + b / N + c / N ** 2, (0.5, 1.0, 0.0)),
        "a + b N^-alpha": (lambda N, a, b, alpha: a + b * N ** (-alpha), (0.5, 1.0, 1.0)),
        "a + b exp(-c N)": (lambda N, a, b, c: a + b * np.exp(-c * N), (0.5, 1.0, 0.1)),
    }


def gap():
    rows = []
    for task in ("reference", "cut-eta3", "cut-g20", "eta0", "kappa", "kappa0", "asymmetric"):
        path = os.path.join(RESULTS, f"{task}.txt")
        if os.path.isfile(path):
            rows.append(load(task))
    data = np.concatenate(rows)
    data = data[np.isfinite(data["gap"])]
    points = sorted({(g, e, k, s) for g, e, k, s in zip(data["g"], data["eta"], data["kappa"], data["Gsym"])})

    print("K6: extrapolated gap G_inf (fit form x fit range); spread = max - min of the nine")
    plt.figure(figsize=(8, 6))
    for g, eta, kappa, gsym in points:
        select = (data["g"] == g) & (data["eta"] == eta) & (data["kappa"] == kappa) & (data["Gsym"] == gsym)
        N, G = data["N"][select], data["gap"][select]
        order = np.argsort(N)
        N, G = N[order], G[order]
        if len(N) < 5:
            continue
        plt.plot(1 / N, G, "o-", ms=3, label=f"g={g:g}, η={eta:g}, κ={kappa:g}")
        estimates = []
        for name, (function, guess) in fit_forms().items():
            for start in (8, 10, 12):
                use = N >= start
                if use.sum() < 4:
                    continue
                try:
                    parameters, _ = curve_fit(function, N[use], G[use], p0=guess, maxfev=20000)
                    estimates.append((name, start, parameters[0]))
                except RuntimeError:
                    pass
        if estimates:
            values = np.array([e[2] for e in estimates])
            print(f"  (g, eta, kappa, Gsym) = ({g:g}, {eta:g}, {kappa:g}, {gsym:.3g}): "
                  f"G_inf = {np.median(values):.4f}, range {values.min():.4f} ... {values.max():.4f}")
            for name, start, value in estimates:
                print(f"      {name:18s} N >= {start:2d}: {value:.4f}")
    plt.xlabel("$1/N$")
    plt.ylabel("gap $\\mathcal{G}$")
    plt.xlim(left=0)
    plt.legend(fontsize=7)
    finish("quantum_gap")


def eta0():
    data = load("eta0")
    Ns = np.unique(data["N"]).astype(int)
    figure, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    for column, label in enumerate(("m1", "m0p", "m0m")):
        for N in Ns:
            select = data["N"] == N
            order = np.argsort(data["g"][select])
            axes[0, column].plot(data["g"][select][order], data[f"{label}_bulk_cos"][select][order],
                                 "o-", ms=2.5, lw=1, label=f"N = {N}")
            axes[1, column].plot(data["g"][select][order], data[f"{label}_real"][select][order],
                                 "o-", ms=2.5, lw=1)
        axes[0, column].axhline(GINIBRE["cos"], color="r", ls="--", lw=1)
        axes[0, column].axhline(0, color="k", ls=":", lw=1)
        axes[0, column].set_title({"m1": "q = 2π/3 (reflection∘† inside)", "m0p": "q = 0, parity +",
                                   "m0m": "q = 0, parity -"}[label])
        axes[1, column].set_xlabel("$g$")
    axes[0, 0].set_ylabel(LABELS["cos"] + " (bulk)")
    axes[1, 0].set_ylabel("fraction of real eigenvalues")
    axes[0, 0].legend(fontsize=7)
    finish("quantum_eta0")


if __name__ == "__main__":
    commands = {"plane": plane, "cuts": cuts, "gap": gap, "eta0": eta0}
    command = sys.argv[1] if len(sys.argv) > 1 else ""
    if command not in commands:
        print(__doc__)
        sys.exit(1)
    commands[command]()
