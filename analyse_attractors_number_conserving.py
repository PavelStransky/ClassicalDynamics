"""Figures for the classical calculations of BHNumberConservingTransient.jl and
BHNumberConservingAttractors.jl (items C13 - C17 of number-conserving-BH-paper-TODO.md).

    transient      C13  t_crit/tau (mean and spread) against g for every distance d, sigma_eps
                        against eps (linear = regular transient, flat = chaotic transient), and
                        the mean FTLE over the transient, for the tasks eta0 and eta3
    bifurcation    C15  local maxima of the bond coherence and of n_1 against g (eta = 3), random
                        initial conditions in grey, the two continuation sweeps in colour
    uncertainty    C16  f(eps) on log-log axes with the fitted exponent
    simplex        C17  the attractors on the population simplex (ternary coordinates), the visit
                        density of the strange attractor and the Poincare section
    correlations   C14  |C(s)| of the two observables on a log scale, with the fitted decay

Usage:  python analyse_attractors_number_conserving.py <command> [--results <dir>] [--no-show]
Default results: ~/results/bh/number-conserving/classical/3/{transient,attractors}
"""

import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt


def option(name, default=None):
    if name in sys.argv:
        return sys.argv[sys.argv.index(name) + 1]
    return default


BASE = os.path.join(os.path.expanduser("~"), "results", "bh", "number-conserving", "classical", "3")
TRANSIENT = option("--results", os.path.join(BASE, "transient"))
ATTRACTORS = option("--results", os.path.join(BASE, "attractors"))
SHOW = "--no-show" not in sys.argv


def finish(name):
    plt.tight_layout()
    plt.savefig(f"{name}.png", dpi=150)
    plt.savefig(f"{name}.pdf")
    if SHOW:
        plt.show()
    plt.close()


def ternary(n1, n2, n3):
    """Ternary coordinates of the population simplex: corners n1 = 1 at (0, 0), n2 = 1 at (1, 0),
    n3 = 1 at (1/2, sqrt(3)/2)."""
    return n2 + 0.5 * n3, np.sqrt(3) / 2 * n3


def simplex_frame(ax):
    ax.plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3) / 2, 0], "k-", lw=0.8)
    ax.plot(*ternary(1 / 3, 1 / 3, 1 / 3), "k+", ms=8)
    ax.text(-0.04, -0.04, "$n_1$")
    ax.text(1.0, -0.04, "$n_2$")
    ax.text(0.48, 0.89, "$n_3$")
    ax.set_aspect("equal")
    ax.axis("off")


def transient():
    for task in ("eta0", "eta3"):
        path = os.path.join(TRANSIENT, f"summary_{task}.txt")
        if not os.path.isfile(path):
            print("missing", path)
            continue
        data = np.atleast_1d(np.genfromtxt(path, names=True))
        names = data.dtype.names
        distances = sorted({n.split("_")[0] for n in names if n.startswith("d") and n.endswith("_mean")})
        epsilons = sorted({n.split("_")[0] for n in names if n.startswith("eps")},
                          key=lambda e: float(e[3:]))
        g = data["g"]

        figure, axes = plt.subplots(1, 3, figsize=(16, 4.8))
        for d in distances:
            axes[0].errorbar(g, data[f"{d}_mean"], yerr=data[f"{d}_sigma"], fmt="o-", ms=3, capsize=2,
                             label=d.replace("d", "d = "))
            axes[2].plot(g, data[f"{d}_ftle"], "o-", ms=3, label=d.replace("d", "d = "))
        axes[0].set_ylabel("$t_{crit}/\\tau$ (mean ± σ_d)")
        axes[2].axhline(0, color="k", lw=0.8)
        axes[2].set_ylabel("mean FTLE over the transient")
        for e in epsilons:
            axes[1].semilogy(g, data[f"{e}_sigma"], "o-", ms=3, label=f"ε = {float(e[3:]):.0e}")
        axes[1].set_ylabel("$\\sigma_\\varepsilon$ of $t_{crit}/\\tau$")
        for ax in axes:
            ax.set_xlabel("$g$")
            ax.legend(fontsize=7)
        axes[0].set_title(f"{task}: relaxation time over the linear time τ")
        axes[1].set_title("sensitivity to the initial condition")
        finish(f"transient_{task}")


def bifurcation():
    for path in glob.glob(os.path.join(ATTRACTORS, "bifurcation_e*.txt")):
        data = np.loadtxt(path)
        figure, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        for ax, observable, label in zip(axes, (1, 2), ("maxima of $C$", "maxima of $n_1$")):
            select = data[:, 2] == observable
            random = select & (data[:, 1] > 0)
            ax.plot(data[random, 0], data[random, 3], ",", color="0.4", alpha=0.6)
            for source, color, name in ((-1, "C0", "continuation, g decreasing"),
                                        (-2, "C3", "continuation, g increasing")):
                sweep = select & (data[:, 1] == source)
                ax.plot(data[sweep, 0], data[sweep, 3], ".", ms=1.5, color=color, label=name)
            ax.set_ylabel(label)
        axes[0].legend(fontsize=8, markerscale=6)
        axes[1].set_xlabel("$g$")
        axes[0].set_title(os.path.basename(path))
        finish(os.path.splitext(os.path.basename(path))[0])


def uncertainty():
    for path in glob.glob(os.path.join(ATTRACTORS, "uncertainty_*.txt")):
        data = np.loadtxt(path)
        use = data[:, 1] > 0
        slope, intercept = np.polyfit(np.log10(data[use, 0]), np.log10(data[use, 1]), 1)
        plt.figure(figsize=(6, 5))
        plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], fmt="o", capsize=2)
        x = np.logspace(np.log10(data[:, 0].min()), np.log10(data[:, 0].max()), 50)
        plt.loglog(x, 10 ** intercept * x ** slope, "r-", label=f"α = {slope:.3f}")
        plt.xlabel("ε (Fubini–Study)")
        plt.ylabel("f(ε), uncertain fraction")
        plt.legend()
        plt.title(os.path.basename(path))
        finish(os.path.splitext(os.path.basename(path))[0])


def simplex():
    directory = os.path.join(ATTRACTORS, "simplex")
    figure, axes = plt.subplots(2, 2, figsize=(12, 11))

    for ax, pattern, title in ((axes[0, 0], "attractor_g-6.00_e3.00_*.txt", "g = -6, η = 3: limit cycles"),
                               (axes[0, 1], "attractor_g-20.00_e1.00_*.txt", "g = -20, η = 1: coexisting attractors")):
        simplex_frame(ax)
        for k, path in enumerate(sorted(glob.glob(os.path.join(directory, pattern)))):
            d = np.loadtxt(path)
            x, y = ternary(d[:, 1], d[:, 2], d[:, 3])
            if np.ptp(x) + np.ptp(y) < 1e-4:        # a fixed point: a marker, not a line
                ax.plot(x[-1], y[-1], "o", ms=8, label=f"attractor {k + 1} (fixed point)")
            else:
                ax.plot(x, y, lw=0.8, label=f"attractor {k + 1}")
        ax.legend(fontsize=7, loc="upper right")
        ax.set_title(title)

    ax = axes[1, 0]
    simplex_frame(ax)
    density = os.path.join(directory, "density_g-20.00_e3.00.txt")
    if os.path.isfile(density):
        d = np.loadtxt(density)
        bins = int(d[:, 0].max())
        grid = np.zeros((bins, bins))
        grid[d[:, 0].astype(int) - 1, d[:, 1].astype(int) - 1] = d[:, 2]
        grid[grid == 0] = np.nan
        ax.imshow(np.log10(grid.T), origin="lower", extent=(0, 1, 0, np.sqrt(3) / 2), cmap="magma",
                  aspect="equal")
    ax.set_title("g = -20, η = 3: visit density (log) of the strange attractor")

    ax = axes[1, 1]
    section = os.path.join(directory, "section_g-20.00_e3.00.txt")
    if os.path.isfile(section):
        d = np.loadtxt(section)
        ax.plot(d[:, 1], d[:, 0], ",", color="k", alpha=0.5)
        ax.set_xlabel("bond coherence $C$")
        ax.set_ylabel("$n_3$")
    ax.set_title("Poincaré section $n_1 = n_2$ (upward), g = -20, η = 3")
    finish("attractors_simplex")


def correlations():
    for path in glob.glob(os.path.join(ATTRACTORS, "correlations_*.txt")):
        d = np.loadtxt(path)
        plt.figure(figsize=(8, 5))
        plt.semilogy(d[:, 0], np.hypot(d[:, 1], d[:, 2]), label=r"$|C_{\tilde n_q}(s)|$, $q = 2\pi/3$")
        plt.semilogy(d[:, 0], np.abs(d[:, 3]), label=r"$|C_{C}(s)|$ (bond coherence)")
        plt.axhline(0.02, color="0.5", lw=0.8)
        plt.xlabel("lag $s$")
        plt.ylabel("autocorrelation")
        plt.legend()
        plt.title(os.path.basename(path) + "  (compare the slopes with the gap)")
        finish(os.path.splitext(os.path.basename(path))[0])


if __name__ == "__main__":
    commands = {"transient": transient, "bifurcation": bifurcation, "uncertainty": uncertainty,
                "simplex": simplex, "correlations": correlations}
    command = sys.argv[1] if len(sys.argv) > 1 else ""
    if command not in commands:
        print(__doc__)
        sys.exit(1)
    commands[command]()
