import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

PATH = 'd:/results/bhx/3/'
PATH_SP = 'd:/results/bh/3/'

THRESHOLD_STABLE = 0.0025
THRESHOLD_METASTABLE_FACTOR = 0.5

U = 1
Js = np.linspace(-0.5, 0.5, 101)
energies = U * np.linspace(0, 1, 101)

energies_mean = np.zeros((101, 101))
lyapunov_mean = np.zeros((101, 101))
lyapunov_var = np.zeros((101, 101))

freg_stable = np.zeros((101, 101))
freg = np.zeros((101, 101))

for ji, J in enumerate(Js):
    for ei, energy in enumerate(energies):
        try:
            data = np.loadtxt(f'{PATH}{J:.2f}_{energy:.2f}.txt')

        except Exception as e:
            print(f"Error loading {energy:.2f}: {e}")
            continue

        data = np.flip(np.sort(data))

        stable = data[data < THRESHOLD_STABLE]
        unstable = data[data >= THRESHOLD_STABLE]

        threshold_mestastable = THRESHOLD_STABLE

        if len(unstable) > 5:
            for i in range(10):
                threshold_mestastable = np.std(unstable) * THRESHOLD_METASTABLE_FACTOR + THRESHOLD_STABLE
                unstable = data[data >= threshold_mestastable]

        metastable = data[data >= THRESHOLD_STABLE]
        metastable = metastable[metastable < threshold_mestastable]

        # Calculate mean and standard deviation
        mean = 0
        std = 0

        freg[ji, ei] = len(metastable) / len(data)
        freg_stable[ji, ei] = len(stable) / len(data)

        if len(unstable) > 3:
            lyapunov_mean[ji, ei] = np.mean(unstable)
            lyapunov_var[ji, ei] = np.std(unstable)

def PlotESQPTs():
    for type in ['hsaddle', 'hstable', 'hunstable']:
        for i in range(9):
            fname = f'{PATH_SP}{type}_{i}.txt'

            data = []
            try:
                data = np.loadtxt(fname, delimiter=',')

            except Exception as e:
                print(f"Error loading {fname}.")
                continue

            if len(data) == 0:
                print(f"File {fname} contains no data.")
                continue

            if len(data.shape) == 1:
                print(f"File {fname} contains only one row.")
                continue

            plt.scatter(data[:, 0], data[:, 1], color='white', s=5)


plt.pcolormesh(Js, energies, np.transpose(lyapunov_mean), cmap="viridis", shading="auto")
PlotESQPTs()
plt.colorbar(label="Lyapunov")
plt.title("Lyapunov exponent")
plt.xlabel("J")
plt.ylabel("E")
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.show()

plt.pcolormesh(Js, energies, np.transpose(lyapunov_var), cmap="viridis", shading="auto")
PlotESQPTs()
plt.colorbar(label="Lyapunov variance")
plt.title("Lyapunov variance")
plt.xlabel("J")
plt.ylabel("E")
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.show()

plt.pcolormesh(Js, energies, np.transpose(freg_stable), cmap="viridis", shading="auto")
PlotESQPTs()
plt.colorbar(label="freg")
plt.title("Fraction of regularity")
plt.xlabel("J")
plt.ylabel("E")
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.show()