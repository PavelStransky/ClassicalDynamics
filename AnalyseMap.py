import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

PATH = 'd:/results/bhx/4/'
PATH_SP = 'd:/results/bh/4/'

THRESHOLD_STABLE = 0.0025
THRESHOLD_METASTABLE_FACTOR = 0.5

U = 1

minJ = 0
maxJ = 0.5
numJ = 101

minE = 0
maxE = 1
numE = 201

Js = np.linspace(minJ, maxJ, numJ)
energies = U * np.linspace(minE, maxE, numE)

energies_mean = np.zeros((numJ, numE))
lyapunov_mean = np.zeros((numJ, numE))
lyapunov_var = np.zeros((numJ, numE))

freg_stable = np.zeros((numJ, numE))
freg = np.zeros((numJ, numE))

for ji, J in enumerate(Js):
    for ei, energy in enumerate(energies):
        try:
            data = np.loadtxt(f'{PATH}{J:.3f}_{U:.3f}_{energy:.3f}.txt')
            # data = np.loadtxt(f'{PATH}{J:.2f}_{energy:.2f}.txt')

        except Exception as e:
            print(f"Error loading {energy:.3f}: {e}")
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
plt.xlim(minJ, maxJ)
plt.ylim(minE, maxE)
plt.show()

plt.pcolormesh(Js, energies, np.transpose(lyapunov_var), cmap="viridis", shading="auto")
PlotESQPTs()
plt.colorbar(label="Lyapunov variance")
plt.title("Lyapunov variance")
plt.xlabel("J")
plt.ylabel("E")
plt.xlim(minJ, maxJ)
plt.ylim(minE, maxE)
plt.show()

plt.pcolormesh(Js, energies, np.transpose(freg_stable), cmap="viridis", shading="auto")
PlotESQPTs()
plt.colorbar(label="freg")
plt.title("Fraction of regularity")
plt.xlabel("J")
plt.ylabel("E")
plt.xlim(minJ, maxJ)
plt.ylim(minE, maxE)
plt.show()