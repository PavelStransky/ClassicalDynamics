import numpy as np
import matplotlib.pyplot as plt

# Read numbers from file
# Replace 'data.txt' with your filename

PATH = 'd:/results/bh/lyapunov/angel/'

NL = 5
THRESHOLD = 0.005

U = 1
J = 0.1

# Multiplying factor
F = 1

energies = [[] for _ in range(NL)]
lyapunovs = [[] for _ in range(NL)]

# STATIONARY_POINTS = [-10.100376908326862, -6.1, -5, -5.4, -4.824109297559017, -4.5, -4.1, -3.283333333333334, -2.500000000000001, -0.5]
STATIONARY_POINTS = U * np.array([0.05, 0.25, 0.328333333333333, 0.41, 0.45, 0.48241092975590166, 0.5, 0.54, 0.61, 1.0100376908326862])
# STATIONARY_POINTS = [-sp for sp in STATIONARY_POINTS]

for energy in U * np.linspace(0, 1, 201):
    try:
        data = np.loadtxt(f'{PATH}{J:.3f}_{U:.3f}_{energy:.3f}.txt')

    except Exception as e:
        print(f"Error loading {energy:.2f}: {e}")
        continue

    stable = data[data < THRESHOLD]
    unstable = np.flip(np.sort(data[data >= THRESHOLD]))

    # Plot histogram
    plt.figure()
    counts, bins, patches = plt.hist(F * data, bins=30, range=(0, F * 0.15 * np.abs(U)), edgecolor='black')

    # Calculate mean and standard deviation
    if len(unstable) > 1:
        mean = F * np.mean(unstable)
        std = F * np.std(unstable)

        y_pos = 1000  # Fixed y position for error bar

        # Add the point with vertical error bar
        plt.errorbar(mean, y_pos, xerr=std, fmt='o', color='red', capsize=5, label='Mean ± Std')

    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(f'Energy = {F * energy:.2f}, Stable = {len(stable)} ({len(data)})' + (f', L = {F * np.mean(unstable):.3f} ± {F * np.std(unstable):.3f}, L5 = {F * unstable[4]:.3f}' if len(unstable) > 5 else ''))
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f'{PATH}histogram_{J:.3f}_{U:.3f}_{energy:.3f}.png')
    plt.close()

    for i in range(NL):
        if len(unstable) > i:
            energies[i].append(energy)
            lyapunovs[i].append(unstable[i])

plt.figure(figsize=(15, 10))
for i in range(NL):
    plt.plot(energies[i], lyapunovs[i], marker='o')

for sp in STATIONARY_POINTS:
    plt.axvline(x=sp, color='green', linestyle='--')

plt.xlabel('Energy')
plt.ylabel('Lyapunov Exponent')
plt.title('Lyapunov Exponents vs Energy')
plt.grid(True)
plt.savefig(f'{PATH}lyapunov_vs_energy_{J:.3f}_{U:.3f}.pdf')
plt.savefig(f'{PATH}lyapunov_vs_energy_{J:.3f}_{U:.3f}.png')
plt.show()