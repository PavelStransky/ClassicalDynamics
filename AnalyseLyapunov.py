import numpy as np
import matplotlib.pyplot as plt

# Read numbers from file
# Replace 'data.txt' with your filename

PATH = 'd:/results/bh/lyapunov/angel/'

THRESHOLD_STABLE = 0.0025
THRESHOLD_METASTABLE_FACTOR = 0.5

U = 1
J = 0.1

# Multiplying factor
F = 1

# NL maximum Lyapunov exponents
NL = 5
energies_ml = [[] for _ in range(NL)]
lyapunovs_ml = [[] for _ in range(NL)]

energies_mean = []
lyapunov_mean = []
lyapunov_var = []

freg_stable = []
freg = []

STATIONARY_POINTS = U * np.array([0.05, 0.25, 0.328333333333333, 0.41, 0.45, 0.48241092975590166, 0.5, 0.54, 0.61, 1.0100376908326862])

energies = U * np.linspace(0, 1, 201)

for energy in energies:
    try:
        data = np.loadtxt(f'{PATH}{J:.3f}_{U:.3f}_{energy:.3f}.txt')

    except Exception as e:
        print(f"Error loading {energy:.2f}: {e}")
        continue

    data = np.flip(np.sort(data))
    data = data[data > 0]

    if len(data) == 0:
        print(f"No valid data found for {energy:.2f}")
        continue

    stable = data[data < THRESHOLD_STABLE]
    unstable = data[data >= THRESHOLD_STABLE]

    threshold_mestastable = THRESHOLD_STABLE

    if len(unstable) > 5:
        for i in range(10):
            threshold_mestastable = np.std(unstable) * THRESHOLD_METASTABLE_FACTOR + THRESHOLD_STABLE
            unstable = data[data >= threshold_mestastable]

    metastable = data[data >= THRESHOLD_STABLE]
    metastable = metastable[metastable < threshold_mestastable]

    # Plot histogram
    plt.figure()
    counts, bins, patches = plt.hist(F * data, bins=30, range=(0, F * 0.15 * np.abs(U)) if F > 0 else (F * 0.15 * np.abs(U), 0), edgecolor='black')

    # Calculate mean and standard deviation
    mean = 0
    std = 0

    if len(unstable) > 1:
        mean = np.abs(F) * np.mean(unstable)
        std = np.abs(F) * np.std(unstable)

        y_pos = 1000  # Fixed y position for error bar

        # Add the point with vertical error bar
        plt.errorbar(mean, y_pos, xerr=std, fmt='o', color='red', capsize=5, label='Mean ± Std')

    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(f'Energy = {F * energy:.2f}, Stable = {len(stable)} ({len(data)})' + (f', L = {mean:.3f} ± {std:.3f}, L5 = {np.abs(F) * unstable[4]:.3f}' if len(unstable) > 5 else ''))
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f'{PATH}histogram_{J:.3f}_{U:.3f}_{energy:.3f}.png')
    plt.close()

    energies_mean.append(energy)
    freg.append(len(metastable) / len(data))
    freg_stable.append(len(stable) / len(data))

    if len(unstable) > 3:
        lyapunov_mean.append(np.mean(unstable))
        lyapunov_var.append(np.std(unstable))
    else:
        lyapunov_mean.append(0)
        lyapunov_var.append(0)

    for i in range(NL):
        if len(unstable) > i:
            energies_ml[i].append(energy)
            lyapunovs_ml[i].append(unstable[i])

plt.figure(figsize=(15, 10))
for i in range(NL):
    plt.plot(F * np.array(energies_ml[i]), np.abs(F) * np.array(lyapunovs_ml[i]), marker='o')

for sp in STATIONARY_POINTS:
    plt.axvline(x=F*sp, color='green', linestyle='--')

plt.xlabel('Energy')
plt.ylabel('Lyapunov Exponent')
plt.title('Lyapunov Exponents vs Energy')
plt.grid(True)
plt.savefig(f'{PATH}lyapunov_vs_energy_{J:.3f}_{U:.3f}.pdf')
plt.savefig(f'{PATH}lyapunov_vs_energy_{J:.3f}_{U:.3f}.png')
plt.show()

fig, ax1 = plt.subplots(figsize=(15, 10))

color = 'tab:blue'
ax1.plot(F * np.array(energies_mean), np.abs(F) * np.array(lyapunov_mean), marker='o', label='Mean', color=color)
ax1.fill_between(F * np.array(energies_mean), np.abs(F) * (np.array(lyapunov_mean) - np.array(lyapunov_var)), np.abs(F) * (np.array(lyapunov_mean) + np.array(lyapunov_var)), alpha=0.2, color=color)

ax1.set_xlabel('Energy')
ax1.set_ylabel('Lyapunov Exponent', color=color)

ax1.tick_params(axis="y", labelcolor=color)

color = 'tab:orange'
ax2 = ax1.twinx()
# ax2.plot(energies_mean, np.array(freg), marker='o', label='freg all', color=color)
ax2.plot(F * np.array(energies_mean), np.abs(F) * np.array(freg_stable), marker='s', label='freg', color=color)

for sp in STATIONARY_POINTS:
    ax1.axvline(x=F * sp, color='green', linestyle='--')

ax2.set_ylabel("freg", color=color)
ax2.tick_params(axis="y", labelcolor=color)

ax1.set_title('Regularity')
ax1.grid()
ax2.grid()
fig.tight_layout()

plt.show()