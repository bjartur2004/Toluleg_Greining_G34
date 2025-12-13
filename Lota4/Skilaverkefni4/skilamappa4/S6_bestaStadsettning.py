import HitadreifingsPlot
import HitaJofnuhneppi_Sparce
import numpy as np
import matplotlib.pyplot as plt



def LeysaCase(case):
    global n, m
    Power_L = case[0]
    power_b = case[1]

    P = 5
    L = 4 #cm
    K = 1.68
    H = 0.005
    delta = 0.1 #cm
    ambient_tempeature = 20 #°C

    arguments = [P,L,K,H,delta,Power_L,ambient_tempeature]

    u = HitaJofnuhneppi_Sparce.solve_u(n,m,arguments, PowerAreaBottom=power_b)
    return u

n = 100
m = 100

lengths = [1]
poses   = [0, 1 , 1.5, 2, 3]
cases = []
for l in lengths:
    for p in poses:
        cases.append([l,p])

results = []
for case in cases:
    results.append(LeysaCase(case).reshape(n,m))

n_rows = len(lengths)
n_cols = len(poses)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))

for i, l in enumerate(lengths):
    for j, p in enumerate(poses):
        idx = i * n_cols + j
        ax = axes[i, j] if n_rows > 1 and n_cols > 1 else axes[max(i,j)]
        im = ax.imshow(results[idx], origin='lower', aspect='auto', cmap='hot')
        ax.set_title(f"Lengd={l}, Fjarlægð frá x-ás={p}")
        ax.axis('off')
        max_temp = np.max(results[idx])
        ax.text(0.5, -0.05, f"Hámarks hitastig = {max_temp:.2f} °C",
                transform=ax.transAxes, ha='center', va='top', fontsize=13)

fig.subplots_adjust(right=0.85, top=0.8)  

cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.8]) 
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("Hitastig [°C}")

fig.suptitle("Hitadreifing fyrir mismunandi lengdir og staðsetningar", fontsize=20)
plt.show()
