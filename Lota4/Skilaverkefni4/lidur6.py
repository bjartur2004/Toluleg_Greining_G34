import numpy as np
from math import ceil
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


P = 5
L = 2.0       # cm
K = 1.68
H = 0.005
delta = 0.1   # cm

def A(n, m):
    A = np.zeros((n*m, n*m))
    h = L / (n - 1)

    # neðri hluti
    j = 0
    for i in range(1, m-1):
        k = i + j*m
        A[k, k]     = -3 + (2*h*H)/K
        A[k, k+m]   = 4
        A[k, k+2*m] = -1

    # efri hluti
    j = n-1
    for i in range(1, m-1):
        k = i + j*m
        A[k, k]     =  3 - (2*h*H)/K
        A[k, k-m]   = -4
        A[k, k-2*m] =  1

    # hægri hluti
    i = m-1
    for j in range(n):
        k = i + j*m
        A[k, k]   = 3 - (2*h*H)/K
        A[k, k-1] = -4
        A[k, k-2] =  1

    # efri helmingur vinstra megin
    i = 0
    for j in range(ceil(n/2), n):
        k = i + j*m
        A[k, k]   = 3 - (2*h*H)/K
        A[k, k+1] = -4
        A[k, k+2] =  1

    # neðri helmingur vinstra megin 
    i = 0
    for j in range(0, ceil(n/2)):
        k = i + j*m
        A[k, k]   = -3
        A[k, k+1] =  4
        A[k, k+2] = -1

    # innri punktar 
    for i in range(1, m-1):
        for j in range(1, n-1):
            k = i + j*m
            A[k, k]     = -4 - (2*H*h*h)/(K*delta)
            A[k, k-1]   = 1
            A[k, k+1]   = 1
            A[k, k-m]   = 1
            A[k, k+m]   = 1

    return A


def solve_chip_location(n, m, y0_frac, y1_frac):
    h = L/(n-1)

    A_mat = A(n, m)

    b_vec = np.zeros(n*m)

    j0 = int(y0_frac * (n-1))
    j1 = int(y1_frac * (n-1))

    for j in range(j0, j1 + 1):
        k = 0 + j*m
        b_vec[k] = -2*h*P/(L * delta * K)

    u = np.linalg.solve(A_mat, b_vec)
    u = u + 20 

    return u

def plot3d(u, n, m, title, filename):
    U = u.reshape((n, m))

    x = np.linspace(0, L, m)
    y = np.linspace(0, L, n)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, U, cmap='hot', edgecolor='none')

    ax.set_title(title)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("Hitastig (°C)")

    fig.colorbar(surf, shrink=0.6)

    fig.savefig(filename, dpi=300)

    plt.close(fig)   
    plt.show()
    #plt.show()
# mismunandi stillingar
n = m = 80
configs = [
    (0.00, 0.25,"[0.00,0.25]"),
    (0.25, 0.50,"[0.25,0,50]"),
    (0.40, 0.60,"[0.40,0.60]"),
    (0.50, 0.75,"[0.50,0.75]"),
    (0.75, 1.00,"[0.75,1.00]")
]

results = []
for idx, (y0, y1, name) in enumerate(configs):
    u = solve_chip_location(n, m, y0, y1)
    Tmax = np.max(u)
    results.append((name, Tmax))

    print(f"{name:20s} → max T = {Tmax:.3f}°C")

    filename = f"heatmap_{idx}.png"   
    plot3d(u, n, m, name, filename)

best = min(results, key=lambda x: x[1])
print("\nBesta stilling til að lágmarka hámarks hitastig:")
print(f"{best[0]}  with  Tmax = {best[1]:.3f}°C")
