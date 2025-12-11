import numpy as np
from math import ceil
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Fastar
P = 5
Lx = 4.0      # breidd
Ly = 4.0      # hæð
Lchip = 2.0   # lengd örgjörva (2 cm FAST)
K = 1.68
H = 0.005
delta = 0.1


def A(n, m):
    A = np.zeros((n*m, n*m))
    h = Lx / (m - 1)     

    # neðri hluti
    j = 0
    for i in range(1, m-1):
        k = i + j*m
        A[k,k]     = -3 + (2*h*H)/K
        A[k,k+m]   = 4
        A[k,k+2*m] = -1

    # efri hluti
    j = n-1
    for i in range(1, m-1):
        k = i + j*m
        A[k,k]     =  3 - (2*h*H)/K
        A[k,k-m]   = -4
        A[k,k-2*m] =  1

    # hægri hluti
    i = m-1
    for j in range(n):
        k = i + j*m
        A[k,k]   = 3 - (2*h*H)/K
        A[k,k-1] = -4
        A[k,k-2] =  1

    # vinstra megin 
    i = 0
    for j in range(n):
        k = i + j*m
        A[k,k]   = -3
        A[k,k+1] =  4
        A[k,k+2] = -1

    # innri punktar
    for i in range(1, m-1):
        for j in range(1, n-1):
            k = i + j*m
            A[k,k]   = -4 - (2*H*h*h)/(K*delta)
            A[k,k-1] = 1
            A[k,k+1] = 1
            A[k,k-m] = 1
            A[k,k+m] = 1

    return A


def solve_chip_location(n, m, y_center_frac):
    A_mat = A(n, m)
    b = np.zeros(n*m)

    h = Lx / (m - 1)

    chip_points = int((Lchip / Ly) * (n - 1))

    # Miðja örgjörvans í j-hnitum
    j_center = int(y_center_frac * (n - 1))


    j0 = max(0, j_center - chip_points // 2)
    j1 = min(n - 1, j_center + chip_points // 2)

    for j in range(j0, j1+1):
        k = 0 + j*m
        b[k] = -2 * h * P / (Lchip * delta * K)

    # Leysa jöfnuhneppið
    u = np.linalg.solve(A_mat, b)
    return u + 20   # bæta umhverfishita við



def plot3d(u, n, m, title, filename):
    U = u.reshape((n, m))
    x = np.linspace(0, Lx, m)
    y = np.linspace(0, Ly, n)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, U, cmap='hot', edgecolor='none')
    ax.set_title(title)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("Hitastig (°C)")

    fig.colorbar(surf, shrink=0.6)
    fig.savefig(filename, dpi=300)

    plt.close(fig)



if __name__ == "__main__":
    n = m = 80

    # Mismunandi y-miðjur örgjörvans (0 → neðst, 1 → efst)
    centers = [0.10, 0.30, 0.50, 0.70, 0.90]

    results = []

    for idx, yc in enumerate(centers):
        u = solve_chip_location(n, m, yc)
        Tmax = np.max(u)
        results.append((yc, Tmax))

        print(f"Miðja = {yc:.2f}  →  T_max = {Tmax:.3f} °C")
        plot3d(u, n, m, f"Liður 6 – miðja {yc:.2f}", f"lidur6_{idx}.png")

    # Finnum bestu staðsetninguna
    best = min(results, key=lambda x: x[1])
    print("\nBesta staðsetningin")
    print(f"y_center = {best[0]:.2f}  →  T_max = {best[1]:.3f} °C")
