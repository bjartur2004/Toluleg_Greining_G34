import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


Lx = 4.0     # cm
Ly = 4.0     # cm
Lchip = 2.0  # cm – örgjörvinn er fast 2 cm að lengd
H = 0.005
delta = 0.1


def A_sparse(n, m, K):

    N = n * m
    A = lil_matrix((N, N))
    h = Ly / (n - 1)

    # Neðri hlið 
    j = 0
    for i in range(1, m - 1):
        k = i + j * m
        A[k, k]       = -3 + (2*h*H)/K
        A[k, k + m]   =  4
        A[k, k + 2*m] = -1

    # Efri hlið 
    j = n - 1
    for i in range(1, m - 1):
        k = i + j*m
        A[k, k]       =  3 - (2*h*H)/K
        A[k, k - m]   = -4
        A[k, k - 2*m] =  1

    # Hægri hlið 
    i = m - 1
    for j in range(n):
        k = i + j*m
        A[k,k]     =  3 - (2*h*H)/K
        A[k,k - 1] = -4
        A[k,k - 2] =  1

    # Vinstri hlið 
    i = 0
    for j in range(n):
        k = i + j*m
        A[k,k]     = -3
        A[k,k+1]   =  4
        A[k,k+2]   = -1

    # Innri punktar 
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            k = i + j*m
            A[k,k]     = -4 - (2*H*h*h)/(K*delta)
            A[k,k - 1] =  1
            A[k,k + 1] =  1
            A[k,k - m] =  1
            A[k,k + m] =  1

    return A.tocsr()


def solve_chip_with_params(n, m, y_center_frac, K, P):
    A_mat = A_sparse(n, m, K)
    b = np.zeros(n*m)
    h = Lx / (m - 1)


    chip_points = int((Lchip / Ly) * (n - 1))

    j_center = int(y_center_frac * (n - 1))
    j0 = max(0, j_center - chip_points//2)
    j1 = min(n - 1, j_center + chip_points//2)

    for j in range(j0, j1 + 1):
        k = 0 + j*m
        b[k] = -2 * h * P / (Lchip * delta * K)

    # Notum sparse leysara
    u = spsolve(A_mat, b)
    return u + 20


def find_max_P(n, m, K, y_center, Tmax_allowed=100):
    P_low = 0
    P_high = 40
    for _ in range(40):
        P_mid = 0.5*(P_low + P_high)
        u = solve_chip_with_params(n, m, y_center, K, P_mid)
        Tmax = np.max(u)

        if Tmax > Tmax_allowed:
            P_high = P_mid
        else:
            P_low = P_mid

    return P_low


if __name__ == "__main__":
    n = m = 80

    y_center_best = 0.50

    Ks = np.linspace(1, 5, 15)
    Pmax_list = []

    for Kval in Ks:
        Pmax = find_max_P(n, m, Kval, y_center_best)
        Pmax_list.append(Pmax)
        print(f"K = {Kval:.2f}   P_max = {Pmax:.3f} W")

    plt.figure(figsize=(8,5))
    plt.plot(Ks, Pmax_list, marker='o')
    plt.xlabel("K (W/cm°C)")
    plt.ylabel("Hámarks leyfilegt afl P_max (W)")
    plt.title("Leyfilegt hámarksafl sem fall af hitaleiðni K")
    plt.grid()
    plt.savefig("lidur8_plot.png", dpi=300)
    plt.show()
