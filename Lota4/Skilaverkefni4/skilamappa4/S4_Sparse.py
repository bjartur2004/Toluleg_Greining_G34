import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


P = 5.0                 # W
Lx = 2.0                # cm 
Ly = 2.0                # cm 
K = 1.68                # W/(cm·°C)
H = 0.005               # W/(cm²·°C)
delta = 0.1             # cm 
T_úti = 20.0        # °C


def A_sparse(n, m):

    N = n * m
    A = lil_matrix((N, N))
    h = Ly / (n - 1)  

    # --- Neðri hlið  ---
    j = 0
    for i in range(1, m - 1):
        k = i + j * m
        A[k, k]       = -3.0 + (2.0 * h * H) / K
        A[k, k + m]   =  4.0
        A[k, k + 2*m] = -1.0

    # --- Efri hlið  ---
    j = n - 1
    for i in range(1, m - 1):
        k = i + j * m
        A[k, k]       =  3.0 - (2.0 * h * H) / K
        A[k, k - m]   = -4.0
        A[k, k - 2*m] =  1.0

    # --- Hægri hlið  ---
    i = m - 1
    for j in range(n):
        k = i + j * m
        A[k, k]     =  3.0 - (2.0 * h * H) / K
        A[k, k - 1] = -4.0
        A[k, k - 2] =  1.0

    # --- Vinstri hlið  ---
    i = 0
    for j in range(n):
        k = i + j * m
        A[k, k]     = -3.0
        A[k, k + 1] =  4.0
        A[k, k + 2] = -1.0

    # --- Innri punktar ---
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            k = i + j * m
            A[k, k]     = -4.0 - (2.0 * H * h**2) / (K * delta)
            A[k, k - 1] =  1.0
            A[k, k + 1] =  1.0
            A[k, k - m] =  1.0
            A[k, k + m] =  1.0

    return A.tocsr()


def b_vigur(n, m):

    N = n * m
    b = np.zeros(N)
    h = Ly / (n - 1)

    i = 0
    for j in range(n):
        k = i + j * m
        b[k] = -2.0 * h * P / (Ly * delta * K)

    return b


def solve_heatsink_sparse(n=200, m=200):
    A = A_sparse(n, m)
    b = b_vigur(n, m)

   
    u = spsolve(A, b)            
    T = u + T_úti             

   
    Tmax = T.max()
    kmax = T.argmax()
    i_max = kmax % m
    j_max = kmax // m

    x_gildi = np.linspace(0, Lx, m)
    y_gildi = np.linspace(0, Ly, n)
    x_max = x_gildi[i_max]
    y_max = y_gildi[j_max]

    print(f"Hæsta hitastig T_max ≈ {Tmax:.4f} °C")
    print(f"Gerist í punkti (i,j) = ({i_max+1}, {j_max+1}) ≈ (x,y) = ({x_max:.3f}, {y_max:.3f})")


    T_grid = T.reshape((n, m))

    return T_grid, x_gildi, y_gildi


def plot_mesh(T_grid, x_gildi, y_gildi, title="Hitadreifing (mesh)"):
    X, Y = np.meshgrid(x_gildi, y_gildi)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, T_grid, cmap='hot', edgecolor='k', linewidth=0.2)
    fig.colorbar(surf, ax=ax, label='Hitastig (°C)')
    ax.set_title(title)
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('y (cm)')
    ax.set_zlabel('T (°C)')
    plt.tight_layout()
    plt.show()




print("n = m = 200 ")
n = m = 200
T200, x200, y200 = solve_heatsink_sparse(n, m)


plot_mesh(T200, x200, y200, title="Hitadreifing fyrir m = n = 200 (sparse)")
