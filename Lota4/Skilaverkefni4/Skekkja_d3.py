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



def Tmax_n(n):
    m = n
    A = A_sparse(n, m)
    b = b_vigur(n, m)
    u = spsolve(A, b)        
    T = u + T_úti          
    return T.max()



if __name__ == "__main__":
    
    n_viðmið = 200      
    print(f"Reikna viðmiðunarlausn fyrir n = 200 = {n_viðmið}")
    T_max_viðmið = Tmax_n(n_viðmið)
    print(f"T_max (viðmið) ≈ {T_max_viðmið:.6f} °C\n")

   
    n_prufur = [10, 20, 30, 40, 50, 60, 70, 80]

    hs = []
    villa = []

    for n in n_prufur:
        print(f"Reikna fyrir n = m = {n} ...")
        T_max_n = Tmax_n(n)
        h = Lx / (n - 1)               
        skekkja = abs(T_max_n - T_max_viðmið) 

        hs.append(h)
        villa.append(skekkja)

        print(f"  h = {h:.5f}, T_max = {T_max_n:.6f}, skekkja = {skekkja:.6e}")

    hs = np.array(hs)
    villa = np.array(villa)

    
    plt.figure(figsize=(8,6))


plt.loglog(hs, villa, 'o-', linewidth=2, markersize=8,)


p, C = np.polyfit(np.log(hs), np.log(villa), 1)
villa_fit = np.exp(C) * hs**p


plt.loglog(hs, villa_fit, '--', linewidth=2,
           label=fr'Best-fit lína  (p ≈ {p:.2f})')

plt.grid(True, which="both", ls="--", alpha=0.6)

plt.xlabel(r'$h = hx =  hy$', fontsize=14)
plt.ylabel(r'Skekkja', fontsize=14)

plt.title(r'Skekkjumat', fontsize=16)


plt.xticks(hs, [f"{val:.3f}" for val in hs], fontsize=12, rotation=45)
plt.yticks(fontsize=12)

plt.legend(fontsize=13)
plt.tight_layout()
plt.show()

print(f"\nHallatala p ≈ {p:.4f}")

