import numpy as np
from math import ceil
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
P = 5
L = 4    # cm
K = 1.68
H = 0.005
delta = 0.1  # cm

def A(n,m):
    A = np.zeros([n*m, n*m])
    h = L/(n-1)

    j = 0
    for i in range(1,m-1):
        k = i + j*m
        A[k, k]     = -3 + (2*h*H)/K
        A[k, k+m]   = 4
        A[k, k+2*m] = -1

    # efri hluti 
    j = n-1
    for i in range(1,m-1):
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

    # vinstra megin - efri helmingur 
    i = 0
    for j in range(ceil(n/2), n):
        k = i + j*m
        A[k, k]   = 3 - (2*h*H)/K
        A[k, k+1] = -4
        A[k, k+2] =  1

    # vinstra megin - neðri helmingur 
    i = 0
    for j in range(0, ceil(n/2)):
        k = i + j*m
        A[k, k]   = -3
        A[k, k+1] =  4
        A[k, k+2] = -1

    # innri nóður 
    for i in range(1,m-1):
        for j in range(1,n-1):
            k = i + j*m
            A[k, k]     = -4 - (2*H*h*h)/(K*delta)
            A[k, k-1]   = 1
            A[k, k+1]   = 1
            A[k, k-m]   = 1
            A[k, k+m]   = 1

    return A


def b(n,m):
    b = np.zeros([n*m])
    h = L/(n-1)

    # powerið (neðri helmingurinn vinstra megin)
    i = 0
    for j in range(0, ceil(n/2)):
        k = i + j*m
        b[k] = -2*h*P/(L * delta * K)

    return b

n = 100
m = 100
A_mat = A(n,m)
b_vec = b(n,m)

u = np.linalg.solve(A_mat, b_vec)
u_corrected = u + 20   

print(f"\nHámarks hitastig: {max(u_corrected):.4f} °C")


def plot_surface_3d(u, n, m):
    U = u.reshape((n, m))

    x = np.linspace(0, L, m)
    y = np.linspace(0, L, n)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, U, cmap='hot', edgecolor='none')

    ax.set_title("3D hitastigsmynd")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("Hitastig (°C)")

    fig.colorbar(surf, shrink=0.6, pad=0.1)
    plt.show()

plot_surface_3d(u_corrected, n, m)
