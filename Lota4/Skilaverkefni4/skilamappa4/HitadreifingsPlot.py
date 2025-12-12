import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def plotHitadreifingu(u, n, m,arguments, titill):
    P, L, K, H, delta, Power_L, ambient_tempeature = arguments

    U = u.reshape((n, m))

    x = np.linspace(0, L, m)
    y = np.linspace(0, L, n)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, U, cmap='hot', edgecolor='none')

    ax.set_title(titill)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("Hitastig (°C)")

    fig.colorbar(surf, shrink=0.6, pad=0.1)
    plt.show()

