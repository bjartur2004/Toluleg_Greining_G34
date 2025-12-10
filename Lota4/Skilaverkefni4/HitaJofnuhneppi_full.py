import numpy as np
from math import ceil

# Consts
P = 5
L = 2 #cm
K = 1.68
H = 0.005
delta = 0.1 #cm

Power_L = 2 #cm

def A(n,m):
    A = np.zeros([n*m, n*m])
    h = L/(n-1)


    # neðri hlið
    j = 0
    for i in range(1,m-1):
        k = i+j*m
        A[k, k] =     -3+(2*h*H)/K
        A[k, k+m] =   4
        A[k, k+2*m] = -1

    # efri hlið
    j = n-1
    for i in range(1,m-1):
        k = i+j*m
        A[k, k] =     3-(2*h*H)/K
        A[k, k-m] = -4
        A[k, k-2*m] = 1
    
    # hægri hlið
    i = m-1
    for j in range(n):
        k = i+j*m
        A[k,k] =     3-(2*h*H)/K
        A[k, k-1] = -4
        A[k, k-2] =   1

    # vinstri efri hlið
    i = 0
    for j in range(ceil(n*(Power_L/L)), n):
        k = i+j*m
        A[k,k] =     3-(2*h*H)/K
        A[k,k+1] = -4
        A[k,k+2] =   1

    # vinstri neðri hlið
    i = 0
    for j in range(0, ceil(n*(Power_L/L))):
        k = i+j*m
        A[k,k] =    -3
        A[k,k+1] =  4
        A[k,k+2] =  -1

    for i in range(1,m-1):
        for j in range(1,n-1):
            k = i+j*m
            A[k,k] = -4-(2*H*h**2)/(K*delta)
            A[k,k-1] = 1
            A[k,k+1] = 1
            A[k,k-m] = 1
            A[k,k+m] = 1
            
    return A


def b(n,m):
    b = np.zeros([n*m])
    h = L/(n-1)

    i = 0
    for j in range(ceil(n*(Power_L/L))):
        k = i+j*m
        b[k] = -2*h*P/(L*delta*K)
    return b

n = 100
m = 100
u = np.linalg.solve(A(n,m), b(n,m))
u_corrected = u + 20 

print("\n")
print(u_corrected)

print(f"u max: {max(u_corrected)}")

import matplotlib.pyplot as plt

def plot_heatmap(u, n, m):
    # Reshape into a 2D field (rows = n, cols = m)
    U = u.reshape((n, m))

    plt.figure(figsize=(5,4))
    plt.imshow(U, cmap='hot', origin='lower')
    plt.colorbar(label='Temperature (°C)')
    plt.title('Temperature Distribution')
    plt.xlabel('x direction')
    plt.ylabel('y direction')
    plt.show()

plot_heatmap(u_corrected,n,m)