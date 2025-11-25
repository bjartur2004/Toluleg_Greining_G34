import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
import math



c2 = 1
def C(c1):
    return np.array([[c1+c2, -c2],\
                     [-c2,   c2]])

M = np.array([[1,0],\
              [0,2]])

k1 = 10
k2 = 5
K = np.array([[k1+k2, -k2],\
              [-k2,   k2]])

# system matrix af kerfinu af dZ = AZ, Z = [x, dx]T 
def A(M,C,K):
    Minv = np.linalg.inv(M)

    A = np.block([
        [np.zeros((2,2)), np.eye(2)],
        [-Minv @ K,      -Minv @ C]
    ])

    #print(A)
    return A

def calc_damping(c1):
    # find system eigen values 
    eig = la.eigvals(A(M,C(c1),K))
    #print(f"\n{eig}")

    # sort eigenvalues by real part (largest to smallest)
    eig_sorted = sorted(eig, key=lambda z: z.real, reverse=True)

    # The oscillatory mode is the second one
    mode = eig_sorted[1]
    damping = mode.real
    return damping


a=-30
b=30
tol = 10**-3

dampings = []
c1s = np.linspace(0,200,1000)
for c1 in c1s:
    dampings.append(calc_damping(c1))

plt.plot(c1s, dampings)
plt.show()
