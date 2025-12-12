import numpy as np
import scipy as sp
from math import ceil

def A(n,m,arguments, PowerAreaBottom = 0):
    P, L, K, H, delta, Power_L, ambient_tempeature = arguments
    h = L/(n-1)

    entries = [[],[],[]] # shape: [[v0,v1,v2], [I0,I1,I2], [J0,J1,J2]]
    def addEntry(value,I,J):
        entries[0].append(value)
        entries[1].append(I)
        entries[2].append(J)

    # neðri hlið
    j = 0
    for i in range(1, m-1):
        k = i + j*m
        addEntry(-3 + (2*h*H)/K, k, k)
        addEntry(4,              k, k+m)
        addEntry(-1,             k, k+2*m)


    # efri hlið
    j = n-1
    for i in range(1, m-1):
        k = i + j*m
        addEntry( 3 - (2*h*H)/K, k, k)
        addEntry(-4,             k, k-m)
        addEntry( 1,             k, k-2*m)

    
    # hægri hlið
    i = m-1
    for j in range(n):
        k = i + j*m
        addEntry( 3 - (2*h*H)/K, k, k)
        addEntry(-4,             k, k-1)
        addEntry( 1,             k, k-2)




    start_j = int(np.floor(PowerAreaBottom / h))          
    end_j   = int(np.ceil((PowerAreaBottom + Power_L) / h))  

    # vinstri efri hlið
    i = 0
    for j in range(end_j, n):
        k = i + j*m
        addEntry( 3 - (2*h*H)/K, k, k)
        addEntry(-4,             k, k+1)
        addEntry( 1,             k, k+2)

    i = 0
    for j in range(start_j, end_j):
        k = i + j*m
        addEntry(-3, k, k)
        addEntry( 4, k, k+1)
        addEntry(-1, k, k+2)

    # fylla í botnin ef power svæðið byrjar ekki á botninum
    i = 0
    for j in range(0, start_j):
        k = i + j*m
        addEntry( 3 - (2*h*H)/K, k, k)
        addEntry(-4,             k, k+1)
        addEntry( 1,             k, k+2)


    # innri stök
    for i in range(1, m-1):
        for j in range(1, n-1):
            k = i + j*m
            addEntry(-4 - (2*H*h**2)/(K*delta), k, k)
            addEntry(1, k, k-1)
            addEntry(1, k, k+1)
            addEntry(1, k, k-m)
            addEntry(1, k, k+m)
    
    A = sp.sparse.csr_matrix((entries[0], (entries[1], entries[2])), shape=(n*m, n*m))
    return A


def b(n,m, arguments, PowerAreaBottom=0):
    P, L, K, H, delta, Power_L, ambient_tempeature = arguments

    b = np.zeros([n*m])
    h = L/(n-1)

    start_j = int(np.floor(n * PowerAreaBottom / L))  # start from nearest lower index
    end_j   = int(np.ceil(n * (PowerAreaBottom + Power_L) / L))  # cover entire fractional length

    i = 0  # left column
    for j in range(start_j, end_j):
        k = i + j*m
        b[k] = -2 * h * P / (L * delta * K)

    return b



def solve_u(n,m,arguments, PowerAreaBottom = 0): # argumens should be [P,L,K,H,Delta,Power_L, ambient_tempeature]
    P, L, K, H, delta, Power_L, ambient_tempeature = arguments
    u = sp.sparse.linalg.spsolve(A(n,m, arguments, PowerAreaBottom=PowerAreaBottom), b(n,m, arguments, PowerAreaBottom=PowerAreaBottom))
    u_corrected = np.array(u) + 20 

    return u_corrected