import numpy as np
import scipy as sp
from math import ceil

def A(n, m, p, arguments):

    P, h, K, H, delta, Power_Lx, Power_Ly, ambient_temperature = arguments
    entries = [[], [], []]

    def add(val, I, J):
        entries[0].append(val)
        entries[1].append(I)
        entries[2].append(J)

    def I(i, j, k):
        return i + j*m + k*m*n

    # Bottom face (z = 0) 
    # takes all its edges
    k = 0
    q_flux = P / (Power_Lx * Power_Ly)

    x0 = (n*h - Power_Lx) / 2
    x1 = (n*h + Power_Lx) / 2
    y0 = (m*h - Power_Ly) / 2
    y1 = (m*h + Power_Ly) / 2
    for j in range(n):
        y = j * h
        for i in range(m):
            x = i * h
            row = I(i, j, k)

            u1 = I(i, j, 1)
            u2 = I(i, j, 2)

            in_power = (x0 <= x <= x1) and (y0 <= y <= y1)

            if in_power:
                add( 3*K/(2*h), row, row)
                add(-2*K/h,      row, u1)
                add( K/(2*h),    row, u2)

            else:
                add( 3*K/(2*h) + H, row, row)
                add(-2*K/h,          row, u1)
                add( K/(2*h),        row, u2)

    # Top face (z = p-1) 
    # takes ALL its edges
    k = p-1
    for j in range(n):
        for i in range(m):
            row = I(i, j, k)
            u1 = I(i, j, p-2)
            u2 = I(i, j, p-3)

            add( 3*K/(2*h) + H, row, row)
            add(-2*K/h,          row, u1)
            add( K/(2*h),        row, u2)

    # Left face (x = 0)
    # skip z = 0 and z = p-1
    i = 0
    for k in range(1, p-1):
        for j in range(n):
            row = I(i, j, k)
            u1 = I(1, j, k)
            u2 = I(2, j, k)

            add( 3*K/(2*h) + H, row, row)
            add(-2*K/h,          row, u1)
            add( K/(2*h),        row, u2)

    # Right face (x = m-1) 
    # skip z = 0 and z = p-1
    i = m-1
    for k in range(1, p-1):
        for j in range(n):
            row = I(i, j, k)
            u1 = I(m-2, j, k)
            u2 = I(m-3, j, k)

            add( 3*K/(2*h) + H, row, row)
            add(-2*K/h,          row, u1)
            add( K/(2*h),        row, u2)

    # Front face (y = 0)
    # skip z = 0 and z = p-1
    j = 0
    for k in range(1, p-1):
        for i in range(1, m-1):
            row = I(i, j, k)
            u1 = I(i, 1, k)
            u2 = I(i, 2, k)

            add( 3*K/(2*h) + H, row, row)
            add(-2*K/h,          row, u1)
            add( K/(2*h),        row, u2)

    # Back face (y = n-1) 
    # skip z = 0 and z = p-1
    j = n-1
    for k in range(1, p-1):
        for i in range(1, m-1):
            row = I(i, j, k)
            u1 = I(i, n-2, k)
            u2 = I(i, n-3, k)

            add( 3*K/(2*h) + H, row, row)
            add(-2*K/h,          row, u1)
            add( K/(2*h),        row, u2)

    for k in range(1, p-1):
        for j in range(1, n-1):
            for i in range(1, m-1):
                row = I(i, j, k)

                add(-6*K/h**2, row, row)

                add( K/h**2, row, I(i+1, j,   k))
                add( K/h**2, row, I(i-1, j,   k))
                add( K/h**2, row, I(i,   j+1, k))
                add( K/h**2, row, I(i,   j-1, k))
                add( K/h**2, row, I(i,   j,   k+1))
                add( K/h**2, row, I(i,   j,   k-1))


    A = sp.sparse.csr_matrix((entries[0], (entries[1], entries[2])),shape=(n*m*p, n*m*p))

    return A


import numpy as np

def b(n, m, p, arguments):
    P, h, K, H, delta, Power_Lx, Power_Ly, ambient_temperature = arguments

    q_flux = P / (Power_Lx * Power_Ly)

    # power boundry
    x0 = (n*h - Power_Lx) / 2
    x1 = (n*h + Power_Lx) / 2
    y0 = (m*h - Power_Ly) / 2
    y1 = (m*h + Power_Ly) / 2

    def I(i, j, k):
        return i + j*m + k*m*n


    b = np.zeros(n * m * p)

    # Bottom face (z = 0)
    k = 0
    for j in range(n):
        y = j * h
        for i in range(m):
            x = i * h

            if x0 <= x <= x1 and y0 <= y <= y1:
                row = I(i, j, k)
                b[row] = q_flux
    return b




def solve_u(n, m, p, arguments):
    P, L, K, H, delta, Power_Lx, Power_Ly, ambient_temperature = arguments
    
    A_ = A(n, m, p, arguments)
    b_ = b(n, m, p, arguments)

    u = sp.sparse.linalg.spsolve(A_, b_)
    return u + ambient_temperature
