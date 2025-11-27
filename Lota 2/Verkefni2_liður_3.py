import numpy as np
import matplotlib.pyplot as plt

QB = 7.0
nu = 1e-3
r = 0.05
L = 100
G = (np.pi * r**4) / (8 * nu * L)

# A, B, C úr lið 2:
A = 1.3209447623073836       
B = 4.015670949654586  
C = 9.45       

omega = 2*np.pi/24 

A_fylki = np.array([
    [ 3, -1, -1,  0,  0],
    [ 1, -2,  0,  1,  0],
    [ 1,  0, -4,  1,  2],
    [ 0,  3,  3, -8,  2],
    [ 0,  0,  6,  2, -11]
])
