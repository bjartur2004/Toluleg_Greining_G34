import numpy as np
import matplotlib.pyplot as plt
from eulerstep import eulerstep

def pendulum_euler(theta0, dtheta0, T, n):
    
    h = T / n
    t = np.linspace(0.0, T, n + 1)

    Y = np.zeros((n + 1, 2))
    Y[0, 0] = theta0
    Y[0, 1] = dtheta0

    for k in range(n):
        Y[k + 1, :] = eulerstep(t[k], Y[k, :], h)
        
    return t, Y


