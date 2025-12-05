import numpy as np
from RK4 import rk4

def double_pendulum_rk4(theta1_0, omega1_0, theta2_0, omega2_0, T, n):
    
    h = T / n
    t = np.linspace(0.0, T, n + 1)
    Y = np.zeros((n + 1, 4))
    Y[0, :] = [theta1_0, omega1_0, theta2_0, omega2_0]
    for k in range(n):
        Y[k + 1, :] = rk4(t[k], Y[k, :], h)

    return t, Y

