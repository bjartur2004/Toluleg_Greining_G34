import numpy as np
import math

def ydot2(t, y):
    # y = [þeta1, omega1, þeta2, omega2]
    theta1, omega1, theta2, omega2 = y

    # fastar
    L1 = L2 = 2.0
    m1 = m2 = 1.0
    g = 9.81

    delta = theta2 - theta1   # ∆

    # nefnari
    den1 = (m1 + m2) * L1 - m2 * L1 * math.cos(delta)**2
    den2 = (m1 + m2) * L2 - m2 * L2 * math.cos(delta)**2

    dydt = np.zeros(4)

    # þeta1' = omega1
    dydt[0] = omega1

    # omega1' = d²(theta1)/dt²
    dydt[1] = ( m2*L1*omega1**2*math.sin(delta)*math.cos(delta)
                + m2*g*math.sin(theta2)*math.cos(delta)
                + m2*L2*omega2**2*math.sin(delta)
                - (m1+m2)*g*math.sin(theta1) ) / den1

    # þeta2' = omega2
    dydt[2] = omega2

    # omega2' = d²(theta2)/dt²
    dydt[3] = ( -m2*L2*omega2**2*math.sin(delta)*math.cos(delta)
                + (m1+m2)*( g*math.sin(theta1)*math.cos(delta)
                - L1*omega1**2*math.sin(delta) - g*math.sin(theta2) ) ) / den2

    return dydt
