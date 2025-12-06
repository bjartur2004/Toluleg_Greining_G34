import numpy as np
import math
import matplotlib.pyplot as plt

def ydot2(t, y):
    theta1, omega1, theta2, omega2 = y

    L1 = L2 = 2.0
    m1 = m2 = 1.0
    g = 9.81

    delta = theta2 - theta1

    den1 = (m1 + m2) * L1 - m2 * L1 * math.cos(delta)**2
    den2 = (m1 + m2) * L2 - m2 * L2 * math.cos(delta)**2

    dydt = np.zeros(4)

    dydt[0] = omega1

    dydt[1] = ( m2*L1*omega1**2*math.sin(delta)*math.cos(delta)
                + m2*g*math.sin(theta2)*math.cos(delta)
                + m2*L2*omega2**2*math.sin(delta)
                - (m1+m2)*g*math.sin(theta1) ) / den1

    dydt[2] = omega2

    dydt[3] = ( -m2*L2*omega2**2*math.sin(delta)*math.cos(delta)
                + (m1+m2)*( g*math.sin(theta1)*math.cos(delta)
                - L1*omega1**2*math.sin(delta) - g*math.sin(theta2) ) ) / den2

    return dydt



def rk4_step(t, x, h):
    k1 = ydot2(t, x)
    k2 = ydot2(t + 0.5*h, x + 0.5*h*k1)
    k3 = ydot2(t + 0.5*h, x + 0.5*h*k2)
    k4 = ydot2(t + h,     x + h*k3)
    return x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)


def rk4_solve(x0, T, n):
    t = np.linspace(0, T, n+1)
    h = T / n

    w = np.zeros((n+1, len(x0)))
    w[0] = x0

    for i in range(n):
        w[i+1] = rk4_step(t[i], w[i], h)

    return t, w


theta1_0 = np.pi/3
omega1_0 = 0
theta2_0 = np.pi/6
omega2_0 = 0

x0 = [theta1_0, omega1_0, theta2_0, omega2_0]

T = 100
n = 10000

t, w = rk4_solve(x0, T, n)

theta1 = w[:, 0]
theta2 = w[:, 2]

plt.figure(figsize=(6,6))
plt.plot(theta1, theta2, linewidth=1)
plt.xlabel(r"$\theta_1(t)$")
plt.ylabel(r"$\theta_2(t)$")
plt.title("Stikaður ferill í R²: (θ1(t), θ2(t))")
plt.grid(True)
plt.axis('equal')
plt.show()
