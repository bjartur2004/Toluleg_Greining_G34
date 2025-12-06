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



theta1_vals = []
theta2_vals = []

theta1_0 = np.pi/3
omega1_0 = 0
theta2_0 = np.pi/6
omega2_0 = 0

x = np.array([theta1_0, omega1_0, theta2_0, omega2_0])

T = 100
n = 20000
h = T / n
t = 0.0

plt.ion()
fig, ax = plt.subplots(figsize=(6,6))

line, = ax.plot([], [], lw=1)
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_xlabel(r"$\theta_1(t)$")
ax.set_ylabel(r"$\theta_2(t)$")
ax.set_title("Rauntíma-teikning á ferlinum (θ1, θ2)")
ax.grid(True)

for i in range(n):

    
    x = rk4_step(t, x, h)
    t += h

    theta1_vals.append(x[0])
    theta2_vals.append(x[2])

    
    line.set_data(theta1_vals, theta2_vals)

    
    plt.pause(0.0001)

plt.ioff()
plt.show()
