import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Constants
m1 = 1.0   # Mass of first pendulum
m2 = 1.0   # Mass of second pendulum
L1 = 1.0   # Length of first rod
L2 = 1.0   # Length of second rod
g  = 9.81  # Gravity

# Initial state: [theta1, theta2, omega1, omega2]
y0 = np.array([np.pi*0.25, np.pi*0.25, 0.0, 0.0])

def ydt(y, t):
    theta1, theta2, omega1, omega2 = y

    delta = theta2 - theta1

    denom1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    denom2 = (L2 / L1) * denom1

    domega1 = (m2 * L1 * omega1**2 * np.sin(delta) * np.cos(delta) +
               m2 * g * np.sin(theta2) * np.cos(delta) +
               m2 * L2 * omega2**2 * np.sin(delta) -
               (m1 + m2) * g * np.sin(theta1)) / denom1

    domega2 = (- m2 * L2 * omega2**2 * np.sin(delta) * np.cos(delta) +
               (m1 + m2) * (g * np.sin(theta1) * np.cos(delta) -
                            L1 * omega1**2 * np.sin(delta) -
                            g * np.sin(theta2))) / denom2

    return np.array([omega1, omega2, domega1, domega2])


def rungeKutta(dy, y0, T, h):
    y = [np.array(y0, dtype=float)] 
    t0 = T[0] if type(T) is list else 0
    T = T[1] if type(T) is list else T
    for t in np.arange(t0, T, h):
        half_h = h/2
        yi = y[-1]
        k1 = dy(yi, t)
        k2 = dy(yi + half_h * k1, t + half_h)
        k3 = dy(yi + half_h * k2, t + half_h)
        k4 = dy(yi + h * k3, t + h)
        y.append(yi+h/6*(k1+2*k2+2*k3+k4))

    return y

T = 20
h = 0.01
t_vals = np.arange(0, T+h, h)
y = np.array(rungeKutta(ydt, y0, T, h))


# Compute (x, y) positions of both pendulums
x1 = L1 * np.sin(y[:,0])
y1 = -L1 * np.cos(y[:,0])

x2 = x1 + L2 * np.sin(y[:,1])
y2 = y1 - L2 * np.cos(y[:,1])

fig, ax = plt.subplots()
ax.set_xlim(-(L1+L2)*1.2, (L1+L2)*1.2)
ax.set_ylim(-(L1+L2)*1.2, (L1+L2)*0.2)
ax.set_aspect('equal')
ax.grid(True)

line, = ax.plot([], [], 'o-', lw=2, color='blue')
trace, = ax.plot([], [], '-', lw=1, color='red', alpha=0.5)  # trace of second bob

# Store previous positions for trace
x2_trace = []
y2_trace = []

def init():
    line.set_data([], [])
    trace.set_data([], [])
    return line, trace

def update(frame):
    x = [0, x1[frame], x2[frame]]
    y_vals = [0, y1[frame], y2[frame]]
    line.set_data(x, y_vals)

    # Update trace
    x2_trace.append(x2[frame])
    y2_trace.append(y2[frame])
    trace.set_data(x2_trace, y2_trace)

    return line, trace

ani = FuncAnimation(fig, update, frames=len(t_vals),
                    init_func=init, blit=True, interval=10)

plt.show()
