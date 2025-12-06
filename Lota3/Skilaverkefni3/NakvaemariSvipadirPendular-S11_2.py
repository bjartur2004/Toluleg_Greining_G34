import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
mpl.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\bin\ffmpeg.exe"


# Constants
m1 = 1.0   
m2 = 1.0   
L1 = 1.0   
L2 = 1.0   
g  = 9.81  
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

T = 40
h = 0.02

N = 5 # number of pengelums
variation = 1e-10
# Initial state: [theta1, theta2, omega1, omega2], [...]
Y0 = np.array([[np.pi/2+variation*i, np.pi/2+variation*i, 0.0, 0.0] for i in range(N)])

t_vals = np.arange(0, T+h, h)
Y_all = []
for y0 in Y0:
    Y_all.append(np.array(rungeKutta(ydt, y0, T, h)))



# Compute pendulum positions
X1_all = [L1*np.sin(Y[:,0]) for Y in Y_all]
Y1_all = [-L1*np.cos(Y[:,0]) for Y in Y_all]
X2_all = [X1 + L2*np.sin(Y[:,1]) for X1,Y in zip(X1_all, Y_all)]
Y2_all = [Y1 - L2*np.cos(Y[:,1]) for Y1,Y in zip(Y1_all, Y_all)]

# Plot setup
fig, ax = plt.subplots()
ax.set_xlim(-(L1+L2)*1.2, (L1+L2)*1.2)
ax.set_ylim(-(L1+L2)*1.2, (L1+L2)*0.2)
ax.set_aspect('equal')
ax.grid(True)

# Create line objects for each pendulum
lines = [ax.plot([], [], 'o-', lw=2)[0] for _ in range(N)]
traces = [ax.plot([], [], '-', lw=1, alpha=0.5)[0] for _ in range(N)]
x2_traces = [[] for _ in range(N)]
y2_traces = [[] for _ in range(N)]

def init():
    for line, trace in zip(lines, traces):
        line.set_data([], [])
        trace.set_data([], [])
    return lines + traces

def update(frame):
    for i in range(N):
        x = [0, X1_all[i][frame], X2_all[i][frame]]
        y_vals = [0, Y1_all[i][frame], Y2_all[i][frame]]
        lines[i].set_data(x, y_vals)

        # Update trace
        x2_traces[i].append(X2_all[i][frame])
        y2_traces[i].append(Y2_all[i][frame])
        traces[i].set_data(x2_traces[i], y2_traces[i])

    return lines + traces

ani = FuncAnimation(fig, update, frames=len(t_vals),
                    init_func=init, blit=True, interval=10)

ani.save("nakvaemariPendull.mp4", writer=FFMpegWriter(fps=60))
