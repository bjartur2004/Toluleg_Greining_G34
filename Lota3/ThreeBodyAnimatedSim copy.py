import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Constants for 3-body problem
m1 = 1.0
m2 = 1.0
m3 = 1.0
G  = 1.0  # Gravitational constant

x1, y1 = 0.97000436, -0.24308753
x2, y2 = -0.97000436, 0.24308753
x3, y3 = 0.0, 0.0

vx1, vy1 = 0.4662036850, 0.4323657300
vx2, vy2 = 0.4662036850, 0.4323657300
vx3, vy3 = -0.93240737, -0.86473146

y0 = np.array([x1, y1, x2, y2, x3, y3,
               vx1, vy1, vx2, vy2, vx3, vy3])



def ydt(y, t):
    # Unpack positions and velocities
    x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3 = y
    
    # Compute pairwise distances
    r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    r13 = np.sqrt((x3-x1)**2 + (y3-y1)**2)
    r23 = np.sqrt((x3-x2)**2 + (y3-y2)**2)
    
    # Compute accelerations using Newton's law of gravitation
    ax1 = G * (m2*(x2-x1)/r12**3 + m3*(x3-x1)/r13**3)
    ay1 = G * (m2*(y2-y1)/r12**3 + m3*(y3-y1)/r13**3)
    
    ax2 = G * (m1*(x1-x2)/r12**3 + m3*(x3-x2)/r23**3)
    ay2 = G * (m1*(y1-y2)/r12**3 + m3*(y3-y2)/r23**3)
    
    ax3 = G * (m1*(x1-x3)/r13**3 + m2*(x2-x3)/r23**3)
    ay3 = G * (m1*(y1-y3)/r13**3 + m2*(y2-y3)/r23**3)
    
    return np.array([vx1, vy1, vx2, vy2, vx3, vy3,
                     ax1, ay1, ax2, ay2, ax3, ay3])


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
h = 0.02
t_vals = np.arange(0, T+h, h)
y = np.array(rungeKutta(ydt, y0, T, h))

def plot_three_body(y, t_vals, trail=True):

    n_bodies = 3
    x = [y[:, i*2]     for i in range(n_bodies)]  # x1, x2, x3
    y_pos = [y[:, i*2+1] for i in range(n_bodies)]  # y1, y2, y3
    colors = ['blue', 'green', 'red']
    traces = [[] for _ in range(n_bodies)] if trail else None

    fig, ax = plt.subplots()
    #margin = 1.5 * max(np.max(np.abs(y[:,0:6])), 1)
    margin = 2
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)
    ax.set_aspect('equal')
    ax.grid(True)

    lines = [ax.plot([], [], 'o', lw=2, color=colors[i])[0] for i in range(n_bodies)]
    trace_lines = [ax.plot([], [], '-', lw=1, color=colors[i], alpha=0.5)[0] for i in range(n_bodies)] if trail else None

    def init():
        for line in lines:
            line.set_data([], [])
        if trail:
            for tline in trace_lines:
                tline.set_data([], [])
        return lines + (trace_lines if trail else [])

    def update(frame):
        for i in range(n_bodies):
            # wrap in list so set_data gets a sequence
            lines[i].set_data([x[i][frame]], [y_pos[i][frame]])
            if trail:
                traces[i].append((x[i][frame], y_pos[i][frame]))
                tx, ty = zip(*traces[i])
                trace_lines[i].set_data(tx, ty)
        return lines + (trace_lines if trail else [])

    ani = FuncAnimation(fig, update, frames=len(t_vals),
                        init_func=init, blit=True, interval=10)
    plt.show()


plot_three_body(y, t_vals)