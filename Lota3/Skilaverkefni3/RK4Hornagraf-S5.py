import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

g = 9.81
l = 2

def f(y):
    y1 = y[0]         
    y2 = y[1]         
    y1_diffrað = y2
    y2_diffrað = -(g/l) * np.sin(y1)
    return np.array([y1_diffrað, y2_diffrað])

def rungeKutta(y, h):
    k1 = f(y)
    k2 = f(y + 0.5 * h * k1)
    k3 = f(y + 0.5 * h * k2)
    k4 = f(y + h * k3)
    return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

def rk_lausn(y0, T, n):
    t = np.zeros(n+1)
    w = np.zeros((n+1, len(y0)))

    w[0, :] = y0
    t[0] = 0
    h = T / n

    for i in range(n):
        w[i+1, :] = rungeKutta(w[i, :], h)
        t[i+1] = t[i] + h

    return t, w

#Fyrra skilyrði
#theta0 = np.pi/12
theta0 = np.pi/2    
thetadiffrað0 = 0          
T = 20
n = 500


tgildi, wgildi = rk_lausn([theta0, thetadiffrað0], T, n)

theta = wgildi[:, 0]
x_data = l * np.sin(theta)
y_data = -l * np.cos(theta)


plt.close("all")
fig = plt.figure(figsize=(10,5))

pad = 0
ax1 = fig.add_subplot(
    121,
    autoscale_on=False,
    xlim=(-l-pad, l+pad),
    ylim=(-l-pad, l+pad)
)

plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Rauntíma hreyfing pendúls (RK4)")

fig.tight_layout()

line_1, = ax1.plot([], [], 'o-', lw=2, ms=8)

def init():
    line_1.set_data([], [])
    return line_1,

def animate(i):
    x = [0, x_data[i]]
    y = [0, y_data[i]]
    line_1.set_data(x, y)
    return line_1,

ax2 = fig.add_subplot(122)
ax2.set_title("Graf af θ(t)")
ax2.set_xlabel("t")
ax2.set_ylabel("θ(t)")
ax2.grid(True)

ax2.plot(tgildi, theta, 'b')
anim = animation.FuncAnimation(
    fig, animate,
    frames=n,
    interval=1,
    blit=True,
    repeat=False,
    init_func=init
)
plt.tight_layout()
plt.show()

Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
anim.save("pendulum_rk4.mp4", writer=writer)
