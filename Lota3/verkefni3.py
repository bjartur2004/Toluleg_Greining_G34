import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


g = 9.81
l = 2

# 1. stigs diffurjöfnuhneppið 
def f(y):
    y1 = y[0]         
    y2 = y[1]         
    y1_diffrað = y2
    y2_diffrað = -(g/l) * np.sin(y1)
    return np.array([y1_diffrað, y2_diffrað])

def eulerstep(y, h):
    return y + h * f(y)

# aðferð eulers 
def eulersolver(y0, T, n):
    t = np.zeros(n+1)
    w = np.zeros((n+1, len(y0)))

    w[0, :] = y0
    t[0] = 0
    h = T / n

    for i in range(n):
        w[i+1, :] = eulerstep(w[i, :], h)
        t[i+1] = t[i] + h

    return t, w


theta0 = np.pi/2    
thetadiffrað0 = 0          
T = 20
n = 500


tgildi, wgildi = eulersolver([theta0, thetadiffrað0], T, n)


theta = wgildi[:, 0]
x_data = l * np.sin(theta)
y_data = -l * np.cos(theta)



plt.close("all")
fig = plt.figure()

pad = 0.5
ax1 = fig.add_subplot(
    111,
    autoscale_on=False,
    xlim=(-l-pad, l+pad),
    ylim=(-l-pad, l+pad)
)

plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Rauntíma hreyfing pendúls (Euler)")

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

anim = animation.FuncAnimation(
    fig, animate,
    frames=n,
    interval=1,
    blit=True,
    repeat=False,
    init_func=init
)

plt.show()


Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
anim.save("pendulum_euler.mp4", writer=writer)

