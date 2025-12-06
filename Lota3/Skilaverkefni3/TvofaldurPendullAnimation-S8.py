import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation #, FFMpegWriter
from doublependul7 import double_pendulum_rk4

L1 = 2.0
L2 = 2.0

def animate_double_pendulum(theta1_0, omega1_0, theta2_0, omega2_0, T, n):

    t, Y = double_pendulum_rk4(theta1_0, omega1_0, theta2_0, omega2_0, T, n)

    theta1 = Y[:, 0]
    theta2 = Y[:, 2]

    # reiknum staðsetningu
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)

    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)

    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-4.5, 1.0)
    ax.set_title("Double Pendulum Animation")

    # pendúl línur
    line1, = ax.plot([], [], 'o-', lw=2, color='blue')
    line2, = ax.plot([], [], 'o-', lw=2, color='red')

    # slóð massans
    path1, = ax.plot([], [], lw=1, color='green')
    path2, = ax.plot([], [], lw=1, color='purple')

    def update(frame):
        line1.set_data([0, x1[frame]], [0, y1[frame]])
        line2.set_data([x1[frame], x2[frame]], [y1[frame], y2[frame]])

        path1.set_data(x1[:frame], y1[:frame])
        path2.set_data(x2[:frame], y2[:frame])
        return line1, line2, path1, path2

    ani = FuncAnimation(fig, update, frames=len(t), interval=1000*T/n)

    # --- SAVE AS MP4 ---
    #writer = FFMpegWriter(fps=30)
    #ani.save("upphafsgildi_3_8.mp4", writer=writer)

    plt.show()



# Hér breytum við handvirkt um upphafsgildi og keyrum svo.
# þetta er t.d fyrsta upphafsgildið, hér breytum við svo bara gildum á thetum og keyrum fyrir mismunandi 5 gildi.
# Hin gildin haldast óbreytt.
if __name__ == "__main__":
    animate_double_pendulum(
        theta1_0=np.pi/3,
        omega1_0=0,
        theta2_0=np.pi/6,
        omega2_0=0,
        T=20,
        n=200
    )
