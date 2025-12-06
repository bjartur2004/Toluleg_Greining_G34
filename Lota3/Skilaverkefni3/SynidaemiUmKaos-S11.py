import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
mpl.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\bin\ffmpeg.exe"

m1 = 1.0; m2 = 1.0
L1 = 1.0; L2 = 1.0
g  = 9.81

# diffur jöfnuhneppið fyrir tvöfaldan pendúl
def ydt(y, t):
    theta1, theta2, omega1, omega2 = y
    delta = theta2 - theta1

    denom1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    denom2 = (L2 / L1) * denom1

    domega1 = (m2 * L1 * omega1**2 * np.sin(delta) * np.cos(delta) +
               m2 * g * np.sin(theta2) * np.cos(delta) +
               m2 * L2 * omega2**2 * np.sin(delta) -
               (m1 + m2) * g * np.sin(theta1)) / denom1

    domega2 = (-m2 * L2 * omega2**2 * np.sin(delta) * np.cos(delta) +
               (m1 + m2) * (g * np.sin(theta1) * np.cos(delta) -
                            L1 * omega1**2 * np.sin(delta) -
                            g * np.sin(theta2))) / denom2

    return np.array([omega1, omega2, domega1, domega2])

# RK4 algorithmi
def rungeKutta(dy, y0, T, h):
    y = [np.array(y0)]
    for t in np.arange(0, T, h):
        yi = y[-1]
        k1 = dy(yi, t)
        k2 = dy(yi + h/2 * k1, t + h/2)
        k3 = dy(yi + h/2 * k2, t + h/2)
        k4 = dy(yi + h   * k3, t + h)
        y.append(yi + h/6*(k1 + 2*k2 + 2*k3 + k4))
    return np.array(y)


# búa til 3 mysmunandi plot
def main():
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    pendulums = []

    # Shared baseline pendulum (the same for all plots)
    y0_base = [np.pi/2, np.pi/2, 0.0, 0.0]
    Y_base = rungeKutta(ydt, y0_base, T=30, h=0.05)

    # Variations for each subplot
    ks = [3,5,8]
    nudges = [10**(-k) for k in ks]
    for e in nudges:
        y0_var = [np.pi/2 + e, np.pi/2 + e, 0.0, 0.0]
        Y_var = rungeKutta(ydt, y0_var, T=30, h=0.05)
        pendulums.append((Y_base, Y_var))  # tuple of (baseline, varied)

    lines = []
    traces = []
    traces_x = []
    traces_y = []

    colors = ['tab:blue', 'tab:orange']  # baseline and varied

    # Setup plots
    for ax in axes:
        ax.set_xlim(-(L1+L2)*1.2, (L1+L2)*1.2)
        ax.set_ylim(-(L1+L2)*1.2, (L1+L2)*0.2)
        ax.set_aspect('equal')
        ax.grid(True)
        line_pair = []
        trace_pair = []
        for c in colors:
            line, = ax.plot([], [], 'o-', lw=2, color=c)
            trace, = ax.plot([], [], '-', lw=1, alpha=0.5, color=c)
            line_pair.append(line)
            trace_pair.append(trace)
        lines.append(line_pair)
        traces.append(trace_pair)
        traces_x.append([[], []])
        traces_y.append([[], []])
        ax.set_title(f'Upphafstöðu Skekkja 10^-{ks[list(axes).index(ax)]}')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
    fig.suptitle('Tveir Tvöfaldir Pendúlar Með Smáan Mun Á Upphafsstöðu', fontsize=16)

    t_len = len(pendulums[0][0][:,0])  # number of frames

    def init():
        for line_pair, trace_pair in zip(lines, traces):
            for line, trace in zip(line_pair, trace_pair):
                line.set_data([], [])
                trace.set_data([], [])
        return [l for pair in lines for l in pair] + [t for pair in traces for t in pair]

    def update(frame):
        for i, (Y_base, Y_var) in enumerate(pendulums):
            for j, Y in enumerate([Y_base, Y_var]):
                X1 = L1*np.sin(Y[frame,0])
                Y1 = -L1*np.cos(Y[frame,0])
                X2 = X1 + L2*np.sin(Y[frame,1])
                Y2 = Y1 - L2*np.cos(Y[frame,1])

                x = [0, X1, X2]
                y_vals = [0, Y1, Y2]
                lines[i][j].set_data(x, y_vals)

                traces_x[i][j].append(X2)
                traces_y[i][j].append(Y2)
                traces[i][j].set_data(traces_x[i][j], traces_y[i][j])

        return [l for pair in lines for l in pair] + [t for pair in traces for t in pair]

    ani = FuncAnimation(fig, update, frames=t_len,
                        init_func=init, blit=True, interval=10)

    # Save video
    ani.save("double_pendulums.mp4", writer=FFMpegWriter(fps=24))
    print("Done")




if __name__ == "__main__":
    main()
