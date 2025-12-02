import numpy as np
import matplotlib.pyplot as plt

run = 1
# --------- Single Variable Runge Kutta ---------
if run == 0:
    y0 = 1
    def ydt(y,t):
        return 2*y*t

    # for non vector problems
    def rungeKutta(dy, y0, T, h): 
        y = [y0]
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

    T = 3
    h = 0.01
    y = rungeKutta(ydt, y0, T, h)

    t_vals = np.arange(0, T+h, h)
    plt.plot(t_vals, y)
    plt.show()


# --------- Vector Runge Kutta ---------
if run == 1:
    y0 = np.array([1,1])
    def ydt(y, t):
        y1, y2 = y
        dy1dt = y1 + y2
        dy2dt = -y1 + 2*y2
        return np.array([dy1dt, dy2dt])


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

    T = 3
    h = 0.01
    y = rungeKutta(ydt, y0, T, h)

    t_vals = np.arange(0, T+h, h)
    y = np.array(rungeKutta(ydt, y0, T, h))

    plt.plot(t_vals, y[:,0], label='y1')
    plt.plot(t_vals, y[:,1], label='y2')
    plt.legend()
    plt.show()