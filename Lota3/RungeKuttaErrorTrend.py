# finna skekkju trendið þegar n breitist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

def ydot2(y,t):
    # y = [þeta1, omega1, þeta2, omega2]
    theta1, omega1, theta2, omega2 = y

    # fastar
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

def rungeKuttaLastElement(dy, y0, T, n):
    h = T/n
    yi = np.array(y0, dtype=float)
    t0 = T[0] if type(T) is list else 0
    T = T[1] if type(T) is list else T
    for i in range(n):
        t = i*h
        k1 = dy(yi, t)
        k2 = dy(yi + h/2 * k1, t + h/2)
        k3 = dy(yi + h/2 * k2, t + h/2)
        k4 = dy(yi + h * k3, t + h)
        yi = yi + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return yi

def randInitialPosition():
    angle = (np.random.rand(2)-0.5)*2*np.pi * 0.25
    speed = (np.random.rand(2)-0.5)*2*np.pi * 0
    return np.array([angle[0],speed[0],angle[1],speed[1]])

T = 20
n = np.array([200, 400, 800, 1600, 3200, 6400])
N = 10*max(n)
N_OF_INITIAL_POSITIONS = 10


Y0 = [randInitialPosition() for _ in range(N_OF_INITIAL_POSITIONS)]

YT = [rungeKuttaLastElement(ydot2, y0, T, N) for y0 in Y0]

error = []
for pidx in range(N_OF_INITIAL_POSITIONS):
    error.append([])
    for ni in n:
        error[pidx].append(np.linalg.norm(rungeKuttaLastElement(ydot2, Y0[pidx], T, ni) - YT[pidx]))

logError = np.log10(error)
logh = np.log10(T/n)
colors = cm.rainbow(np.linspace(0, 1, N_OF_INITIAL_POSITIONS))
for logerror_for_initial_pos, color in zip(logError, colors):
    plt.scatter(logh, logerror_for_initial_pos, color=color)


loghTiled = np.tile(logh, N_OF_INITIAL_POSITIONS)
logErrorFlat = logError.flatten()
exponent, offset = np.polyfit(loghTiled, logErrorFlat, 1)

print("Exponent:", exponent)

logFitH = np.array([min(logh), max(logh)])
logFitErr = exponent * logFitH + offset

plt.plot(logFitH, logFitErr)

plt.text(-1.5,-6,f"Halli: {exponent:.3f}")

plt.show()
