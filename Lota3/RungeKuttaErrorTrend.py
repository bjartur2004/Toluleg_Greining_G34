# finna skekkju trendið þegar n breitist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def ydt(y, t):
    y1, y2 = y
    dy1dt = y1 + y2
    dy2dt = -y1 + 2*y2
    return np.array([dy1dt, dy2dt])


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
    return (np.random.rand(2)-0.5)*4*np.pi
T = 20
n = np.array([200, 400, 800, 1600, 3200, 6400])
N = 10*max(n)
N_OF_INITIAL_POSITIONS = 3

Y0 = [randInitialPosition() for _ in range(N_OF_INITIAL_POSITIONS)]

YT = [rungeKuttaLastElement(ydt, y0, T, N) for y0 in Y0]

error = []
for pidx in range(N_OF_INITIAL_POSITIONS):
    error.append([])
    for ni in n:
        error[pidx].append(np.linalg.norm(rungeKuttaLastElement(ydt, Y0[pidx], T, ni) - YT[pidx]))

logError = np.log10(error)
logh = np.log10(T/n)
colors = cm.rainbow(np.linspace(0, 1, len(n)))
for logerror_for_initial_pos, color in zip(error, colors):
    plt.scatter(logh, logerror_for_initial_pos, color=color)


loghTiled = np.tile(logh, N_OF_INITIAL_POSITIONS)
logErrorFlat = logError.flatten()
print(np.shape(loghTiled), np.shape(logErrorFlat))
exponent, offset = np.polyfit(loghTiled, logErrorFlat, 1)
print("Exponent:", exponent)

logFitH = np.array([min(logh), max(logh)])
logFitErr = exponent * logFitH + offset

plt.plot(logFitH, logFitErr)

plt.text(-1.5,6,f"Halli: {exponent:.3f}")

plt.show()
