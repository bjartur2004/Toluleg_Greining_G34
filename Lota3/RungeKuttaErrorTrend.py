# finna skekkju trendið þegar n breitist
import numpy as np
import matplotlib.pyplot as plt

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

y0 = np.array([1,1])
T = 10

n = np.array([200, 400, 800, 1600, 3200, 6400])
N = 10*max(n)

YT = rungeKuttaLastElement(ydt, y0, T, N)

error = []
for ni in n:
    error.append(np.linalg.norm(rungeKuttaLastElement(ydt, y0, T, ni) - YT))

logError = np.log10(error)
logh = np.log10(T/n)

plt.scatter(logh, logError)

exponent, offset = np.polyfit(logh, logError, 1)
print("Exponent:", exponent)

logFitH = np.array([min(logh), max(logh)])
logFitErr = exponent * logFitH + offset

plt.plot(logFitH, logFitErr)

plt.text(-2,-2,f"Halli: {exponent:.3f}")

plt.show()
