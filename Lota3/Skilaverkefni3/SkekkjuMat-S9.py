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

# við búum til mismunandi gildi með því að búa til tveggja laga tré
# af valmöguleikum fyrir angles og para þá saman þangað til að við
# erum komin með nógu mörg pör.
# Við höfum hraðan alltaf núll í upphafsstöðu
def initialPositionsTree(n_positions, angle_range):
    levels = int(np.ceil(np.log2(n_positions)))
    leaves = 2**levels
    angles = np.linspace(-angle_range, angle_range, leaves)

    Y0 = []
    for i in range(leaves):
        for j in range(leaves):
            Y0.append(np.array([angles[i], 0.0, angles[j], 0.0]))
            if len(Y0) >= n_positions:
                return Y0
    
    return Y0

T = 20
n = np.array([200, 400, 800, 1600, 3200, 6400])
N = 10 * max(n)
N_OF_INITIAL_POSITIONS = 30

Y0 = initialPositionsTree(N_OF_INITIAL_POSITIONS, 0.7*np.pi)

YT = [rungeKuttaLastElement(ydot2, y0, T, N) for y0 in Y0]

error = []
for pidx in range(N_OF_INITIAL_POSITIONS):
    error.append([])
    for ni in n:
        error[pidx].append(np.linalg.norm(rungeKuttaLastElement(ydot2, Y0[pidx], T, ni) - YT[pidx]))

logError = np.log10(error)
logh = np.log10(T/n)

last_errors = logError[:, 0] # we swich from n to h, so order flips          
sorted_indices = np.argsort(last_errors)  

cmap = cm.viridis
norm = plt.Normalize(vmin=last_errors.min(), vmax=last_errors.max())
mapped_colors = cmap(norm(last_errors))     

print("sorted indices:", sorted_indices)
print("sorted final errors:", last_errors[sorted_indices])

plt.figure(figsize=(9,5))
for z, idx in enumerate(sorted_indices):   
    c = mapped_colors[idx]
    plt.scatter(logh, logError[idx], color=c, s=40, edgecolor='k', linewidth=0.3, zorder=z)


loghTiled = np.tile(logh, N_OF_INITIAL_POSITIONS)
logErrorFlat = logError.flatten()
exponent, offset = np.polyfit(loghTiled, logErrorFlat, 1)

print("Exponent:", exponent)

logFitH = np.array([min(logh), max(logh)])
logFitErr = exponent * logFitH + offset

#plt.plot(logFitH, logFitErr)

#plt.text(-1.5,-6,f"Halli: {exponent:.3f}")

plt.title("Endapunkts skekksja með tilliti til skref stærðar")
plt.xlabel("h (skrefa stærð) [log10]")
plt.ylabel("Endapunkts Skekkja [log10]")

plt.show()


