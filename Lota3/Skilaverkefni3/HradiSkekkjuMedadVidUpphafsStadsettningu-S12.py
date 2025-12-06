import numpy as np
import matplotlib.pyplot as plt

# cosnt
m1 = 1.0
m2 = 1.0  
L1 = 1.0   
L2 = 1.0   
g  = 9.81  


def ydt(y, t):
    theta1, theta2, omega1, omega2 = y

    delta = theta2 - theta1

    denom1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    denom2 = (L2 / L1) * denom1

    domega1 = (m2 * L1 * omega1**2 * np.sin(delta) * np.cos(delta) +
               m2 * g * np.sin(theta2) * np.cos(delta) +
               m2 * L2 * omega2**2 * np.sin(delta) -
               (m1 + m2) * g * np.sin(theta1)) / denom1

    domega2 = (- m2 * L2 * omega2**2 * np.sin(delta) * np.cos(delta) +
               (m1 + m2) * (g * np.sin(theta1) * np.cos(delta) -
                            L1 * omega1**2 * np.sin(delta) -
                            g * np.sin(theta2))) / denom2

    return np.array([omega1, omega2, domega1, domega2])


    return y

def SingleStepRungeKutta(dy, y_current, t_current, h):
        half_h = h/2
        k1 = dy(y_current, t_current)
        k2 = dy(y_current + half_h * k1, t_current + half_h)
        k3 = dy(y_current + half_h * k2, t_current + half_h)
        k4 = dy(y_current + h * k3, t_current + h)
        return (y_current+h/6*(k1+2*k2+2*k3+k4))

def findEndPointPosition(y):
    theta1, theta2, _, _ = y
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)
    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)

    return (x2,y2)

def upphafsgildi(i):
    return np.pi/500 * i

def simulate_divergence(args):
    i, j, e, t_max, h, diverganceTreshold = args
    o1 = upphafsgildi(i)
    o2 = upphafsgildi(j)

    yacc = np.array([o1, o2, 0.0, 0.0])
    yoff = np.array([o1+e, o2+e, 0.0, 0.0])

    for t in np.arange(0, t_max, h):
        yacc = SingleStepRungeKutta(ydt, yacc, t, h)
        yoff = SingleStepRungeKutta(ydt, yoff, t, h)

        accpos = findEndPointPosition(yacc)
        offpos = findEndPointPosition(yoff)
        dist = np.linalg.norm(np.array(accpos) - np.array(offpos))

        if dist >= diverganceTreshold:
            #print(f"[{i},{j}] > {t:.5f}")
            return (i, j, t)

    #print(f"[{i},{j}] > Did Not diverge")
    return (i, j, np.inf)



# multi thread 
if __name__ == "__main__":
    from multiprocessing import Pool

    N = 500 # fjöldi o1 upphafs gilda
    M = 500 # fjöldi o2 upphafs gilda
    diverganceTreshold = 0.1
    e = 1e-8
    t_max = 200
    h = 0.05

    args_list = []
    for i in range(N):
        for j in range(M):
            args_list.append([i,j,e,t_max,h,diverganceTreshold])

    with Pool(processes=16) as pool:
        results = pool.map(simulate_divergence, args_list)

    diveranceTime = np.full((N, M), np.inf)
    for i, j, t in results:
        diveranceTime[i, j] = t

    # plotta
    theta_vals = [upphafsgildi(i) for i in range(N)]
    masked = np.ma.masked_where(np.isinf(diveranceTime), diveranceTime)

    cmap = plt.colormaps['turbo_r'].copy()
    cmap.set_bad(color='black')

    plt.figure(figsize=(6,5))
    plt.imshow(masked, origin='lower', cmap=cmap,
               extent=[theta_vals[0], theta_vals[-1], theta_vals[0], theta_vals[-1]],
               aspect='auto')
    plt.colorbar(label='Fráviks Tími (s)')
    plt.xlabel("Upphafs θ2 (rad)")
    plt.ylabel("Upphafs θ1 (rad)")
    plt.title("Fráviks tími tvöfaldra pendúla")
    plt.show()
