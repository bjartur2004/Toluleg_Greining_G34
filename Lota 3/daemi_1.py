import numpy as np

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


def eulerstep(y, h):
    return y + h * f(y)


def f(y):
    l = 2
    g =  9.81   
    y1= y[0]
    y2 = y[1]
    y1_diffrad = y2
    y2_diffrad = -g/l * np.sin(y1)

    return np.array(y1_diffrad,y2_diffrad)
    
