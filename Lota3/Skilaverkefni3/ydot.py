import numpy as np
import math

def ydot(t, y):
    g = 9.81
    L = 2.0

    dydt = np.zeros(2)
    dydt[0] = y[1]                        
    dydt[1] = -(g / L) * math.sin(y[0])   
    return dydt
