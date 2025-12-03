from ydot2 import ydot2
import numpy as np

def rk4(t, x, h):
    k1 = ydot2(t, x)
    k2 = ydot2(t + 0.5*h, x + 0.5*h*k1)
    k3 = ydot2(t + 0.5*h, x + 0.5*h*k2)
    k4 = ydot2(t + h,     x + h*k3)

    return x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)


