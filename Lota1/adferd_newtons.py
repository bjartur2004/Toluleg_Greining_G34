import matplotlib.pyplot as mpl
import numpy as np

def f(x):
    return x**2-2

def df(x):
    return 2*x

# plot
x = np.linspace(-10,10,20)
f_x = f(x)

mpl.plot(x,f_x)
mpl.show()

I_GUESS = 10
TOLERANCE = 0.001

def adferd_newtons(I_guess, tolerance):
    x = I_guess
    x_old = 99999 # eða bara I_guess+2*tolerance

    while abs(x-x_old)>tolerance:
        x_old = x
        x = x_old - f(x_old) / df(x_old)

    return x


print(adferd_newtons(I_GUESS, TOLERANCE))