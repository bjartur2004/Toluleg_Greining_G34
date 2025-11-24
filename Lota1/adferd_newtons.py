import matplotlib.pyplot as mpl
import numpy as np

def f(x):
    return 2*x**4 + 23*x**3 + 45*x**2 + 60*x + 50


def df(x):
    return 8*x**3 + 69*x**2 + 80*x + 60

# plot
x = np.linspace(-10,10,20)
f_x = f(x)

mpl.plot(x,f_x)
#mpl.show()

I_GUESS = -9.41
TOLERANCE = 0.00001

def adferd_newtons(I_guess, tolerance):
    x = I_guess
    x_old = 99999 # eða bara I_guess+2*tolerance

    while abs(x-x_old)>tolerance:
        x_old = x
        x = x_old - f(x_old) / df(x_old)
    return x


print(adferd_newtons(I_GUESS, TOLERANCE))