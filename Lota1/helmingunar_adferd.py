import matplotlib.pyplot as mpl
import numpy as np

def f(x):
    return 2*x**4 + 23*x**3 + 45*x**2 + 60*x + 50


# plot
x = np.linspace(-10,10,20)
f_x = f(x)
# huh
mpl.plot(x,f_x)
mpl.show()

MIN = -10
MAX = -5
TOLERANCE = 1e-6


# helmingunar aðferð sem finnur rót e. bysection = binery search
# í staðin fyrir minna en eða meira en notum við margföldun til að tékka hvort fallið hafi mismunandi formerki
def root_binary_search(I_MIN, I_MAX, tolerance):
    min_bound = I_MIN
    max_bound = I_MAX

    f_min = f(min_bound)
    f_max = f(max_bound)

    while (max_bound - min_bound)/2 > tolerance:
        m = (max_bound+min_bound)/2
        f_m = f(m)

        if f_m * f_max > 0: # ef lausnin er ekki á milli m og max
            max_bound = m
            f_max = f(max_bound)

        else:
            min_bound = m
            f_min = f(min_bound)

    # solution found
    return (min_bound+max_bound)/2

print(root_binary_search(MIN, MAX, TOLERANCE))

