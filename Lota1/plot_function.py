import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return 2*x**4 + 23*x**3 + 45*x**2 + 60*x + 50

xRange = [-20,20]
samples = 1000

# plot
x = np.linspace(xRange[0], xRange[1], samples)
f_x = f(x)
# huh
plt.plot(x,f_x)
plt.xlim(xRange)

plt.show()