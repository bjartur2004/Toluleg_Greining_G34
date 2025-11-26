import numpy as np
import cmath, math
import matplotlib.pyplot as plt

allar_roots = np.roots([2, 23, 45, 60, 50])

def dp(s):
    return 8*s**3 + 69*s**2 + 90*s + 60

def p(s):
    return 2*s**4 + 23*s**3 + 45*s**2 + 60*s + 50


def newton_converge_index(z0, roots, tol=1e-8, maxit=100):
    z = complex(z0)
    for z in range(maxit):
        dz = dp(z)
        if dz == 0:
            return -1
        z = z - p(z)/dz
        for i, r in enumerate(roots):
            if abs(z - r) < tol:
                return i
    return -1


nx = ny = 400
re_vals = np.linspace(-8, -6, nx)
im_vals = np.linspace(-1, 1, ny)
X = []
Y = []
labels = []
for i, re0 in enumerate(re_vals):
    for j, im0 in enumerate(im_vals):
        z0 = complex(re0, im0)
        idx = newton_converge_index(z0, allar_roots, tol=1e-8, maxit=80)
        X.append(re0); Y.append(im0); labels.append(idx if idx>=0 else 4)


plt.figure(figsize=(6,6))
plt.scatter(X, Y, c=labels, s=1)
plt.title("Newton fractal")
plt.xlabel("Raunhluti")
plt.ylabel("Þverhluti")
plt.show()


