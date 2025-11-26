import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


litir = ["b", "g", "k", "y", "r"]  

allar_roots = np.roots([2, 23, 45, 60, 50])

def p(s):
    return 2*s**4 + 23*s**3 + 45*s**2 + 60*s + 50

def dp(s):
    return 8*s**3 + 69*s**2 + 90*s + 60
t = (0)
def Newton_Fraktalar(z0, roots, tol=1e-10, maxit=1000):
    z = complex(z0)
    for t in range(maxit):
        dz = dp(z)
        if dz == 0:
            return -1
        z = z - p(z)/dz
        for i, r in enumerate(roots):
            if abs(z - r) < tol:
                return i
    return -1 
    
z0 = 5 + 7j
result = Newton_Fraktalar(z0, allar_roots)
print("Rót nr", result)



nx = ny = 400
raunhluti = np.linspace(-8, -6, nx)
þverhluti = np.linspace(-1, 1, ny)


m = np.zeros((ny, nx), dtype=int)


for ix, re0 in enumerate(raunhluti):
    for iy, im0 in enumerate(þverhluti):
        z0 = complex(re0, im0)
        idx = Newton_Fraktalar(z0, allar_roots)
        if idx == -1:
            idx = 4   
        m[iy, ix] = idx


cmap = ListedColormap(litir)

plt.figure(figsize=(6, 6))
plt.imshow(m, cmap=cmap, origin="lower",
           extent=[raunhluti[0], raunhluti[-1], þverhluti[0], þverhluti[-1]])
plt.scatter([], [], color="b", label="Rót 0")
plt.scatter([], [], color="g", label="Rót 1")
plt.scatter([], [], color="k", label="Rót 2")
plt.scatter([], [], color="y", label="Rót 3")
plt.scatter([], [], color="r", label="Engin samleitni")
plt.legend(loc="upper right", fontsize=9)
plt.title("Newton-fraktalar")
plt.xlabel("Raunhluti")
plt.ylabel("Þverhluti")
plt.show()