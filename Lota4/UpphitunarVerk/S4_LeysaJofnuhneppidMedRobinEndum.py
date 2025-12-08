import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin

a = 0
b = np.pi/2
n = 50
h = (b-a)/(n-1)

b = np.concatenate([[0], np.zeros(n-2), [-2*h]])

alpha = -2+h**2
c = -3-2*h
d = -3-4*h
A = np.concatenate([
    [np.concatenate([[c,4,-1], np.zeros(n-3)])],
    [np.concatenate([np.zeros(i),[1,alpha,1], np.zeros(n-i-3)]) for i in range(n-2)],
    [np.concatenate([np.zeros(n-3), [-1,4,d]])]
])

u = np.linalg.solve(A, b)

#print(b)
#print(A)

print(f"u = {u}")
x = [a+h*i for i in range(n)]
y = [cos(xi)+sin(xi) for xi in x]

plt.plot(x,u, ".")
plt.plot(x,y)
plt.legend(["Nálgun [u]", "Nálvæm lausn [y]"])
plt.title("Lausn Við Dæmi 4")
plt.xlabel("x")
plt.show()

