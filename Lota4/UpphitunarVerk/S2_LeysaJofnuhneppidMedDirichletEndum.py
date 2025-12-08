import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin

a = 0
b = np.pi/2
n = 50
h = (b-a)/(n-1)

b = np.concatenate([[1], np.zeros(n-2), [-5]])

alpha = -2+h**2
A = np.concatenate([
    [np.concatenate([[1], np.zeros(n-1)])],
    [np.concatenate([np.zeros(i),[1,alpha,1], np.zeros(n-i-3)]) for i in range(n-2)],
    [np.concatenate([np.zeros(n-1), [1]])]
])

u = np.linalg.solve(A, b)

print(A)
#print(f"u = {u}")
x = [a+h*i for i in range(n)]
y = [cos(xi)-5*sin(xi) for xi in x]

plt.plot(x,u, ".")
plt.plot(x,y)
plt.legend(["Nálgun [u]", "Nálvæm lausn [y]"])
plt.title("Lausn Við Dæmi 2")
plt.xlabel("x")
plt.show()

