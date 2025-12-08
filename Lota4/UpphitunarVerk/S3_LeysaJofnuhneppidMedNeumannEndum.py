import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin

a = 0
b = np.pi/2
n = 50
h = (b-a)/(n-1)

b = np.concatenate([[2*h], np.zeros(n-2), [-2*h*3]])

alpha = -2+h**2
A = np.concatenate([
    [np.concatenate([[-3,4,-1], np.zeros(n-3)])],
    [np.concatenate([np.zeros(i),[1,alpha,1], np.zeros(n-i-3)]) for i in range(n-2)],
    [np.concatenate([np.zeros(n-3), [-1,4,-3]])]
])

u = np.linalg.solve(A, b)

#print(b)
#print(A)

print(f"u = {u}")
x = [a+h*i for i in range(n)]
y = [-3*cos(xi)+sin(xi) for xi in x]

e = [np.linalg.norm(y[i]-u[i]) for i in range(n)]

# plot
fig, ax1 = plt.subplots()

ax1.plot(x, u, ".", label="Nálgun [u]")
ax1.plot(x, y, label="Nákvæm lausn [y]")
ax1.set_xlabel("x")
ax1.set_ylabel("u, y")

ax2 = ax1.twinx()
ax2.plot(x, e, "r--", alpha=0.5, label="Skekkja [e]")
ax2.set_ylabel("Skekkja |y−u|", color="r")

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc="upper right")

plt.title("Lausn Við Dæmi 3")
plt.show()

