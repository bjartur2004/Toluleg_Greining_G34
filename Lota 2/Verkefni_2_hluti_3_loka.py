import numpy as np
import matplotlib.pyplot as plt


t_m = np.array([1,3,5,7,9,11,13,15,17,19,21,23], float)
q_m = np.array([12.0,14.0,13.9,12.3,11.7,8.9,6.3,6.5,6.1,5.8,6.6,9.3], float)

w = 2*np.pi/24

M = np.column_stack([np.cos(w*t_m), np.sin(w*t_m), np.ones_like(t_m)])
A, B, C = np.linalg.lstsq(M, q_m, rcond=None)[0]


def qE0(t):
    return A*np.cos(w*t) + B*np.sin(w*t) + C




r = 0.05
nu = 1e-3
L = 100
QB = 7.0

G = np.pi*r**4/(8*nu*L)

A_p = np.array([
    [3,-1,-1,0,0],
    [1,-2,0,1,0],
    [1,0,-4,1,2],
    [0,3,3,-8,2],
    [0,0,6,2,-11]
], float)


def þrystingur(p1):
    b = np.array([p1, -QB/G, 0, 0, 0], float)
    return np.linalg.solve(A_p, b)


def q(p1):
    p = þrystingur(p1)
    return G*p[4]


def find_p1(q_target, a=0, b=1e7, tol=1e-10):
    fa = q(a) - q_target
    fb = q(b) - q_target

    for _ in range(100):
        m = 0.5*(a+b)
        fm = q(m) - q_target
        if abs(fm) < tol:
            return m
        if fa*fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5*(a+b)



t_100 = np.linspace(0, 24, 100)
q_allt = qE0(t_100)

p1_allt = np.zeros_like(t_100)
for i,t in enumerate(t_100):
    p1_allt[i] = find_p1(q_allt[i])




fig, ax1 = plt.subplots(figsize=(10,5))

ax1.plot(t_100, q_allt, label="qE0(t)", color="blue")
ax1.set_xlabel("Tími [klst]")
ax1.set_ylabel("qE0 [m³/s]", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")

ax2 = ax1.twinx()
ax2.plot(t_100, p1_allt, label="p1(t)", color="red")
ax2.set_ylabel("p1 [Pa]", color="red")
ax2.tick_params(axis="y", labelcolor="red")

plt.title("qE0(t) og samsvarandi p1(t)")
plt.tight_layout()
plt.show()
midgildi = len(t_100)//2

print("\n--- Niðurstöður ---")
print(f"Fyrsta gildi :  t = {t_100[0]:.2f}   q = {q_allt[0]:.4f}   p1 = {p1_allt[0]:.4f}")
print(f"Miðju gildi  :  t = {t_100[midgildi]:.2f}   q = {q_allt[midgildi]:.4f}   p1 = {p1_allt[midgildi]:.4f}")
print(f"Síðasta gildi:  t = {t_100[-1]:.2f}   q = {q_allt[-1]:.4f}   p1 = {p1_allt[-1]:.4f}")
