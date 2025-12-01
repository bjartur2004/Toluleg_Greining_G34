import numpy as np
import matplotlib.pyplot as plt
# Gefin gögn
t = np.array([1,3,5,7,9,11,13,15,17,19,21,23], dtype=float)
q = np.array([12.0, 14.0, 13.9, 12.3, 11.7, 8.9, 6.3, 6.5, 6.1, 5.8, 6.6, 9.3], dtype=float)

# Tíðni
w = 2 * np.pi / 24

# Fylkið M (12x3)
M = np.column_stack([
    np.cos(w * t),
    np.sin(w * t),
    np.ones_like(t)
])

# komum ofákvörðuðu fylki yfir í ferhyrnt kerfi
x, _, _, _= np.linalg.lstsq(M, q, rcond=None)
A, B, C = x

print("A",A)
print("B",B)
print("C",C)


t_fine = np.linspace(0, 25, 300)
q_fit = A * np.cos(w * t_fine) + B * np.sin(w * t_fine) + C


plt.figure(figsize=(8,5))
plt.plot(t, q, "o")
plt.plot(t_fine, q_fit)
plt.xlabel("Tími (klst)")
plt.ylabel("Flæði $q_{E0}(t)$")
plt.title("Nálgun á mæligögnum með aðferð minnstu kvaðrata")
plt.show()