import numpy as np r = 0.05
L = 100
nu = 1.0e-3
p1 = 4.2e6
QB = 7.0

# Compute G
G = np.pi * r**4 / (8 * nu * L)

# Matrix A (5×5)
A = np.array([
    [3, -1, -1,  0,  0],       # Eq (8)
    [1, -2,  0,  1,  0],       # Eq (9)
    [1,  0, -4,  1,  2],       # Eq (10)
    [0,  3,  3, -8,  2],       # Eq (11)
    [0,  0,  6,  2, -11]       # Eq (12)
], dtype=float)

# RHS vector b
b = np.array([
    p1,
    -QB / G,
    0,
    0,
    0
], dtype=float)

# Solve the linear system
p = np.linalg.solve(A, b)

pA, pB, pC, pD, pE = p

print("pA =", pA)
print("pB =", pB)
print("pC =", pC)
print("pD =", pD)
print("pE =", pE)

# Once pressures are known, compute flows using q = G_ij*(pi - pj)
# GCE = 2G, GDE = 2G/3
def q(pi, pj, Gij):
    return Gij * (pi - pj)

q1A = G*(p1 - pA)
qAB = G*(pA - pB)
qAC = G*(pA - pC)
qBD = G*(pB - pD)
qCD = G*(pC - pD)
qCE = 2*G*(pC - pE)
qDE = (2/3)*G*(pD - pE)
qE0 = G*(pE - 0)

print("\nFlows q:")
print("q1A =", q1A)
print("qAB =", qAB)
print("qAC =", qAC)
print("qBD =", qBD)
print("qCD =", qCD)
print("qCE =", qCE)
print("qDE =", qDE)


print("qE0 =", qE0)
