import numpy as np

# hér erum við með fastana sem voru gefnir í dæminu 
p1 = 4.2e6      # Pa
QB = 7.0        # m^3/s
V = 1e-3       # Pa·s
r = 0.05        # m
L = 100         # m

#skilgreinum G 
G = (np.pi * r**4) / (8 * V * L)

# fylki A sem eru stuðlar af þr
A = np.array([
    [ 3, -1, -1,  0,  0],
    [ 1, -2,  0,  1,  0],
    [ 1,  0, -4,  1,  2],
    [ 0,  3,  3, -8,  2],
    [ 0,  0,  6,  2, -11]
])

b = np.array([
    p1,
    -QB / G,
    0,
    0,
    0
])
 
p = np.linalg.solve(A, b) # leysir fylkið

pA, pB, pC, pD, pE = p
print("pA =", pA)
print("pB =", pB)
print("pC =", pC)
print("pD =", pD)
print("pE =", pE)

G1A = GAB = GAC = GBD = GCD = GE0 = G #þar sem G er 2,54 * 10^-5
GCE = 2 * G 
GDE = (2/3) * G

p0 = 0.0  

# Flæði milli nóða

q1A = G1A * (p1 - pA)
qAB = GAB * (pA - pB)
qAC = GAC * (pA - pC)
qBD = GBD * (pB - pD)
qCD = GCD * (pC - pD)
qCE = GCE * (pC - pE)
qDE = GDE * (pD - pE)
qE0 = GE0 * (pE - p0)

flows = {
    "q1A": q1A,
    "qAB": qAB,
    "qAC": qAC,
    "qBD": qBD,
    "qCD": qCD,
    "qCE": qCE,
    "qDE": qDE,
    "qE0": qE0
}


print("\nRúmmálsflæði (m^3/s):")
for name, q in flows.items():
    print(name, "=", q)
