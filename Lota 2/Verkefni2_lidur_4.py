import numpy as np
import math

draw_graph = True
if draw_graph: # fyrir yfirfarara ef þið viljið keyra án þess að hafa auka graf kóðann:)
    from pipeSystemGraf import plot_pressure_flow_graph

QB = 7
p1 = 4.2*10**6
p0 = 0

p = 998
L = 100
f = 0.02
r = 0.05
K = (8*f*p*L)/(32*math.pi**2 * r**5)



def F(X):
    global K, p1, p0, QB
    q1A, qAB, qAC, qBD, qCD, qCE, qDE, qE0 = X
    F = np.array([
        q1A - qAB - qAC,
        QB  + qAB - qBD,
        qAC - qCD - qCE,
        qBD + qCD - qDE,
        qCE + qDE - qE0,
        qAB*abs(qAB) - qAC*abs(qAC) + qBD*abs(qBD) - qCD*abs(qCD),
        qCD*abs(qCD) - 0.5*qCE*abs(qCE) + 1.5*qDE*abs(qDE),
        (p0-p1)/K + q1A*abs(q1A) + qAC*abs(qAC) + 0.5*qCE*abs(qCE) + qE0*abs(qE0)
    ])
    return F



def J(X):
    q1A, qAB, qAC, qBD, qCD, qCE, qDE, qE0 = X

    A1 = np.array([
        [1,-1,-1,0,0,0,0,0],
        [0,1,0,-1,0,0,0,0],
        [0,0,1,0,-1,-1,0,0],
        [0,0,0,1,1,0,-1,0],
        [0,0,0,0,0,1,1,-1]
    ])

    A2 = np.array([
        [0,2*abs(qAB), -2*abs(qAC), 2*abs(qBD), -2*abs(qCD), 0,0,0],
        [0,0,0,0,2*abs(qCD), -1*abs(qCE), 3*abs(qDE), 0],
        [2*abs(q1A), 0,2*abs(qAC), 0,0, abs(qCE), 0, 2*abs(qE0)]
    ])

    J = np.block([[A1], [A2]])

    return J

def newtonSearch(F, J, x0, tolerance):
    x = np.array(x0, dtype=float)
    x_old = np.array([999999 for i in range(8)], dtype=float)

    while np.linalg.norm(x-x_old)>tolerance:
        x_old = x
        x = x_old - np.linalg.solve(J(x), F(x))
    return x


x0 = [1 for i in range(8)]
tolerance = 1e-10

q_lausn = newtonSearch(F,J,x0, tolerance)
q1A, qAB, qAC, qBD, qCD, qCE, qDE, qE0 = q_lausn

pA = p1 - K*q1A*abs(q1A)
pB = pA - K*qAB*abs(qAB)
pC = pA - K*qAC*abs(qAC)
pD = pC - K*qCD*abs(qCD)
pE = pC - K*0.5*qCE*abs(qCE)


# OUTPUT:
flows = [
    ("1", "A", q1A),
    ("A", "B", qAB),
    ("A", "C", qAC),
    ("B", "D", qBD),
    ("C", "D", qCD),
    ("C", "E", qCE),
    ("D", "E", qDE),
    ("E", "0", qE0)
]

# Pressures at nodes in MPa
pressures = {
    "1": p1*1e-6,
    "A": pA*1e-6,
    "B": pB*1e-6,
    "C": pC*1e-6,
    "D": pD*1e-6,
    "E": pE*1e-6,
    "0": p0*1e-6
}

print("\n=== FLOW SOLUTION (m³/s) ===")
print(f"{'Edge':<7} | {'Flow':>12}")
print("-"*22)
for u, v, q in flows:
    print(f"{u}->{v:<4} | {q:12.6f}")

print("\n=== PRESSURE SOLUTION (MPa) ===")
print(f"{'Node':<5} | {'Pressure':>12}")
print("-"*22)
for node, p in pressures.items():
    print(f"{node:<5} | {p:12.6f}")

if draw_graph:
    positions = {
        "1": (0, 2),
        "A": (1, 2),
        "B": (2, 2),
        "C": (1, 1),
        "D": (2, 1),
        "E": (1, 0.5),
        "0": (0, 0.5)
    }
    plot_pressure_flow_graph(flows, pressures, positions, node_cmap="bright_plasma")

