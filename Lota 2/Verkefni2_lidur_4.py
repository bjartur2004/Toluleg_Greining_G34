import numpy as np
import math

draw_graph = True
if draw_graph: # fyrir yfirfarara ef þið viljið keyra án þess að hafa auka graf kóðann:)
    from pipeSystemGraf import plot_pressure_flow_graph

def F(X):
    q1A, qAB, qAC, qBD, qCD, qCE, qDE, qE0 = X
    QB = 7
    p1 = 4.2*10**6
    p0 = 0

    p = 998
    L = 100
    f = 0.02
    r = 0.05
    K = (8*f*p*L)/(32*math.pi**2 * r**5)

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
        #print("  ".join(f"{xi:.3f}" for xi in x))   # prentar hvert skref
    return x


x0 = [1 for i in range(8)]
tolerance = 1e-10

x_lausn = newtonSearch(F,J,x0, tolerance)
F_lausn = F(x_lausn)


print("Lausnin er:", "  ".join(f"{xi:.3f}" for xi in x_lausn))

if draw_graph:
    # Extract solved flows
    q1A, qAB, qAC, qBD, qCD, qCE, qDE, qE0 = x_lausn

    # Define edges (direction = positive direction)
    edges = [
        ("1", "A", q1A),
        ("A", "B", qAB),
        ("A", "C", qAC),
        ("B", "D", qBD),
        ("C", "D", qCD),
        ("C", "E", qCE),
        ("D", "E", qDE),
        ("E", "0", qE0)
    ]

    # Compute pressures using the energy balance equation
    # p1 and p0 must match your F() function
    p1 = 4.2 * 10**6
    p0 = 0

    pressures = {
        "1": p1,
        "A": p1 - q1A*abs(q1A),
        "B": p1 - q1A*abs(q1A) - qAB*abs(qAB),
        "C": p1 - q1A*abs(q1A) - qAC*abs(qAC),
        "D": p1 - q1A*abs(q1A) - qAC*abs(qAC) - qCD*abs(qCD),
        "E": p1 - q1A*abs(q1A) - qAC*abs(qAC) - 0.5*qCE*abs(qCE),
        "0": p0
    }

    # Node layout for plotting
    positions = {
        "1": (0, 2),
        "A": (1, 2),
        "B": (2, 2),
        "C": (1, 1),
        "D": (2, 1),
        "E": (1, 0),
        "0": (0, 0)
    }

    # Call your plotting function
    plot_pressure_flow_graph(edges, pressures, positions)
