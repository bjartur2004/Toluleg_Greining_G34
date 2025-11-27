import numpy as np




def F(X):
    q1A, qAB, qAC, qBD, qCD, qCE, qDE, qE0 = X
    QB = 1

    F = np.array([
        [q1A - qAB - qAC],
        [QB + qAB - qBD],
        [qAC - qCD - qCE],
        [qBD + qCE - qDE],
        [qCE + qDE - qE0],
        [qAB * abs(qAB) - ]
    ])



def J(X):
    q1A, qAB, qAC, qBD, qCD, qCE, qDE, qE0 = X

    A1 = np.array([
        [1,-1,-1,0,0,0,0,0],
        [0,1,0,-1,0,0,0,0],
        [0,0,1,0,-1,-1,0,0],
        [0,0,0,1,1,-1,0,0],
        [0,0,0,0,0,1,1,-1]
    ])

    A2 = np.array([
        [0,2*abs(qAB), -2*abs(qAC), 2*abs(qBD), -2*abs(qCD), 0,0,0],
        [0,0,0,0,2*abs(qCD), -1*abs(qCE), 3*abs(qDE), 0],
        [2*abs(q1A), 0,2*abs(qAC), 0,0, abs(qCE), 0, 2*abs(qE0)]
    ])

    J = np.block([[A1], [A2]])

    return J

def newtonSearch(F, J, x0, treshold):
    x = x0
    x_old = 99999

    while np.linalg.norm(x-x_old)>tolerance:
        x_old = x
        x = x_old - np.linalg.solve(J(x), F(x))
        print(f"{x:.3f}")   # prentar hvert skref
    return x





