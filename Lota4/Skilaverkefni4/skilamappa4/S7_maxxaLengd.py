import HitadreifingsPlot
import HitaJofnuhneppi_Sparce
import numpy as np
import matplotlib.pyplot as plt


def LeysaCase(case):
    global n, m
    P = case
    Power_L = 2          
    L = 4                
    K = 1.68
    H = 0.005
    delta = 0.1         
    ambient_tempeature = 20  

    arguments = [P, L, K, H, delta, Power_L, ambient_tempeature]

    u = HitaJofnuhneppi_Sparce.solve_u(
        n, m, arguments,
        PowerAreaBottom=(L-Power_L)/2,
    )
    return u


def FinnaMaxTemp(case):
    u = LeysaCase(case)
    return max(u)


def root_binary_search(f, I_MIN, I_MAX, tolerance):
    min_bound = I_MIN
    max_bound = I_MAX

    f_min = f(min_bound)
    f_max = f(max_bound)

    while (max_bound - min_bound) / 2 > tolerance:
        m = (max_bound + min_bound) / 2
        f_m = f(m)

        if f_m * f_max > 0:
            max_bound = m
            f_max = f(max_bound)
        else:
            min_bound = m
            f_min = f(min_bound)

    return (min_bound + max_bound) / 2


n = 50
m = 50
maxTemp = 100

f = lambda P: FinnaMaxTemp(P) - maxTemp

maxP = root_binary_search(f, 0, 20, 1e-2)

print(maxP)
