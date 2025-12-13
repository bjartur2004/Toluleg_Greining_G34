import HitadreifingsPlot
import HitaJofnuhneppi_Sparce
import numpy as np
import matplotlib.pyplot as plt



def LeysaCase(case):
    global n, m
    Power_L = case[0]
    power_b = case[1]

    P = 5
    L = 4 #cm
    K = 1.68
    H = 0.005
    delta = 0.1 #cm
    ambient_tempeature = 20 #°C

    arguments = [P,L,K,H,delta,Power_L,ambient_tempeature]

    u = HitaJofnuhneppi_Sparce.solve_u(n,m,arguments, PowerAreaBottom=power_b)
    return u

def FinnaMaxTemp(case):
    u = LeysaCase(case)
    return max(u)


def root_binary_search(f, I_MIN, I_MAX, tolerance):
    min_bound = I_MIN
    max_bound = I_MAX

    f_min = f(min_bound)
    f_max = f(max_bound)

    while (max_bound - min_bound)/2 > tolerance:
        m = (max_bound+min_bound)/2
        f_m = f(m)

        if f_m * f_max > 0: # ef lausnin er ekki á milli m og max
            max_bound = m
            f_max = f(max_bound)

        else:
            min_bound = m
            f_min = f(min_bound)

    # solution found
    return (min_bound+max_bound)/2

n = 50
m = 50
maxTemp = 50
f = lambda L: FinnaMaxTemp([L,0])-maxTemp

maxL = root_binary_search(f, 0, 4, 1e-2)

print(maxL)