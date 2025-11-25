import matplotlib.pyplot as plt
import numpy as np
import math

class DidNotConverge(Exception):
    pass


# einkennismargliða kerfisins er:
def characteristic_polinomial(s,m1,m2,k1,k2,c1,c2):
    t1 = m1*m2*s**4
    t2 = c2*(m1+m2*(c1+c2))*s**3
    t3 = (m1*k2+m2*(k1+k2)+c1*c2)*s**2
    t4 = (c1*k2+c2*k1)*s
    t5 = k1*k2

    return t1+t2+t3+t4+t5

def d_characteristic_polinomial(s,m1,m2,k1,k2,c1,c2):
    t1 = 4*m1*m2*s**3
    t2 = 3*c2*(m1+m2*(c1+c2))*s**2
    t3 = 2*(m1*k2+m2*(k1+k2)+c1*c2)*s**1
    t4 = (c1*k2+c2*k1)

    return t1+t2+t3+t4


# við finnum rót með aðferð newtons
def adferd_newtons(f, df, I_guess, tolerance, itterLimit=1000):
    x = I_guess
    x_old = 999999 # eða bara I_guess+2*tolerance

    ittr = 0
    while abs(x-x_old)>tolerance:
        if ittr>=itterLimit: raise DidNotConverge("newtons aðferð fann ekki lausn innan itterLimit")
        x_old = x
        x = x_old - f(x_old) / df(x_old)
        ittr +=1
    return x


# finnum raungildið af tvinntölu rót af caracteristic margliðunni okkar
def find_damping_rate(c1):

    # skilgreinum föstu inntökin okkar
    m1 = 1
    m2 = 2
    k1 = 10
    k2 = 5
    c2 = 1

    char_polinomial = lambda s: characteristic_polinomial(s, m1,m2,k1,k2,c1,c2)
    d_char_polinomial = lambda s: d_characteristic_polinomial(s, m1,m2,k1,k2,c1,c2)

    tolerance = 1e-10
    # þar sem þetta er vel þekkt línulegt kerfi vitum við að það mun
    # einungis hafa eina tvinntölu pars rót
    innitial_value = 0+1j # reinum að finna tvinntölu rótina
    root = adferd_newtons(char_polinomial, d_char_polinomial, innitial_value, tolerance)
    return np.real(root)


# basic gullsníðar leit
def goldensearch(f,a,b,tol):
    phi=(math.sqrt(5)-1)/2
    x1=a+(1-phi)*(b-a)
    x2=a+phi*(b-a)
    f1=f(x1)
    f2=f(x2)
    while (b-a)/2>tol:
        if f1<f2:
            b=x2
            x2=x1
            x1=a+(1-phi)*(b-a)
            f2=f1
            f1=f(x1)
        else:
            a=x1
            x1=x2
            x2=a+phi*(b-a)
            f1=f2
            f2=f(x2)
    return((a+b)/2)


# finnum lágildi á damping 
c1_min = 5
c1_max = 20
c1_tolerance = 1e-10
c1_minimised = goldensearch(find_damping_rate, c1_min,c1_max, c1_tolerance)


print(f"C1 for minimised damping rate: {c1_minimised:.3f}")

# búum til graf af damping rate vs c1
c1_range_min = 0
c1_range_max = 30
samples = 300

c1_values = np.linspace(c1_range_min, c1_range_max, samples)
damping_values = np.array([find_damping_rate(c1) for c1 in c1_values])

damping_min = find_damping_rate(c1_minimised)

plt.figure(figsize=(10, 6))
plt.plot(c1_values, damping_values, linewidth=2)

# setjum inn lágildið sem við fundum
plt.scatter([c1_minimised], [damping_min], color="red", s=80, zorder=5)

plt.text(
    c1_minimised-4, damping_min+0.05,
    f"Lágildis punkturinn er: [{c1_minimised:.3f}, {damping_min:.4f}]",
    fontsize=11,
    verticalalignment='bottom'
)

# og látum það lúkka næs
plt.title("Damping Rate sem fall af C1", fontsize=16)
plt.xlabel("C1", fontsize=14)
plt.ylabel("Damping Rate", fontsize=14)

plt.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout()
plt.show()





