import matplotlib.pyplot as mpl
import numpy as np

def f(x):
    return 2*x**4 + 23*x**3 + 45*x**2 + 60*x + 50

def Df(x):
    return 8*x**3 + 69*x**2 + 90*x + 60


xn = [] # býr til töflu fyrir X_n gildunum. 
def adferd_newtons(Gisk, Nákvæmni):
    x1 = Gisk
    x2 = x1 + 10 * Nákvæmni   
    while abs(x1 - x2) > Nákvæmni:
        xn.append(x1)
        x2 = x1
        x1 = x2 - f(x2) / Df(x2)
        #print(x1) notað til að sjá gildin í hverju skrefi fyrir sig.
    return x1


Gisk = 5 + 7j 
Nákvæmni = 1e-10

root = adferd_newtons(Gisk, Nákvæmni)
print("Niðurstaða Newtons:", root) 

en = [np.linalg.norm(xi-root) for xi in xn] #tafla fyrir en 
print(en)
print(xn)
