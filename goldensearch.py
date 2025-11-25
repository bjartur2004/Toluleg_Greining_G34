import math
import matplotlib.pyplot as mpl

def goldensearch(a,b,tol):
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

a=0
b=30
tol = 10**-3

def f(s, c1):
    c1 = input("") 
    return (
        2*s**4
        + 3*s**3
        + 2*c1*s**3
        + 35*s**(c1 + 2)
        + 15*s
        + 50
    )
print(goldensearch)



                
