def f(x):
    return 2*x**4 + 23*x**3 + 45*x**2 + 60*x + 50

def bisection(a, b, tol):
    while (b - a)/2 > tol:
        c = (a + b)/2
        print(f"c = {c:.6f}")   # prentar aðeins nálgun rótarinnar

        if f(a) * f(c) < 0:
            b = c
        else:
            a = c

    return c

# rót ~ -9.42
root = bisection(-10, -5, 1e-6)
print("Rót fundin:", root)

