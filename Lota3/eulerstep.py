from ydot import ydot
def eulerstep(t, x, h):
    return x + h * ydot(t, x)