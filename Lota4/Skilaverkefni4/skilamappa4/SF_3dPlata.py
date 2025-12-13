import numpy as np
import pyvista as pv
import HitajofnuHneppi_3d

P = 50
Lx = 10  # cm
Ly = 10  # cm
Lz = 1  # cm
h = 0.25

K = 1.68
H = 0.005
delta = 0.1  # cm
Power_Lx = 5  # cm
Power_Ly = 2  # cm
ambient_temperature = 20  # °C
arguments = [P, h, K, H, delta, Power_Lx, Power_Ly, ambient_temperature]


n = int(round(Lx/h))
m = int(round(Ly/h))
p = int(round(Lz/h))

u = HitajofnuHneppi_3d.solve_u(n, m, p, arguments)
U = u.reshape((p, n, m))

U = np.swapaxes(U, 0, 2)

# Setja í pyvista grid
x = np.linspace(0, Lx, m)
y = np.linspace(0, Ly, n)
z = np.linspace(0, Lz, p)

grid = pv.StructuredGrid(*np.meshgrid(x, y, z, indexing='ij'))

grid["Hitastig [°C]"] = U.flatten(order="F")

# Plot
t_min = U.min()
t_max = U.max()

plotter = pv.Plotter()


def s_curve(x, k=2.0, x0=0.0):
    x = np.asarray(x, dtype=float)
    t = np.clip(x - x0, 0.0, 1.0)
    return t**k / (t**k + (1 - t)**k)

k = 2.2
x0 = 0.2
c=0.05
opacity = np.array([min(s_curve(x,k,x0)+c, 1.0) for x in np.linspace(0,1,10)])

plotter.add_volume(grid, scalars="Hitastig [°C]",
                   cmap="plasma", opacity=opacity,
                   shade=True,
                   clim=[t_min, t_max])

plotter.add_mesh(grid.outline(), color="black", line_width=2)
plotter.add_axes()
plotter.set_background("white")
plotter.add_axes()
plotter.add_text("Hitadreifing í þunnum ferstrendingi með sama rúmmál", position='upper_edge', font_size=15, color='black')
plotter.show()
