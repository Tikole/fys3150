import numpy as np
import matplotlib.pyplot as plt

V0 = 1
d = 1
def V(x,y,z):
    return V0/(2*d**2)*(z**2 - x**2 - y**2)

X = np.linspace(-1, 1, 101)
Y = np.linspace(-1, 1, 101)
Z = np.linspace(-1, 1, 101)
xg, zg = np.meshgrid(X,Z)
Vg = V(xg, 0, zg)


plt.contourf(xg, zg, Vg)
plt.axis('scaled')
plt.xlabel("X [d]")
plt.ylabel("Z [d]")
plt.savefig("Vcont.pdf")

xg, yg = np.meshgrid(X, Y)
Vg= V(xg, yg, 0)

plt.contourf(xg, yg, Vg)
plt.axis('scaled')
plt.xlabel("X [d]")
plt.ylabel("Y [d]")
plt.savefig("Vcont_XY.pdf")