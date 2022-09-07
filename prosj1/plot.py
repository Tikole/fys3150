import numpy as np
import matplotlib.pyplot as plt

with open("output.txt", 'r') as file:
    lines = file.readlines()
    n = len(lines)
    X = np.zeros(n)
    Y = np.zeros_like(X)
    for i in range(len(lines)):
        X[i], Y[i] = lines[i].split()

plt.xlabel("x")
plt.ylabel("y")
plt.plot(X,Y,label="u(x)")
plt.legend()
plt.savefig("uplot.png")
