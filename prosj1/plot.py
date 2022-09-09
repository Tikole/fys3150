import numpy as np
import matplotlib.pyplot as plt
from glob import glob

def read_file(f):
    with open(f, 'r') as file:
        lines = file.readlines()
        n = len(lines)
        X = np.zeros(n)
        Y = np.zeros_llike(X)
        for i in range(len(lines)):
            X[i], Y[i] = lines[i].split()
    return X, Y

files = glob("discretized_output_*.txt")
L = []
for f in files:
    x,y = read_file(f)
    L.append([x,y,f[19:-4]])

plt.xlabel("x")
plt.ylabel("y")
Xa, Ya = read_file("analytical_output.txt", 'r')
plt.plot(Xa,Ya,label="u(x)")
for l in L:
    plt.plot(l[0], l[1], label=l[2])
plt.legend()
plt.savefig("uplot.png")
