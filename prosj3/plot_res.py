import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
with open(filename, 'r') as file:
    Aline = file.readline()
    A = [float(w) for w in Aline.split()]
    fline = file.readline()
    f = [float(f) for f in fline.split()]
    D = np.zeros(shape=(len(f), len(A)))
    lines = file.readlines()
    for i in range(len(f)):
        D[i,:] = np.array([int(w) for w in lines[i].split()])

for i in range(D.shape[1]):
    plt.plot(f, D[:,i], label="A=%g" % A[i])
plt.legend()
plt.xlabel("frekvens [rad/us]")
plt.ylabel("Partikkler fanget")
plt.savefig("res.pdf")