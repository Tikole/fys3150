from glob import glob

import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

try:
    N_tridiag = []
    it_tridiag = []
    t_tridiag = []
    with open("scaling_tridiag.txt", 'r') as file:
        lines = file.readlines()
        for l in lines:
            N, it, t = l.split()
            N, it, t = int(N), int(it), float(t)
            N_tridiag.append(N)
            it_tridiag.append(it)
            t_tridiag.append(t)
    N_tridiag = np.array(N_tridiag)
    it_tridiag = np.array(it_tridiag)
    t_tridiag = np.array(t_tridiag)
except:
    print("Couldn't find/open 'scaling_tridiag.txt', or it is ill formed.")

try:
    N_dense = []
    it_dense = []
    t_dense = []
    with open("scaling_dense.txt", 'r') as file:
        lines = file.readlines()
        for l in lines:
            N, it, t = l.split()
            N, it, t = int(N), int(it), float(t)
            N_dense.append(N)
            it_dense.append(it)
            t_dense.append(t)
    N_dense = np.array(N_dense)
    it_dense = np.array(it_dense)
    t_dense = np.array(t_dense)
except:
    print("Couldn't find/open 'scaling_dense.txt', or it is ill formed.")

# Find all files matching 'solution_<n>.txt these are solutions with n steps
# Gather all data in a list of 2-dimensional arrays with one entry per file.
try:
    files = glob("solution_*.txt")
    if not files:
        raise "error"
    Data = []
    for f in files:
        data = []
        with open(f, 'r') as file:
            lines = file.readlines()
            for l in lines:
                x, j1, j2, j3, a1, a2, a3 = l.split();
                data.append(
                    [
                        float(x),
                        float(j1),
                        float(j2),
                        float(j3),
                        float(a1),
                        float(a2),
                        float(a3)
                    ]
                )
            data = np.array(data)
            Data.append(data)
except:
    print("Couldn't find/open any 'solution_*.txt', or some are ill formed.")

# Plot results for scaling with a tridiagonal matrix and a linear regression to quantize
# scaling behavior
if len(N_tridiag):
    plt.figure()
    plt.xlabel("Matrisestørrelse NxN")
    plt.ylabel("Transformasjoner")
    plt.scatter(N_tridiag, it_tridiag, label="Tridiagonal")
    plt.legend()
    plt.savefig("tridiag_scaling.pdf")
    plt.close()

    plt.figure()
    plt.xlabel("Matrisestørrelse NxN")
    plt.ylabel("Tid [s]")
    plt.scatter(N_tridiag, t_tridiag, label="Tridiagonal")
    plt.legend()
    plt.savefig("tridiag_time.pdf")
    plt.close()

    print("Lineærtilpasning av log(transformasjoner) = f(log(N)) for tridiagonal:")
    R = linregress(np.log(N_tridiag), np.log(it_tridiag))
    plt.figure()
    plt.xlabel("log(N)")
    plt.ylabel("log(Transformasjoner)")
    plt.title("linreærtilpasning: stigning = %.2f" % R.slope)
    plt.scatter(np.log(N_tridiag), np.log(it_tridiag))
    plt.axline((0.0, R.intercept), slope=R.slope)
    plt.savefig("linregress.pdf")
    plt.close()

# Plot results for scaling with a dense matrix
if len(N_dense):
    plt.figure()
    plt.xlabel("Matrisestørrelse NxN")
    plt.ylabel("Transformasjoner")
    plt.scatter(N_dense, it_dense, label="Dense")
    plt.legend()
    plt.savefig("dense_scaling.pdf")
    plt.close()

    plt.figure()
    plt.xlabel("Matrisestørrelse NxN")
    plt.ylabel("Tid [s]")
    plt.scatter(N_dense, t_dense, label="Dense")
    plt.legend()
    plt.savefig("dense_time.pdf")
    plt.close()

# Plot results for comparing scaling of dense and sparse matrices
if len(N_tridiag) and len(N_dense):
    plt.figure()
    plt.xlabel("Matrisestørrelse NxN")
    plt.ylabel("Transformasjoner")
    plt.scatter(N_tridiag, it_tridiag, label="Tridiagonal")
    plt.scatter(N_dense, it_dense, label="Dense")
    plt.legend()
    plt.savefig("tridiag_and_dense.pdf")
    plt.close()

    plt.figure()
    plt.xlabel("Matrisestørrelse NxN")
    plt.ylabel("Tid [s]")
    plt.scatter(N_tridiag, t_tridiag, label="Tridiagonal")
    plt.scatter(N_dense, t_dense, label="Dense")
    plt.legend()
    plt.savefig("tridiag_dense_time.pdf")
    plt.close()

# Since a flipped eigenvector is an equally good solution, we must sometimes
# flip one of the vectors, multiply by -1, to compare eigenvectors produced
# by different methods. This method is not very efficient, but leaves nothing
# to chance. If flipping either A or B makes a better match between them it
# returns -1, otherwise 1.
def flipped(A,B):
    a1 = np.mean(np.abs(A-B))
    a2 = np.mean(np.abs(A+B))
    if a1 > a2:
        return -1
    else:
        return 1

if len(Data):
    for D in Data:
        n = len(D[:, 0]) - 1
        plt.figure()
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.plot(D[:, 0], D[:, 4], color='black', label="Analytiske løsninger")
        plt.plot(D[:, 0], D[:, 5], color='black')
        plt.plot(D[:, 0], D[:, 6], color='black')
        i = flipped(D[:, 1], D[:, 4])
        plt.plot(D[:, 0], i*D[:, 1], linestyle='dotted', label="1. Jacobiløsning")
        i = flipped(D[:, 2], D[:, 5])
        plt.plot(D[:, 0], i*D[:, 2], linestyle='dotted', label="2.")
        i = flipped(D[:, 3], D[:, 6])
        plt.plot(D[:, 0], i*D[:, 3], linestyle='dotted', label="3.")
        plt.legend()
        plt.savefig("solution_%d.pdf" % n)
        plt.close()