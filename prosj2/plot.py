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
    print("Couldn't open 'scaling_tridiag.txt'")

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
    print("Couldn't open 'scaling_dense.txt'")

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