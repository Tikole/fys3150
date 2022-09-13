import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt

file_ending = "pdf"

def print_usage():
    print(
        """Plots the data contained in 'output_<eval_points>.txt' files.\n
        Arguments must be provided. Avaiable:\n
        * 'exact' : use the highest resolution available to plot the exact
          solution.
        * 'general' : plots all available solutions calculated with the general
          solution, and the highest resolution data available for the exact
          solution.
        * 'abs_error' : plots the absolute error in the numerical solution
          against x for all available number of steps.
        * 'rel_error' : plots the error the numerical solution  divided by the
          analytical solution for all available number of steps. Also prints
          max relative error for all number of steps.
        """
    )

def read_file(filename):
    """
    Reads output files into 4 vectors containing X-values and corresponding
    Y-values for analytic and numeric solutions.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
        n = len(lines)
        X = np.zeros(n)
        Ya = np.zeros_like(X)
        Yg = np.zeros_like(X)
        Ys = np.zeros_like(X)
        for i in range(len(lines)):
            X[i], Ya[i], Yg[i], Ys[i] = lines[i].split()
    return X, Ya, Yg, Ys

def find_highest_resolution(L):
    """
    Finds the index of the data set in L with the highest number of data points.
    """
    N = 0 # Highest resolution seen
    j = None # Index of highest resolution
    for i in range(len(L)):
        n = L[i][0] # Resolution
        if n > N:
            j = i
            N = n
    return j

def plot_exact(L):
    """
    Use higest resolution data in L to make a plot of the exact solution.
    """
    j = find_highest_resolution(L)
    plt.figure()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.plot(L[j][1], L[j][2], label="u(x)")
    plt.legend()
    plt.savefig("exact(%d).%s" % (file_ending))
    plt.close()

def plot_general(L):
    """
    Plots all the available solutions calculated with the general algorithm
    and the highest resolution available for the exact solution.
    """
    j = find_highest_resolution(L)
    plt.figure()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.plot(L[j][1], L[j][2], label="Exact", linestyle='dashed')
    for l in L:
        n = l[0]
        plt.plot(l[1], l[3], label="n=%d" % n)
    plt.legend()
    plt.savefig("general_solver.%s" % (file_ending))

def plot_abs_error(L):
    """
    Plots the absolute error of the general algorithm when compared to the
    analytical solution, taken to be exact.
    """
    plt.figure()
    plt.xlabel("X")
    plt.ylabel("Log10 (absolute error)")
    for i in range(len(L)):
        X = L[i][1]
        u = L[i][2]
        v = L[i][3]
        n = L[i][0]
        Y = np.log10(np.absolute(u[1:-1] - v[1:-1]))
        plt.plot(X, Y, label="n=%d" % n)
    plt.legend()
    plt.savefig("abs_error.%s" % file_ending)

def plot_rel_error(L):
    """
    Plots the relative error of the general algorithm when compared to the
    analytical solution, taken to be exact.
    """
    plt.figure()
    plt.xlabel("X")
    plt.ylabel("Log10 of relative error")
    for i in range(len(L)):
        X = L[i][1]
        u = L[i][2]
        v = L[i][3]
        n = L[i][0]
        Y = np.log10(np.absolute((u[1:-1] - v[1:-1])/u[1:-1]))
        max = np.max(Y)
        print("n=%d : %.2f" % (n, max))
        plt.plot(X[1:-1],Y, label="n=%d" % n)
    plt.legend()
    plt.savefig("rel_error.%s" % file_ending)

# Read all files 
files = glob("output_*.txt")
L = []
for f in files:
    X,Ya,Yg,Ys = read_file(f)
    L.append([int(f[7:-4]), X,Ya,Yg,Ys])

if len(sys.argv) == 1:
    print_usage()


# Which plots are wanted
else:
    if len(L) == 0:
        print("No data found.")
    else:
        if "exact" in sys.argv:
            plot_exact(L)
        if "general" in sys.argv:
            plot_general(L)
        if "abs_error" in sys.argv:
            plot_abs_error(L)
        if "rel_error" in sys.argv:
            plot_rel_error(L)
        

