from functools import partial
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np


F = h5py.File("slitex.hdf5", 'r')
N = F['N'][0][0]
M = F['M'][0][0]
T = F['T'][0][0]
Dt = F['Dt'][0][0]
X = np.linspace(0.0, 1.0, M)
S_real = np.reshape(np.array(F['S_real']), newshape=(N,M,M))
S_imag = np.reshape(np.array(F['S_imag']), newshape=(N,M,M))
P = S_real**2 + S_imag**2
V = F['V']

def plot_potential(ax=None):
    """ Plots the potential in the provided axes. If None plot in a standalone
    plot.
    """
    b = False
    if not ax:
        b = True
        fig, ax = plt.subplots()
    ax.contour(X,X, V[:,:], cmap='Reds')
    if b:
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        fig.savefig("potential.pdf")

# class Animation:
#     def __init__(self):
#         S_real = np.array(F['S_real'])

def plot_p(t=0):
    """ Plots probability plot at time t"""
    i = round(t/Dt)
    if i >= N:
        i = N-1
    p_max = np.max(P)
    fig, ax = plt.subplots()
    plot_potential(ax)
    ax.contourf(X,X, P, levels=100, cmap="inferno")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    fig.savefig("P_%d.png" % i)

def plot_prob_t():
    """ Plots total probability as function of time."""
    Psum = np.sum(P, axis=1)
    Psum = 1 - Psum
    X = np.arange(0, len(Psum))
    plt.figure()
    plt.plot(X, Psum)
    plt.xlabel("Tid")
    plt.ylabel("Total sannsynlighet")
    plt.savefig("Pt.png")

def plot_times(n):
    """Make probalitiy plots at n evenly spaced times including endpoints"""
    times = np.linspace(0.0, T, n)
    for t in times:
        plot_p(t)

def ani_func(frame, ax, data):
    print("\r               \r", end='')
    print("Frame: %d/%d" % (frame, N), end='', flush=True)
    data.collections = []
    data = ax.contourf(X,X,P[frame],levels=100, cmap="inferno")

def animate():
    """Creates an animation of the whole simulation"""
    fig, ax = plt.subplots()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.contour(X,X,V,cmap="Reds")
    data = ax.contourf(X,X,P[0],levels=100,cmap="inferno")
    print("Animating...")
    ani = FuncAnimation(fig, partial(ani_func, ax=ax, data=data), frames=N, blit=False)
    ani.save("animation.mp4", writer="ffmpeg", bitrate=-1, fps=30)

animate()

# plot_prob_t()
# plot_times(20)