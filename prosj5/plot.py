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
        plt.close(fig)

def plot_detection(t=0.002):
    """ Plot probability along the y-axis at x=0.8 at time t"""
    i = round(t/Dt)
    j = int(0.8*(M-1))
    if i >= N:
        i = N-1
    fig, ax = plt.subplots()
    Y = np.linspace(-0.5, 0.5, M)
    Pdet = P[i,:,j]
    Pdet /= np.sqrt(np.sum(Pdet**2))
    ax.plot(Y, Pdet)
    ax.set_xlabel("Y")
    ax.set_ylabel("P")
    fig.savefig("DetP.pdf")
    

def plot_p(t=0):
    """ Plots probability plot at time t"""
    i = round(t/Dt)
    if i >= N:
        i = N-1
    fig, ax = plt.subplots()
    plot_potential(ax)
    ax.contourf(X,X, P[i], levels=10, cmap="inferno")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    fig.savefig("P_%d.pdf" % i)
    plt.close(fig)

def plot_RnI(t=0):
    """ Plots the magnitude of the real and imaginary components of wavefunction
    at time t"""
    i = round(t/Dt)
    if i >= N:
        i = N-1
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    plot_potential(ax1)
    ax1.contourf(X,X, S_real[i], levels=11, cmap="plasma")
    ax1.set_xlabel("X")
    ax1.set_xlabel("Y")
    fig1.savefig("R_%d.pdf" % i)
    plt.close(fig1)

    plot_potential(ax2)
    ax2.contourf(X,X, S_imag[i], levels=11, cmap="seismic")
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    fig2.savefig("I_%d.pdf" %i)
    plt.close(fig2)

def plot_prob_t():
    """ Plots deviation of total probability from 1 as function of time."""
    Psum = np.sum(P, axis=(1,2))
    Psum = 1 - Psum
    X = np.arange(0, len(Psum))
    plt.figure()
    plt.plot(X, Psum)
    plt.xlabel("Steg")
    plt.ylabel("1 - P(t)")
    plt.savefig("P_dev.pdf")
    plt.close()

def plot_times(n):
    """Make probalitiy plots and real and imaginary at n evenly spaced times
    including endpoints"""
    times = np.linspace(0.0, T, n)
    for t in times:
        plot_p(t)
        plot_RnI(t)

def ani_func(frame, ax, data):
    print("\r               \r", end='')
    print("Frame: %d/%d" % (frame, N), end='', flush=True)
    data.collections = []
    data = ax.contourf(X,X,P[frame],levels=10, cmap="inferno")

def animate():
    """Creates an animation of the whole simulation"""
    fig, ax = plt.subplots()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.contour(X,X,V,cmap="Reds")
    data = ax.contourf(X,X,P[0],levels=10,cmap="inferno")
    print("Animating...")
    ani = FuncAnimation(fig, partial(ani_func, ax=ax, data=data), frames=N, blit=False)
    ani.save("animation.mp4", writer="ffmpeg", bitrate=-1, fps=30)
    print()

plot_prob_t()
plot_times(16)
plot_detection()
# animate()