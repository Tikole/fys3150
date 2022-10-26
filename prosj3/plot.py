import matplotlib.pyplot as plt
import numpy as np
import sys
from glob import glob

"""Open and read all files on form "PT_*.txt". Extract data into a list of 
dictionaries called 'data'.
"""
data = []
filenames = glob("PT*.txt")
if not filenames:
    print("ERROR: No data files found.")
    sys.exit()
for fn in filenames:
    # Extract data from filename
    fn_words = fn[3:-3].split('_')
    # Coulomb interactions?
    coulomb = bool(int(fn_words[3]))
    # Integration method
    method = fn_words[4]
    with open(fn, 'r') as file:
        desc = file.readline()[:-1]
        file.readline() # Discard second line
        lines = file.readlines()
        steps = len(lines) - 1
        # Count number of particles
        words = lines[0].split()
        # 4 fields describe trap, 6 fields for each particle
        n_particles = (len(words) - 4)//6
        # Process data lines
        t = np.zeros(steps + 1)
        B = np.zeros(steps + 1)
        Vo = np.zeros(steps + 1)
        d = np.zeros(steps + 1)
        r = np.zeros(shape=(steps+1, n_particles, 3))
        v = np.zeros(shape=(steps+1, n_particles, 3))
        for i in range(len(lines)):
            w = lines[i].split()
            t[i] = float(w[0])
            B[i] = float(w[1])
            Vo[i] = float(w[2])
            d[i] = float(w[3])
            for j in range(n_particles):
                r[i,j] = np.array([float(w[6*j + 4]), float(w[6*j + 5]), float(w[6*j + 6])])
                v[i,j] = np.array([float(w[6*j + 7]), float(w[6*j + 8]), float(w[6*j + 9])])
        data.append(
            {
                "coulomb": coulomb,
                "method": method,
                "desc": desc,
                "t": t,
                "B": B,
                "Vo": Vo,
                "d": d,
                "r": r,
                "v": v
            }
        )

def plot_zt(D):
    """
    Plot z(t) for data object D. 
    """
    for i in range(D['r'].shape[1]):
        plt.plot(D['t'], D['r'][:,i,2])
    plt.xlabel("tid [us]")
    plt.ylabel("z [um]")
    plt.savefig("z(t).pdf")

def plot_xy(D):
    """
    Plot trajectory in xy-plane of data objects in D.
    """
    for i in range(D['r'].shape[1]):
        rx = D['r'][:,i,0]
        ry = D['r'][:,i,1]
        l = plt.plot(rx, ry)
        c = l[0].get_color()
        plt.scatter(rx[0], ry[0], marker='X', color=c)
        plt.scatter(rx[-1], ry[-1], marker='o', color=c)
    plt.xlabel("x [um]")
    plt.ylabel("y [um]")
    plt.axis('equal')
    plt.savefig("xy-plane.pdf")

def plot_phase(D, ax):
    """
    Plot movement in phase space of data object D in direction given by ax
    ax=0 for x-axis, ax=1 for y-axis, ax=2 for z-axis.
    """
    for i in range(D['r'].shape[1]):
        rx = D['r'][:,i,ax]
        vx = D['v'][:,i,ax]
        l = plt.plot(rx,vx)
        c = l[0].get_color()
    plt.xlabel("X [um]")
    plt.ylabel("Vx [um/us")
    plt.savefig("phase%d.pdf" % ax)

def plot_3d_trajectory(D):
    """
    Plots 3d trajectories of particles of data object D
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i in range(D['r'].shape[1]):
        rx = D['r'][:,i,0]
        ry = D['r'][:,i,1]
        rz = D['r'][:,i,2]
        ax.plot(rx,ry,rz)
    fig.savefig("threeD.pdf")    

class Analytic:
    """
    Function object for solving analytical case
    """
    def __init__(self, x0,z0,v0,q,V,B,m,d):
        self.z0 = z0
        self.om_z = np.sqrt(2*q*V/(m*d*d))
        self.om_0 = q*B/m
        self.om_m = 1/2 * (self.om_0 - np.sqrt(self.om_0**2 - 2*self.om_z**2))
        self.om_p = 1/2 * (self.om_0 + np.sqrt(self.om_0**2 - 2*self.om_z**2))
        self.A_p = (v0 + self.om_m*x0)/(self.om_m - self.om_p)
        self.A_m = -(v0 + self.om_p*x0)/(self.om_m - self.om_p)

    def __call__(self, t):
        xy = self.A_p*np.exp(-1j * self.om_p*t) + self.A_m*np.exp(-1j * self.om_m*t)
        z = self.z0*np.cos(self.om_z*t)
        res = np.array([xy.real, xy.imag, z])
        if type(t) == np.ndarray:
            res = np.array([xy.real, xy.imag, z]).transpose()
        return res

ana = Analytic(20.0, 20.0, 25.0, 1.0, 25*9.64852558e4, 1.0*9.64852558e1, 40.078, 500)    

def plot_analytic_rel_error(*D):
    """
    Plots relative errors of objects in compared to the exact solution of the analytical case.
    """
    D = sorted(D, key=lambda D: len(D['t']) - 1)
    max_err = np.zeros(len(D))
    h = np.zeros(len(D))
    for i in range(len(D)):
        d = D[i]
        t = d['t']
        r = d['r']
        exact = ana(t)
        err = r[:,0,:]-exact
        err_mag = np.sqrt(err[:,0]**2 + err[:,1]**2 + err[:,2]**2)
        max_err[i] = np.max(err_mag)
        h[i] = t[-1]/(len(t)-1)
        rel_err = err_mag/np.sqrt(exact[:,0]**2 + exact[:,1]**2 + exact[:,2]**2)
        plt.plot(t, rel_err, label="%s n=%d" % (d['method'], len(t)-1))
    s = 0.0
    for i in range(1, len(max_err)):
        s += 1/3 * np.log(max_err[i]/max_err[i-1])/np.log(h[i]/h[i-1])
    print ("Convergence rate of error: %g" % s)
    plt.xlabel("t [us]")
    plt.ylabel("relativ feil")
    plt.legend()
    plt.title("Estimert onvergenshastighet: %g" % s)
    plt.savefig("rel_err.pdf")

def plot_compare_analytic_case(d):
    """
    Plots trajectory of d together with analytical case.
    """
    t = d['t']
    r = d['r']
    exact = ana(t)
    plt.plot(t, exact[:,2], label="Exact")
    lab = "%s, n=%d" %(d['method'], len(t)-1)
    plt.plot(t, r[:,0,2], linestyle="dotted", label=lab)
    plt.xlabel("t [us]")
    plt.ylabel("Z [um]")
    plt.legend()
    plt.savefig("cmpana_zt.pdf")
    
    plt.figure()
    plt.plot(exact[:,0], exact[:,1], label="Exact")
    plt.plot(r[:,0,0], r[:,0,1], linestyle='dotted', label=lab)
    plt.xlabel("X [um]")
    plt.ylabel("Y [um]")
    plt.legend()
    plt.savefig("cmpana_xy.pdf")

def plot_resfreq():
    """
    
    """

def menu():
    print("Plot types are identified by string before ':'")
    print("Data objects are identified by numbers (1,2,...)")
    print("ex. Plot type 'B' on data object 3: >>> B 3")
    print("---")
    print("*** Plot types ***")
    print("A: - Plots graph of z(t) for a single data object")
    print("B: - Plots movement in xy-plane for a single data object")
    print("Cx: - Plots movement in phase space for x-axis of a single data object")
    print("Cy: - ... along y-axis")
    print("Cz: - ... along z-axis")
    print("D: - Plots trajectory in 3D")
    print("E: - Plots relative error of data objects to analytic case")
    print("F: - Plots tracjectory in xy-plane, and z(t) compared to analytic case")
    print("     for a single data object.")
    print("*** Data objects ***")
    for i in range(len(data)):
        print("%d. %s" % (i, data[i]['desc']))
    print("***")
    c = input(">>> ")
    words = c.split()
    if (words[0] == "A") :
        d = int(words[1])
        plot_zt(data[d])
    elif (words[0] == "B"):
        d = int(words[1])
        plot_xy(data[d])
    elif (words[0] == "Cx"):
        d = int(words[1])
        plot_phase(data[d],0)
    elif (words[0] == "Cy"):
        d = int(words[1])
        plot_phase(data[d],1)
    elif (words[0] == "Cz"):
        d = int(words[1])
        plot_phase(data[d],2)
    elif (words[0] == "D"):
        d = int(words[1])
        plot_3d_trajectory(data[d])
    elif (words[0] == "E"):
        D = [data[int(w)] for w in words[1:]]
        plot_analytic_rel_error(*D)
    elif (words[0] == "F"):
        d = int(words[1])
        plot_compare_analytic_case(data[d])
    else:
        print("Not a valid plot type.")

menu()

# def compare_analytic_case(*D):
#     """
#     Plots graphs comparing RK4 and FE-solutions to the analytical solution.
#     """
#     q = 1
#     m = 40.078
#     B = 9.64852558
#     V = 25*9.64852558e4
#     d = 500.0
#     om_z = om_z_func(q,V,m,d)
#     z0 = 20.0
#     T = 50.0
#     time = np.linspace(0, T, 1000)
#     plt.plot(time, analytic_z(z0, om_z, time), label="analytical")
#     for d in D:
#         plt.plot(d['t'], d['r'][:,0,2], linestyle='dotted', label=d['method'])
#     plt.xlabel("tid [us]")
#     plt.ylabel("z [um]")
#     plt.legend()
#     plt.savefig("com_ana_z(t).png")
