from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import sys

"""Open and read all files on form "data/output_*.txt". Extract data into a list of 
dictionaries called 'data'.
"""
data = []
filenames = glob("data/output*.txt")
if not filenames:
    print("ERROR: No data files found")
    sys.exit()
for i in range(len(filenames)):
    with open(filenames[i], 'r') as file:
        file.readline() # discard
        fn = filenames[i][5:]
        data.append({"filename": fn})
        words = file.readline().split()
        data[i]["L"] = int(words[0])
        data[i]["T"] = float(words[1])
        data[i]["random"] = int(words[2])
        data[i]["max_cycles"] = float(words[3])
        data[i]["target_cycles"] = float(words[4])
        data[i]["chunk_sz"] = int(words[5])
        data[i]["sample_freq"] = int(words[6])
        data[i]["epsilon"] = float(words[7])
        data[i]["max_chunks"] = int(words[8])
        data[i]["target_chunks"] = int(words[9])
        file.readline() # discard
        words = file.readline().split()
        data[i]["E_chunk_avg"] = np.array([float(w) for w in words])
        words = file.readline().split()
        data[i]["E_chunk_index"] = np.array([float(w) for w in words])
        words = file.readline().split()
        data[i]["E"] = np.array([float(w) for w in words])
        words = file.readline().split()
        data[i]["M"] = np.array([float(w) for w in words])
        file.readline() # discard
        words = file.readline().split()
        data[i]["equilibrium_index"] = int(words[0])

def plot_E():
    D = data
    for i in range(len(D)):
        plt.figure()
        L = data[i]["L"]
        spc = L*L/data[i]["sample_freq"]
        l = plt.plot(np.arange(len(D[i]["E"]))/spc, D[i]["E"]/(L*L), zorder=0.0)
        c = l[0].get_color()
        c = "blue"
        # plt.scatter(D[i]["E_chunk_index"]/spc, D[i]["E_chunk_avg"]/(L*L), marker='o', color=c, zorder=1.0)
        # eqi = D[i]["equilibrium_index"]
        # if not eqi == -1:
        #     plt.axvline(D[i]["equilibrium_index"]/spc, ls="dashed", color=c)
        plt.xlabel("Sykluser")
        plt.ylabel("$\epsilon [J]$")
        plt.savefig("Energy%d.pdf" % i)

def plot_M():
    D = data
    for i in range(len(D)):
        L = data[i]["L"]
        spc = L*L/data[i]["sample_freq"]
        l = plt.plot(np.arange(len(D[i]["M"]))/spc, np.abs(D[i]["M"])/(L*L), zorder=0.0)
        c = l[0].get_color()
        # eqi = D[i]["equilibrium_index"]
        # if not eqi == -1:
        #     plt.axvline(D[i]["equilibrium_index"]/spc, ls="dashed", color=c)
    plt.xlabel("Sykluser")
    plt.ylabel("$m$")
    plt.savefig("Mag.pdf")

def plot_E_and_M_avg():
    D = data
    for i in range(len(D)):
        L = D[i]["L"]
        T = D[i]["T"]
        if D[i]["random"]:
            s = "tilfeldig"
        else:
            s = "ordnet"
        spc = L*L/D[i]["sample_freq"]
        M = np.abs(D[i]["M"]/(L*L))
        M_cummean = np.cumsum(M)/np.arange(1, len(M)+1)
        plt.plot(np.arange(len(M))/spc, M_cummean, label="T=%.1f, %s" % (T,s))
    plt.xlabel("Sykluser")
    plt.ylabel("$m$ kumulativt snitt")
    plt.legend()
    plt.savefig("M_cummean.pdf")

    plt.figure()
    for i in range(len(D)):
        T = D[i]["T"]
        L = D[i]["L"]
        if D[i]["random"]:
            s = "tilfeldig"
        else:
            s = "ordnet"
        spc = L*L/D[i]["sample_freq"]
        E = data[i]["E"]/(L*L)
        E_cummean = np.cumsum(E)/np.arange(1, len(E)+1)
        plt.plot(np.arange(len(M))/spc, E_cummean, label="T=%.1f, %s" % (T,s))
    plt.xlabel("Sykluser")
    plt.ylabel("$\epsilon$ kumulativt snitt")
    plt.legend()
    plt.savefig("E_cummean.pdf")

def plot_dist():
    D = data
    L = D[0]["L"]
    cut = int(100 * D[0]["L"]**2 / D[0]["sample_freq"])
    mvecs = np.concatenate([d["M"][cut:]/L**2 for d in D])
    plt.hist(mvecs, bins=30, density=True)
    plt.savefig("mdist.pdf")

    plt.figure()
    evecs = np.concatenate([d["E"][cut:]/L**2 for d in D])
    plt.hist(evecs, bins=26, density=True)
    plt.xlabel("$\epsilon$ [J]")
    plt.savefig("edist.pdf")

def retrieve(L, T):
    """List of data entries with L=L and T=T"""
    li = []
    for d in data:
        if (d["L"] == L and d["T"] == T):
            li.append(d)
    return li

def plot_ofT():
    D = data
    L = set()
    T = set()
    for d in D:
        L.add(d["L"])
        T.add(d["T"])
    L = sorted(L)
    T = sorted(T)

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    for i in range(len(L)):
        l = L[i] 
        e = np.zeros(len(T))
        m = np.zeros_like(e)
        Cv = np.zeros_like(e)
        X = np.zeros_like(e)
        for j in range(len(T)):
            t = T[j]
            objects = retrieve(l, t)
            ej = np.concatenate([ob["E"]/(l*l) for ob in objects])
            e[j] = np.mean(ej)
            mj = np.concatenate([np.abs(ob["M"]/(l*l)) for ob in objects])
            m[j] = np.mean(mj)
            e_sq_j = ej**2
            m_sq_j = mj**2
            Cv[j] = 1/(t**2) * (np.mean(e_sq_j) - e[j]**2)
            X[j] = 1/t * (np.mean(m_sq_j) - m[j]**2)
        ax1.plot(T, e, label="L=%d" % l)
        ax2.plot(T, m, label="L=%d" % l)
        ax3.plot(T, Cv, label="L=%d" % l)
        ax4.plot(T, X, label="L=%d" % l)
    ax1.set_xlabel("$T [J/k_B]$")
    ax1.set_ylabel("$\epsilon [J]$")
    ax2.set_xlabel("$T [J/k_B]$")
    ax2.set_ylabel("$m$")
    ax3.set_xlabel("$T [J/k_B]$")
    ax3.set_ylabel("$C_V$")
    ax4.set_xlabel("$T [J/k_B]$")
    ax4.set_ylabel("$\chi$")
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    fig1.savefig("energy_T.pdf")
    fig2.savefig("mag_T.pdf")
    fig3.savefig("Cv_T.pdf")
    fig4.savefig("chi_T.pdf")

if sys.argv[1] == "A":
    plot_E()
elif sys.argv[1] == "B":
    plot_E_and_M_avg()
elif sys.argv[1] == "C":
    plot_M()
elif sys.argv[1] == "D":
    plot_dist()
elif sys.argv[1] == "E":
    plot_ofT()
