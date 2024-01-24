import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olib import *
# a=fig4, b=fig5, ..
exec = "c"


# fig, ax = plt.subplots(4, 4)
def a():
    FILE_PATH = "data/data_fig4.csv"
    Monte_Carlo_Iterations = 100
    data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(Monte_Carlo_Iterations), np.arange(90, 100, 1))).to_numpy()
    data = data[-1]
    
    data = data.flatten()

    fig, ax = plt.subplots()
    radius = 4
    ax.set_xlim((-radius, radius))
    ax.set_title(f"fig4: Trajectory")
    ax.plot(data[:100], np.arange(data.shape[0])[:100], label=f"Trajectory", linewidth=0.5)
    fig.set_size_inches(11.69, 8.27)
    fig.savefig("fig4.pdf", dpi=500)

def b(): 
    fig, ax = plt.subplots()
    FILE_PATH = "data/data_fig4.csv"
    ax.set_title(f"fig5: Probability Density")
    Monte_Carlo_Iterations = 1360
    radius = 4 
    resolution = radius*0.05
    tot_counts = np.zeros(np.arange(2*radius/resolution+1).shape[0]-1)
    MCPerHistogram = 10
    for i in range(int(Monte_Carlo_Iterations/MCPerHistogram)):
        data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(Monte_Carlo_Iterations), np.arange(i*MCPerHistogram, (i+1)*MCPerHistogram, 1))).to_numpy()
        counts, bins = np.histogram(data.flatten(), np.arange(2*radius/resolution+1)*resolution-radius)
        tot_counts = tot_counts + counts
    psi = tot_counts/tot_counts.sum()/(bins[1]-bins[0])
    ax.scatter(bins[:-1]+resolution*0.5, psi, linewidths=1, marker="x")
    x = np.linspace(-radius, radius, 100)
    ax.plot(x, 3.14**(-1/2)*np.exp(-x**2), linestyle="--", color="r", linewidth=0.5)
    ax.plot(x, 0.59*np.exp(-1.1*x**2), linestyle="--", color="b", linewidth=0.5)
    fig.set_size_inches(11.69, 8.27)
    fig.savefig("fig5.pdf", dpi=500)

def c():
    #FILE_PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/CPU/Simulation/Simulation/data_fig6.csv"
    FILE_PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/GPU/data/data_fig6.csv"
    totalConfigurations = 2720
    configurations_to_use = 2000
    MCPerCALC = 10
    Iterations = int(configurations_to_use/MCPerCALC)
    m = 0.5
    Lambda = 0
    a = 0.5
    deltatao = 1
    result = {"zero_tao": [], "zero_tao_err": [], "E0": [], "E1": [], "E1-E0": []}
    err = {"zero_tao": []}
    tao = np.arange(0, 4)

    fig, ax = plt.subplots()
    for i in range(Iterations):
        data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(totalConfigurations), np.arange(i*MCPerCALC, (i+1)*MCPerCALC, 1))).to_numpy()

        zero_tao_bin = []
        zero_tao_bin_err = []
        # zero_tao.shape = (bins, tao)
        for i in tao:
            zero_tao_bin.append((data * np.roll(data, i, axis=1)).mean(axis=(0,1)))
        result["zero_tao"].append(zero_tao_bin)
        zero_tao_bin = np.array(zero_tao_bin)
        #result["E0"].append(m**2*zero_tao_bin+Lambda*zero_tao_bin)
        #result["E1-E0"].append(-np.log(np.roll(zero_tao_bin, -deltatao)/zero_tao_bin)/deltatao)

    for key in result:
        result[key].append(np.array(result[key]).mean(axis=0))
        print(f"{key}: {result[key][-1]}")

    #print("test2", np.array(result["zero_tao"]))
    #f = np.log(result["zero_tao"])
    f = np.array(result["zero_tao"])
    #print("test ", f[:,-1])
    for i in tao:
        err["zero_tao"].append((((f[-1, i]-f)[:-1, i])**2).mean()**(1/2)/np.abs(f[-1, i]))
    

    print("Uncertainties:")
    for key in err:
        print(f"{key}: {err[key]}")

    X = tao
    Xerr = np.zeros_like(tao)
    #Y = np.array(result["zero_tao"][-1])
    Y = np.log10(result["zero_tao"][-1])
    #Yerr = np.abs(1e-2/np.array(result["zero_tao"][-1]))
    #Yerr = np.array(err["zero_tao"])
    Yerr = np.abs(np.array(err["zero_tao"])/np.array(result["zero_tao"])[-1]/np.log(10))
    #Yerr = np.abs(np.log(err["zero_tao"]))
    ax, model = plotData(ax, X, Xerr, Y, Yerr, polyfit=1, fmt="x", label="Daten")
    model.printParameter()
    #ax = setSpace(ax, X, Y)
    
    #ax.set_yscale("log")
    print("plotting")
    plt.show()


def d():
    N = 50
    s = 20
    ff = np.array([0.5, 1, 2, 4, 8, 16])
    name = ["a", "b", "c", "d", "e", "f"]
    radius = 5

    fig, ax = plt.subplots(1, 6)
    for i in range(len(name)):
        FILE_PATH = f"data_fig7{name[i]}.csv"
        Monte_Carlo_Iterations = 100
        data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(Monte_Carlo_Iterations), np.arange(90, 100, 1))).to_numpy()
        data = data[-1]
        
        data = data.flatten()

        ax[i].set_xlim((-radius, radius))
        ax[i].set_title(f"{name[i]}: (f²={ff[i]})")
        ax[i].plot(data[:N], np.arange(data.shape[0])[:N], label=f"Trajectory", linewidth=0.5)
        ax[i].scatter(data[:N], np.arange(data.shape[0])[:N], marker="x", linewidth=0.5, s=s)
        ax[i].plot([-ff[i]**(1/2), -ff[i]**(1/2)], [0, N], linewidth=0.5, linestyle="--", color="#333333")
        ax[i].plot([ff[i]**(1/2), ff[i]**(1/2)], [0, N], linewidth=0.5, linestyle="--", color="#333333")
        fig.set_size_inches(11.69, 8.27)
        fig.savefig("fig7.pdf", dpi=500)

def e():

    fig, ax = plt.subplots()


    FILE_PATH = "data_fig8.csv"
    Monte_Carlo_Iterations = 2000
    radius = 2.4
    resolution = radius*0.1
    tot_counts = np.zeros(np.arange(2*radius/resolution+1).shape[0]-1)
    MCPerHistogram = 10
    xj_count = 0
    for i in range(int(Monte_Carlo_Iterations/MCPerHistogram)):
        data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(Monte_Carlo_Iterations), np.arange((i)*MCPerHistogram, (i+1)*MCPerHistogram, 1))).to_numpy()
        xj_count += data.size
        counts, bins = np.histogram(data.flatten(), np.arange(2*radius/resolution+1)*resolution-radius)
        tot_counts = tot_counts + counts
    print(xj_count)
    psi = tot_counts/tot_counts.sum()/(bins[1]-bins[0])

    ax.set_title(f"fig8: Probability Density")
    ax.scatter(bins[:-1]+resolution*0.5, psi, linewidths=1, marker="x")

    x = np.linspace(-radius, radius, 100)
    ax.set_xlim((-2.4, 2.4))
    ax.set_ylim((0,0.7))
    #ax.plot([-2**(1/2), -2**(1/2)], [0, 100], linewidth=0.5, linestyle="--", color="#333333")
    #ax.plot([2**(1/2), 2**(1/2)], [0, 100], linewidth=0.5, linestyle="--", color="#333333")
    # ax.plot(x, 0.59*np.exp(-1.1*x**2), linestyle="--", color="b", linewidth=0.5)
    fig.set_size_inches(11.69, 8.27)
    fig.savefig("fig8.pdf", dpi=500)


        




for char in exec:
    locals()[char]()

#plt.show()

