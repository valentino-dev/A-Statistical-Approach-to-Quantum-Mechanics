import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#from olib import *
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
    FILE_PATH = r"C:\Users\angel\Documents\git\A-Statistical-Approach-to-Quantum-Mechanics\CPU\Simulation\Simulation\data_fig6.csv"
    #FILE_PATH = "data_fig6.csv"
    Monte_Carlo_Iterations = 100
    MCPerCALC = 10
    Iterations = int(Monte_Carlo_Iterations/MCPerCALC)
    x_sq_mean = np.zeros((Iterations))
    x_quad_mean = x_sq_mean
    m = 0.5
    Lambda = 0
    a = 0.5
    deltatao = 1
    tao = np.arange(5)
    E0 = np.zeros(tao.shape[0])
    zero_deltatao_mean = np.zeros((Iterations, tao.shape[0]))
    zero_tao_mean = np.zeros((Iterations, tao.shape[0]))


    fig, ax = plt.subplots()
    for i in range(Iterations):
        data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(Monte_Carlo_Iterations), np.arange(i*MCPerCALC, (i+1)*MCPerCALC, 1))).to_numpy()
        for k in range(tao.shape[0]):
            zero_deltatao_mean[i, k] = (data[:,0]*data[:, int(tao[k] + deltatao)]).mean()
            print("0: ", data[:, 0])
            print("0: ", data[:, tao[k]])
            zero_tao_mean[i, k] = (data[:,0]*data[:, tao[k]]).mean()
        x_sq_mean[i] = (data.flatten()**2).mean()
        x_quad_mean[i] = (data.flatten()**4).mean()



    E0 = m**2*x_sq_mean+3*Lambda*x_quad_mean
    E0 = E0.mean()
    E1 = -1/deltatao*np.log(zero_deltatao_mean/zero_tao_mean)+E0
    E1 = E1.mean(axis=0)
    print("zero_deltatao_mean ", zero_deltatao_mean)
    print(zero_tao_mean)
    print("x mean ", x_sq_mean.mean())
    print("E0 ", E0)
    print("E1 ", E1)
    print("E1 - E0 ", E1 - E0)

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

