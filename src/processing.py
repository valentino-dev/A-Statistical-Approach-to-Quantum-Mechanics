import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.odr
plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

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

#def func(beta, x):
    #return beta[0]+beta[1]*np.exp(beta[2]*x)

def func(beta, x):
    return beta[0]*np.exp(beta[1]*x)

def c():
    PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/GPU/"
    FILE_PATH = PATH + "data/data_fig6.csv"
    total_number_of_f_sq = 1
    totalConfigurations = 2720
    configurations_to_use = 2000
    MCPerCALC = 100
    a = 0.5
    Iterations = int(configurations_to_use/MCPerCALC)
    ts = np.arange(6) # tau divided by a
    zero_tau = np.zeros((total_number_of_f_sq, Iterations, MCPerCALC, ts.shape[0]))
    zero_tau_err = np.zeros_like(zero_tau)
    zero_tau_mean = np.zeros((total_number_of_f_sq, ts.shape[0]))
    zero_tau_mean_err = np.zeros_like(zero_tau_mean)
    delta_E = np.zeros((total_number_of_f_sq))

    fig, ax = plt.subplots()
    for n_f_sq in range(total_number_of_f_sq):
        for i in range(Iterations):
            print(f"({0.01*int(10000*(n_f_sq*Iterations + i)/total_number_of_f_sq/Iterations)}%)")
            data = pd.read_csv(FILE_PATH, sep=";", header=None, skiprows=np.delete(np.arange(totalConfigurations), np.arange(i*MCPerCALC, (i+1)*MCPerCALC, 1))).to_numpy()

            for t in ts:
                tempzt = data * np.roll(data, -t, axis=1) 
                tempztmean = tempzt.mean(axis=1)
                tempzterr = ((tempzt - np.expand_dims(tempztmean, axis=0).T)**2).mean(axis=1)**(1/2)/np.abs(tempztmean)
                zero_tau[n_f_sq, i, :, t] = tempztmean[:]
                zero_tau_err[n_f_sq, i, :, t] = tempzterr[:]

    zero_tau_mean = zero_tau.mean(axis=2).mean(axis=1)
    zero_tau_mean_err[:, :] = np.linalg.norm(np.linalg.norm(zero_tau_err, axis=2)/zero_tau_err.shape[2], axis=1)/zero_tau_err.shape[1]
    X = ts*a
    Xerr = np.zeros_like(X)
    Y = zero_tau_mean
    Yerr = zero_tau_mean_err

    for n_f_sq in range(total_number_of_f_sq):
        ax, model = plotData(ax, X, Xerr, Y[n_f_sq], Yerr[n_f_sq], polyfit=0, fmt="x")
        # ax.scatter(X, Y[n_f_sq], marker="x", linewidths=0.7)
        ax.set_yscale("log")
        model = scipy.odr.Model(func)
        odrFit = scipy.odr.ODR(scipy.odr.RealData(X, Y[n_f_sq], Xerr, Yerr[n_f_sq]), model, [-1, 1, 0])
        odrFit.set_job(fit_type = 2)
        output = odrFit.run()
        delta_E[n_f_sq] = -output.beta[0]
        x = X
        y = func(output.beta, x)
        ax.plot(x, y, linewidth = 0.7)
        
    print(delta_E)


    ax.set_title(r"Fig. 6: $\Delta E="+f"{delta_E[0]}"+r"$")
    print("plotting")
    fig.savefig("fig6.pdf", dpi=500)


def d():
    N = 50
    s = 20
    ff = np.array([0.0, 0.5, 1, 1.5, 2])
    name = ["a", "b", "c", "d", "e"]
    radius = 5

    fig, ax = plt.subplots(1, 6)
    for i in range(len(name)):
        FILE_PATH = f"data/data_fig7{name[i]}.csv"
        Monte_Carlo_Iterations = 2720
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
        fig.savefig("fig7with10data.pdf", dpi=500)

def e():

    #FILE_PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/CPU/Simulation/Simulation/data_fig6.csv"
    #PATH = "~/Git"
    PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/GPU/"
    FILE_PATH = PATH + "data/data_fig8.csv"
    totalConfigurations = 10880
    configurations_to_use = 1000
    MCPerCALC = 10
    Iterations = int(configurations_to_use/MCPerCALC)
    m = 0.5
    Lambda = 1
    a = 0.25
    deltatao = 0.5 
    radius = 4 
    resolution = radius*0.05
    tot_counts = np.zeros(np.arange(2*radius/resolution+1).shape[0]-1)

    fig, ax = plt.subplots()
    for i in range(Iterations):
        data = pd.read_csv(FILE_PATH, sep=";", header=None, skiprows=np.delete(np.arange(totalConfigurations), np.arange(i*MCPerCALC, (i+1)*MCPerCALC, 1))).to_numpy()
        counts, bins = np.histogram(data.flatten(), np.arange(2*radius/resolution+1)*resolution-radius)
        tot_counts = tot_counts + counts
    psi = tot_counts/tot_counts.sum()/(bins[1]-bins[0])
    ax.scatter(bins[:-1]+resolution*0.5, psi, linewidths=1, marker="x")
    x = np.linspace(-radius, radius, 100)
    fig.set_size_inches(11.69, 8.27)
    fig.savefig("fig8.pdf", dpi=500)

def f():
    PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/GPU/"
    FILE_PATH = PATH + "data/data_fig9"
    totalConfigurations = 10880
    configurations_to_use = 20 
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

    Lambda = 0
    a = 0.5
    deltatao = 0.5 
    fig, ax = plt.subplots()
    n_bars = [5, 10, 15, 20]
    n_bar_label = ["a", "b", "c", "d"]
    for n_idx in range(len(n_bars)):
        result = {"zero_tao": [], "zero_tao_err": [], "E0": [], "E1": [], "E1-E0": []}
        err = {"zero_tao": []}
        tao_div_a = np.arange(0, 20)

        for i in range(Iterations):
            data = pd.read_csv(FILE_PATH+n_bar_label[n_idx]+".csv", sep=";", skiprows=np.delete(np.arange(totalConfigurations), np.arange(i*MCPerCALC, (i+1)*MCPerCALC, 1))).to_numpy()

            zero_tao_bin = []
            for i in tao_div_a:
                zero_tao_bin.append((data * np.roll(data, -i, axis=1)).mean(axis=(0,1)))
            result["zero_tao"].append(zero_tao_bin)
            zero_tao_bin = np.array(zero_tao_bin)

        for key in result:
            if result[key] != []:
                result[key].append(np.array(result[key]).mean(axis=0))
                print(f"{key}: {result[key][-1]}")

        f = np.array(result["zero_tao"])
        for i in tao_div_a:
            err["zero_tao"].append((((f[-1, i]-f)[:-1, i])**2).mean()**(1/2)/np.abs(f[-1, i]))


        X = tao_div_a*a
        Y = np.log10(result["zero_tao"][-1])

        ax.scatter(X, Y, marker=["s", "x", "s", "v"][n_idx], label=r"$\bar{n}="+str(n_bars[n_idx])+r"$")
    ax.set_title(r"Anharmoic Oscillator: $f^{2}=2.0$")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$\log_{10} (<x(0)x(\tau)>)$")
    ax.legend()
    print("plotting")
    fig.savefig("fig9.pdf", dpi=500)


def func(beta, x):
    return beta[2] + beta[1] * np.exp(beta[0]*x)

def g():
    PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/GPU/"
    FILE_PATH = PATH + "data/data_fig10_"
    total_number_of_f_sq = 5
    totalConfigurations = 5440
    configurations_to_use = 5000
    MCPerCALC = 100
    a = 0.25
    Iterations = int(configurations_to_use/MCPerCALC)
    ts = np.arange(8) # tau divided by a
    zero_tau = np.zeros((total_number_of_f_sq, Iterations, MCPerCALC, ts.shape[0]))
    zero_tau_err = np.zeros_like(zero_tau)
    zero_tau_mean = np.zeros((total_number_of_f_sq, ts.shape[0]))
    zero_tau_mean_err = np.zeros_like(zero_tau_mean)
    delta_E = np.zeros((total_number_of_f_sq))

    fig, ax = plt.subplots(1, 2)
    for n_f_sq in range(total_number_of_f_sq):
        for i in range(Iterations):
            print(f"({0.01*int(10000*(n_f_sq*Iterations + i)/total_number_of_f_sq/Iterations)}%)")
            data = pd.read_csv(FILE_PATH + str(n_f_sq) + ".csv", sep=";", header=None, skiprows=np.delete(np.arange(totalConfigurations), np.arange(i*MCPerCALC, (i+1)*MCPerCALC, 1))).to_numpy()

            for t in ts:
                tempzt = data * np.roll(data, -t, axis=1) 
                tempztmean = tempzt.mean(axis=1)
                tempzterr = ((tempzt - np.expand_dims(tempztmean, axis=0).T)**2).mean(axis=1)**(1/2)/np.abs(tempztmean)
                zero_tau[n_f_sq, i, :, t] = tempztmean[:]
                zero_tau_err[n_f_sq, i, :, t] = tempzterr[:]

    zero_tau_mean = zero_tau.mean(axis=2).mean(axis=1)
    zero_tau_mean_err[:, :] = np.linalg.norm(np.linalg.norm(zero_tau_err, axis=2)/zero_tau_err.shape[2], axis=1)/zero_tau_err.shape[1]
    X = ts*a
    Xerr = np.zeros_like(X)
    Y = zero_tau_mean
    Yerr = zero_tau_mean_err

    for n_f_sq in range(total_number_of_f_sq):
        # ax[1].errorbar(X, Y[n_f_sq], Yerr[n_f_sq], fmt="x", label=str(n_f_sq*0.5))
        ax[1].scatter(X, Y[n_f_sq], marker="x", label=r"$f^2=" + str(n_f_sq*0.5) + r"$", linewidths=0.7)
        ax[1].set_yscale("log")
        model = scipy.odr.Model(func)
        odrFit = scipy.odr.ODR(scipy.odr.RealData(X, Y[n_f_sq], Xerr, Yerr[n_f_sq]), model, [-1, 1, 0])
        odrFit.set_job(fit_type = 2)
        output = odrFit.run()
        delta_E[n_f_sq] = -output.beta[0]
        x = X
        y = func(output.beta, x)
        ax[1].plot(x, y, linewidth = 0.7)
    print(delta_E)
    ax[0].scatter(np.arange(total_number_of_f_sq)/2, delta_E, marker="x", linewidths=0.7, label="simulated")
    ax[0].scatter(np.arange(total_number_of_f_sq)/2, [4.3, 3.4, 2.4, 1.4, 0.7], marker="x", linewidths=0.7, label="paper")

    ax[1].legend()
    ax[0].legend()
    ax[1].set_xlabel(r"$\tau$")
    ax[1].set_ylabel(r"$\langle x(0)x(\tau)\rangle$")
    ax[0].set_xlabel(r"$f^2$")
    ax[0].set_ylabel(r"$\Delta E(f^2)$")
    ax[0].set_title(r"Anharmonic Oscillator")
    ax[1].grid(visible=True, which="minor", linewidth=0.2, linestyle="--")
    ax[1].grid(visible=True, which="major", linestyle="--")
    ax[0].grid(visible=True, which="minor", linewidth=0.2, linestyle="--")
    ax[0].grid(visible=True, which="major", linestyle="--")
    fig.tight_layout()
    fig.savefig("fig10.pdf", dpi=500)


def printAction():
    PATH = "~/Git/A-Statistical-Approach-to-Quantum-Mechanics/GPU/"
    FILE_PATH = PATH + "data/data_fig6.csv"
    totalConfigurations = 500
    configurations_to_use = 500
    MCPerCALC = 10
    a = 0.5
    m0 = 0.5
    Lambda = 0.0
    mu = 2.0
    f = 2.0


    Iterations = int(configurations_to_use/MCPerCALC)
    action = np.zeros((configurations_to_use))

    for i in range(Iterations):
        data = pd.read_csv(FILE_PATH, sep=";", header=None, skiprows=np.delete(np.arange(totalConfigurations), np.arange(i*MCPerCALC, (i+1)*MCPerCALC, 1))).to_numpy()
        for k in range(data.shape[0]):
            action[i*data.shape[0]+k] = (m0 /a/ 2 * (data[k] - np.roll(data[k], 1))**2 + a*mu/2*data[k]**2+Lambda*data[k]**4).sum()
    
    X = np.arange(configurations_to_use)
    Y = action

    fig, ax = plt.subplots()
    ax.plot(X, Y)

    fig.savefig("action.pdf", dpi=500)


# printAction()


for char in exec:
    locals()[char]()

#plt.show()

