import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olib import *


if True:

    FILE_PATH = "data_fig5.csv"
    radius = 2
    resolution = 0.1
    tot_counts = np.zeros(int(2*radius/resolution))
    for i in range(10):
        data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(100), np.arange(0+i*10, 10+i*10))).to_numpy()
        data = data[:, :-1]
        counts, bins = np.histogram(data.flatten(), np.arange(2*radius/resolution+1)*resolution-radius)
        tot_counts = tot_counts + counts
    fig, ax = plt.subplots()
    print(tot_counts.sum())
    ax.stairs(tot_counts, bins)
    ax.scatter(bins[:-1]+0.05, tot_counts, marker="x", linewidths=0.5)
    plt.show()
    exit()

# data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(100), np.arange(0, 100, 10)+9)).to_numpy()

FILE_PATH = "data_fig4.csv"
data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(100), np.arange(0, 10, 1))).to_numpy()
# data = pd.read_csv(FILE_PATH, sep=";").to_numpy()
# print(data.shape)
# print(data[:,-1])
data = data[:,:-1]
# print(data.shape)
# print(data[-1,:])
equilibrium_ensambles = data
# data = np.array_split(data, 1)

fig, ax = plt.subplots(1, equilibrium_ensambles.shape[0])
# fig, ax = plt.subplots()

for i in range(equilibrium_ensambles.shape[0]):
   if True:
       ax[i].set_xlim((-1, 1))
       ax[i].set_title(f"{i}")
       ax[i].plot(data[i,:100], np.arange(equilibrium_ensambles.shape[1])[:100], label=f"Path {i}", linewidth=0.5)
       #ax[0, i].legend()

# counts, bins = np.histogram(data.flatten(), 20)
# ax.stairs(counts, bins)
#for i in range(data.shape[0]):
#    if i == data.shape[0]-1:
#        ax[0, 1].plot(data[i], np.arange(data.shape[1]), label=f"Markov Iter. {i}", linewidth=0.5)
#        counts, bins = np.histogram(data[i], 20)
#        ax[1, 1].stairs(counts, bins)
        
#ax[0, 0].set_xlim((-4, 4))
#ax[0, 1].set_xlim((-4, 4))
    
#ax[0, 0].legend()
#ax[0, 1].legend()
plt.show()


