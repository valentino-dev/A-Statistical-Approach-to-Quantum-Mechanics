import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olib import *


if True: 

    FILE_PATH = "data_test.csv"
    radius = 2 
    resolution = radius*0.02
    tot_counts = np.zeros(np.arange(2*radius/resolution+1).shape[0]-1)
    for i in range(100):
        data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(1000), np.arange(0+i*10, 10+i*10))).to_numpy()
        data = data[:, :-1]
        counts, bins = np.histogram(data.flatten(), np.arange(2*radius/resolution+1)*resolution-radius)
        tot_counts = tot_counts + counts
    fig, ax = plt.subplots()
    print(tot_counts.sum())
    psi = tot_counts/tot_counts.sum()/(bins[1]-bins[0])/2
    ax.stairs(psi, bins)
    ax.plot(bins[:-1]+resolution*0.5, psi, linewidth=0.5)
    plt.show()
    exit()

# data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(100), np.arange(0, 100, 10)+9)).to_numpy()

FILE_PATH = "data_test.csv"
data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(1000), np.arange(0, 10, 1))).to_numpy()
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

radius = 1
for i in range(equilibrium_ensambles.shape[0]):
   if True:
        
       ax[i].set_xlim((-radius, radius))
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


