import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olib import *

FILE_PATH = "data.csv"



data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(100), np.arange(0, 100, 10)+9)).to_numpy()
# data = pd.read_csv(FILE_PATH, sep=";", skiprows=np.delete(np.arange(100), np.arange(0, 10, 1))).to_numpy()

data = data[:,:-1]
equilibrium_ensambles = data
# data = np.array_split(data, 1)

fig, ax = plt.subplots(1, equilibrium_ensambles.shape[0])

for i in range(equilibrium_ensambles.shape[0]):
    if True:
        ax[i].set_title(f"Ensamble {i}")
        ax[i].plot(data[i,:100], np.arange(equilibrium_ensambles.shape[1])[:100], label=f"Path {i}", linewidth=0.5)
        ax[i].legend()

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


