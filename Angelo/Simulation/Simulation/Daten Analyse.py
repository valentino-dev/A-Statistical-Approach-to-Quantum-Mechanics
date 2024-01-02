import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olib import *

FILE_PATH = "data.csv"

data = pd.read_csv(FILE_PATH, sep=";").to_numpy()
data = data[:,:-1]

fig, ax = plt.subplots(2, 2)
for i in range(data.shape[0]):
    if i == 0:
        ax[0, 0].plot(data[i], np.arange(data.shape[1]), label=f"Markov Iter. {i}", linewidth=0.5)
        counts, bins = np.histogram(data[i], 20)
        ax[1, 0].stairs(counts, bins)

for i in range(data.shape[0]):
    if i == data.shape[0]-1:
        ax[0, 1].plot(data[i], np.arange(data.shape[1]), label=f"Markov Iter. {i}", linewidth=0.5)
        counts, bins = np.histogram(data[i], 20)
        ax[1, 1].stairs(counts, bins)


        
        
ax[0, 0].set_xlim((-4, 4))
ax[0, 1].set_xlim((-4, 4))
    
ax[0, 0].legend()
ax[0, 1].legend()
plt.show()


