



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dask.dataframe as dd
import dask.array as da
import scipy.odr
import sys
from scipy.optimize import curve_fit

# sys.path.append("/home/angelo/Git/olib")
# from olib import *

plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

# path = "data/data_fig5_1k_RS_100.csv"

# Recording Start = 100, 8192 Configurations, Statistical independend iterations (every SII one mesurement is taken) = 1
path = "data/data_RS100_CONF8192_SII1.csv"
# path = "data/data_fig5_100k.csv"


def autoCorrelationErrorPerBinSize():
    # read and crop datau

    '''
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 2 ** (da.floor(da.log2(data.shape[0])))
    print(length)
    data = data[:length, :]
    print(data.shape)

    # calc correltaion
    correlator = (data * da.roll(data, 10, 1)).mean(1)
    '''
    path = "data/correlation.csv"
    correlator = da.reshape(dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True), (5000, 26))[:, 0]
    length = 2 ** (da.floor(da.log2(correlator.shape[0])))
    correlator = correlator[:length]

    # rebin and calc error
    bin_sizes = 2 ** np.arange(1,  int(np.log2(int(correlator.shape[0]))))
    bin_count = correlator.shape[0] / bin_sizes
    print(bin_count)
    deltas = np.zeros(bin_sizes.shape[0])
    for i in range(deltas.shape[0]):
        print("i: ", i)
        print(bin_sizes[i], bin_count[i])

        bined_correlator = da.reshape(correlator, (bin_count[i], bin_sizes[i])).mean(1)
        delta = da.sqrt(
            ((bined_correlator - bined_correlator.mean()) ** 2).sum()
            / (bin_count[i] - 1) / bin_count[i]
        )
        deltas[i] = delta.compute()
        
    fig, ax = plt.subplots()
    ax.scatter(bin_sizes, deltas)
    plt.show()




def f(x, beta):
    return beta[0] + beta[1] * np.exp(beta[2] * x + beta[3])

def calcCorrelation():
    bin_size = 8
    time_slice = 0.5

    # read and crop data
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 2 ** (da.floor(da.log2(data.shape[0])))
    bin_count = length / bin_size
    print(f"{length} sites are going to be used")

    # cutting and binning data
    data = data[:length, :]
    data_bin = da.reshape(data, (bin_count, bin_size, data.shape[1]))
    lattice_time = np.arange(0, 5)

    correlator_bin = da.zeros(
        (lattice_time.shape[0], data_bin.shape[0], data_bin.shape[1])
    )

    for i in range(lattice_time.shape[0]):
        correlator_bin[i] = (data_bin * da.roll(data_bin, i, 2)).mean(2)

    correlator_estimate = correlator_bin.mean(2)
    correlator_estimate_error = (
        (((correlator_bin - da.expand_dims(correlator_estimate, 2)) ** 2).sum(2))
        / (bin_size - 1)
        / bin_size
    ) ** (1 / 2)

    delta_energies = da.zeros(10)
    resample_size = correlator_estimate.size

    for i in range(0, delta_energies.shape[0]):

        # choosing random correlation times
        time_indecies = np.array(
            da.floor(
                da.random.random_sample((resample_size)) * correlator_estimate.shape[0]
            ).compute(),
            dtype=int,
        )

        # choosing random uncorrelated correlators
        correlator_indecies = np.array(
            da.floor(
                da.random.random_sample((resample_size)) * correlator_estimate.shape[1]
            ).compute(),
            dtype=int,
        )
        correlator_estimate_sample_estimate = correlator_estimate.compute()[
            time_indecies, correlator_indecies
        ]
        correlator_estimate_sample_estimate_error = correlator_estimate_error.compute()[
            time_indecies, correlator_indecies
        ]

        popt, pcov = curve_fit(
            exponentialFunction,
            time_indecies * time_slice,
            correlator_estimate_sample_estimate,
            p0=np.array([0, 1, 1, 0]),
            sigma=correlator_estimate_sample_estimate_error,
        )

        delta_energies[i] = popt[2]
        print(popt)

    delta_E = delta_energies.mean().compute()
    delta_E_error = (
        (((delta_energies - delta_E) ** 2).sum() / (delta_energies.shape[0] - 1))
        ** (1 / 2)
    ).compute()
    print(f"dE={delta_E}+-{delta_E_error}")


autoCorrelationErrorPerBinSize()
# calcCorrelation()
