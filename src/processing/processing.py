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

def delta(obs: np.ndarray, ax=0):
    return np.std(obs, ddof=1, axis=ax)/np.sqrt(obs.shape[ax])
    
def autoCorrelationErrorPerBinSize():
    # read and crop datau

    """
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    correlator = (data * da.roll(data, 1, 1)).mean(1).compute()
    """
    path = "data/correlation.csv"
    correlator = da.reshape(
        dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True), (5000, 26)
    )[:, 0]

    # Hier die größe der bins angeben :

    # bin_sizes = 2 ** np.arange(1,  int(np.log2(int(correlator.shape[0]))))
    bin_sizes = 2 ** np.arange(1, 1 + int(np.log2(256)))
    # bin_sizes = np.arange(1, 100, 1)

    bin_count = np.floor(correlator.shape[0] / bin_sizes)
    lengths = np.array(bin_count * bin_sizes, dtype=np.uint)
    print(lengths)
    deltas = np.zeros(bin_sizes.shape[0])
    for i in range(bin_sizes.shape[0]):
        print("i: ", i)
        
        # cropping to resize
        cropped_correlator = correlator[: int(lengths[i])]

        # resize
        bined_correlator = np.reshape(
            cropped_correlator, (bin_count[i], bin_sizes[i])
        ).mean(1)

        # Fehler berechnen
        deltas[i] = delta(bined_correlator)

    fig, ax = plt.subplots()
    ax.scatter(bin_sizes, deltas)
    plt.show()


def f(x, beta):
    return beta[0] + beta[1] * np.exp(beta[2] * x + beta[3])


def func(x, b1, b2):
    return b1 * np.exp(-b2 * x)


def pickIndecies(shape, resample_size):
    # choosing random correlation times
    time_indecies = np.array(
        np.floor(np.random.random_sample((resample_size)) * shape[0]),
        dtype=np.uint,
    )

    # choosing random uncorrelated correlators
    correlator_indecies = np.array(
        np.floor(np.random.random_sample((resample_size)) * shape[1]),
        dtype=np.uint,
    )

    return time_indecies, correlator_indecies


def pickSample(obs: np.ndarray, obs_error: np.ndarray, size):
    time_idx, obs_idx = pickIndecies(obs.shape, size)
    sample = obs[time_idx, obs_idx]
    error_sample = obs_error[time_idx, obs_idx]
    return sample, error_sample, time_idx


def calcCorrelation():
    time_slice = 0.5
    lattice_time = np.arange(0, 5)

    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True).compute()
    correlator = np.zeros((lattice_time.shape[0], data.shape[0]))

    for i in range(lattice_time.shape[0]):
        correlator[i] = (data * np.roll(data, i, 1)).mean(1)

    """
    path = "data/correlation.csv"
    correlator = (
        da.reshape(
            dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True),
            (5000, 26),
        )
        .compute()
        .T[lattice_time, :]
    )
    """

    print(type(correlator))
    print(correlator.shape)

    bin_size = 16
    bin_count = int(correlator.shape[1] / bin_size)
    length = bin_count * bin_size

    # cutting and binning data
    correlator = correlator[:, :length]
    correlator = np.reshape(correlator, (lattice_time.shape[0], bin_count, bin_size))
    uncorrelated_correlator = correlator.mean(2)
    uncorrelated_correlator_error = delta(correlator, axis=2)

    delta_energies = np.zeros(1000)
    resample_size = correlator.size

    # Bootstraping
    for i in range(0, delta_energies.shape[0]):
        print(i, delta_energies.shape[0])

        # chose random sample
        correlator_sample, correlator_error_sample, time_idx = pickSample(
            uncorrelated_correlator, uncorrelated_correlator_error, resample_size
        )
        
        # fit to sample
        popt, pcov = curve_fit(
            func,
            time_idx * time_slice,
            correlator_sample,
            p0=np.array([1, 1]),
            sigma=correlator_error_sample,
        )

        delta_energies[i] = popt[1]

    delta_E = delta_energies.mean()
    delta_E_error = np.std(delta_energies, ddof=1)
    print(f"dE={delta_E}+-{delta_E_error}")


autoCorrelationErrorPerBinSize()
# calcCorrelation()
