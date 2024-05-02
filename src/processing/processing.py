import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dask.dataframe as dd
import dask.array as da
import scipy.odr
import sys
import scipy.optimize import curve_fit

# sys.path.append("/home/angelo/Git/olib")
# from olib import *

plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

# path = "data/data_fig5_1k_RS_100.csv"
path = "data/data_RS100_CONF8192_SII1.csv"
# path = "data/data_fig5_100k.csv"

"""
def autoCorrelationErrorPerBinStep():
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    correlation = data * np.roll(data, 1, 1)

    std_of_bins = np.zeros(int(np.log10(correlation.shape[0])))
    correlation_bined_mean = correlation
    for i in range(std_of_bins.shape[0]):
        print("i:", i)
        print(std_of_bins.shape[0]-i)
        correlation_bined = correlation_bined_mean.reshape(10**(std_of_bins.shape[0]-i-1), 10, correlation_bined_mean.shape[1])
        print(correlation_bined.shape)
        correlation_bined_mean = correlation_bined.mean(1)
        print(correlation_bined_mean.shape)
        std = np.std(correlation_bined_mean, 0)
        std_of_bins[i] = std.mean()

    X = np.arange(0, std_of_bins.shape[0])
    Xerr = np.zeros_like(X)
    Y = std_of_bins
    Yerr = np.zeros_like(Y)
    table = Table(X, Xerr, Y, Yerr, "Correlation Error", "bin step", r"$\sigma_x$")

    fig, ax = plt.subplots()
    ax = plotData(ax, table, polyfit=0, fmt="x")
    ax = setSpace(ax, table)
    plt.savefig("plots/CorrelationError.pdf", dpi=500)
"""


def autoCorrelationErrorPerBinSize():
    # read and crop data
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 2 ** (da.floor(da.log2(data.shape[0])))
    print(length)
    data = data[:length, :]
    print(data.shape)

    # calc correltaion
    correlation = data * da.roll(data, 1, 1)
    print(correlation.shape)

    # rebin and calc error
    bin_sizes = 2 ** np.arange(1, int(np.log2(int(data.shape[0]))))
    deltas = np.zeros(bin_sizes.shape[0])
    for i in range(deltas.shape[0]):
        print("i: ", i)
        corr_shape = (
            correlation.shape[0] / bin_sizes[i],
            bin_sizes[i],
            correlation.shape[1],
        )
        print("Corr Shape: ", corr_shape)
        correlation_bined = da.reshape(correlation, corr_shape)
        cb_mean = da.mean(correlation_bined, 1)
        cbm_mean = da.mean(cb_mean)
        delta = (
            da.sum((cb_mean - cbm_mean) ** 2)
            / (cb_mean.shape[0] * (cb_mean.shape[0] - 1))
        ) ** (1 / 2)

        print(deltas.shape[0] - i)
        deltas[i, :] = delta.compute()
        std = (
            da.sum(
                (correlation_bined_mean - da.mean(correlation_bined_mean, 0)) ** 2, 0
            )
            / ((correlation_bined_mean.shape[0] - 1) * correlation_bined_mean.shape[0])
        ) ** (1 / 2)
        # std = da.std(correlation_bined_mean, 0)/correlation_bined_mean.shape[0]**(1/2)
        # print(std.dask)
        print(std_of_bins.shape[0] - i)
        std_of_bins[i, :] = std.compute()[:]
        # std_of_bins[i] = std.mean().compute()
    # plot
    X = bin_sizes[1:-3]
    Xerr = np.zeros(X.shape)
    Y = deltas.mean(1)[1:-3]
    Yerr = np.zeros(Y.shape)
    # table = Table(X, Xerr, Y, Yerr, "Correlation Error per bin size", "bin size", r"$\sigma_x$")

    fig, ax = plt.subplots()
    # for i in [10, 13, 15, 110, 123, 400]:
    # ax.plot(bin_sizes[:-3], std_of_bins[:-3, i])
    # ax = plotData(ax, table, polyfit=0, fmt="x")
    # ax = setSpace(ax, table)
    plt.savefig("plots/CorrelationErrorPerBinSizeMean.pdf", dpi=500)
    # ax.plot(X, Y)


def calculateCorrelation(observables, tao):
    return data * da.roll(observables, tao)


def binData(data, bin_count, bin_size):
    return da.reshape(data, (bin_count, bin_size))


def exponentialFunction(beta, x):
    return beta[0] + beta[1] * np.exp(beta[2] * x + beta[3])


def f(beta, x):
    return beta[0] + beta[1] * np.exp(beta[2] * x + beta[3])


def binCorrelator(correlator, bin_size):
    bin_count = correlator.shape[0] / bin_size
    bined_correlator = da.reshape(data, (bin_count, bin_size, correlator.shape[-1]))
    bin_estimate = bined_correlator.mean(1)
    observable_errors = (
        da.sum((bin_estimate - bin_estmiate.mean(0) ** 2)) / (bin_count - 1) / bin_count
    ) ** (1 / 2)

    error = correlator_errors.mean()

    return bined_correlator.compute(), error.compute()

def correlator(sites: np.ndarray, tao: float) -> np.ndarray:
    return (sites*da.roll(sites, tao)).mean(1)

def binObservables(observables: np.ndarray, bin_size: int) -> np.ndarray:
    return da.reshape(observables, observables.shape[0]/bin_size, bin_size).mean(1)

def varianz(observables: np.ndarray) -> float:
    return ((observables - observables.mean())**2).sum()/(observables.shape[0]-1)

def error(observables: np.ndarray) -> float:
    return (varianz(observables)/observables.shape[0])**(1/2)

def delta_Energy():
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 2 ** (da.floor(da.log2(data.shape[0])))
    bin_count = length / bin_size
    print(f"{length} sites are going to be used")

    lattice_times = np.arange(0, 5)

    data = data[:length, :]

    correlated_correlator = np.zeros((lattice_times.shape[0], data.shape[0]))
    uncorrelated_correlator = np.zeros((lattice_times.shape[0], data.shape[0]))
    for i in lattice_times:
        correlated_correlator[i,:] = correlator(data, i)
        uncorrelated_correlator[i,:] = binObservables(correlated_correlator[i,:], 256)
    
    
    energies = da.zeros(20000)
    for i in range(energies.shape[0]):
        chosen_lattice_times_index_mask = np.choice(np.arange(lattice_times.shape[0]))
        chosen_lattice_times_index_mask = np.choice(np.arange(uncorrelated_correlator.shape[0]))
        chosen_uncorrelated_correlator = 
        curve_fit(exponentialFunction, lattice_times, np.choice(uncorrelated_correlator))
    



def calcCorrelation():
    bin_size = 256
    time_slice = 0.5

    # read and crop data
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 2 ** (da.floor(da.log2(data.shape[0])))
    bin_count = length / bin_size
    print(f"{length} sites are going to be used")

    data = data[:length, :]
    data = da.reshape(data, (bin_count, bin_size, data.shape[1]))

    lattice_time = np.arange(0, 5)
    time = lattice_time * time_slice
    correlator_estimate = da.zeros(lattice_time.shape[0], bin_count)
    correlator_estimate_error = da.zeros(lattice_time.shape[0])

    for i in lattice_time:
        # calc correltaion
        correlator = data * da.roll(data, i, 2)
        correlator_estimate[i], correlator_estimate_error[i] = binCorrelator(correlator, bin_size)


    delta_energys = da.zeros(20000)

    resample_size = correlator_estimate.shape[1]

    for i in range(0, delta_energys.shape[0]):
        correlator_estimate_sample_estimate = da.zeros(lattice_time.shape[0])
        for i in lattice_time:
            correlator_estimate_sample_estimate[i] = da.choice(
                correlator_estimate[i], size=resample_size, replace=True
            ).mean()

        curve_fit(exponentialFunction, lattice_time, correlator_estimate_sample_estimate)

    d_E = delta_energys.mean()
    print(d_E)

    # fig, ax = plt.subplots()

    # ax.errorbar(time, correlator_estimate, correlator_estimate_error)
    # ax.set_yscale("log")
    # print(correlator_estimate_error)

    # plt.show()


# autoCorrelationErrorPerBinSize()
calcCorrelation()
