import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dask.dataframe as dd
import scipy.odr
from olib import *

plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

#path = "data/data_fig5_1k_RS_100.csv"
path = "data/data_fig5_10k_RS_1k.csv"

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

def autoCorrelationErrorPerBinSize():
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    correlation = data * np.roll(data, 1, 1)

    std_of_bins = np.zeros(int(np.log10(correlation.shape[0]))+1)
    for i in range(std_of_bins.shape[0]):
        print("i:", i)
        print(std_of_bins.shape[0]-i)
        correlation_bined = correlation.reshape(10**(std_of_bins.shape[0]-i-1), 10**i, correlation.shape[1])
        print(correlation_bined.shape)
        correlation_bined_mean = correlation_bined.mean(1)
        print(correlation_bined_mean.shape)
        std = np.std(correlation_bined_mean, 0)
        std_of_bins[i] = std.mean()

    X = np.arange(0, std_of_bins.shape[0])
    Xerr = np.zeros_like(X)
    Y = std_of_bins
    Yerr = np.zeros_like(Y)
    table = Table(X, Xerr, Y, Yerr, "Correlation Error per bin size", "bin size", r"$\sigma_x$")

    fig, ax = plt.subplots()
    ax = plotData(ax, table, polyfit=0, fmt="x")
    ax = setSpace(ax, table)
    plt.savefig("plots/CorrelationErrorPerBinSize.pdf", dpi=500)


autoCorrelationErrorPerBinSize()
