import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dask.dataframe as dd
import dask.array as da
import scipy.odr
from olib import *

plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

#path = "data/data_fig5_1k_RS_100.csv"
path = "data/data_fig5_10k_RS_1k.csv"
#path = "data/data_fig5_100k.csv"

'''
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
'''

def autoCorrelationErrorPerBinSize():
    # read and crop data
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 2**(da.floor(da.log2(data.shape[0])))
    print(length)
    data = data[:length, :]
    print(data.shape)

    # calc correltaion
    correlation = data * da.roll(data, 1, 1)
    print(correlation)

    # rebin and calc error
    bin_sizes = 2**np.arange(1, int(np.log2(int(correlation.shape[0]))))
    #bin_sizes = np.arange(1, correlation.shape[0], 10)
    print(type(bin_sizes))
    std_of_bins = np.zeros((bin_sizes.shape[0], correlation.shape[1]))
    print(std_of_bins)
    for i in range(std_of_bins.shape[0]):
        print("i:", i)
        corr_shape = (correlation.shape[0]/bin_sizes[i], bin_sizes[i], correlation.shape[1])
        correlation_bined = da.reshape(correlation, corr_shape)
        correlation_bined_mean = correlation_bined.mean(1)
        std = da.std(correlation_bined_mean, 0)/correlation_bined.shape[0]
        print(std.dask)
        std_of_bins[i, :] = std.compute()[:]
        #std_of_bins[i] = std.mean().compute()


    #plot
    X = bin_sizes[:-3]
    Xerr = np.zeros(X.shape)
    Y = std_of_bins.mean(1)[:-3]
    Yerr = np.zeros(Y.shape)
    table = Table(X, Xerr, Y, Yerr, "Correlation Error per bin size", "bin size", r"$\sigma_x$")

    fig, ax = plt.subplots()
    #for i in [10, 13, 15, 110, 123, 400]:
        #ax.plot(bin_sizes[:-3], std_of_bins[:-3, i])
    ax = plotData(ax, table, polyfit=0, fmt="x")
    ax = setSpace(ax, table)
    plt.savefig("plots/CorrelationErrorPerBinSizeMean.pdf", dpi=500)

def calcCorrelation():
    # read and crop data
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 10**(da.floor(da.log10(data.shape[0])))
    print(length)
    data = data[:length, :]
    print(data.shape)

    fig, ax = plt.subplots()
    tao = np.arange(1, 5)
    dt = 0.5
    t = toa*dt
    x_zero = np.zeros(tao.shape[0])

    for i in range(1, 5):
        # calc correltaion
        correlation = data * da.roll(data, i, 1)
        x_zero = correlation.mean()

        # rebin and calc error
        bin_sizes = 10**np.arange(1, int(np.log10(int(correlation.shape[0]))))
        std_of_bins = np.zeros((bin_sizes.shape[0]))
        corr_shape = (correlation.shape[0]/bin_sizes[i], bin_sizes[i], correlation.shape[1])
        correlation_bined = da.reshape(correlation, corr_shape)
        correlation_bined_mean = correlation_bined.mean(1)
        std = da.std(correlation_bined_mean, 0)/correlation_bined.shape[0]
        err = std.mean().compute()
        ax.errorbar(t, correlation[])



    #plot
    #X = bin_sizes[:-3]
    #Xerr = np.zeros(X.shape)
    #Y = std_of_bins.mean(1)[:-3]
    #Yerr = np.zeros(Y.shape)
    #table = Table(X, Xerr, Y, Yerr, "Correlation Error per bin size", "bin size", r"$\sigma_x$")

    #fig, ax = plt.subplots()
    #for i in [10, 13, 15, 110, 123, 400]:
        #ax.plot(bin_sizes[:-3], std_of_bins[:-3, i])
    #ax = plotData(ax, table, polyfit=0, fmt="x")
    #ax = setSpace(ax, table)
    plt.savefig("plots/CorrelationErrorPerBinSizeMean.pdf", dpi=500)



autoCorrelationErrorPerBinSize()
