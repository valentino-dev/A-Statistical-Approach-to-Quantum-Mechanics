import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dask.dataframe as dd
import dask.array as da
import scipy.odr
import sys

sys.path.append("/home/angelo/Git/olib")
from olib import *

plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

# path = "data/data_fig5_1k_RS_100.csv"
path = "data/data_fig5_10k_RS_1k.csv"
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
    correlation = da.mean(data * da.roll(data, 1, 1), 1).compute()
    print("Configurations: ", correlation.shape[0])

    # rebin and calc error
    bin_sizes = 2 ** np.arange(1, int(np.log2(int(data.shape[0]))))
    deltas = np.zeros(bin_sizes.shape[0])
    for i in range(deltas.shape[0]):
        print("i: ", i)
        corr_shape = (
            int(correlation.shape[0] / bin_sizes[i]),
            bin_sizes[i]
        )
        print("Corr Shape: ", corr_shape)
        correlation_bined = da.reshape(correlation, corr_shape)
        cb_mean = da.mean(correlation_bined, 1)
        cbm_mean = da.mean(cb_mean)
        delta = (da.sum((cb_mean-cbm_mean)**2)/(cb_mean.shape[0]*(cb_mean.shape[0]-1)))**(1/2)

        print(deltas.shape[0] - i)
        deltas[i, :] = delta.compute()

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
    # plt.savefig("plots/CorrelationErrorPerBinSizeMean.pdf", dpi=500)
    ax.plot(X, Y)
    plt.show()

def f(beta, x):
    return beta[0]+beta[1]*np.exp(beta[2]*x+beta[3])


def calcCorrelation():
    # read and crop data
    data = dd.read_csv(path, sep=";", header=None).to_dask_array(lengths=True)
    length = 10 ** (da.floor(da.log10(data.shape[0])))
    print(length)
    data = data[:length, :]
    print(data.shape)

    tao = np.arange(0, 5)
    dt = 0.5
    t = tao * dt
    x_zero = np.zeros(tao.shape[0])
    x_zero_autocorr = np.zeros(tao.shape[0])

    for i in range(0, 5):
        # calc correltaion
        correlation = data * da.roll(data, i, 1)
        x_zero[i] = correlation.mean()
        x_zero_autocorr[i] = correlation.std()/length**(1/2)

    X = t
    Xerr = np.zeros(X.shape)+1
    Y = x_zero
    Yerr = x_zero_autocorr
    print(X, Y, Xerr, Yerr)

    scipymodel = scipy.odr.Model(f)
    fitData = scipy.odr.RealData(X, Y, sy=Yerr)
    odr = scipy.odr.ODR(fitData, scipymodel, beta0=[0, 1, -1, 0])
    output = odr.run()
    #output.pprint()
    

    print(output.beta)

    

    fig, ax = plt.subplots()
    ax.errorbar(X, Y, Yerr, fmt=",")
    ax.plot(X, f(output.beta, X))
    ax.set_yscale("log")
    plt.show()


autoCorrelationErrorPerBinSize()
#calcCorrelation()
