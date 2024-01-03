import pdfkit
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class model:
    def __init__(self) -> None:
        pass

    def fit(self, X, Y):
        pass

    def predict(self, x) -> np.ndarray:
        pass


def roundToError(x, xerr=[""], digits=2):
    if xerr == [""]:
        xerr = x
    print(xerr.shape)
    error_digits = np.zeros_like(xerr).astype("int")
    print("esr", xerr, xerr!=0)
    print(error_digits[xerr!=0])
    error_digits[xerr!=0] = np.floor(np.log10(xerr[xerr!=0])).astype("int")
    #error_digits[error_digits**2==0] = 0
    roundedx = []
    roundedxerr = []
    for i in range(len(x)):
        roundedx.append(np.round(x[i], digits-error_digits[i]-1))
    for i in range(len(xerr)):
        roundedxerr.append(np.round(xerr[i], digits-error_digits[i]-1))
    return np.array(roundedx), np.array(roundedxerr)


class LinearRegression(model):
    def __init__(self):
        self.V_m = None
        self.V_n = None
        self.m = None
        self.n = None

        return None

    def weighted_average(self, Z, Zerr):
        return (Z/Zerr**2).sum()/(1/Zerr**2).sum()

    def fit(self, X, Y, Yerr):
        self.y_wa = self.weighted_average(Y, Yerr)
        self.xy_wa = self.weighted_average(X*Y, Yerr)
        self.x_wa = self.weighted_average(X, Yerr)
        self.x_sq_wa = self.weighted_average(X**2, Yerr)

        self.sigma_wa = Yerr.size/(1/Yerr**2).sum()
        self.V_m = self.sigma_wa/(X.size*(self.x_sq_wa-self.x_wa**2))
        self.V_n = self.sigma_wa / \
            (X.size*(self.x_sq_wa-self.x_wa**2))*self.x_sq_wa
        self.m = (self.xy_wa-self.x_wa*self.y_wa)/(self.x_sq_wa-self.x_wa**2)
        self.n = (self.x_sq_wa*self.y_wa-self.x_wa*self.xy_wa) / \
            (self.x_sq_wa-self.x_wa**2)

        return self

    def setParameters(self, m, V_m, n, V_n):
        self.m = m
        self.n = n
        self.V_m = V_m
        self.V_n = V_n

    def getParameters(self):
        return {"m": self.m, "V_m": self.V_m, "n": self.n, "V_n": self.V_n}
    
    def printParameter(self):
        m_error_digits = math.floor(math.log10(self.V_m**(1/2)))
        n_error_digits = math.floor(math.log10(self.V_n**(1/2)))
        print(f"m: {np.round(self.m, 1-m_error_digits)}+-{np.round(self.V_m**(1/2), 1-m_error_digits)}; n: {np.round(self.n, 1-n_error_digits)}+-{np.round(self.V_n**(1/2), 1-n_error_digits)}")
        print(f"m: {self.m}+-{self.V_m**(1/2)}; n: {self.n}+-{self.V_n**(1/2)}")

    def predict(self, x):
        return x*self.m+self.n

    def predict_most_and_least(self, x):
        return [x*(self.m+self.V_m**(1/2))+self.n+self.V_n**(1/2), x*(self.m-self.V_m**(1/2))+self.n-self.V_n**(1/2)]


def CSV_to_PDF(pandas_table, output_file_name):
    html_table = pandas_table.to_html()

    options = {'page-size': 'Letter', 'margin-top': '0mm',
               'margin-right': '0mm', 'margin-bottom': '0mm', 'margin-left': '0mm'}

    pdfkit.from_string(html_table, output_file_name, options=options)


def setSpace(axis, allX, allY, title="plot", xlabel="x", ylabel="y", border=0.1, resolution=0.1):
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_title(title)

    xRange = np.max(allX)-np.min(allX)
    yRange = np.max(allY)-np.min(allY)
    print(yRange)

    xstep = np.round((xRange * resolution)*2, -1)/2
    print(yRange)
    ystep = np.round((yRange * resolution)*2, -1)/2
    print(ystep)

    xLimits = [math.floor((np.min(allX)-xRange*border)/10)*10,
               math.ceil((np.max(allX)+xRange*border+xstep)/10)*10]
    yLimits = [math.floor((np.min(allY)-yRange*border)/10)*10,
               math.ceil((np.max(allY)+yRange*border+ystep)/10)*10]
    print(yLimits[0], yLimits[1], ystep)
    axis.set_xticks(np.arange(xLimits[0], xLimits[1], xstep))
    axis.set_xticks(np.arange(xLimits[0], xLimits[1], xstep * 0.2), minor=True)
    axis.set_yticks(np.arange(yLimits[0], yLimits[1], ystep))
    axis.set_yticks(np.arange(yLimits[0], yLimits[1], ystep * 0.2), minor=True)

    axis.set_xlim(xLimits)
    axis.set_ylim(yLimits)

    axis.grid(which="both")
    axis.grid(which="minor", alpha=0.4)
    axis.grid(which="major", alpha=0.7)

    return axis


def plotLine(axis, x, y, label=None, linewidth=0.7, linestyle="solid", color="r"):
    axis.plot(x, y, linewidth=linewidth, label=label,
              linestyle=linestyle, color=color)

    return axis


def plotData(axis, X, Xerr, Y, Yerr, label=None, capsize=1, elinewidth=0.5, fmt=",", polyfit=1, color="r", xscaleing=1, yscaleing=1, errorbar=True):
    X = X*xscaleing
    Xerr = Xerr*xscaleing
    Y = Y*yscaleing
    Yerr = Yerr*yscaleing
    if errorbar:
        axis.errorbar(X, Y, Yerr, Xerr, label=label, capsize=capsize,
                      elinewidth=elinewidth, fmt=fmt, color=color)
    else:
        axis.scatter(X, Y, label=label, marker="x", color=color, linewidth=0.6)

    model = None

    if polyfit != 0:
        xrange = np.max(X)-np.min(X)
        x = np.linspace(np.min(X)-xrange*0.1, np.max(X)+xrange*0.1, 100)

        if polyfit == 1:
            model = LinearRegression().fit(X, Y, Yerr)
            y_most_and_least = model.predict_most_and_least(x)
            plotLine(axis, x, y_most_and_least[0],
                     linestyle="dashdot", color=color, linewidth=1)
            plotLine(axis, x, y_most_and_least[1],
                     linestyle="dashdot", color=color, linewidth=1)
            y = model.predict(x)
        elif polyfit > 1:
            y = np.poly1d(np.polyfit(X, Y, polyfit))(x)

        plotLine(axis, x, y, label=label, color=color, linewidth=1)

    return axis, model
