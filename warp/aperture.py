#!/usr/bin/env python
# -*- coding:utf-8 -*-
import sys

import numpy as np
import time
from scipy import interpolate
import scipy.optimize
import copy


class apertureSet:
    def __init__(self, apref, arrayLength=2048, reduceFullData=True, selectedOrders=[]):
        if apref.find(".fits") != -1:
            apref = apref.rstrip("fits").rstrip(".")

        # open file

        apfile = open("./database/ap" + apref, "r")
        aplines = apfile.readlines()
        apfile.close()

        # search the last definition in the aperture file

        apnum = []
        beginline = []
        for i in range(len(aplines)):
            if aplines[i].find("begin") != -1:
                apnum.append(int(aplines[i + 2].split()[1]))
                beginline.append(i)
        beginline.append(len(aplines))

        # read parameters from files

        coefficients = [[] for i in range(len(apnum))]
        curve = 0
        centerx = []
        centery = []
        aplow = []
        aphigh = []
        for i in range(len(beginline) - 1):
            for j in range(beginline[i], beginline[i + 1]):
                if aplines[j].find("center") != -1:
                    centerx.append(float(aplines[j].split()[1]))
                    centery.append(float(aplines[j].split()[2]))
                    aplow.append(float(aplines[j + 1].split()[1]))
                    aphigh.append(float(aplines[j + 2].split()[1]))
                if aplines[j].find("curve") != -1:
                    curve = int(aplines[j].split()[1])
                    for k in range(1, curve + 1):
                        if k <= 2:
                            coefficients[i].append(int(float(aplines[j + k].split()[0])))
                        else:
                            coefficients[i].append(float(aplines[j + k].split()[0]))

        ftype = [coefficients[i][0] for i in range(len(apnum))]
        yorder = [coefficients[i][1] for i in range(len(apnum))]
        ymin = [coefficients[i][2] for i in range(len(apnum))]
        ymax = [coefficients[i][3] for i in range(len(apnum))]
        cm = [coefficients[i][4:4 + yorder[i]] for i in range(len(apnum))]

        self.echelleOrders = apnum
        self.arrayLength = arrayLength
        self.apertures = {}
        if reduceFullData == True:
            for i in range(len(self.echelleOrders)):
                self.apertures[apnum[i]] = self.aperture(apnum[i], centerx[i], centery[i], aplow[i], aphigh[i],
                                                         ftype[i], yorder[i], ymin[i], ymax[i], cm[i],
                                                         arrayLength=self.arrayLength)
        else:
            resetEchelleOrders = []
            for i in range(len(self.echelleOrders)):
                if self.echelleOrders[i] in selectedOrders:
                    self.apertures[apnum[i]] = self.aperture(apnum[i], centerx[i], centery[i], aplow[i], aphigh[i],
                                                             ftype[i], yorder[i], ymin[i], ymax[i], cm[i],
                                                             arrayLength=self.arrayLength)
                    resetEchelleOrders.append(self.echelleOrders[i])
            self.echelleOrders = resetEchelleOrders
            if len(self.echelleOrders) == 0:
                print("Error: No echelle orders were selected. Please change conf.reduceFullData to True, or "
                      "input effective echelle order numbers into conf.selectedOrders.")
                sys.exit()

    def apmaskArray(self, lowlim="INDEF", upplim="INDEF", margin=10):
        x = np.array([[i + 1. for i in range(self.arrayLength)] for j in range(self.arrayLength)])
        y = np.array([[j + 1. for i in range(self.arrayLength)] for j in range(self.arrayLength)])
        mask = np.zeros((self.arrayLength, self.arrayLength))

        xlim1, xlim2 = margin, self.arrayLength - margin
        ylim1, ylim2 = margin, self.arrayLength - margin

        aplow = {}
        if lowlim != "INDEF":
            for m in self.echelleOrders:
                aplow[m] = lowlim
        else:
            for m in self.echelleOrders:
                aplow[m] = self.apertures[m].apLow

        aphigh = {}
        if lowlim != "INDEF":
            for m in self.echelleOrders:
                aphigh[m] = upplim
        else:
            for m in self.echelleOrders:
                aphigh[m] = self.apertures[m].apHigh

        for m in self.echelleOrders:
            apxarray = np.tile(self.apertures[m].tracex, (self.arrayLength, 1)).T
            residue = np.absolute(x - apxarray)
            mask[(aplow[m] < residue) & (residue < aphigh[m])] = self.apertures[m].order

        mask[x < xlim1] = 0
        mask[x > xlim2] = 0
        mask[y < ylim1] = 0
        mask[y > ylim2] = 0

        return mask

    def slitcoordArray(self, lowlim="INDEF", upplim="INDEF", interOrderValue=-10000):
        x = np.array([[i + 1. for i in range(self.arrayLength)] for j in range(self.arrayLength)])
        y = np.array([[j + 1. for i in range(self.arrayLength)] for j in range(self.arrayLength)])
        slitcoord = np.zeros((self.arrayLength, self.arrayLength)) + interOrderValue

        aplow = {}
        if lowlim != "INDEF":
            for m in self.echelleOrders:
                aplow[m] = lowlim
        else:
            for m in self.echelleOrders:
                aplow[m] = self.apertures[m].apLow

        aphigh = {}
        if lowlim != "INDEF":
            for m in self.echelleOrders:
                aphigh[m] = upplim
        else:
            for m in self.echelleOrders:
                aphigh[m] = self.apertures[m].apHigh

        for m in self.echelleOrders:
            mask = np.zeros((self.arrayLength, self.arrayLength))
            apxarray = np.tile(self.apertures[m].tracex, (self.arrayLength, 1)).T
            residue = x - apxarray
            residue_abs = np.absolute(x - apxarray)
            mask[(aplow[m] < residue_abs) & (residue_abs < aphigh[m])] = 1
            slitcoord[mask == 1] += residue[mask == 1] - interOrderValue

        return slitcoord

    def write_apdatabase(self, apfile, bgsample):
        # What time is it now?

        timeline = time.strftime("# %a %H:%M:%S %d-%b-%Y\n", time.gmtime())

        # if apfile is string, open the file here.

        if apfile.find(".fits") != -1:
            apfile = apfile.rstrip("fits").rstrip(".")

        if isinstance(apfile, str):
            dbwf = open("database/ap%s" % apfile, "w")

        for i in range(len(self.echelleOrders)):
            m = self.echelleOrders[i]

            # if apfile is list, open the files here.
            if isinstance(apfile, list):
                dbwf = open("database/ap%s" % apfile[i], "w")

            # write the current time in the file first
            dbwf.write(timeline)

            # determine xmin and xmax
            apcoord = [self.apertures[m].apLow, self.apertures[m].apHigh]
            if bgsample[i].find("INDEF") == -1:
                for j in range(len(bgsample[i].split(","))):
                    for k in range(len(bgsample[i].split(",")[j].split(":"))):
                        apcoord.append(float(bgsample[i].split(",")[j].split(":")[k]))

            xmin = min(apcoord)
            xmax = max(apcoord)

            # write the aperture parameters

            if isinstance(apfile, list):
                dbwf.write("begin\taperture\t%s %d %.2f %.2f\n" % (
                    apfile[m], m, self.apertures[m].apertureCenterX, self.apertures[m].apertureCenterY))
                dbwf.write("\timage\t%s\n" % (apfile[i]))

            if isinstance(apfile, str):
                dbwf.write("begin\taperture\t%s %d %.2f %.2f\n" % (
                    apfile, m, self.apertures[m].apertureCenterX, self.apertures[m].apertureCenterY))
                dbwf.write("\timage\t%s\n" % (apfile))
                dbwf.write("\taperture\t%d\n" % m)
                dbwf.write("\tbeam\t%d\n" % (m))
                dbwf.write(
                    "\tcenter\t%.2f %f\n" % (self.apertures[m].apertureCenterX, self.apertures[m].apertureCenterY))
                dbwf.write("\tlow\t%.2f %f\n" % (self.apertures[m].apLow, -self.apertures[m].apertureCenterY + 1.))
                dbwf.write("\thigh\t%.2f %f\n" % (self.apertures[m].apHigh, self.apertures[m].apertureCenterY + 1.))
                dbwf.write("\tbackground\n")
                dbwf.write("\t\txmin %.2f\n" % (xmin))
                dbwf.write("\t\txmax %.2f\n" % (xmax))
                dbwf.write("\t\tfunction chebyshev\n")
                dbwf.write("\t\torder 1\n")
                dbwf.write("\t\tsample %s\n" % (bgsample[i]))
                dbwf.write("\t\tnaverage -3\n")
                dbwf.write("\t\tniterate 0\n")
                dbwf.write("\t\tlow_reject 3.\n")
                dbwf.write("\t\thigh_reject 3.\n")
                dbwf.write("\t\tgrow 0.\n")
                dbwf.write("\taxis 1\n")
                dbwf.write("\tcurve %d\n" % (4 + self.apertures[m].yorder))
                dbwf.write("\t\t%f\n" % self.apertures[m].functionTypeID)
                dbwf.write("\t\t%f\n" % self.apertures[m].yorder)
                dbwf.write("\t\t%f\n" % self.apertures[m].ymin)
                dbwf.write("\t\t%f\n" % self.apertures[m].ymax)
                for j in range(self.apertures[m].yorder):
                    dbwf.write("\t\t%f\n" % self.apertures[m].coefficient[j])
                dbwf.write("\n")

                # close the file

            if isinstance(apfile, list):
                dbwf.close()

        if isinstance(apfile, str):
            dbwf.close()

    def renewOrders(self, newaplist):
        newEchelleOrders = []
        for m in self.echelleOrders:
            if m in newaplist:
                newEchelleOrders.append(m)
            else:
                del self.apertures[m]
                print("m={} aperture was deleted from the apertureSet object.".format(m))
        self.echelleOrders = newEchelleOrders

    class aperture:
        def __init__(self, m: int, centerx: float, centery: float, aplow: float, aphigh: float, ftype: int, yorder: int,
                     ymin: float, ymax: float, cm: list, arrayLength=2048):
            self.order = m
            self.apertureCenterX = centerx
            self.apertureCenterY = centery
            self.apLow = aplow
            self.apHigh = aphigh
            self.functionTypeID = ftype
            if ftype == 1:
                self.function = self.chebyshev
                self.functionDeriv = self.chebyshevDeriv
            elif ftype == 2:
                self.function = self.legendre
                self.functionDeriv = self.legendreDeriv
            self.yorder = yorder
            self.ymin = ymin
            self.ymax = ymax
            self.coefficient = cm
            self.arrayLength = arrayLength
            self.calculateTrace()

        def calculateTrace(self, changeParameter=False, changedParam="INDEF", changedCenter="INDEF"):
            # define the array, in which the aperture function will be stored

            self.tracey = np.arange(self.arrayLength) + 1.
            if changeParameter:
                self.tracex = np.zeros(self.arrayLength) + changedCenter
            else:
                self.tracex = np.zeros(self.arrayLength) + self.apertureCenterX

            # calculate the functions
            self.traceNormy = (2 * self.tracey - self.ymin - self.ymax) / (self.ymax - self.ymin)
            if changeParameter:
                for j in range(self.yorder):
                    self.tracex += self.function(j) * changedParam[j]
            else:
                for j in range(self.yorder):
                    self.tracex += self.function(j) * self.coefficient[j]

            self.tracexfunc = interpolate.interp1d(self.tracey, self.tracex, kind="cubic")

        def chebyshev(self, n):
            p = [0. for i in range(n + 1)]
            p[0] = 1.
            if n >= 1:
                p[1] = self.traceNormy
            if n >= 2:
                for i in range(2, n + 1):
                    p[i] = 2 * self.traceNormy * p[i - 1] - p[i - 2]
            return p[n]

        def legendre(self, n):
            p = [0. for i in range(n + 1)]
            p[0] = 1.
            if n >= 1:
                p[1] = self.traceNormy
            if n >= 2:
                for i in range(2, n + 1):
                    p[i] = ((2 * i - 1) * self.traceNormy * p[i - 1] - (i - 1) * p[i - 2]) / i
            return p[n]

        def make_apnumtext(self):
            if self.order < 10:
                ap = "000%d" % self.order
            elif self.order < 100:
                ap = "00%d" % self.order
            elif self.order < 1000:
                ap = "0%d" % self.order
            else:
                ap = "%d" % self.order
            return ap

        def adjustParameter(self, x, y, iteration=3, clipsig=5, clip=True, printDetail=True):
            def residue_trace(p, xdata, ydata):
                self.calculateTrace(changeParameter=True, changedParam=p, changedCenter=self.apertureCenterX)
                res = (xdata - self.tracexfunc(ydata))
                return (res)

            clipbool = np.array([True for _ in range(len(x))])
            p0 = copy.copy(self.coefficient)
            self.calculateTrace()
            for n in range(iteration):
                param_output1 = scipy.optimize.leastsq(residue_trace, p0, args=(x, y), full_output=True)
                p0 = param_output1[0]
                self.coefficient = copy.copy(p0)
                self.calculateTrace()
                residueFit = x - self.tracexfunc(y)
                residueAve = np.average(residueFit[clipbool])
                residueStd = np.std(residueFit[clipbool])
                if clip:
                    clipbool[np.absolute(residueFit - residueAve) > clipsig * residueStd] = False
                    if printDetail:
                        print("Iteration {}: {} points were clipped. (std={}, {} sig)".format(
                            n + 1, len(x) - np.sum(clipbool), residueStd, clipsig))

            self.calculateTrace()
            if clip:
                return clipbool

        def traceDerivative(self):
            self.traceNormy = (2 * self.tracey - self.ymin - self.ymax) / (self.ymax - self.ymin)
            self.tracexDeriv = np.zeros(self.arrayLength)
            for i in range(self.yorder):
                self.tracexDeriv += self.functionDeriv(i) * self.coefficient[i] * 2. / (self.ymax - self.ymin)

        def chebyshevDeriv(self, n):
            p = [0. for i in range(n + 1)]
            p[0] = 1.
            if n >= 1:
                p[1] = self.traceNormy
            if n >= 2:
                for i in range(2, n + 1):
                    p[i] = 2 * self.traceNormy * p[i - 1] - p[i - 2]

            p_deriv = [0. for i in range(n + 1)]
            p_deriv[0] = 0.
            if n >= 1:
                p_deriv[1] = 1.
            if n >= 2:
                for i in range(2, n + 1):
                    p_deriv[i] = 2 * p[i - 1] + 2 * self.traceNormy * p_deriv[i - 1] - p_deriv[i - 2]

            return p_deriv[n]

        def legendreDeriv(self, n):
            p = [0. for i in range(n + 1)]
            p[0] = 1.
            if n >= 1:
                p[1] = self.traceNormy
            if n >= 2:
                for i in range(2, n + 1):
                    p[i] = ((2 * i - 1) * self.traceNormy * p[i - 1] - (i - 1) * p[i - 2]) / i

            p_deriv = [0. for i in range(n + 1)]
            p_deriv[0] = 0.
            if n >= 1:
                p_deriv[1] = 1.
            if n >= 2:
                for i in range(2, n + 1):
                    p_deriv[i] = 2 * p[i - 1] + 2 * self.traceNormy * p_deriv[i - 1] - p_deriv[i - 2]

            return p_deriv[n]
