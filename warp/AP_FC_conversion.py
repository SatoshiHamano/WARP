# -*- coding:utf-8 -*-


import sys, os, datetime
import numpy as np
import math
import scipy.optimize
from warp.aperture import apertureSet

__version__ = "1.2"


# made by hamano (2014-02-20)
# revised by hamano (2014-2-24): ver 2
#   add a "hedit" command to define the dispersion axis before "transform".
# parameters are adjusted by hamano (2014-08-13)
#   using the comparison frame obtained in 2014-08-12


# usage:
#   python AP_FC_conversion.py <file1>

# <file1>
#   FITS file with which aperture trace functions are obtained with "aptrace".
#   This script will read "database/ap<file1>".
#

###################################################################
### parameters for the angle of tilted slit image as a function of
### detector coordinate (x,y).
### These parameters should be edited if they are changed.
###################################################################
### theta[rad] = a*ynorm + b
### ynorm = (2*y - (ymin + ymax))/(ymax-ymin)
#############################################################
########################################################


def chebyshev(n, x):
    p = [0. for i in range(n + 1)]
    p[0] = 1.
    if n >= 1:
        p[1] = x
    if n >= 2:
        for i in range(2, n + 1):
            p[i] = 2 * x * p[i - 1] - p[i - 2]
    return p[n]


def legendre(n, x):
    p = [0. for i in range(n + 1)]
    p[0] = 1.
    if n >= 1:
        p[1] = x
    if n >= 2:
        for i in range(2, n + 1):
            p[i] = ((2 * i - 1) * x * p[i - 1] - (i - 1) * p[i - 2]) / i
    return p[n]


def residue2(pk, normyarray, func, p, jj, ftype):
    tmpp = p
    tmpp[jj] = pk[0]
    ymax = len(normyarray)
    x = np.zeros(ymax)
    if ftype == 1:
        for k in range(len(p)):
            x += chebyshev(k, normyarray) * tmpp[k]
    elif ftype == 2:
        for k in range(len(p)):
            x += legendre(k, normyarray) * tmpp[k]

    res = (func - x)
    return (res)


def AP_to_FC(apfname, p0_ynorm, p0_order, iterate=3):
    [t0, t1, t2], [ap0, ap1, ap2] = p0_ynorm, p0_order

    # read parameters defined with aptrace

    print("Reading the aperture parameters...\n")
    apset = apertureSet(apfname)

    # reproduce the function of aptrace using read parameters

    rymin = 1
    rymax = 2048
    dry = rymax - rymin + 1
    print("Reset aperture function...\n")

    # param = []
    x0 = {}
    x1 = {}
    for i in apset.echelleOrders:
        ap = apset.apertures[i]
        ap.ymin = 1.
        ap.ymax = 2048.
        ap.adjustParameter(ap.tracex, ap.tracey)
        ap.traceDerivative()
        x0[i] = ap.tracex[0]
        x1[i] = ap.tracex[-1]

    # cut input data with imcopy

    xmargin = 100
    xmin = []
    xmax = []
    for i in apset.echelleOrders:
        xmin.append(1)
        if int(x0[i] - xmargin) < 1:
            ltv1 = 0
        else:
            ltv1 = 1 - int(x0[i] - xmargin)
        apset.apertures[i].apertureCenterX += ltv1
        apset.apertures[i].calculateTrace()
        if int(x1[i] + xmargin) > dry:
            xmax.append(dry + ltv1)
        else:
            xmax.append(int(x1[i] + xmargin) + ltv1)

    # calcurate a matrix for the coordinate transformation

    print("Calculating transform matrix...\n")

    addn = 2
    iterate2 = 3

    fctransmatrix = [
        [[1. for k in range(apset.apertures[m].yorder + addn)], [-1. for k in range(apset.apertures[m].yorder + addn)]]
        for m in apset.echelleOrders]

    for i in range(len(apset.echelleOrders)):
        m = apset.echelleOrders[i]
        ap = apset.apertures[m]
        normyarray = ap.traceNormy
        tantheta = np.tan((t2 * normyarray ** 2 + t1 * normyarray + t0) * (ap2 * m ** 2 + ap1 * m + ap0) / 180. * math.pi)

        print(
            "m=%d" % apset.echelleOrders[i], xmin[i], xmax[i], min(ap.tracex), max(ap.tracex), (min(ap.tracex) + max(ap.tracex)) / 2.,
            ap.apertureCenterX, ap.tracexDeriv.shape, np.average(ap.tracexDeriv), np.average(tantheta))

        p1func = -0.5 * (xmax[i] - xmin[i]) * tantheta / (1.0 - np.array(ap.tracexDeriv) * tantheta)
        param0 = [fctransmatrix[i][1][k] for k in range(ap.yorder + addn)]
        for ite in range(iterate2):
            for k in range(len(param0)):
                pk = np.array(param0[k])
                param_outputx1 = scipy.optimize.leastsq(residue2, pk, args=(normyarray, p1func, param0, k, ap.functionTypeID),
                                                        full_output=True)
                param0[k] = param_outputx1[0][0]
        for k in range(ap.yorder + addn):
            fctransmatrix[i][1][k] = param0[k]

        p0func = (normyarray * (rymax - rymin) + rymin + rymax) * 0.5 + (
                (xmin[i] + xmax[i]) * (-0.5) + ap.tracex) * tantheta / (
                         1.0 - ap.tracexDeriv * tantheta)
        param0 = [fctransmatrix[i][0][k] for k in range(ap.yorder + addn)]
        for ite in range(iterate2):
            for k in range(len(param0)):
                pk = np.array(param0[k])
                param_outputx0 = scipy.optimize.leastsq(residue2, pk, args=(normyarray, p0func, param0, k, ap.functionTypeID),
                                                        full_output=True)
                param0[k] = param_outputx0[0][0]

        for k in range(ap.yorder + addn):
            fctransmatrix[i][0][k] = param0[k]

        fcfile = open("database/fc%s_%d" % (apfname.rstrip("fits").rstrip("."), m), "w")

        fcfile.write("# %s\n" % datetime.datetime.today())
        fcfile.write("begin\t%s_%d\n" % (apfname.rstrip("fits").rstrip("."), m))
        fcfile.write("\ttask\tfitcoords\n")
        fcfile.write("\taxis\t2\n")
        fcfile.write("\tunits\tangstroms\n")
        fcfile.write("\tsurface\t%.0f\n" % (8. + 2 * (ap.yorder + addn)))
        fcfile.write("\t\t%.0f\n" % ap.functionTypeID)
        fcfile.write("\t\t2.\n")  # xorder
        fcfile.write("\t\t%.0f\n" % (ap.yorder + addn))
        fcfile.write("\t\t1.\n")
        fcfile.write("\t\t%.0f.\n" % xmin[i])
        fcfile.write("\t\t%.0f.\n" % xmax[i])
        fcfile.write("\t\t%.0f.\n" % rymin)
        fcfile.write("\t\t%.0f.\n" % rymax)

        for k in range(ap.yorder + addn):
            for j in range(2):
                fcfile.write("\t\t%.8f\n" % fctransmatrix[i][j][k])

        fcfile.close()


if __name__ == "__main__":
    filename = sys.argv[1:]
    [t0, t1, t2] = [3.28119e+01, -3.19816e+00, 5.50413e-01]
    [ap0, ap1, ap2] = [7.26861e-01, 1.02516e-02, -9.48522e-05]

    AP_to_FC(filename[0], [t0, t1, t2], [ap0, ap1, ap2])
