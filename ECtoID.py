# -*- coding:utf-8 -*-

__version__ = "1.4"

from pyraf import iraf
import math
import sys
import numpy as np

iraf.imred()
iraf.echelle()


# usage:
# python ECtoID7.py <eclist> <comp>

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

def make_apnumtext(apnum):
    if apnum < 10:
        ap = "000%d" % apnum
    elif apnum < 100:
        ap = "00%d" % apnum
    elif apnum < 1000:
        ap = "0%d" % apnum
    else:
        ap = "%d" % apnum
    return ap


def readparamEC(ecfname):
    ecfile = open("database/ec" + ecfname.rstrip("fits").rstrip("."), "r")
    ecdata = ecfile.readlines()
    ecfile.close()

    pbegin = []
    for i in range(len(ecdata)):
        if ecdata[i].find("begin") != -1:
            pbegin.append(i)

    begin = pbegin[-1]

    slope = 1
    offset_flag = False
    for i in range(begin, len(ecdata)):
        l = ecdata[i].split()
        if ecdata[i].find("offset") != -1:
            last_fea = ecdata[i - 1].split()
            offset = int(l[1])
            offset_flag = True
        if ecdata[i].find("slope") != -1:
            slope = int(l[1])
        if ecdata[i].find("coefficients") != -1:
            ncoeff = int(l[1])
            coefficients = [0. for k in range(ncoeff)]
            for j in range(ncoeff):
                coefficients[j] = float(ecdata[i + j + 1])
            break
        if ecdata[i].find("image") != -1:
            image = l[1]

    if not offset_flag:
        offset = 0

    ftype = round(coefficients[0])
    xpow = round(coefficients[1])
    opow = round(coefficients[2])
    xterms = round(coefficients[3])
    xmin = round(coefficients[4])
    xmax = round(coefficients[5])
    omin = round(coefficients[6])
    omax = round(coefficients[7])

    cmn = [[coefficients[8 + j + i * xpow] for i in range(opow)] for j in range(xpow)]

    apn = list(range(omin, omax + 1))
    apn_num = len(apn)

    features = [[] for i in range(apn_num)]
    nfap = [0 for i in range(apn_num)]
    for i in range(begin, len(ecdata)):
        if ecdata[i].find("features") != -1:
            l = ecdata[i].split()
            nf = int(l[1])
            for j in range(i + 1, i + 1 + nf):
                if ecdata[j].find("INDEF") == -1:
                    ll = ecdata[j].split()
                    features[int(ll[0]) - (offset + slope * apn[0])].append(
                        [float(ll[2]), float(ll[3]), float(ll[4]), float(ll[5]), int(ll[6]), int(ll[7])])
                    nfap[int(ll[0]) - (offset + slope * apn[0])] += 1

    return offset, slope, apn, ftype, xpow, opow, xterms, xmin, xmax, omin, omax, cmn, features, nfap


def calcucmnid(offset, slope, apn, ftype, xpow, xmin, xmax, opow, omin, omax, cmn):
    onorm = [((offset + slope * i) * 2 - (omin + omax)) / float(omax - omin) for i in apn]
    cmnid = [[0. for i in range(len(apn))] for j in range(xpow)]

    if ftype == 1:
        for n in range(len(apn)):
            for i in range(xpow):
                for j in range(opow):
                    cmnid[i][n] += cmn[i][j] * chebyshev(j, onorm[n]) / float(offset + slope * apn[n])
    elif ftype == 2:
        for n in range(len(apn)):
            for i in range(xpow):
                for j in range(opow):
                    cmnid[i][n] += cmn[i][j] * legendre(j, onorm[n]) / float(offset + slope * apn[n])

    return cmnid


def ECtoID(ecfname):
    offset, slope, apn, ftype, xpow, opow, xterms, xmin, xmax, omin, omax, cmn, features, nfap = readparamEC(ecfname)

    cmnid = calcucmnid(offset, slope, apn, ftype, xpow, xmin, xmax, opow, omin, omax, cmn)

    apnum = []
    for i in range(len(apn)):
        ap = make_apnumtext(apn[i])
        apnum.append(ap)

        fname = ecfname.rstrip("fits").rstrip(".") + ".%s" % ap
        idfile = open("database/id" + fname, "a+")
        idfile.write("begin\tidentify %s\n\tid \t %s\n" % (fname, fname))
        idfile.write("\ttask\tidentify\n\timage\t%s\n\tunits\tAngstroms\n" % (fname))
        idfile.write("\tfeatures\t%d\n" % (nfap[i]))
        for j in range(nfap[i]):
            idfile.write("\t\t%.2f\t%.4f\t%.4f\t%.1f\t%d\t%d\n" % (
                features[i][j][0], features[i][j][1], features[i][j][2], features[i][j][3], features[i][j][4],
                features[i][j][5]))
        idfile.write("\tcoefficients\t%d\n" % (xpow + 4))
        idfile.write("\t\t%d.\n\t\t%d.\n\t\t%f\n\t\t%f\n" % (ftype, xpow, xmin, xmax))
        for j in range(xpow):
            idfile.write("\t\t%f\n" % cmnid[j][i])
        idfile.write("\n")
        idfile.close()

    return features, nfap, apnum


if __name__ == "__main__":
    filename = sys.argv[1:]

    iraf.ecidentify(filename[0])

    features, nfap, apnum = ECtoID(filename[0])

    fitlam_all = []
    reflam_all = []
    for i in range(len(nfap)):
        fitlam = []
        reflam = []
        for j in range(len(features[i])):
            fitlam.append(features[i][j][1])
            fitlam_all.append(features[i][j][1])
            reflam.append(features[i][j][2])
            reflam_all.append(features[i][j][2])
        fitlam = np.array(fitlam)
        reflam = np.array(reflam)

        iraf.hedit(filename[0], "DISP%s" % apnum[i],
                   "%d %.4f %.5f" % (nfap[i], np.average(fitlam), np.std(fitlam - reflam)), add="yes",
                   verify="no",
                   update="yes")

    fitlam_all = np.array(fitlam_all)
    reflam_all = np.array(reflam_all)
    iraf.hedit(filename[0], "DISPALL",
               "%d %.4f %.5f" % (np.sum(nfap), np.average(fitlam_all), np.std(fitlam_all - reflam_all)),
               add="yes", verify="no", update="yes")
