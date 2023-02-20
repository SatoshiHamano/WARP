# -*- coding:utf-8 -*-

import numpy as np
import math
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.io.fits as fits

from warp.Spec1Dtools import open_ec_specfiles


def gaussianfunc(p,x,c):
    y = p[0] * np.exp(-(x-c)**2 / (2 * p[1] **2) ) + p[2]
    return y

def residue(p,y,x,c):
    res = (y-gaussianfunc(p,x,c))
    return (res)

def resolution_measure(fitsf, linef, output, mode, slit):


    resolution_mode = {"WIDE": 28000, "HIRES-Y": 68000, "HIRES-J": 68000}
    center_lam = {"WIDE": 11000., "HIRES-Y":10300., "HIRES-J": 12300.}
    resolution_factor = {100: 1, 140: 1.2, 200: 1.5, 400: 2.5}


    base_res = resolution_mode[mode] / resolution_factor[slit]


    linefile = open(linef, "r")
    width = center_lam[mode] / base_res * 5.
    print(width)
    linelines = linefile.readlines()
    lineorders = []
    linewave = []
    linepos = []
    for i in range(len(linelines)):
        lineorders.append(int(linelines[i].split()[0]))
        linewave.append(float(linelines[i].split()[3]))
        linepos.append(float(linelines[i].split()[2]))

    orders = [i for i in range(lineorders[0],lineorders[-1]+1)]
    wave = [[] for i in range(len(orders))]
    pos = [[] for i in range(len(orders))]

    for i in range(len(linelines)):
        wave[orders.index(lineorders[i])].append(linewave[i])
        pos[orders.index(lineorders[i])].append(linepos[i])

    rfile = open(fitsf,"r")
    rlines = rfile.readlines()
    rfile.close()

    pp = PdfPages(output+".pdf")
    fig = plt.figure()

    wave_ch = []
    peak = []
    centerlams = []
    gwidth = []
    gwidth_err = []
    for l in range(len(rlines)):
        apfits = fits.open(rlines[l].split()[0])
        apdata = apfits[0].data
        aplength = int(apfits[0].header["NAXIS1"])
        apcrval1 = float(apfits[0].header["CRVAL1"])
        apcdelt1 = float(apfits[0].header["CDELT1"])
        apcrpix1 = float(apfits[0].header["CRPIX1"])

        alam = [apcrval1 + apcdelt1 * (a - apcrpix1 + 1.) for a in range(aplength)]
        alamarray = np.array([apcrval1 + apcdelt1 * (a - apcrpix1 + 1.) for a in range(aplength)])
        apfits.close()

        wave_ch.append([])
        peak.append([])
        centerlams.append([])
        gwidth.append([])
        gwidth_err.append([])
        for i in range(len(wave[l])):
            minlpix = int((((wave[l][i]-width) - apcrval1) / apcdelt1)+apcrpix1-1.)
            maxlpix = int((((wave[l][i]+width) - apcrval1) / apcdelt1)+apcrpix1-1.)
            if minlpix > 1 and maxlpix < aplength:
                p0 = [1.e+3, base_res/wave[l][i], np.median(apdata[minlpix:maxlpix])]
                offset = np.median(apdata[minlpix:maxlpix])
                param_output = scipy.optimize.leastsq(residue, p0, args=(apdata[minlpix:maxlpix], alamarray[minlpix:maxlpix], wave[l][i]), full_output=True)

                if wave[l][i] / math.fabs(param_output[0][1]*2.35) > base_res / 3. and wave[l][i] / math.fabs(param_output[0][1]*2.35) < base_res * 1.5:
                    wave_ch[l].append(wave[l][i])
                    peak[l].append(param_output[0][0])
                    # centerlams[l].append(param_output[0][1])
                    gwidth[l].append(param_output[0][1] * 2.35)

                    plt.scatter(alamarray[minlpix:maxlpix], apdata[minlpix:maxlpix])
                    plt.plot(alamarray[minlpix:maxlpix], gaussianfunc(param_output[0],alamarray[minlpix:maxlpix],wave[l][i]))
                    plt.title("%s: %.2f" % (rlines[l].split()[0],wave[l][i]))
                    plt.savefig(pp, format="pdf", dpi=200)
                    plt.clf()


    pp.close()



    wfile = open(output + ".dat", "w")
    for i in range(len(orders)):
        for j in range(len(wave_ch[i])):
            if orders[i] % 2 == 0:
                plt.scatter(wave[i][j], wave[i][j]/math.fabs(gwidth[i][j]), c="b", s=6.)
            if orders[i] % 2 == 1:
                plt.scatter(wave[i][j], wave[i][j]/math.fabs(gwidth[i][j]), c="g", s=6.)
            wfile.write("%d\t%.8f\t%.8f\t%.8f\n" % (orders[i], pos[i][j], wave[i][j], math.fabs(gwidth[i][j])))#, gwidth_err[i][j]))
    #    wfile.write("\n####################\n\n")

    plt.xlabel("Wavelength")
    plt.ylabel("R")
    plt.ylim(base_res/2.5,base_res*1.2)
    plt.grid()
    plt.savefig(output + ".png", format="png")

    wfile.close()


