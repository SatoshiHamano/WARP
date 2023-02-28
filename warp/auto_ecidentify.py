
# -*- coding:utf-8 -*-

import numpy
import astropy.io.fits as pyfits
from pyraf import iraf
import scipy.optimize
from scipy.signal import find_peaks_cwt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, os, datetime
import math
from warp.Spec1Dtools import open_ec_specfiles

__version__ = "1.4"

#Description:
#
#   Execute IRAF/ecidentify automatically.
#
#
#Usage:
#
#   $ python Auto_ecidentify.py ecfile slit linelist
#
#
#Output:
#
#   Database file of ecidentify
#
#Updates:
#
#   ver1.2 updated by Hamano in 2019.10.30
#       - open_ec_specfiles was moved to Spec1Dtools.py.
#       - 140 um slit was added.
#       - refactoring.
#
#   ver1.3 updated by Hamano in 2020.01.13
#       - bug fix (caused by the refrenece line whose position get out from the array range.)
#
#   ver1.4 updated by Hamano in 2020.02.16
#       - algorithm to identify registered emission lines are improved.
#
##


def peak_center_single_line(spx, spy, width, c_index, lowlim=0.25, upplim=1.):
    # input: (x,y) of spectrum and roughly estimated peak positions of emission lines.
    # output: peak positions

    high_bg = numpy.arange(max(10,width*2), 30)
    low_bg = numpy.arange(-30, min(-10,-width*2))
    high_bg_c = numpy.average(high_bg)
    low_bg_c = numpy.average(low_bg)


    sregion = numpy.arange(c_index - 2 * width, c_index + 2 * width + 1)
    high_bg_med = numpy.median(spy[c_index + high_bg])
    low_bg_med = numpy.median(spy[c_index + low_bg])
    inc_bg = (high_bg_med - low_bg_med) / (high_bg_c - low_bg_c)
    int_bg = high_bg_med - inc_bg * (high_bg_c + spx[c_index])
    bg_count = inc_bg * spx[sregion] + int_bg
    grav_c = numpy.sum(spx[sregion] * numpy.absolute(spy[sregion] - bg_count)) / numpy.sum(
        numpy.absolute(spy[sregion] - bg_count))
    stddev = math.fabs((numpy.sum(
        (spx[sregion] - grav_c) ** 2 * numpy.absolute(spy[sregion] - bg_count)) / numpy.sum(
        numpy.absolute(spy[sregion] - bg_count)))) ** 0.5

    if width * lowlim < stddev < width * upplim:
        return grav_c
    else:
        return "N/A"


def auto_ecidentify(ecfile, slit, linelist):
    xpix, y, m, _, delx, _ = open_ec_specfiles(ecfile)
    x = [numpy.array([i+1 for i in range(y[j].size)]) for j in range(len(m))]
    pp = PdfPages("%s.pdf" % ecfile.rstrip("fits").rstrip("."))
    plt.figure(figsize=(20,6))

    low_thres = 50
    high_thres = 2000

    gc = []
    snr = 8.
    width_max = 15.
    upplimfactor = 1.
    if slit == 100:
        linewidth = 2
    elif slit == 140:
        linewidth = 3
    elif slit == 200:
        linewidth = 4
    elif slit == 400:
        linewidth = 8
        upplimfactor = 1.5
    else:
        linewidth = 2

    linewidth = int(linewidth/delx)

    for i in range(len(m)):
        gc_tmp = []
        indexes = find_peaks_cwt(y[i], numpy.arange(1, width_max), min_snr=snr)
        for j in range(len(indexes)):
            if low_thres < xpix[i][indexes[j]] < high_thres:
                gc_sgl = peak_center_single_line(x[i], y[i], linewidth, indexes[j], upplim=upplimfactor)
                if gc_sgl != "N/A":
                    gc_tmp.append(gc_sgl)

        gc.append(numpy.array(gc_tmp))


    gc_arr = numpy.array(gc)

    linef = open(linelist, "r")
    linel = linef.readlines()
    linef.close()

    orders_ref = numpy.array([int(i.split()[0]) for i in linel])
    position_ref = numpy.array([float(i.split()[2]) for i in linel])
    wavelength_ref = numpy.array([float(i.split()[4]) for i in linel])

    position_ref_m = numpy.array([position_ref[orders_ref == i] for i in m])
    wavelength_ref_m = numpy.array([wavelength_ref[orders_ref == i] for i in m])
    lineid_ref_m = numpy.array([[j for j in range(len(position_ref_m[i]))] for i in range(len(m))])

    lineid_gc_m = numpy.array([[j for j in range(len(gc_arr[i]))] for i in range(len(m))])

    minsepid = numpy.array(
        [numpy.array(
            [numpy.argmin(numpy.absolute(gc_arr[i] - position_ref_m[i][j])) for j in range(len(position_ref_m[i]))])
            for i in range(len(m))])
    minsep = numpy.array([gc_arr[i][minsepid[i]] - position_ref_m[i] for i in range(len(m))])

    minsep_1d = []
    for i in range(len(m)):
        minsep_1d.extend(minsep[i])
    minsep_1darr = numpy.array(minsep_1d)

    sep_a, sep_b = numpy.polyfit(position_ref, minsep_1darr, 1)
    sep_std = numpy.std(minsep_1darr - sep_a * position_ref - sep_b)
    nite = 10
    for k in range(nite):
        sep_bool = numpy.absolute(minsep_1darr - sep_a * position_ref - sep_b) / sep_std < 2.

        sep_a, sep_b = numpy.polyfit(position_ref[sep_bool], minsep_1darr[sep_bool], 1)
        sep_std = numpy.std(minsep_1darr[sep_bool] - sep_a * position_ref[sep_bool] - sep_b)

        minsepid = numpy.array(
            [numpy.array(
                [numpy.argmin(numpy.absolute(gc_arr[i] - position_ref_m[i][j] - position_ref_m[i][j] * sep_a - sep_b)) for j in range(len(position_ref_m[i]))])
             for i in range(len(m))])
        minsep = numpy.array([gc_arr[i][minsepid[i]] - position_ref_m[i] for i in range(len(m))])

        minsep_1d = []
        for i in range(len(m)):
            minsep_1d.extend(minsep[i])
        minsep_1darr = numpy.array(minsep_1d)

        plt.scatter(position_ref, minsep_1darr - sep_a * position_ref - sep_b, label="Clipped")
        plt.scatter(position_ref[sep_bool], minsep_1darr[sep_bool] - sep_a * position_ref[sep_bool] - sep_b, lable="Used for calculation")
        plt.grid()
        plt.ylabel("Residual (pix)")
        plt.xlabel("Line center (pix)")
        plt.legend()
        plt.title("%d lines (Iteration %d / %d)" % (numpy.sum(sep_bool), (k+1), nite))
        plt.savefig(pp, format="pdf")
        plt.clf()

    minsep_sig = [[(minsep[i][j] - sep_a * position_ref_m[i][j] - sep_b) / sep_std for j in range(len(minsep[i]))] for i in range(len(m))]

    id_id = [[] for i in range(len(m))]
    id_center = [[] for i in range(len(m))]
    id_wavelength = [[] for i in range(len(m))]
    noid_center_estimate = [[] for i in range(len(m))]
    noid_wavelength = [[] for i in range(len(m))]

    for i in range(len(m)):
        for j in range(len(position_ref_m[i])):
            if math.fabs(minsep_sig[i][j]) < 3.:
                id_id[i].append(lineid_ref_m[i][j])
                id_center[i].append(gc_arr[i][minsepid[i][j]])
                id_wavelength[i].append(wavelength_ref_m[i][j])



    count_idline = 0
    count_refline = 0
    print("Identification of comparison emission lines.")
    for i in range(len(m)):
        plt.step(x[i],y[i],where="mid",color="k")
        ymax = max(y[i])
        for j in range(len(id_center[i])):
            if j == 0:
                plt.plot([id_center[i][j], id_center[i][j]], [ymax*1.15, ymax*1.25], "b", label="Identified", lw=2.)
            else:
                plt.plot([id_center[i][j], id_center[i][j]], [ymax*1.15, ymax*1.25], "b", lw=2.)

        for j in range(len(position_ref_m[i])):
            if j == 0:
                if lineid_ref_m[i][j] in id_id[i]:
                    plt.plot([position_ref_m[i][j], position_ref_m[i][j]], [ymax*1.05, ymax*1.15], "r", label="Reference", lw=2.)
                else:
                    plt.plot([position_ref_m[i][j], position_ref_m[i][j]], [ymax*1.05, ymax*1.15], "r", label="Reference")
            else:
                if lineid_ref_m[i][j] in id_id[i]:
                    plt.plot([position_ref_m[i][j], position_ref_m[i][j]], [ymax*1.05, ymax*1.15], "r", lw=2.)
                else:
                    plt.plot([position_ref_m[i][j], position_ref_m[i][j]], [ymax*1.05, ymax*1.15], "r")

        for j in range(len(lineid_ref_m[i])):
            if not lineid_ref_m[i][j] in id_id[i]:
                noid_center = int(round(position_ref_m[i][j] + sep_a * position_ref_m[i][j] + sep_b))
                if len(xpix[i]) > noid_center and noid_center > 0:
                    noid_center_estimate[i].append(int(round(position_ref_m[i][j] + sep_a * position_ref_m[i][j] + sep_b)))
                    noid_wavelength[i].append(wavelength_ref_m[i][j])

        for j in range(len(noid_center_estimate[i])):
            if low_thres < xpix[i][noid_center_estimate[i][j]] < high_thres:
                gc_sgl = peak_center_single_line(x[i], y[i], linewidth, noid_center_estimate[i][j], upplim=upplimfactor)
                if gc_sgl != "N/A":
                    plt.plot([gc_sgl, gc_sgl], [ymax*1.15, ymax*1.25], "b", lw=1.)
                    id_center[i].append(gc_sgl)
                    id_wavelength[i].append(noid_wavelength[i][j])

        print("m=%d: %d / %d (%.1f percent)" % (m[i], len(id_center[i]), len(position_ref_m[i]), float(len(id_center[i]))/len(position_ref_m[i]) * 100))
        count_idline += len(id_center[i])
        count_refline += len(position_ref_m[i])

        plt.legend()
        plt.title("m=%d" % m[i])
        plt.ylim(0., ymax*1.3)
        plt.ylabel("Flux (arbitrary unit)")
        plt.xlabel("$Y$ (pix)")

        plt.savefig(pp, format="pdf")
        plt.clf()

    print("Total: %d / %d (%.1f percent)" % (count_idline, count_refline, float(count_idline)/count_refline * 100))

    wf = open("database/ec" + ecfile.rstrip("fits").rstrip("."), "w")
    wf.write("# %s\n" % datetime.datetime.today())
    wf.write("begin   ecidentify %s\n" % ecfile.rstrip("fits").rstrip("."))
    wf.write("\tid\t%s\n" % ecfile.rstrip("fits").rstrip("."))
    wf.write("\ttask\tecidentify\n")
    wf.write("\timage\t%s\n" % ecfile.rstrip("fits").rstrip("."))
    wf.write("\tunits\tangstroms\n")
    wf.write("\tfeatures\t%d\n" % count_idline)
    for i in range(len(m)):
        for j in range(len(id_center[i])):
            wf.write("\t\t%d\t%d\t%.2f\t%.5f\t%.5f   4.0 1  1\n" % (
            m[i], m[i], id_center[i][j], id_wavelength[i][j], id_wavelength[i][j]))

    wf.close()
    pp.close()

if __name__ == "__main__":
    filename = sys.argv[1:]

    auto_ecidentify(filename[0], filename[1], filename[2])
