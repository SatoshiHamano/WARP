# -*- coding:utf-8 -*-

import time, shutil, os, sys
from pyraf import iraf
import numpy as np
import math
import astropy.io.fits as fits
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import find_peaks_cwt
import matplotlib.cm as cm
from astropy.visualization import ZScaleInterval

from warp.Spec1Dtools import pyapall
from warp.Spec2Dtools import header_key_read
from warp.aperture import apertureSet
from warp.auto_ecidentify import open_ec_specfiles, peak_center_single_line

__version__ = "1.3"

# ver 1.3 updated by S.Hamano (2020.01.12)
#   - bug fix (caused by nearlineid including groupedids.
#

# def apx_estimate(apy, shift, center, ftype, ymin, ymax, yorder, cm):
#     # define the array, in which the aperture function will be stored
#
#     apx = center[0] + shift
#     normy = (2 * apy - ymin - ymax) / (ymax - ymin)
#
#     # calculate the functions
#
#     if ftype == 1:
#         for j in range(yorder):
#             apx += chebyshev(j, normy) * cm[j]
#     elif ftype == 2:
#         for j in range(yorder):
#             apx += legendre(j, normy) * cm[j]
#
#     return apx


def comparison_extract(compfile, apref, shift, apfs):
    apset = apertureSet(apref)
    bgregion = ["0:1" for j in range(len(apset.echelleOrders))]

    for i in range(len(shift)):
        for m in apset.echelleOrders:
            apset.apertures[m].apLow = -0.5 + shift[i]
            apset.apertures[m].apHigh = 0.5 + shift[i]
        apset.write_apdatabase(apfs[i], bgregion)
        pyapall(compfile, apfs[i], apfs[i], "none", "echelle")

def residue_linear(p, xarray, yarray):
    res = (yarray - xarray * p[0] - p[1])
    return (res)

def quadraticFunction(p, xarray):
    return p[0] + p[1] * xarray + p[2] * xarray ** 2

def residue_niji(p, yarray, xarray):
    nijifunc = quadraticFunction(p, xarray)
    res = (yarray - nijifunc)
    return (res)

def auto_angle_measurement(compfname, shift, apfs, apnum, paramnpz):
    if compfname.find(".fits") == -1:
        compfname += ".fits"

    print("\nMeasuring the angles of emission line images in the comparison frame...\n")

    compdata = fits.open(compfname)
    comphdr = compdata[0].header
    slit = header_key_read(comphdr, "SLIT")
    interval = ZScaleInterval()

    if slit == 100:
        linewidth = 2
    elif slit == 140:
        linewidth = 3
    elif slit == 200:
        linewidth = 4
    elif slit == 400:
        linewidth = 8
    else:
        linewidth = 2

    params = np.load(paramnpz)
    p0 = params["p0_yn"]
    p0_m = params["p0_m"]

    apset = apertureSet(apfs[0])
    aplength = len(apset.echelleOrders)

    snr = 4.
    width_max = 15.
    nite = 10
    ns = len(shift)
    low_thres = 50
    high_thres = 2000

    [t0, t1, t2] = p0
    [ap0, ap1, ap2] = p0_m

    spx = []
    spy = []
    ids = []
    c_grav = []
    for i in range(ns):
        x, y, _, _, _, _ = open_ec_specfiles(apfs[i])
        ids.append([])
        c_grav.append([])
        spx.append(x)
        spy.append(y)
        for j in range(aplength):
            peakids = find_peaks_cwt(y[j], np.arange(1, width_max), min_snr=snr)
            ids[i].append([])
            c_grav[i].append([])
            for k in range(len(peakids)):
                if low_thres < x[j][peakids[k]] < high_thres:
                    c_single = peak_center_single_line(x[j], y[j], linewidth, peakids[k])
                    if c_single != "N/A":
                        c_grav[i][j].append(c_single)
                        ids[i][j].append(peakids[k])

    orders_fit = []
    angles_fit = []
    centers_fit = []

    plt.figure(figsize=(10, 10))
    colors = ["g", "y", "r", "k", "c", "m", "b"]
    pp = PdfPages("%s.pdf" % compfname.rstrip("fits").rstrip("."))

    plt.imshow(interval(compdata[0].data), origin="lower", interpolation="none", cmap=cm.gray)

    for j in range(aplength):
        m = apset.echelleOrders[j]
        peaks = np.array([np.array([c_grav[i][j][k] for k in range(len(c_grav[i][j]))]) for i in range(ns)])
        peaks_onearray = []
        for i in range(len(peaks)):
            for k in range(len(peaks[i])):
                peaks_onearray.append(peaks[i][k])
        peaks_onearray = np.array(peaks_onearray)

        peaks_norm = [(2 * peaks[i] - 1. - apset.arrayLength) / (apset.arrayLength - 1.) for i in range(ns)]

        shifts = np.array([np.array([shift[i] for k in range(len(c_grav[i][j]))]) for i in range(ns)])
        shifts_onearray = []
        for i in range(len(shifts)):
            for k in range(len(shifts[i])):
                shifts_onearray.append(shifts[i][k])
        shifts_onearray = np.array(shifts_onearray)

        angles_ref = [np.tan(quadraticFunction(p0, peaks_norm[i]) * quadraticFunction(p0_m, m) / 180. * math.pi) for i
                      in range(ns)]
        angles_onearray = []
        for i in range(len(angles_ref)):
            for k in range(len(angles_ref[i])):
                angles_onearray.append(angles_ref[i][k])
        angles_onearray = np.array(angles_onearray)
        angles_ave = np.average(angles_onearray)

        positionarray = peaks_onearray - shifts_onearray * angles_onearray
        ids = np.array([a for a in range(len(peaks_onearray))])

        groupedids = []
        groups = []
        for i in range(len(peaks_onearray)):
            if not i in groupedids:
                nearlinesid = ids[np.absolute(positionarray - positionarray[i]) < np.absolute(
                    shifts_onearray - shifts_onearray[i]) * 0.176 + 2.]  # 0.176 = tan(10deg)
                groupedids_set = set(groupedids)
                nearlinesid_set = set(nearlinesid)
                nearlinesid = list(nearlinesid_set - groupedids_set)
                nearlinesid.sort()
                nearlinesid = np.array(nearlinesid)
                if len(nearlinesid) >= 3:
                    positionNL = positionarray[nearlinesid]
                    shiftsNL = shifts_onearray[nearlinesid]
                    peaksNL = peaks_onearray[nearlinesid]
                    id_group = []

                    for k in range(len(shift)):
                        if shift[k] in shiftsNL:
                            candids = nearlinesid[shiftsNL == shift[k]]
                            candpositions = positionNL[shiftsNL == shift[k]]
                            minid = candids[np.argmin(np.absolute(candpositions - positionarray[i]))]
                            if not minid in groupedids:
                                id_group.append(minid)
                                groupedids.append(minid)
                    id_group = np.array(id_group)
                    groups.append(id_group)

        # for i in range(ns):
        #     plt.scatter(peaks[i], shifts[i], color="0.5")

        peaks_group = []
        shifts_group = []
        apxs_group = []
        for i in range(len(groups)):
            peaks_group.append(peaks_onearray[groups[i]])
            shifts_group.append(shifts_onearray[groups[i]])
            apxs_group.append(apset.apertures[m].tracexfunc(peaks_group[i]) + shifts_group[i])

            # plt.scatter(peaks_group[i], shifts_group[i], color=colors[i % len(colors)])
            plt.scatter(apxs_group[i], peaks_group[i], s=1., color=colors[i % len(colors)])

        # plt.grid()
        # plt.title("m=%d" % m)
        # plt.xlabel("Y on array (pix)")
        # plt.ylabel("Shift from aperture center (pix)")
        # plt.savefig(pp, format="pdf")
        # plt.clf()

        for k in range(len(peaks_group)):
            if len(shifts_group[k]) > 2:
                plinear0 = [angles_ave, 1000]
                param_output = scipy.optimize.leastsq(residue_linear, plinear0,
                                                      args=(np.array(apxs_group[k]), np.array(peaks_group[k])),
                                                      full_output=True)
                orders_fit.append(m)
                angles_fit.append(math.atan(param_output[0][0]) * 180. / math.pi)
                centers_fit.append((2 * np.average(
                    np.array(peaks_group[k]) - param_output[0][0] * np.array(shifts_group[k])) - 1. - apset.arrayLength) / (apset.arrayLength - 1.))

        print("m=%d Finished" % m)

    plt.xlabel("$X$ (pix)")
    plt.ylabel("$Y$ (pix)")
    plt.savefig(pp, format="pdf")
    plt.clf()

    orders_fit = np.array(orders_fit)
    angles_fit = np.array(angles_fit)
    centers_fit = np.array(centers_fit)

    clip_flag = np.array([True for i in range(len(angles_fit))])

    sigma_thr = 3.

    plt.figure()
    for n in range(nite):
        p0 = [t0, t1, t2]

        plt.scatter(centers_fit[clip_flag], angles_fit[clip_flag] / quadraticFunction(p0_m, orders_fit[clip_flag]),
                    c="r", s=5, label="Clipped")

        param_output_y = scipy.optimize.leastsq(residue_niji, p0, args=(angles_fit[clip_flag], centers_fit[clip_flag]),
                                                full_output=True)
        res = angles_fit - quadraticFunction(param_output_y[0], centers_fit)
        res_std = angles_fit[clip_flag] - quadraticFunction(param_output_y[0], centers_fit[clip_flag])
        std = np.std(res_std)
        res_abs = np.absolute(res)

        clip_flag[res_abs > std * sigma_thr] = False

        plt.scatter(centers_fit[clip_flag], angles_fit[clip_flag] / quadraticFunction(p0_m, orders_fit[clip_flag]),
                    c="b", label="Comp. lines", s=5)
        plt.plot(np.arange(-1, 1, 0.01), quadraticFunction(param_output_y[0], np.arange(-1, 1, 0.01)), "k--",
                 label="Fitting function")
        plt.title("Iteration %d / %d" % (n + 1, nite))
        plt.ylabel("Angle (degree)")
        plt.xlabel("Normalized Y")
        plt.legend()

        plt.savefig(pp, format="pdf")
        plt.clf()

        [t0, t1, t2] = param_output_y[0]

        p0_m = [ap0, ap1, ap2]

        plt.scatter(orders_fit[clip_flag],
                    angles_fit[clip_flag] / quadraticFunction(param_output_y[0], centers_fit[clip_flag]), c="r", s=5,
                    label="Clipped")

        param_output_m = scipy.optimize.leastsq(residue_niji, p0_m, args=(
            angles_fit[clip_flag] / quadraticFunction(param_output_y[0], centers_fit[clip_flag]), orders_fit[clip_flag]),
                                                full_output=True)

        res = angles_fit / quadraticFunction(param_output_y[0], centers_fit) - quadraticFunction(param_output_m[0], orders_fit)
        res_std = angles_fit[clip_flag] / quadraticFunction(param_output_y[0], centers_fit[clip_flag]) - quadraticFunction(
            param_output_m[0], orders_fit[clip_flag])
        std = np.std(res_std)
        res_abs = np.absolute(res)

        clip_flag[res_abs > std * sigma_thr] = False

        plt.scatter(orders_fit[clip_flag],
                    angles_fit[clip_flag] / quadraticFunction(param_output_y[0], centers_fit[clip_flag]), c="b",
                    label="Comp. lines", s=5)
        plt.plot(np.array(apnum), quadraticFunction(param_output_m[0], np.array(apnum)), "k--",
                 label="Fitting function")

        [ap0, ap1, ap2] = param_output_m[0]
        plt.title("Iteration %d / %d" % (n + 1, nite))
        plt.ylabel("Normalized angle")
        plt.xlabel("Echelle order m")
        plt.legend()

        plt.savefig(pp, format="pdf")
        plt.clf()

    pp.close()

    return [t0, t1, t2], [ap0, ap1, ap2]


if __name__ == "__main__":
    filename = sys.argv[1:]

    ns = 15

    shift = [3 * int((i + 1) / 2) * (-1) ** i for i in range(ns)]
    apfs = ["shift%d.fits" % shift[i] if shift[i] <= 0 else "shift+%d.fits" % shift[i] for i in range(ns)]

    # apnum, center, aplow, aphigh, ftype, yorder, ymin, ymax, cm = read_apdatabase(
    #     filename[1].rstrip("fits").rstrip("."))
    #
    comparison_extract(filename[0], filename[1], shift, apfs)

    [t0, t1, t2], [ap0, ap1, ap2] = auto_angle_measurement(filename[0], shift, apfs, apnum, filename[2])
    print([t0, t1, t2], [ap0, ap1, ap2])
