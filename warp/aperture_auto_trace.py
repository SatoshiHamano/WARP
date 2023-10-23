# -*- coding:utf-8 -*-

__version__ = "1.1"

import numpy as np
import astropy.io.fits as pyfits
import scipy.optimize
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, os, copy

from warp.aperture import apertureSet
from warp.Spec2Dtools import header_key_read


def open_tracefiles(tracefile):
    if tracefile.find(".fits") == -1:
        tracefile += ".fits"

    spf = pyfits.open(tracefile)
    spdata = spf[0].data
    hdulist = spf[0].header

    crval1 = header_key_read(hdulist, "CRVAL1")
    cdelt1 = header_key_read(hdulist, "CDELT1")
    crpix1 = header_key_read(hdulist, "CRPIX1")
    crval2 = header_key_read(hdulist, "CRVAL2")
    cdelt2 = header_key_read(hdulist, "CDELT2")
    crpix2 = header_key_read(hdulist, "CRPIX2")
    naxis1 = spf[0].header["NAXIS1"]
    naxis2 = spf[0].header["NAXIS2"]

    if "N/A" in [crval1, cdelt1, crpix1]:
        crval1, cdelt1, crpix1 = 1, 1, 1
    else:
        crval1 = float(crval1)
        cdelt1 = float(cdelt1)
        crpix1 = float(crpix1)

    if "N/A" in [crval2, cdelt2, crpix2]:
        crval2, cdelt2, crpix2 = 1, 1, 1
    else:
        crval2 = float(crval2)
        cdelt2 = float(cdelt2)
        crpix2 = float(crpix2)

    apx = np.array([crval1 + cdelt1 * (l - crpix1 + 1.) for l in range(naxis1)])
    apy = np.array([crval2 + cdelt2 * (l - crpix2 + 1.) for l in range(naxis2)])

    spf.close()

    return apx, apy, spdata, naxis2, cdelt1


def peak_single_line_trace(spx, spy, width, c_index):
    # input: (x,y) of spectrum and roughly estimated peak positions of emission lines.
    # output: peak positions

    sregion = np.arange(max(c_index - 2 * width, 0), min(c_index + 2 * width + 1, len(spx) - 1))
    grav_c = np.sum(spx[sregion] * np.absolute(spy[sregion])) / np.sum(np.absolute(spy[sregion]))

    return grav_c


def residue_trace(pk, ap: apertureSet.aperture, apxarray, apy, p0, jj):
    temp = copy.copy(p0)
    temp[jj] = pk[0]
    x = trace_func(ap.ymin, ap.ymax, apy, ap.yorder, ap.apertureCenterX, temp, ap.functionTypeID)

    res = (apxarray - x)
    return (res)


def trace_func(ymin, ymax, apy, n, center, p, ftype):
    x = np.array([0. for j in range(len(apy))])
    normy = []
    for j in range(len(apy)):
        normy.append((2 * apy[j] - ymin - ymax) / (ymax - ymin))
    if ftype == 1:
        for j in range(len(apy)):
            x[j] += center[0]
            for k in range(n):
                x[j] += chebyshev(k, normy[j]) * p[k]
    elif ftype == 2:
        for j in range(len(apy)):
            x[j] += center[0]
            for k in range(n):
                x[j] += legendre(k, normy[j]) * p[k]

    return x


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


def auto_aptrace(inputdata, apfile_ref, refnpz, shiftmax):
    params = np.load(refnpz)
    orders_ref = params["orders"]
    center_g_ref = params["centers"]
    center_index = params["traceid"]

    apx, apy, spdata, naxis2, cdelt2 = open_tracefiles(inputdata)
    firstid = int(naxis2 / 2 - 1)
    spd = np.average(spdata[firstid - 3:firstid + 2], axis=0)

    apset = apertureSet(apfile_ref, arrayLength=naxis2)

    indexes_peaks = argrelmax(spd, order=5)

    indexes = []
    center_g = []

    for i in list(indexes_peaks[0]):
        if spd[i] > 100 and i < 2000 and i > 40:
            indexes.append(i)
            center_single = peak_single_line_trace(apx, spd, 2, i)
            center_g.append(center_single)

    shift = list(range(-shiftmax, shiftmax + 1))
    chi2 = []

    for j in range(len(shift)):
        minsepid_shift = np.array(
            [np.argmin(np.absolute(center_g - center_g_ref[i] + shift[j])) for i in range(len(center_g_ref))])
        minsep_shift = np.array(
            [center_g[minsepid_shift[i]] - center_g_ref[i] + shift[j] for i in range(len(center_g_ref))])
        chi2.append(np.sum(minsep_shift ** 2))

    shift_min_id = np.argmin(chi2)

    minsepid = np.array(
        [np.argmin(np.absolute(center_g - center_g_ref[i] + shift[shift_min_id])) for i in
         range(len(center_g_ref))])
    minsep = np.array(
        [center_g[minsepid[i]] - center_g_ref[i] + shift[shift_min_id] for i in range(len(center_g_ref))])

    pp = PdfPages("%s.pdf" % os.path.splitext(inputdata)[0])

    peak_center_m = {}
    peaks_m = {}
    for i in range(len(orders_ref)):
        peaks_m[orders_ref[i]] = []
    peak_index_m = {}
    for i in range(len(orders_ref)):
        if center_index[i] == 1:
            peak_center_m[orders_ref[i]] = center_g[minsepid[i]]
            peak_index_m[orders_ref[i]] = indexes[minsepid[i]]
        peaks_m[orders_ref[i]].append(center_g[minsepid[i]])

    searchid = firstid
    apy_search_upp = []
    while apy[searchid] < 2000.:
        apy_search_upp.append(searchid)
        searchid += 10
    searchid = firstid
    apy_search_low = []
    while apy[searchid] > 40.:
        apy_search_low.append(searchid)
        searchid -= 10

    search_width = 30
    search_range = np.arange(-search_width, search_width + 1)
    center_new = []

    plt.figure(figsize=(10, 10))
    for i in range(len(apset.echelleOrders)):
        m = apset.echelleOrders[i]
        ap = apset.apertures[m]
        id_cur = peak_index_m[m]
        center_cur = peak_center_m[m]
        center_new.append([center_cur, firstid + 1])
        apset.apertures[m].apertureCenterX = center_new[i][0]
        apset.apertures[m].apertureCenterY = center_new[i][1]
        shift_prev = 0.
        peak_x = []
        peak_y = []
        for j in range(len(apy_search_low)):
            if id_cur > search_width:
                apx_cut = apx[search_range + id_cur]
                spdata_ave = np.average(spdata[apy_search_low[j] - 3:apy_search_low[j] + 2], axis=0)
                spdata_cut = spdata_ave[search_range + id_cur]
                id_cut = argrelmax(spdata_cut, order=3)[0]
                center_peaks = np.array(
                    [peak_single_line_trace(apx_cut, spdata_cut, 2, id_cut[k]) for k in range(len(id_cut))])
                minsepid_peaks = np.array(
                    [np.argmin(np.absolute(center_peaks - peaks_m[m][k] - shift_prev)) for k in
                     range(len(peaks_m[m]))])
                minsep_peaks = np.array(
                    [center_peaks[minsepid_peaks[k]] - peaks_m[m][k] for k in range(len(peaks_m[m]))])
                shift_prev = np.median(minsep_peaks)
                id_cur = id_cur + int(np.average(id_cut[minsepid_peaks])) - search_width

                peak_x.append(center_peaks[minsepid_peaks[3]])
                peak_y.append(apy_search_low[j])

        id_cur = peak_index_m[m]
        shift_prev = 0.
        for j in range(len(apy_search_upp)):
            if id_cur < apx.size - search_width - 1:
                apx_cut = apx[search_range + id_cur]
                spdata_ave = np.average(spdata[apy_search_upp[j] - 3:apy_search_upp[j] + 2], axis=0)
                spdata_cut = spdata_ave[search_range + id_cur]
                id_cut = argrelmax(spdata_cut, order=3)[0]
                center_peaks = np.array(
                    [peak_single_line_trace(apx_cut, spdata_cut, 2, id_cut[k]) for k in range(len(id_cut))])
                minsepid_peaks = np.array(
                    [np.argmin(np.absolute(center_peaks - peaks_m[m][k] - shift_prev)) for k in
                     range(len(peaks_m[m]))])
                minsep_peaks = np.array(
                    [center_peaks[minsepid_peaks[k]] - peaks_m[m][k] for k in range(len(peaks_m[m]))])
                shift_prev = np.median(minsep_peaks)
                id_cur = id_cur + int(np.average(id_cut[minsepid_peaks])) - search_width

                peak_x.append(center_peaks[minsepid_peaks[3]])
                peak_y.append(apy_search_upp[j])

        peak_x = np.array(peak_x)
        peak_y = np.array(peak_y)
        clipbool = ap.adjustParameter(peak_x, peak_y, clip=True, clipsig=5.)
        ap.calculateTrace()

        plt.scatter(peak_y[np.logical_not(clipbool)], peak_x[np.logical_not(clipbool)], s=5, label="Clipped",
                    marker="x")
        plt.scatter(peak_y[clipbool], peak_x[clipbool], s=5, label="Fitted")
        plt.plot(apset.apertures[m].tracey, apset.apertures[m].tracex, label="Fitting function")
        plt.title("m=%d" % m)
        plt.xlabel("Y (pix)")
        plt.ylabel("X (pix)")
        plt.legend()

        plt.savefig(pp, format="pdf")
        plt.clf()

        print("m=%d Finished. %d points clipped" % (m, len(peak_x) - np.sum(clipbool)))

    pp.close()

    return apset


if __name__ == "__main__":
    filename = sys.argv[1:]
    apset = auto_aptrace(filename[0], filename[1], filename[2], 30)
    bgregion = ["-10:-5,5:10" for i in range(len(apnum))]

    flatapfile = "test_flat"
    flat_low = [-32 for i in range(len(apnum))]
    flat_upp = [32 for i in range(len(apnum))]
    flat_upp[-1] = 500

    apset.write_apdatabase(flatapfile, bgregion)

    # apnum_wide = [42 + i for i in range(20)]
    apnum_wide = [44]
    bgregion_wide = []
    center_wide = []
    ftype_wide = []
    yorder_wide = []
    ymin_wide = []
    ymax_wide = []
    cm_wide = []
    for i in range(len(apnum)):
        if apnum[i] in apnum_wide:
            bgregion_wide.append(bgregion[i])
            center_wide.append(center_new[i])
            ftype_wide.append(ftype[i])
            yorder_wide.append(yorder[i])
            ymin_wide.append(ymin[i])
            ymax_wide.append(ymax[i])
            cm_wide.append(cm_new[i])

    onflatapfile = "test_flaton"
    onflat_low = [-27. for i in range(len(apnum_wide))]
    onflat_upp = [26. for i in range(len(apnum_wide))]
    write_apdatabase(onflatapfile, center_wide, bgregion_wide, apnum_wide, apnum_wide, onflat_low, onflat_upp,
                     ftype_wide, yorder_wide, ymin_wide, ymax_wide, cm_wide)

    apfile = "test_aperture"
    ap_low = [-28. for i in range(len(apnum_wide))]
    ap_high = [28. for i in range(len(apnum_wide))]
    write_apdatabase(apfile, center_wide, bgregion_wide, apnum_wide, apnum_wide, ap_low, ap_high, ftype_wide,
                     yorder_wide, ymin_wide, ymax_wide, cm_wide)
