# -*- coding:utf-8 -*-

__version__ = "2.3"

import sys
from pyraf import iraf
from astropy.io import fits
import numpy as np
import math
import scipy.optimize
from warp.aperture import *
from warp.Spec2Dtools import *

iraf.images()
iraf.imutil()
iraf.imgeom()


#
# Description:
#
#   This script search for the center position and the width of the spatial profile in the aperture from input 2d-image.
#
# Usage:
#
#   python centersearch_fortrans_v1_0.py input reference badpixfile
#
#   input is the fits image after the transformation.
#   reference is the file name, with which the aperture center and range are defined in IRAF format.
#   badpixfile is the output image file name. the map of the detected bad pix is saved.
#
#
# Updates:
#   ver1.0: made by hamano (2016-12-06)
#
#   ver2.0: updated by Hamano (2017-04-14)
#       - By ignoring the pixels having negative values during fitting, the systematic error of measured FWHM is reduced.
#       - centersearch_fortrans function start to produce spatial profiles.
#       - bad pix detection function is separated and stopped.
#
#   ver2.1: updated by Hamano (2017-10-17)
#       - 500 pix from array edge is not used.
#       - The lines with low peak count and/or strange peak positions are not used.
#       - If the peak counts around center lines are lower than 100 counts, all lines are used for measurement.
#       - Stacked spatial profile is made first, and then the Gaussian function is fitted to the profile.
#
#   ver2.2: updated by Hamano (2018-04-26)
#       - In case that the nodding pattern is "ABBA", "O" position is rejected from the center-search area.
#       - The clipping algorithm was much improved.
#       - Speed up! Calculation time was much shortened.
#
#   ver2.3: udpated by Hamano (2019-09-22)
#       - the version number is removed from the name of the script.
#       - Not used function "bpdetection" was removed.
#       - Refactoring.
#
#

def gaussianfunc_cs(p, x, c):
    y = math.sqrt(p[0] ** 2) * np.exp(-(x - p[1] - c) ** 2 / (2 * p[2] ** 2)) + p[3]
    return y


def residue_cs(p, y, x, c):
    res = (y - gaussianfunc_cs(p, x, c))
    return (res)


def gaussianfunc_onep(p, x, c, s):
    y = math.sqrt(p[0] ** 2) * np.exp(-((x - c) ** 2) / (2 * s ** 2)) + p[1]
    return y


def residue_onep(p, y, x, c, s):
    res = (y - gaussianfunc_onep(p, x, c, s))
    return (res)


def clipping(dataarray, fitfunc, clipsig):
    difdata = (dataarray - fitfunc) / fitfunc ** 0.5
    stddata = np.std(difdata)
    clipped_pos = []
    for i in range(len(dataarray)):
        if math.fabs(difdata[i]) > clipsig * stddata:
            clipped_pos.append(i)
    return clipped_pos, stddata


def make_slit_profile(transedimage, apdatabase, datfile):
    # set parameters
    lowlim = 500  # the lower limit of the region to be used for the center search
    upplim = 2048 - lowlim  # the upper limit of the region to be used for the center search
    step_sampling = 5  # the step of the center search in pix
    distthres = 8.
    center_width = 50

    # read transformed image
    if transedimage.find(".fits") == -1:
        transedimage += ".fits"

    trimg_fits = fits.open(transedimage)
    trimg_fits_hdr = trimg_fits[0].header
    trimg_fits_data = trimg_fits[0].data
    trimg_fits.close()

    naxis1 = trimg_fits_hdr["NAXIS1"]
    naxis2 = trimg_fits_hdr["NAXIS2"]
    crpix1 = trimg_fits_hdr["CRPIX1"]
    crpix2 = trimg_fits_hdr["CRPIX2"]
    cdelt1 = trimg_fits_hdr["CDELT1"]
    cdelt2 = trimg_fits_hdr["CDELT2"]
    crval1 = trimg_fits_hdr["CRVAL1"]
    crval2 = trimg_fits_hdr["CRVAL2"]

    # trimg_xaxis = crval1 + cdelt1 * (np.arange(naxis1) - crpix1 + 1.)
    # trimg_yaxis = crval2 + cdelt2 * (np.arange(naxis2) - crpix2 + 1.)
    lowlim_ypix = int((lowlim - crval2) / cdelt2 + crpix2 - 1) - 1
    upplim_ypix = int((upplim - crval2) / cdelt2 + crpix2 - 1)
    center_ypix = int(((lowlim + upplim) / 2. - crval2) / cdelt2 + crpix2 - 1)

    # read parameters defined with aptrace

    apset = apertureSet(apdatabase, arrayLength=naxis2)
    # apnum, center, aplow, aphigh, ftype, yorder, ymin, ymax, cm = read_apdatabase(apdatabase)

    # make aperture function

    # apxtmp0, apy = make_apfunction(apnum, center, naxis2, ftype, yorder, ymin, ymax,cm)
    ap0 = apset.apertures[apset.echelleOrders[0]]
    apx = ap0.tracex # apxtmp0[0]

    apxlow = apx + ap0.apLow - 1
    apxlow_pix = apxlow.astype(np.int32)
    apxlow_pix[apxlow_pix < 0] = 0

    apxupp = apx + ap0.apHigh - 1
    apxupp_pix = apxupp.astype(np.int32)
    apxupp_pix[apxupp_pix > naxis1] = naxis1 - 1
    apcenter = apx + (ap0.apLow + ap0.apHigh) / 2.

    ysample = np.arange(lowlim_ypix, upplim_ypix, step_sampling)

    pf_xdata = np.array([])
    pf_ydata = np.array([])

    for i in range(len(ysample)):
        nocutdata_y = trimg_fits_data[ysample[i]][apxlow_pix[ysample[i]]:apxupp_pix[ysample[i]]]
        nocutdata_x = np.arange(apxlow_pix[ysample[i]] + 1,
                                   apxupp_pix[ysample[i]] + 1)  # +1 is to match the coordinate defined in IRAF
        pf_xdata = np.r_[pf_xdata, nocutdata_x - apx[ysample[i]]]
        pf_ydata = np.r_[pf_ydata, nocutdata_y / np.max(nocutdata_y)]

    min_xdata = np.amin(pf_xdata)
    max_xdata = np.amax(pf_xdata)
    step = 100
    pf_width = (max_xdata - min_xdata) / float(step)
    pf_bool = [(min_xdata + pf_width * i < pf_xdata) & (pf_xdata < min_xdata + pf_width * (i + 1)) for i in range(step)]
    med_x = np.array([np.median(pf_xdata[pf_bool[i]]) for i in range(step)])
    med_y = np.array([np.median(pf_ydata[pf_bool[i]]) for i in range(step)])

    wf = open(datfile, "w")
    for i in range(step):
        wf.write("%.4f\t%.4f\n" % (med_x[i], med_y[i]))
    wf.close()




def centersearch_fortrans(transedimage, apdatabase, datfile, abbaflag=False):
    # set parameters

    lowlim = 500  # the lower limit of the region to be used for the center search
    upplim = 2048 - lowlim  # the upper limit of the region to be used for the center search
    step_sampling = 5  # the step of the center search in pix
    iterate = 3  # iterate number of sigma clipping
    distthres = 8.
    center_width = 50
    oposition = [-8, 8]

    # read transformed image

    if transedimage.find(".fits") == -1:
        transedimage += ".fits"

    trimg_fits = fits.open(transedimage)
    trimg_fits_hdr = trimg_fits[0].header
    trimg_fits_data = trimg_fits[0].data
    trimg_fits.close()

    naxis1 = trimg_fits_hdr["NAXIS1"]
    naxis2 = trimg_fits_hdr["NAXIS2"]
    crpix1 = trimg_fits_hdr["CRPIX1"]
    crpix2 = trimg_fits_hdr["CRPIX2"]
    cdelt1 = trimg_fits_hdr["CDELT1"]
    cdelt2 = trimg_fits_hdr["CDELT2"]
    crval1 = trimg_fits_hdr["CRVAL1"]
    crval2 = trimg_fits_hdr["CRVAL2"]

    lowlim_ypix = int((lowlim - crval2) / cdelt2 + crpix2 - 1) - 1
    upplim_ypix = int((upplim - crval2) / cdelt2 + crpix2 - 1)
    center_ypix = int(((lowlim + upplim) / 2. - crval2) / cdelt2 + crpix2 - 1)

    # read parameters defined with aptrace

    apset = apertureSet(apdatabase, arrayLength=naxis2)
    ap0 = apset.apertures[apset.echelleOrders[0]]
    apx = ap0.tracex # apxtmp0[0]

    apxlow = apx + ap0.apLow - 1
    apxlow_pix = apxlow.astype(np.int32)
    apxlow_pix[apxlow_pix < 0] = 0

    apxupp = apx + ap0.apHigh - 1
    apxupp_pix = apxupp.astype(np.int32)
    apxupp_pix[apxupp_pix > naxis1] = naxis1 - 1
    apcenter = apx + (ap0.apLow + ap0.apHigh) / 2.

    cutdata_x = [np.arange(apxlow_pix[i] + 1, apxupp_pix[i] + 1) for i in range(naxis2)]
    cutdata_y = [trimg_fits_data[i][apxlow_pix[i]:apxupp_pix[i]] for i in range(naxis2)]
    for i in range(naxis2):
        for j in range(len(cutdata_x[i]) - len(cutdata_y[i])):
            cutdata_y[i] = np.append(cutdata_y[i], 0.)

    if abbaflag:
        for i in range(naxis2):
            oposbool = (cutdata_x[i] < oposition[0] + apcenter[i] + 1) | (cutdata_x[i] > oposition[1] + apcenter[i] + 1)
            cutdata_x[i] = cutdata_x[i][oposbool]
            cutdata_y[i] = cutdata_y[i][oposbool]
    cutdata_x = np.array(cutdata_x)
    cutdata_y = np.array(cutdata_y)

    # construct low-noise profile

    # set step_sampling = 1 if even the maximum count is lower than 100.

    centermax = [np.max(cutdata_y[center_ypix - center_width + i]) for i in
                 range(0, 2 * center_width, step_sampling)]
    if np.average(centermax) < 100.:
        step_sampling = 1

    # determine maximum position in each row

    ysample = np.arange(lowlim_ypix, upplim_ypix, step_sampling)
    maxid = np.array([np.argmax(cutdata_y[i]) for i in ysample])
    maxpix = np.array([cutdata_x[ysample[i]][maxid[i]] for i in range(len(ysample))])
    maxflux = np.array([cutdata_y[ysample[i]][maxid[i]] for i in range(len(ysample))])
    maxdist = np.array([maxpix[i] - apx[ysample[i]] for i in range(len(ysample))])

    # rejecting low S/N data with sigma clipping of the distance and flux.

    dist_med = np.median(maxdist)
    dist_scat = np.std(maxdist)
    dist_clip = np.fabs(maxdist - dist_med) < distthres

    flux_clip = (maxflux > np.median(maxflux) * 0.3) & (maxflux < np.median(maxflux) * 2.0)

    for i in range(iterate):
        flux_med = np.median(maxflux[flux_clip])
        flux_scat = np.std(maxflux[flux_clip])
        flux_clip = flux_clip & (maxflux - flux_med > -1.5 * flux_scat) & (maxflux - flux_med < 4. * flux_scat)

    # for counting the data points not clipped.
    pf_xdata = np.array([])
    pf_ydata = np.array([])
    counter = 0
    if np.sum(np.logical_and(dist_clip, flux_clip)) > 0:
        clipArray = np.logical_and(dist_clip, flux_clip)
    elif np.sum(flux_clip) > 0:
        clipArray = flux_clip
    else:
        return math.nan, math.nan

    for i in range(len(ysample)):
        if clipArray[i]:
            nocutdata_y = trimg_fits_data[ysample[i]][apxlow_pix[ysample[i]]:apxupp_pix[ysample[i]]]
            nocutdata_x = np.arange(apxlow_pix[ysample[i]] + 1,
                                       apxupp_pix[ysample[i]] + 1)  # +1 is to match the coordinate defined in IRAF
            pf_xdata = np.r_[pf_xdata, nocutdata_x - apx[ysample[i]]]
            pf_ydata = np.r_[pf_ydata, nocutdata_y / maxflux[i]]
            counter += 1

    if counter == 0:
        print("Skipped %s because of its low flux." % (counter, len(ysample), transedimage))
        return math.nan, math.nan
    else:
        print("%d/%d rows in %s are used." % (counter, len(ysample), transedimage))

    min_xdata = np.amin(pf_xdata)
    max_xdata = np.amax(pf_xdata)
    step = 100
    pf_width = (max_xdata - min_xdata) / float(step)
    pf_bool = [(min_xdata + pf_width * i < pf_xdata) & (pf_xdata < min_xdata + pf_width * (i + 1)) for i in range(step)]
    med_x = np.array([np.median(pf_xdata[pf_bool[i]]) for i in range(step)])
    med_y = np.array([np.median(pf_ydata[pf_bool[i]]) for i in range(step)])

    wf = open(datfile, "w")
    for i in range(step):
        wf.write("%.4f\t%.4f\n" % (med_x[i], med_y[i]))
    wf.close()

    # input parameters for the gaussian fitting

    trimrange = int(8 / pf_width)
    if abbaflag:
        medx_center = np.average(med_x)
        oposbool2 = (med_x < oposition[0] + medx_center) | (med_x > oposition[1] + medx_center)
        maxid = np.argmax(med_y[oposbool2])
        maxpix = med_x[oposbool2][maxid]
        maxid = np.where(med_x == maxpix)[0][0]
        # print(oposbool2, medx_center, maxid, maxpix)
    else:
        maxid = np.argmax(med_y)
        maxpix = med_x[maxid]

    if maxid - trimrange < 0:
        lowlim_trim = 0
    else:
        lowlim_trim = maxid - trimrange
    if maxid + trimrange > len(med_x):
        upplim_trim = len(med_x)
    else:
        upplim_trim = maxid + trimrange

    med_xtrim = med_x[lowlim_trim:upplim_trim]
    med_ytrim = med_y[lowlim_trim:upplim_trim]

    peak = 1.
    centershift = 0.1
    gwidth = 5.
    offset = 0.01

    # fit first guess
    p0 = [peak, centershift, gwidth, offset]
    param_output = scipy.optimize.leastsq(residue_cs, p0, args=(med_ytrim, med_xtrim, maxpix), full_output=True)

    peak = param_output[0][0]
    centershift = param_output[0][1]
    gwidth = math.fabs(param_output[0][2])
    offset = param_output[0][3]

    xshift_fix = centershift + maxpix  # the difference of the peak position and the center position of the aperture
    fwhm_fix = gwidth * 2.3548

    return xshift_fix, fwhm_fix



if __name__ == "__main__":
    inputfile = str(sys.argv[1])
    referencename = str(sys.argv[2])
    profdat = str(sys.argv[3])

    xs, gw = centersearch_fortrans(inputfile, referencename, profdat)

    print("The shift from reference:\t%.3f pix" % xs)
    print("The width of Gaussian fit:\t%.3f pix" % gw)
