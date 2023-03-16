# -*- coding:utf-8 -*-
import copy
import sys, os
import math
from astropy.io import fits
import numpy as np
from pyraf import iraf
from iraf import onedspec
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.ndimage
from warp.Spec2Dtools import savefitsimage, header_key_read
from warp.aperture import *
from warp.config import constant_str_length

__version__ = "1.4"


# Description:
#
#   Make a hot pixel mask by picking up the outlier pixels in the image after the median filter.
#
#
# Usage:
#
#   $ python badpixmask_v1_1.py inputimage outputmask medianfilter_1 medianfilter_2
#
#
# Output:
#
#   The mask of hot pixels and filered images
#
# Updates:
#
#   ver1.1 updated by Hamano in Nov. 5th 2017
#       - 3x3pixel size was used in ver1.0. In this ver., the 1x5 mask and (i,j)=(1 when i+j=6) mask are used.
#       - From this ver., this script can be used from command line.
#
#   ver1.2 updated by Hamano in 2018.05.24
#       - Python 3
#
#   ver1.3 updated by Hamano in 2019.08.24
#       - version was moved from file name to __version__ parameter.
#       - fixpix was integrated.
#       - badpixmask_flatoff, badpixmask_flaton were added from calibration mode pipeline.
#       - savefitsimage was moved to Spec2Dtools and imported.
##


def badpixmask(inputimage, outputmask, medianfilter_1, medianfilter_2):
    fimg = fits.open(inputimage + ".fits")
    pixdata = fimg[0].data
    fimg.close()

    medsize = 5

    unitarray = np.identity(medsize)

    bmask = np.fabs(pixdata) < 1.e+5
    ave_input = np.average(np.fabs(pixdata[bmask]))

    bp_thres = max(10. * ave_input, 200.)

    lowerlim = -1. * bp_thres
    upperlim = bp_thres

    mfilter_1 = scipy.ndimage.filters.median_filter(pixdata, size=(medsize, 1))
    pd_mfilter1 = pixdata - mfilter_1

    mfilter_2 = scipy.ndimage.filters.median_filter(pd_mfilter1, footprint=unitarray)

    pd_mfilter_2 = pd_mfilter1 - mfilter_2

    savefitsimage(pd_mfilter1, medianfilter_1)
    savefitsimage(pd_mfilter_2, medianfilter_2)

    iraf.imcopy(medianfilter_2, outputmask)
    iraf.imreplace(outputmask, "0", lower=lowerlim, upper=upperlim)
    iraf.imreplace(outputmask, "1", lower="INDEF", upper=lowerlim)
    iraf.imreplace(outputmask, "1", lower=upperlim, upper="INDEF")

    iraf.imarith(outputmask, "*", outputmask, outputmask)

    fimg = fits.open(outputmask + ".fits")
    maskdata = fimg[0].data
    fimg.close()

    return np.sum(maskdata, axis=(0, 1)), bp_thres


def NDRreader(hdr):
    ndrvalue = header_key_read(hdr, "NDR")
    if ndrvalue == "N/A":
        noutputs = header_key_read(hdr, "NOUTPUTS")
        exptime = header_key_read(hdr, "EXPTIME")
        if noutputs == 32:
            if exptime <= 6:
                ndr = 1
            elif exptime <= 12:
                ndr = 2
            elif exptime <= 30:
                ndr = 4
            elif exptime <= 300:
                ndr = 8
            else:
                ndr = 16
        else:
            if exptime <= 40:
                ndr = 1
    else:
        ndr = ndrvalue

    return ndr


def cosmicRayMask(diffimg, rawimg1, rawimg2, outputmask, medianfilter_1, medianfilter_2, apfile, bpmaskflat,
                  abbaflag, noisefits="INDEF", xlim1=-30, xlim2=30, gain=2.27, medsize=5, clipsigma=5., threshold=10.,
                  bins=2, ystep=100, iteration=3, sigstep=2., varatio=2., slitposratio=1.5, maxsigma=20, ndr=16,
                  fixsigma=False):
    if threshold > maxsigma:
        threshold = maxsigma

    difff = fits.open(diffimg + ".fits")
    diffdata = difff[0].data
    difff.close()

    rawf1 = fits.open(rawimg1 + ".fits")
    rawdata1 = rawf1[0].data
    ndr1 = NDRreader(rawf1[0].header)
    rawf1.close()

    rawf2 = fits.open(rawimg2 + ".fits")
    rawdata2 = rawf2[0].data
    ndr2 = NDRreader(rawf2[0].header)
    rawf2.close()

    bpmaskf = fits.open(bpmaskflat)
    bpfdata = bpmaskf[0].data
    bpmaskf.close()

    readn = {1: 19.2, 2: 14., 4: 10., 8: 8., 16: 6., 32: 5.3}
    opos = [-4, 4]
    apos = [14, 22]
    bpos = [-22, -14]

    unitarray = np.identity(medsize)
    diffmf1 = scipy.ndimage.filters.median_filter(diffdata, size=(medsize, 1))
    diffmf1Sub = diffdata - diffmf1
    diffmf2 = scipy.ndimage.filters.median_filter(diffmf1Sub, footprint=unitarray)
    diffmf2Sub = diffmf1Sub - diffmf2
    diffsign = diffmf2Sub > 0
    diffmf2SubAbs = np.absolute(diffmf2Sub)
    savefitsimage(diffmf1Sub, medianfilter_1)
    savefitsimage(diffmf2Sub, medianfilter_2)

    raw1mf1 = scipy.ndimage.filters.median_filter(rawdata1, size=(medsize, 1))
    raw1mf1Sub = rawdata1 - raw1mf1
    raw1mf2 = scipy.ndimage.filters.median_filter(raw1mf1Sub, footprint=unitarray)
    raw1mf2Sub = raw1mf1Sub - raw1mf2
    raw1sign = raw1mf2Sub > 0

    raw2mf1 = scipy.ndimage.filters.median_filter(rawdata2, size=(medsize, 1))
    raw2mf1Sub = rawdata2 - raw2mf1
    raw2mf2 = scipy.ndimage.filters.median_filter(raw2mf1Sub, footprint=unitarray)
    raw2mf2Sub = raw2mf1Sub - raw2mf2
    raw2sign = raw2mf2Sub > 0

    unitarray = unitArrayMake(85., 1, 1, windowSize=10)
    mfimage = scipy.ndimage.filters.median_filter(np.absolute(rawdata1 + rawdata2), footprint=unitarray)
    noiseimg = (mfimage * gain + readn[ndr1] ** 2. + readn[ndr2] ** 2.) ** 0.5 / gain
    if noisefits != "INDEF":
        savefitsimage(noiseimg, noisefits)

    apset = apertureSet(apfile)
    maskap = apset.apmaskArray(lowlim=xlim1, upplim=xlim2)
    slitcoord = apset.slitcoordArray(lowlim=xlim1, upplim=xlim2)

    apy = np.arange(apset.arrayLength) + 1.
    Xarray, Yarray = np.meshgrid(apy - 1, apy - 1)

    x = np.arange(apset.arrayLength)
    y = np.arange(apset.arrayLength)

    missdetectionflag = True
    sigThres = threshold
    itenum = 0

    y_s_list, y_e_list = [], []
    factorlist = []
    for i in range(len(apset.echelleOrders)):
        print("Reducing m={}...".format(apset.echelleOrders[i]))
        y_s = 0
        y_e = ystep
        y_s_list.append([])
        y_e_list.append([])
        factorlist.append([])
        while y_e < apset.arrayLength:
            y_s_list[i].append(y_s)
            y_e_list[i].append(y_e)

            r_s = xlim1 + bins
            r_e = xlim1 + bins * 2

            sclist = []
            stdlist = []
            noiselist = []

            while r_s < xlim2 - bins:
                req_sc = (r_s < slitcoord) & (slitcoord <= r_e) & (maskap == apset.echelleOrders[i]) & (
                        Yarray > y_s) & (Yarray <= y_e)
                mf_sc = diffmf2Sub[req_sc]
                noise_sc = noiseimg[req_sc]
                mfstd_sc1 = np.std(mf_sc)
                for k in range(iteration):
                    mfstd_last = mfstd_sc1
                    mf_req = np.absolute(mf_sc) < clipsigma * mfstd_sc1
                    mfstd_sc1 = np.std(mf_sc[mf_req])
                    if mfstd_last == mfstd_sc1:
                        break

                r_s += bins
                r_e += bins

                stdlist.append(mfstd_sc1)
                noiselist.append(np.average(noise_sc))
                sclist.append((r_s + r_e) / 2.)

            stdmax = max(stdlist)
            noisemax = max(noiselist)
            factor = max(1., stdmax / noisemax)
            factorlist[i].append(factor)

            y_s += ystep
            y_e += ystep
            if apset.arrayLength - y_e < ystep:
                y_e = apset.arrayLength

    while (missdetectionflag) and (sigThres <= maxsigma):
        maskarray = np.zeros(diffmf2Sub.shape, dtype="int16")

        for i in range(len(apset.echelleOrders)):
            for j in range(len(factorlist[i])):
                reqcr = diffmf2SubAbs > sigThres * noiseimg * factorlist[i][j]
                req1 = np.logical_and(diffsign, raw1sign)
                req2 = np.logical_and(np.logical_not(diffsign), raw2sign)
                reqy = (Yarray > y_s_list[i][j]) & (Yarray <= y_e_list[i][j])
                reqm = maskap == apset.echelleOrders[i]
                maskarray[(req1 | req2) & reqy & reqm & reqcr] += 1

        reqmask = (maskarray == 1) & (bpfdata == 0)
        slitcoordMask = slitcoord[reqmask]
        pixnum = np.sum(reqmask)

        hist, bin = np.histogram(slitcoordMask, bins=20, range=(xlim1, xlim2))
        var = np.var(hist)
        ave = np.average(hist)
        varave = var / ave

        if abbaflag:
            reqapos = (slitcoordMask > apos[0]) & (slitcoordMask < apos[1])
            reqbpos = (slitcoordMask > bpos[0]) & (slitcoordMask < bpos[1])
            anum = np.sum(reqapos)
            bnum = np.sum(reqbpos)
            spratio = max(anum, bnum) / pixnum * (xlim2 - xlim1) / max(apos[1] - apos[0], bpos[1] - bpos[0])
            p = "A" if anum > bnum else "B"
            if (spratio > slitposratio) & (varave > varatio) & (not fixsigma):
                missdetectionflag = True
                print(
                    "Iteration {} (sigma={}): bad pix={}, var/ave = {:.2f} > {}, {} position/all = {:.2f} > {}".format(
                        itenum, sigThres, pixnum, varave, varatio, p, spratio, slitposratio))
                itenum += 1
                sigThres += sigstep
            else:
                missdetectionflag = False
                print(
                    "Iteration {} (sigma={}): bad pix={}, var/ave = {:.2f}, {} position/all = {:.2f}".format(
                        itenum, sigThres, pixnum, varave, p, spratio))

        else:
            reqopos = (slitcoordMask > opos[0]) & (slitcoordMask < opos[1])
            onum = np.sum(reqopos)
            spratio = onum / pixnum * (xlim2 - xlim1) / (opos[1] - opos[0])
            if (spratio > slitposratio) & (varave > varatio) & (not fixsigma):
                missdetectionflag = True
                print(
                    "Iteration {} (sigma={}): bad pix={}, var/ave = {:.2f} > {}, O position/all = {:.2f} > {}".format(
                        itenum, sigThres, pixnum, varave, varatio, spratio, slitposratio))
                itenum += 1
                sigThres += sigstep
            else:
                print(
                    "Iteration {} (sigma={}): bad pix={}, var/ave = {:.2f}, O position/all = {:.2f}".format(
                        itenum, sigThres, pixnum, varave, spratio))
                missdetectionflag = False

    savefitsimage(maskarray, outputmask)
    return pixnum, sigThres


def unitArrayMake(angle, samplingSize, width, windowSize=5):
    ksize = samplingSize * windowSize

    b = np.zeros((ksize, ksize))
    x = np.array([[i for i in range(ksize)] for j in range(ksize)])
    x1d = np.array([i for i in range(ksize)])
    y = np.array([[j for i in range(ksize)] for j in range(ksize)])

    angletan = math.tan(angle / 180. * math.pi)

    xcen = (ksize - 1) / 2.
    ycen = (ksize - 1) / 2.

    if angle <= 45.:
        f = y - ycen - angletan * (x - xcen)
    else:
        f = x - xcen - (y - ycen) / angletan

    b[(f > -width / 2.) & (f < width / 2.)] = 1

    return b


def badpixmask_flatoff(inputimage, outputmask):
    fimg = fits.open(inputimage + ".fits")
    pixdata = fimg[0].data
    fimg.close()

    medsize = 15

    bmask = np.fabs(pixdata) < 1.e+2

    bp_thres = np.std(pixdata[bmask])
    clipsig = 5

    lowerlim = -1. * bp_thres * clipsig
    upperlim = bp_thres * clipsig

    mfilter_1 = scipy.ndimage.filters.median_filter(pixdata, size=(medsize, medsize))
    pd_mfilter1 = pixdata - mfilter_1

    medianfilter = inputimage + "_medfilter"
    savefitsimage(pd_mfilter1, medianfilter)

    iraf.imcopy(medianfilter, outputmask)
    iraf.imreplace(outputmask, "0", lower=lowerlim, upper=upperlim)
    iraf.imreplace(outputmask, "1", lower="INDEF", upper=lowerlim)
    iraf.imreplace(outputmask, "1", lower=upperlim, upper="INDEF")

    iraf.imarith(outputmask, "*", outputmask, outputmask)

    fimg = fits.open(outputmask + ".fits")
    maskdata = fimg[0].data
    fimg.close()

    fig = plt.figure(facecolor='white', figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    ax1.imshow(maskdata, cmap=cm.gray, vmin=0, vmax=1, origin='lower')
    plt.title(outputmask)

    plt.savefig("%s.png" % outputmask)

    return np.sum(maskdata, axis=(0, 1)), bp_thres


def badpixmask_flaton(inputimage, outputmask, deadpixmask):
    iraf.apnorm1.background = iraf.apnormalize.background
    iraf.apnorm1.skybox = iraf.apnormalize.skybox
    iraf.apnorm1.weights = iraf.apnormalize.weights
    iraf.apnorm1.pfit = iraf.apnormalize.pfit
    iraf.apnorm1.clean = iraf.apnormalize.clean
    iraf.apnorm1.saturation = iraf.apnormalize.saturation
    iraf.apnorm1.readnoise = iraf.apnormalize.readnoise
    iraf.apnorm1.gain = iraf.apnormalize.gain
    iraf.apnorm1.lsigma = iraf.apnormalize.lsigma
    iraf.apnorm1.usigma = iraf.apnormalize.usigma

    iraf.apnormalize(inputimage + ".fits", inputimage + "_norm.fits", interactive="no", find="no",
                     recenter="no", resize="no", edit="no", trace="no", fittrace="no", normalize="yes", fitspec="yes",
                     function="legendre", background="none", order=15)

    deadpixlim = 0.2
    lowerlim = 0.5
    upperlim = 1.5

    iraf.imcopy(inputimage + "_norm.fits", outputmask)
    iraf.imreplace(outputmask, "1", lower=lowerlim, upper=upperlim)
    iraf.imreplace(outputmask, "-10", lower="INDEF", upper=lowerlim)
    iraf.imreplace(outputmask, "-10", lower=upperlim, upper="INDEF")
    iraf.imreplace(outputmask, "0", lower=lowerlim, upper=upperlim)
    iraf.imreplace(outputmask, "1", lower="INDEF", upper="-5")

    iraf.imcopy(inputimage + "_norm.fits", deadpixmask)
    iraf.imreplace(deadpixmask, "-10", lower="INDEF", upper=deadpixlim)
    iraf.imreplace(deadpixmask, "0", lower=deadpixlim, upper="INDEF")
    iraf.imreplace(deadpixmask, "1", lower="INDEF", upper="-5")

    fimg = fits.open(outputmask + ".fits")
    pixdata = fimg[0].data
    fimg.close()

    fig = plt.figure(facecolor='white', figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    ax1.imshow(pixdata, cmap=cm.gray, vmin=0, vmax=1, origin='lower')
    plt.title(outputmask)

    plt.savefig("%s.png" % outputmask)

    num = np.sum(pixdata, axis=(0, 1))

    return num


def deadpixmap_inter(flaton_inter, widemask, deadpixmap):
    lowerlim = 40

    print(widemask)
    maskfile = fits.open(widemask)
    maskdata = maskfile[0].data
    maskfile.close()

    flatfile = fits.open(flaton_inter + ".fits")
    flatdata = flatfile[0].data
    flatfile.close()

    dpmarray = np.zeros(flatdata.shape, dtype="int16")
    dpmarray[(maskdata == 0) & (flatdata < lowerlim)] = 1
    savefitsimage(dpmarray, deadpixmap)

    # iraf.imcopy(flaton_inter, deadpixmap)
    # iraf.imreplace(deadpixmap, "-10", lower="INDEF", upper=lowerlim)
    # iraf.imreplace(deadpixmap, "0", lower=lowerlim, upper="INDEF")
    # iraf.imreplace(deadpixmap, "1", lower="INDEF", upper="-5")


def pyfixpix(input_file, output_file, mask, linterpolate="INDEF"):
    iraf.imcopy(input_file, output_file)
    iraf.fixpix(output_file, mask, linterp=linterpolate, cinterp="INDEF")

    iraf.hedit(output_file, "BP_MASK", mask, add="yes", verify="no")


if __name__ == "__main__":
    filename = sys.argv[1:]

    a, b = badpixmask(filename[0], filename[1], filename[2], filename[3])
    print(a, b)
