# -*- coding:utf-8 -*-

__version__ = "2.3"

import sys, os, datetime, time
import numpy as np
import math
from astropy.io import fits
from pyraf import iraf
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Spec1Dtools import openspecfits

iraf.noao()
iraf.onedspec()


# Description:
#
#   The shift amount between the first input file and the other input files are searched and shift the other spectrum files to align with the first spectrum file.
#
#
# Usage:
#
#   $ python ccwaveshift_v1_3.py inputlist
#
#   inputlist is the list of the lists, in which the spectrum files with the same echelle orders are listed.
#
#
# Output:
#
#   The shifted spectrum files. The name of the output file is input file name + "s".
#
# Updates:
#
#   ver1.0-1.3 made and updated by Hamano.
#   ver1.0-1.2 are for debug. They cannot be used. ver1.3 is the first version.
#
#   ver2.0 updated by Hamano in Nov. 5th 2017
#       - The resampling method is changed. IRAF scombine is used for the resampling.
#       - PySpecshift is updated. Resampling is done after shifting.
#       - The functions which are not used in this and other scripts are left for the memorial of numerous try-and-error.
#
#   ver2.1 updated by Hamano in Apr. 19th 2018
#       - Python 3
#       - Import updated ver. of spectools (ver1.1).
#       - Not used functions are removed. (They can be seen in old version.)
#       - Interpolation method is replaced from IRAF/scombine with scipy.ndimage.interpolation.shift.
#           With this change, the computational time is shortened.
#
#   ver2.2 updated by Hamano in Sep. 22nd 2019
#       - the version number is removed from the name of the script.
#       - Evaluation function for estimating the shift is changed from absolute to squared.
#
#   ver2.3 updated by Hamano in Feb. 2nd 2020
#       - reference id can be changed from 0 (1st frame) in waveshift_oneorder, waveshift_multiorder.
#       - waveshiftClip is slightly changed.
##


def cc_spec_subp_sampling(sp1, sp2, cshift, width, step):
    # This function returns the shift with which the sp2 matches best with sp1.
    # The searching region of shift is from cshift-width to cshift+width.
    # The searching step is set by "step"
    # cshift, width, step are in pix.

    # open the spectra files

    imfits1 = "tmpfits_%s-%d%d%d%d.fits" % (
        datetime.date.today(), datetime.datetime.now().hour, datetime.datetime.now().minute,
        datetime.datetime.now().second,
        datetime.datetime.now().microsecond)
    PySpecshift(sp1, imfits1, 0.)
    sp1x, sp1y, crval1, cdelt1, crpix1 = openspecfits(imfits1)
    os.remove(imfits1)

    sp2x, sp2y, crval2, cdelt2, crpix2 = openspecfits(sp2)

    # normalize them

    normalization_sp2 = np.median(sp2y)
    sp1y_norm = sp1y / np.median(sp1y)
    sp2y_norm = sp2y / np.median(sp2y)

    len_sp2 = len(sp2x)

    # set the parameters

    div = 100  # (1.-2./div)*100 percent of the spectrum will be used
    ec_short = 800  # the region size (in pix) to be removed in the calculation of difference of the spectrum.
    ec_long = 300  # the region size (in pix) to be removed in the calculation of difference of the spectrum.

    # calculate the region of sp1 to be used

    naxis1 = len(sp1x)
    start_sp1 = int(naxis1 / div)
    end_sp1 = int(naxis1 / div * (div - 1))

    # calculate the shift vector for sp2

    intwidth = int(width)
    shiftnum = int(width / step * 2 + 1)
    shift = [cshift - width + float(i) * step for i in range(shiftnum)]  # calculate the shifts from the input param.

    # shift sp2, resample it and calculate their difference from sp1.

    dify = []
    for i in range(shiftnum):
        sp2y_shifted = scipy.ndimage.interpolation.shift(sp2y, shift[i], order=3)
        sp2y_norm_shifted = sp2y_shifted / normalization_sp2
        tmpydif = (sp1y_norm[start_sp1 + ec_short:end_sp1 - ec_long] - sp2y_norm_shifted[
                                                                       start_sp1 + ec_short:end_sp1 - ec_long]) ** 2.
        dify.append(np.average(tmpydif))
    # the shift with which the difference becomes lowest will be returned.

    return np.array(shift), np.array(dify)


def PySpecshift(input, output, shift):
    imfits = "tmpfits_%s-%d%d%d%d.fits" % (
        datetime.date.today(), datetime.datetime.now().hour, datetime.datetime.now().minute,
        datetime.datetime.now().second,
        datetime.datetime.now().microsecond)
    if input.find(".fits") == -1:
        hdulist = fits.open(input + ".fits")
    else:
        hdulist = fits.open(input)
    prihdr = hdulist[0].header
    naxis1 = prihdr["NAXIS1"]
    cdelt1 = prihdr["CDELT1"]
    iraf.scopy(input, imfits)
    iraf.specshift(imfits, shift)
    iraf.scombine(imfits, output, w1=1., dw=cdelt1, nw=naxis1, logfile="null")
    os.remove(imfits)
    hdulist.close()


def waveshift_oneorder(files, refid):
    # the wavelength shift vector
    # the shift of 1st file from 1st file is set as 0.0.
    wshift_vec = [0.0 for i in range(len(files))]
    sp1x, sp1y, crval1, cdelt1, crpix1 = openspecfits(files[refid])
    for i in range(len(files)):
        if i != refid:
            print("Calculating spectra shift between %s and %s" % (files[refid], files[i]))

            # the first guess. center is 0. step is large (0.5pix). width is large (2pix)
            shift_1, dify_1 = cc_spec_subp_sampling(files[refid], files[i], 0.001, 2.0, 0.5)
            cshift_p = shift_1[np.argmin(dify_1)]

            # the second guess. center is from first guess. the width is set as the step in first guess.
            shift_2, dify_2 = cc_spec_subp_sampling(files[refid], files[i], cshift_p, 0.4, 0.1)
            cshift_sp = shift_2[np.argmin(dify_2)]

            # the last guess. center is from 2nd guess. the width is set as the step in second guess.
            shift_3, dify_3 = cc_spec_subp_sampling(files[refid], files[i], cshift_sp, 0.07, 0.01)
            cshift_ssp = shift_3[np.argmin(dify_3)]

            wshift_vec[i] = cshift_ssp * cdelt1

    return wshift_vec


def waveshift_multiorder(files, refid):
    # calculate the wavelength shift for each orders and make matrix of wavelength shift
    # ane axis is input files. the other axis is echelle orders.
    wshift_matrix = []
    for i in range(len(files)):
        wshift_matrix.append(waveshift_oneorder(files[i], refid))

    # calculate the median of wavelength shift.
    shift_median = []
    for j in range(len(files[0])):
        tmplist = []
        for i in range(len(files)):
            tmplist.append(wshift_matrix[i][j])
        shift_median.append(np.median(tmplist))

    return shift_median

def waveshiftClip(wshift_matrix, objnum, aplength, sigma_1st=1., sigma=2., iterate=5, stdthres=0.1):
    shift_average = []
    shift_calcnum = []
    shift_stddev = []
    for i in range(objnum):
        shiftvec = []
        shiftclip = np.array([0 for jj in range(aplength)])
        for j in range(aplength):
            shiftvec.append(wshift_matrix[j][i])
        shiftvecarray = np.array(shiftvec)

        for k in range(iterate):
            shift_average_ite = np.average(shiftvecarray[shiftclip == 0])
            shift_std_ite = np.std(shiftvecarray[shiftclip == 0])
            if k == 0:
                shiftclip[np.absolute(shiftvecarray - shift_average_ite) / shift_std_ite > sigma_1st] += 1
            elif k != iterate - 1:
                shiftclip[np.absolute(shiftvecarray - shift_average_ite) / shift_std_ite > sigma] += 1

        shift_average.append(shift_average_ite)
        shift_stddev.append(shift_std_ite)
        shift_calcnum.append(np.sum(shiftclip == 0))

    shift_stddev_arr = np.array(shift_stddev)
    shift_calcnum_arr = np.array(shift_calcnum)
    wshift_flag = np.zeros(objnum)
    wshift_flag[shift_stddev_arr > stdthres] += 1
    wshift_flag[shift_calcnum_arr < aplength/2] += 1

    return shift_average, shift_calcnum, shift_stddev, wshift_flag

if __name__ == "__main__":

    filename = sys.argv[1:]

    rfile = open(filename[0], "r")
    rlines = rfile.readlines()
    listfiles = [rlines[i].split()[0] for i in range(len(rlines))]
    rfile.close()

    fileend = "s"

    files = []

    for i in range(len(listfiles)):
        rflist = open(listfiles[i], "r")
        rflines = rflist.readlines()
        rflist.close()
        files.append([rflines[j].split()[0] for j in range(len(rflines))])

    shift_median = waveshift_multiorder(files)

    for i in range(len(rlines)):
        for j in range(len(files[i])):
            PySpecshift(files[i][j], files[i][j].rstrip("fits").rstrip(".") + fileend, shift_median[j])
