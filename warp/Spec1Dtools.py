# -*- coding:utf-8 -*-

__version__ = "1.2"

import os, glob
import numpy as np
from pyraf import iraf
from iraf import imred
from iraf import echelle
from iraf import onedspec
import scipy.optimize
import astropy.io.fits as fits
import scipy.constants
from warp.aperture import *
from scipy import interpolate
from warp.Spec2Dtools import savefitsimage

# Description:
#   This script contains some basic and useful functions to manipulate the spectrum files.
#
# Updates:
#   ver1.0 made by Hamano in 2017.01.25
#
#   ver1.1 updated by Hamano in 2018.05.24
#       - pyfits --> astropy.io.fits
#
#   ver1.2 updated by Hamano in 2019.08.24
#       version was moved from file name to __version__ parameter.
#       apall.py was integrated.
#       The name was changed to Spec1Dtools.
#       FSR reading function was integrated.
#       pix2wave was integrated.
#       open_ec_specfiles was integrated from Auto_ecidentify.py.
#


def openspecfits(fitsfile):
    if fitsfile.find("fits") == -1:
        fitsfile += ".fits"

    spfits = fits.open(fitsfile)
    splength = spfits[0].header["NAXIS1"]
    spdata = spfits[0].data

    rcrval1 = float(spfits[0].header["CRVAL1"])
    rcdelt1 = float(spfits[0].header["CDELT1"])
    rcrpix1 = float(spfits[0].header["CRPIX1"])

    lamx = np.array([rcrval1 + rcdelt1 * (l - rcrpix1 + 1.) for l in range(splength)])
    spfits.close()

    return lamx, spdata, rcrval1, rcdelt1, rcrpix1

def open_ec_specfiles(specfile):
    # open echelle format spectra fits file
    # return the list of x-coordinate, y-coordinate and echelle orders

    spf = fits.open(specfile)
    spdata = spf[0].data
    apnum = len(spdata)

    rcrval1 = float(spf[0].header["CRVAL1"])
    rcdelt1 = float(spf[0].header["CDELT1"])
    rcrpix1 = float(spf[0].header["CRPIX1"])
    splength = spf[0].header["NAXIS1"]

    commonx = np.array([[rcrval1 + rcdelt1 * (l - rcrpix1 + 1.) for l in range(splength)] for i in range(apnum)])

    aporders = []
    for i in range(apnum):
        aporders.append(int(spf[0].header["APNUM%d" % (i + 1)].split()[0]))

    spf.close()

    return commonx, spdata, np.array(aporders), rcrval1, rcdelt1, rcrpix1


def binning_spec(lambdax, fluxy, binning_size):
    lambdax_bin = []
    fluxy_bin = []

    for i in range(len(lambdax) // binning_size):
        tmp_x = 0.
        tmp_y = 0.

        for j in range(binning_size):
            tmp_x += lambdax[binning_size * i + j]
            tmp_y += fluxy[binning_size * i + j]

        lambdax_bin.append(tmp_x / binning_size)
        fluxy_bin.append(tmp_y / binning_size)

    return np.array(lambdax_bin), np.array(fluxy_bin)


def PyScombine(inputlist, output):
    flist = ""

    for i in range(len(inputlist)):
        flist += inputlist[i] + ","

    flist = flist.rstrip(",")

    iraf.scombine(flist, output, combine="average", group="all")


def pyapall(inputimage, outputfile, referencename, bgsubs, mode):
    echelle.apall(inputimage, output=outputfile, reference=referencename, interactive="no", find="no", rece="no",
                  resize="no", edit="no", trace="no", fittrac="no", format=mode, background=bgsubs, extras="no")

def resample2Dspec(inputimage, outputfile, outputhdr, ref, interpolation="cubic", finepix=0.01):
    fitsdata = fits.open(inputimage + ".fits")
    dataArray = fitsdata[0].data
    naxis2 = fitsdata[0].header["NAXIS2"]
    apset = apertureSet(ref, arrayLength=naxis2)
    fitsdata.close()
    m = apset.echelleOrders[0]
    lowlim = int(apset.apertures[m].apLow)
    upplim = int(apset.apertures[m].apHigh)
    xnew = list(range(lowlim, upplim+1))
    xsize = len(xnew)
    resampledData = np.zeros((xsize, apset.arrayLength))
    for y in range(apset.arrayLength):
        center = apset.apertures[m].tracex[y]
        centerI = int(center)
        f = interpolate.interp1d(np.arange(max(centerI + lowlim * 2, 1), min(centerI + upplim * 2 + 1, apset.arrayLength)),
                                 dataArray[y, max(centerI + lowlim * 2,1) - 1:min(centerI + upplim * 2, apset.arrayLength)], kind=interpolation)
        xfine = np.arange(max(centerI + lowlim - 3, 1), min(centerI + upplim + 4, apset.arrayLength), finepix)
        datanew = []
        for x in xnew:
            datanew.append(np.average(f(xfine[np.logical_and(xfine > center + x - 0.5, xfine <= center + x + 0.5)])))
        resampledData[:,y] += np.array(datanew)# / np.sum(datanew) * np.sum(dataArray[y,centerI-lowlim:centerI+upplim])

    outputFits = fits.open(outputhdr + ".fits")
    outputFits[0].data = resampledData
    outputFits.writeto(outputfile + ".fits")
    outputFits.close()


def truncate(rawspec, outputfile, p1=1., p2=2048.):
    iraf.scopy(rawspec, outputfile, w1=p1, w2=p2)

    iraf.hedit(outputfile, "CRVAL1", "1.", verify="no")
    iraf.hedit(outputfile, "CRPIX1", "1.", verify="no")
    iraf.hedit(outputfile, "LTV1", "0.", verify="no")


def FSR_angstrom(select_date="latest"):
    datapath = os.path.dirname(os.path.abspath(__file__))
    basename = "FSR/FSR_winered_"
    extention = "txt"
    datalist = glob.glob("%s/%s*%s" % (datapath, basename, extention))
    print(datapath)
    dates = [int(i.split(basename)[-1].rstrip(extention).rstrip(".")) for i in datalist]
    latest_index = dates.index(max(dates))

    selected_data = datalist[latest_index]
    if select_date != "latest":
        try:
            for i in range(len(dates)):
                if int(select_date) == dates[i]:
                    selected_data = datalist[i]
                    break
        except:
            print("Warning: %s cannot be converted to integer." % select_date)
        if selected_data == datalist[latest_index]:
            print("Your input \"%s\" is not found. Existing data list:" % select_date)
            for i in range(len(dates)):
                print("\t%s" % datalist[i])
            print("\nLatest data \"%s\" is shown.\n" % dates[latest_index])
    else:
        selected_data = datalist[latest_index]

    print(selected_data)

    datafile = open(selected_data, "r")
    datalines = datafile.readlines()
    datafile.close()

    for i in range(len(datalines)):
        if datalines[i][0] != "#" or datalines[i] != "":
            if datalines[i].find("WIDE") != -1:
                wideline = i
            elif datalines[i].find("HIRES-Y") != -1:
                hiresyline = i
            elif datalines[i].find("HIRES-J") != -1:
                hiresjline = i

    lines = [wideline, hiresyline, hiresjline]
    lines.append(len(datalines) - 1)
    lines.sort()

    fsr_angstrom = {}
    for i in range(len(lines) - 1):
        for j in range(lines[i], lines[i + 1] + 1):
            if len(datalines[j].split()) >= 3:
                linecomp = datalines[j].split()
                order = int(linecomp[0])
                lowlim = float(linecomp[1])
                upplim = float(linecomp[2])
                fsr_angstrom[order] = [lowlim, upplim]

    return fsr_angstrom


def cut_1dspec(input, output, cutrange, morder, fsr):
    center = (fsr[morder][0] + fsr[morder][1]) / 2.

    cut_lowlim = (center - (center - fsr[morder][0]) * cutrange)
    cut_upplim = (center + (fsr[morder][1] - center) * cutrange)

    onedspec.scopy(input, output, w1=cut_lowlim, w2=cut_upplim)


def dispcor_multi(input_file, output_file, referencename):
    iraf.hedit(input_file, "REFSPEC1", referencename, add="yes", verify="no")

    onedspec.dispcor.dw = "INDEF"
    onedspec.dispcor.flux = "NO"
    onedspec.dispcor.linear = "YES"
    tmpout = "tmpdispcor_dereniw"
    onedspec.dispcor(input_file, output=tmpout)

    tmpdispcor = fits.open("%s.fits" % tmpout)
    tmpdispcor_hdr = tmpdispcor[0].header
    apnum = tmpdispcor_hdr["NAXIS2"]
    tmpdispcor.close()

    dwlist = []
    for i in range(apnum):
        onedspec.scopy("%s[*,%d]" % (tmpout, i + 1), "%s_%d" % (tmpout, i + 1))
        tmpdispcor_sep = fits.open("%s_%d.fits" % (tmpout, i + 1))
        tmpdispcor_sep_hdr = tmpdispcor_sep[0].header
        dwlist.append(float(tmpdispcor_sep_hdr["CDELT1"]))
        tmpdispcor_sep.close()
        os.remove("%s_%d.fits" % (tmpout, i + 1))

    os.remove("%s.fits" % tmpout)

    dwmin = min(dwlist)

    onedspec.dispcor.dw = dwmin
    onedspec.dispcor.flux = "NO"
    onedspec.dispcor.linear = "YES"
    onedspec.dispcor(input_file, output=output_file)


def dispcor_single(input_file, output_file, referencename):
    iraf.hedit(input_file, "REFSPEC1", referencename, add="yes", verify="no")
    onedspec.dispcor.dw = "INDEF"
    onedspec.dispcor.flux = "NO"
    onedspec.dispcor.linear = "YES"
    onedspec.dispcor(input_file, output=output_file)

