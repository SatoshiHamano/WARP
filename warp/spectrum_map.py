# -*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import os
from matplotlib.backends.backend_pdf import PdfPages
import astropy.io.fits as fits
from warp.Spec1Dtools import FSR_angstrom
from warp.aperture import *
from warp.ECtoID import readparamEC, calcucmnid, chebyshev, legendre
from warp.config import config
from warp.Spec2Dtools import header_key_read

__version__ = "1.2"


# ver 1.2 updated by S.Hamano (2019.09.29)
#   - version was moved from file name to __version__ parameter.
#   - readparamEC and calcucmnid are imported from ECtoID.
#   - read_inputfiles is imported from WINERED_reduction_science.
#


def pix_to_wave(xarray, offset, slope, apn, ftype, xpow, xmin, xmax, opow, omin, omax, cmn):
    onorm = [((offset + slope * i) * 2 - (omin + omax)) / float(omax - omin) for i in apn]
    cmnid = [[0. for i in range(len(apn))] for j in range(xpow)]
    wavearray = [np.zeros(xarray.size) for i in range(len(apn))]
    xnorm = (2. * xarray - (xmin + xmax)) / float(xmax - xmin)

    if ftype == 1:
        for n in range(len(apn)):
            for i in range(xpow):
                for j in range(opow):
                    cmnid[i][n] += cmn[i][j] * chebyshev(j, onorm[n]) / float(offset + slope * apn[n])
                wavearray[n] += cmnid[i][n] * chebyshev(i, xnorm)

    elif ftype == 2:
        for n in range(len(apn)):
            for i in range(xpow):
                for j in range(opow):
                    cmnid[i][n] += cmn[i][j] * legendre(j, onorm[n]) / float(offset + slope * apn[n])
                wavearray[n] += cmnid[i][n] * chebyshev(i, xnorm)

    return wavearray


def tanshift(m, x, param):
    xmin = 1.
    xmax = 2048.
    normx = (2. * x - (xmin + xmax)) / (xmax - xmin)

    return (param[0][0] + param[0][1] * normx + param[0][1] * normx ** 2) * (
            param[1][0] + param[1][1] * m + param[1][2] * m ** 2)


def spectrum_mapping(outputpdf, comp_file, apfile, dy, param):
    if outputpdf.split(".")[-1] != "pdf":
        print("Error: %s is not pdf file." % outputpdf)
        sys.exit()

    pp = PdfPages(outputpdf)
    apset = apertureSet(apfile)
    st = 10

    xpixarray = np.array([1. + i for i in range(int(apset.arrayLength / dy))])
    wavestep = [9000 + i for i in range(0, 5000, 10)]

    offset, slope, apn, ftype, xpow, opow, xterms, xmin, xmax, omin, omax, cmn, features, nfap = readparamEC(comp_file)

    wave = pix_to_wave(xpixarray, offset, slope, apn, ftype, xpow, xmin, xmax, opow, omin, omax, cmn)
    wavemin = [np.amin(wave[i]) for i in range(len(apn))]
    wavemax = [np.amax(wave[i]) for i in range(len(apn))]

    fsr = FSR_angstrom()
    fsr130 = {}
    cutrange = 1.3

    for m in apnum:
        center = (fsr[m][0] + fsr[m][1]) / 2.
        cut_lowlim = (center - (center - fsr[m][0]) * cutrange)
        cut_upplim = (center + (fsr[m][1] - center) * cutrange)
        fsr130[m] = [cut_lowlim, cut_upplim]

    plt.figure(figsize=(10, 10))
    plt.plot([1., 1., apset.arrayLength, apset.arrayLength, 1.], [1., apset.arrayLength, apset.arrayLength, 1., 1.], "k")

    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)

    for i in range(len(apset.echelleOrders)):
        m = apset.echelleOrders[i]
        ap = apset.apertures[m]
        plt.fill(np.hstack((ap.tracex[::st] - 28, ap.tracex[::-st] + 28)), np.hstack((ap.tracey[::st], ap.tracey[::-st])),
                 facecolor="0.5", alpha=0.3)
        plt.text(ap.tracex[-1], ap.tracey[-1] + 10., m, ha="center", fontsize=8)
        plt.text(ap.tracex[0], ap.tracey[0] - 10., m, ha="center", va="top", fontsize=8)
        plt.plot(ap.tracex, ap.tracey, "0.5")

        for j in range(len(wavestep)):
            if wavestep[j] > wavemin[i] and wavestep[j] < wavemax[i]:
                pixid = np.argmin(np.absolute(wave[i] - wavestep[j]))
                if pixid > 50 / dy and pixid < 2000 / dy:
                    theta = tanshift(apnum[i], ap.tracey[int(pixid * dy)], param)
                    if wavestep[j] % 100 == 0:
                        plt.plot([ap.tracex[int(pixid * dy)] - 28, ap.tracex[int(pixid * dy)] + 28],
                                 [ap.tracey[int(pixid * dy)] - 28. * np.tan(theta / 180. * math.pi),
                                  ap.tracey[int(pixid * dy)] + 28. * np.tan(theta / 180. * math.pi)], "0.5",
                                 lw=0.5)
                        plt.text(ap.tracex[int(pixid * dy)], ap.tracey[int(pixid * dy)], wavestep[j], color="white", fontsize=5.,
                                 rotation=theta, va="center", ha="center")
                    else:
                        plt.plot([ap.tracex[int(pixid * dy)] - 14, ap.tracex[int(pixid * dy)] + 14],
                                 [ap.tracey[int(pixid * dy)] - 14. * np.tan(theta / 180. * math.pi),
                                  ap.tracey[int(pixid * dy)] + 14. * np.tan(theta / 180. * math.pi)], "0.5",
                                 lw=1.)

        minid = int(np.argmin(np.absolute(wave[i] - fsr[apnum[i]][0])) * dy)
        minid130 = int(np.argmin(np.absolute(wave[i] - fsr130[apnum[i]][0])) * dy)
        theta = tanshift(apnum[i], ap.tracey[minid], param)
        plt.text(ap.tracex[minid], ap.tracey[minid] + 20, fsr[apnum[i]][0], color="b", fontsize=7., rotation=theta, va="center",
                 ha="center")

        maxid = int(np.argmin(np.absolute(wave[i] - fsr[apnum[i]][1])) * dy)
        maxid130 = int(np.argmin(np.absolute(wave[i] - fsr130[apnum[i]][1])) * dy)
        theta = tanshift(m, ap.tracey[maxid], param)
        plt.text(ap.tracex[maxid], ap.tracey[maxid] - 20, fsr[apnum[i]][1], color="b", fontsize=7., rotation=theta, va="center",
                 ha="center")

        theta_minus = tanshift(m, ap.tracey[maxid:minid:st], param)
        theta_plus = tanshift(m, ap.tracey[minid:maxid:-st], param)
        theta_minus_130 = tanshift(m, ap.tracey[maxid130:minid130:st], param)
        theta_plus_130 = tanshift(m, ap.tracey[minid130:maxid130:-st], param)
        plt.fill(np.hstack((ap.tracex[maxid130:minid130:st] - 28, ap.tracex[minid130:maxid130:-st] + 28)),
                 np.hstack((ap.tracey[maxid130:minid130:st] - 28. * np.tan(theta_minus_130 / 180. * math.pi),
                               ap.tracey[minid130:maxid130:-st] + 28. * np.tan(theta_plus_130 / 180. * math.pi))),
                 facecolor="b", alpha=0.3)
        plt.fill(np.hstack((ap.tracex[maxid:minid:st] - 28, ap.tracex[minid:maxid:-st] + 28)),
                 np.hstack((ap.tracey[maxid:minid:st] - 28. * np.tan(theta_minus / 180. * math.pi),
                               ap.tracey[minid:maxid:-st] + 28. * np.tan(theta_plus / 180. * math.pi))), facecolor="y",
                 alpha=0.8)

    plt.ylim(-50, 2080)
    plt.xlim(-50, 2080)
    plt.savefig(pp, format="pdf", dpi=40)
    plt.clf()

    pp.close()


if __name__ == '__main__':
    conf = config()
    conf.readInputCalib("input_files.txt")
    root_path = os.path.dirname(os.path.abspath(__file__))
    winered_mode = ["WIDE", "HIRES-Y", "HIRES-J"]
    reference_files = {winered_mode[0]: "wide_reference_files.txt",
                       winered_mode[1]: "hiresy_reference_files.txt",
                       winered_mode[2]: "hiresj_reference_files.txt"}

    print(conf.comp_file)

    if os.path.exists(conf.comp_file):
        compfits = fits.open(conf.comp_file)
        comphdr = compfits[0].header
        compfits.close()
        inputmode = header_key_read(comphdr, "INSTMODE")
        inputslit = header_key_read(comphdr, "SLIT")
        compdate = header_key_read(comphdr, "DATE-OBS").replace("-", "")
    else:
        sys.exit("ERROR: \"%s\" does not exist." % (conf.comp_file))

    refpath = root_path + "/reference/" + inputmode + "/"

    if not inputmode in winered_mode:
        print("Caution: WINERED has %s, %s and %s modes." % (winered_mode[0], winered_mode[1], winered_mode[2]))
        sys.exit("ERROR: \"%s\" is not in the list of WINERED observational modes" % inputmode)

    if os.path.exists(refpath + reference_files[inputmode]):
        reffile = open(refpath + reference_files[inputmode], "r")
        reflines = reffile.readlines()
        reffile.close()
    else:
        sys.exit("ERROR: \"%s\" does not exist." % (reference_files[inputmode]))

    refkeywords = ["order", "aperture(pinhole)", "aperture(flat)", "linelist(comp)", "aperture_position",
                   "aperture_position(transform)", "angle_function", "transform_dy"]
    refflag = [False for i in range(len(refkeywords))]
    reffiles = {}
    for i in reflines:
        strline = i.split(":")
        if strline[0] in refkeywords:
            reffiles[strline[0]] = strline[1].split()[0]
            refflag[refkeywords.index(strline[0])] = True

    if False in refflag:
        for i in range(len(refkeywords)):
            if not inputflag[refkeywords[i]]:
                print("ERROR: \"%s\" file is not listed in %s." % (refkeywords[i], reference_files[inputmode]))
        sys.exit()

    [order_word, apfile_pinhole, apfile_flat, linelist, ap_npz, ap_trans_npz_core, anglefile_npz, transdy] = [
        reffiles[refkeywords[i]] for i in range(8)]

    linelist = refpath + linelist
    ap_npz = refpath + ap_npz
    ap_trans_npz_core = refpath + ap_trans_npz_core
    anglefile_npz = refpath + anglefile_npz

    for i in range(8):
        print("%s: %s" % (refkeywords[i], reffiles[refkeywords[i]]))

    transdy = float(transdy)
    apnum = range(int(order_word.split("-")[0]), int(order_word.split("-")[1]) + 1)

    params = np.load(anglefile_npz)
    [t0, t1, t2] = params["p0_yn"]
    [ap0, ap1, ap2] = params["p0_m"]

    spectrum_mapping(sys.argv[1], conf.comp_file, conf.ap_file, conf.dyinput, [[t0, t1, t2], [ap0, ap1, ap2]])
