# -*- coding:utf-8 -*-

__version__ = "1.9"

import sys, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
from matplotlib.ticker import ScalarFormatter
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import ZScaleInterval
import scipy.constants
from scipy import interpolate
from PIL import Image
from aperture import *
import datetime

from Spec1Dtools import openspecfits, FSR_angstrom

colornames = ["b", "g", "r", "c", "m", "y"]
maxfnum = 6
fsrs = FSR_angstrom()


def savefig_small(pngim):
    im = Image.open(pngim)
    im2 = im.convert('RGB').convert('P', palette=Image.ADAPTIVE)
    im2.save(pngim)


def grayConverter(pngim):
    im = Image.open(pngim)
    gray = im.convert("L")
    gray.save(pngim)


def spectrum_plot(ax, spx, spy, xrange, yrange, spx_func="NA", spy_func="NA", colors=["k", "b", "r", "g"], linew=[1.],
                  lines=["-"], labels=[""], colors_func=["k", "b", "r", "g"], linew_func=[1.],
                  lines_func=["--"], labels_func=[""], yshift=[0], yfactor=[1.], xlshift=[0], xvshift=[0],
                  xticks_val="NA", yticks_val="NA", xticks_label="NA", yticks_label="NA", xaxis_label="NA",
                  yaxis_label="NA", grid_flag=False, legend_flag=False, legend_loc=0):
    # spx, spy: list of numpy.array
    # xrange, yrange: [xmin, xmax], [ymin, ymax]

    # colors, linew, lines, labels: list
    # yshift, yfactor, xlinear, xvshift: list of floats
    # xticks_val, yticks_val, xticks_label, yticks_label: list
    # grid_flag, legend_flag: boolean
    # legend_loc: string or int

    if type(spx) is not list:
        print("Parameter \"spx\" is not list.")
        sys.exit()
    if type(spy) is not list:
        print("Parameter \"spy\" is not list.")
        sys.exit()

    Nspec = len(spx)

    colors = colors * math.ceil(Nspec / len(colors)) if len(colors) < Nspec else colors
    linew = linew * math.ceil(Nspec / len(linew)) if len(linew) < Nspec else linew
    lines = lines * math.ceil(Nspec / len(lines)) if len(lines) < Nspec else lines
    labels = labels * math.ceil(Nspec / len(labels)) if len(labels) < Nspec else labels
    yshift = yshift + [0] * (Nspec - len(yshift)) if len(yshift) < Nspec else yshift
    yfactor = yfactor + [1] * (Nspec - len(yfactor)) if len(yfactor) < Nspec else yfactor
    xlshift = xlshift + [0] * (Nspec - len(xlshift)) if len(xlshift) < Nspec else xlshift
    xvshift = xvshift + [0] * (Nspec - len(xvshift)) if len(xvshift) < Nspec else xvshift

    for i in range(Nspec):
        spx_plot = (spx[i] + xlshift[i]) * (1. + xvshift[i] / (scipy.constants.c * 1.e-3))
        spy_plot = spy[i] * yfactor[i] + yshift[i]
        ax.step(spx_plot, spy_plot, where="mid", c=colors[i], label=labels[i], lw=linew[i], ls=lines[i])

    if spx_func != "NA" and spy_func != "NA":
        if type(spx_func) is not list:
            print("Parameter \"spx_func\" is not list.")
            sys.exit()
        if type(spy_func) is not list:
            print("Parameter \"spy_func\" is not list.")
            sys.exit()

        Nspec_func = len(spx_func)

        colors_func = colors_func * math.ceil(Nspec_func / len(colors_func)) if len(
            colors_func) < Nspec_func else colors_func
        linew_func = linew_func * math.ceil(Nspec_func / len(linew_func)) if len(
            linew_func) < Nspec_func else linew_func
        lines_func = lines_func * math.ceil(Nspec_func / len(lines_func)) if len(
            lines_func) < Nspec_func else lines_func
        labels_func = labels_func * math.ceil(Nspec_func / len(labels_func)) if len(
            labels_func) < Nspec_func else labels_func

        for i in range(Nspec_func):
            ax.step(spx_func, spy_func, where="mid", c=colors_func[i], label=labels_func[i], lw=linew_func[i],
                    ls=lines_func[i])

    if xticks_val != "NA":
        ax.set_xticks(xticks_val)
        if xticks_label != "NA":
            ax.set_xticklabels(xticks_label)
    if yticks_val != "NA":
        ax.set_yticks(yticks_val)
        if yticks_label != "NA":
            ax.set_yticklabels(yticks_label)

    if xaxis_label != "NA":
        ax.set_xlabel(xaxis_label)
    if yaxis_label != "NA":
        ax.set_ylabel(yaxis_label)

    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(yrange[0], yrange[1])

    if grid_flag:
        ax.grid()
    if legend_flag:
        if legend_loc != "out":
            ax.legend(loc=legend_loc)
        else:
            ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")


def transmittance_resampling(spx, spy, telx, tely):
    telx_subp = np.linspace(min(telx), max(telx), 10 * telx.size)

    it = interpolate.interp1d(telx, tely, kind="cubic")
    tely_subp = it(telx_subp)

    spy_tm = np.copy(spy)
    for i in range(len(spx)):
        spy_tm[i] = tely_subp[np.argmin(np.absolute(telx_subp - spx[i]))]

    return spy_tm


def plot_all_frames_norm(splist, outputf, apnum, fnum):
    fig = plt.figure(figsize=(25, 6))
    ax = plt.axes([0.04, 0.1, 0.88, 0.8])

    fnum = len(splist)
    divnum = fnum / maxfnum

    pp = PdfPages(outputf)

    for i in range(len(apnum)):
        leftfnum = fnum
        repeatnum = 0

        while leftfnum > 0:
            spxlist, spylist, frameno = [], [], []
            for k in range(min(leftfnum, maxfnum)):
                spx, spy, _, _, _ = openspecfits(splist[repeatnum * maxfnum + k][i])
                spxlist.append(spx)
                spylist.append(spy)
                frameno.append("frame_NO%d" % (k + 1 + repeatnum * maxfnum))

            spectrum_plot(ax, spxlist, spylist, fsrs[apnum[i]], [0, 2.], colors=colornames, linew=[0.5], labels=frameno,
                          xticks_label="Wavelength ($\AA$)", yticks_label="Normalized flux", legend_flag=True,
                          legend_loc="out")
            plt.title("m=%d" % apnum[i])
            plt.savefig(pp, format="pdf")
            plt.cla()

            repeatnum += 1
            leftfnum -= maxfnum

    pp.close()
    plt.close()


def plot_all_frames_flux(splist, outputf, apnum, fnum):
    fig = plt.figure(figsize=(25, 6))
    ax = plt.axes([0.04, 0.1, 0.88, 0.8])

    divnum = fnum / maxfnum

    pp = PdfPages(outputf)

    for i in range(len(apnum)):

        leftfnum = fnum
        repeatnum = 0

        while leftfnum > 0:
            spxlist, spylist, frameno, spymed = [], [], [], []
            for k in range(min(leftfnum, maxfnum)):
                spx, spy, _, _, _ = openspecfits(splist[repeatnum * maxfnum + k][i])
                spxlist.append(spx)
                spylist.append(spy)
                spymed.append(np.median(spy))
                frameno.append("frame_NO%d" % (k + 1 + repeatnum * maxfnum))

            spectrum_plot(ax, spxlist, spylist, fsrs[apnum[i]], [0, max(spymed) * 2.], colors=colornames,
                          linew=[0.5], labels=frameno, xticks_label="Wavelength ($\AA$)", yticks_label="Flux",
                          legend_flag=True, legend_loc="out")
            plt.title("m=%d" % apnum[i])

            plt.savefig(pp, format="pdf")
            plt.cla()

            repeatnum += 1
            leftfnum -= maxfnum

    pp.close()
    plt.close()


def plot_all_frames_flux_BG(splist, bgspec, outputf, apnum, fnum):
    fig = plt.figure(figsize=(25, 7))
    ax1 = plt.axes([0.04, 0.1, 0.88, 0.15])
    ax2 = plt.axes([0.04, 0.3, 0.88, 0.65])

    divnum = fnum / maxfnum

    pp = PdfPages(outputf)

    for i in range(len(apnum)):

        leftfnum = fnum
        repeatnum = 0

        while leftfnum > 0:

            spxlist, spylist, frameno, spymed = [], [], [], []
            bgxlist, bgylist, bgymed = [], [], []
            for k in range(min(leftfnum, maxfnum)):
                spx, spy, _, _, _ = openspecfits(splist[repeatnum * maxfnum + k][i])
                bgx, bgy, _, _, _ = openspecfits(bgspec[repeatnum * maxfnum + k][i])
                spxlist.append(spx)
                spylist.append(spy)
                bgxlist.append(bgx)
                bgylist.append(bgy)
                spymed.append(np.median(spy))
                bgymed.append(np.median(bgy))
                frameno.append("frame_NO%d" % (k + 1 + repeatnum * maxfnum))

            spectrum_plot(ax1, bgxlist, bgylist, fsrs[apnum[i]], [min(bgymed) - 5. * (max(bgymed) - min(bgymed)),
                                                                  max(bgymed) + 5. * (max(bgymed) - min(bgymed))],
                          colors=colornames, linew=[0.5], labels=frameno, xticks_label="Wavelength ($\AA$)",
                          yticks_label="Flux", legend_flag=False)

            spectrum_plot(ax2, spxlist, spylist, fsrs[apnum[i]], [0, max(spymed) * 2.], colors=colornames,
                          linew=[0.5], labels=frameno, xticks_label="Wavelength ($\AA$)", yticks_label="Flux",
                          legend_flag=True, legend_loc="out")

            plt.title("m=%d" % apnum[i])

            plt.savefig(pp, format="pdf")
            plt.cla()

            repeatnum += 1
            leftfnum -= maxfnum

    pp.close()
    plt.close()


def plot_combined_norm(splist, outputf, apnum):
    plt.figure(figsize=(18, 4))

    for i in range(len(apnum)):
        plt.xlim(fsrs[apnum[i]][0], fsrs[apnum[i]][1])
        spx, spy, rcrval1, rcdelt1, rcrpix1 = openspecfits(splist[i])
        plt.step(spx, spy, color="k", lw=3., where="mid")
        plt.ylim(0, np.median(spy) * 1.5)
        plt.xlabel("Wavelength")
        plt.ylabel("Normalized flux")
        plt.title("m=%d" % apnum[i])
        plt.savefig(outputf[i], format="png")
        savefig_small(outputf[i])
        plt.clf()

    plt.close()


def plot_2dimages_mask(maskname, outputfig):
    fimg = fits.open(maskname)
    pixdata = fimg[0].data
    fimg.close()

    fig = plt.figure(facecolor='white', figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    ax1.imshow(pixdata, cmap=cm.gray, vmin=0, vmax=1, origin='lower', interpolation="none")
    plt.title(maskname)

    plt.savefig(outputfig, dpi=400)
    savefig_small(outputfig)
    plt.clf()

    plt.close()


def plot_2dimages(inputimage, outputfig):
    fimg = fits.open(inputimage)
    pixdata = fimg[0].data
    fimg.close()

    blarray = pixdata > 10.
    thres = np.median(pixdata[blarray])

    fig = plt.figure(facecolor='white', figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    ax1.imshow(pixdata, cmap=cm.gray, vmin=-3. * thres, vmax=3. * thres, origin='lower', interpolation="none")

    plt.title(inputimage)

    plt.savefig(outputfig)
    savefig_small(outputfig)
    plt.clf()
    plt.close()


def plot_2dimages_sv(inputimage, outputfig, lowlim="NA", upplim="NA", half=True):
    hdu = fits.open(inputimage)[0]
    wcs = WCS(hdu.header)
    object = hdu.header["OBJECT"]
    nodpos = hdu.header["NODPOS"]
    slitpa = hdu.header["SLT-PA"]
    try:
        seeing = hdu.header["SEEING"]
    except:
        try:
            seeing = float(hdu.header["SVSEEING"])
        except:
            seeing = -10.
    ut = hdu.header["UT"]
    date = hdu.header["DATE-OBS"]
    exptime = hdu.header["EXPTIME"]
    slitwidth = hdu.header["SLIT"]
    naxis1 = int(hdu.header["NAXIS1"])
    naxis2 = int(hdu.header["NAXIS2"])
    pixdata = hdu.data

    thres_date3 = datetime.date(2021, 1, 1)
    SVXoffset1 = 0.
    SVXoffset2 = -22.
    SVXoffset = 0.
    if len(date.split("-")) == 3:
        [year_obs, month_obs, day_obs] = date.split("-")
        obsdate_date = datetime.date(int(year_obs), int(month_obs), int(day_obs))
        if obsdate_date < thres_date3:
            SVXoffset = SVXoffset1
        elif thres_date3 < obsdate_date:
            SVXoffset = SVXoffset2

    if half:
        plt.figure(facecolor='white', figsize=(4, 3))
        ax1 = plt.axes([0.23, 0.16, 0.73, 0.73], projection=wcs)
        fontfactor = 0.75
    else:
        plt.figure(facecolor='white', figsize=(8, 6))
        ax1 = plt.axes([0.16, 0.11, 0.8, 0.8], projection=wcs)
        fontfactor = 1.0

    interval = ZScaleInterval()
    ax1.imshow(interval(pixdata), cmap=cm.hot, origin='lower', interpolation="none")

    aperture_length = 50.
    ps_ratio = 3.392
    if slitwidth == 100:
        [x1, x2, y1, y2] = [553.86 + SVXoffset, 751.86 + SVXoffset, 476.78, 484.88]
    elif slitwidth == 200:
        [x1, x2, y1, y2] = [554.58 + SVXoffset, 752.58 + SVXoffset, 478.28, 494.48]
    else:
        [x1, x2, y1, y2] = [553.86 + SVXoffset, 751.86 + SVXoffset, 476.78, 484.88]

    ax1.plot([x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], color="cyan", linewidth=1.)
    if lowlim != "NA" and upplim != "NA":
        ax1.plot([(x1 + x2) / 2. - lowlim * ps_ratio, (x1 + x2) / 2. - lowlim * ps_ratio],
                 [(y1 + y2) / 2. - aperture_length, (y1 + y2) / 2. + aperture_length], color="springgreen",
                 linewidth=2.)
        ax1.plot([(x1 + x2) / 2. - upplim * ps_ratio, (x1 + x2) / 2. - upplim * ps_ratio],
                 [(y1 + y2) / 2. - aperture_length, (y1 + y2) / 2. + aperture_length], color="springgreen",
                 linewidth=2.)

    ax1.coords.grid(color="white", ls="solid")
    ax1.coords[0].set_axislabel("Right ascension (J2000)", fontsize=13 * fontfactor)
    ax1.coords[1].set_axislabel("Declination (J2000)", fontsize=13 * fontfactor)
    if half:
        ax1.set_xlim(naxis1 / 4., naxis1 * 3. / 4.)
        ax1.set_ylim(naxis2 / 4., naxis2 * 3. / 4.)
    ax1.coords[0].set_ticklabel(size=10 * fontfactor)
    ax1.coords[1].set_ticklabel(size=10 * fontfactor)
    if half:
        plt.figtext(0.25, 0.20,
                    " Object: %s \n Date (UT): %s \n Time (UT): %s \n SV Exposure: %.1f sec \n Nod. position: %s \n Slit PA: %.2f deg \n Seeing: %.2f arcsec" % (
                        object, date, ut, exptime, nodpos, slitpa, seeing),
                    fontsize=6, bbox=dict(facecolor='white', alpha=0.8))
    else:
        plt.figtext(0.19, 0.15,
                    " Object: %s \n Date (UT): %s \n Time (UT): %s \n SV Exposure: %.1f sec \n Nod. position: %s \n Slit PA: %.2f deg \n Seeing: %.2f arcsec" % (
                        object, date, ut, exptime, nodpos, slitpa, seeing),
                    fontsize=12 * fontfactor, bbox=dict(facecolor='white', alpha=0.8))

    # overlay = ax1.get_coords_overlay("fk5")
    # overlay.grid(color="white", ls="dotted")
    # overlay[0].set_axislabel("RA")
    # overlay[1].set_axislabel("Dec")

    ax1.set_title("%s" % (inputimage), fontsize=20 * fontfactor)

    plt.savefig(outputfig, dpi=100)
    savefig_small(outputfig)
    # grayConverter(outputfig)
    plt.clf()
    plt.close()


def snr_plots(lams_sn, snr_val, spfiles, aplength, outputpng):
    fig, ax1 = plt.subplots(figsize=(6, 4))
    twocolors_1 = ["orangered", "r"]
    twocolors_2 = ["dodgerblue", "b"]

    ax1.set_xlabel(r"Wavelength in air ($\AA$)")
    ax1.set_ylabel("Signal-to-noise ratio")
    ax1.set_ylim(0., np.max(snr_val) * 1.2)

    ax2 = ax1.twinx()
    spmax = []
    lammin = []
    lammax = []

    for j in range(aplength):
        spx, spy, a, b, c = openspecfits(spfiles[j])
        ax2.plot(spx, spy, color=twocolors_1[j % 2], alpha=0.5)
        spmax.append(np.median(spy))
        lammin.append(np.min(spx))
        lammax.append(np.max(spx))

    ax2.set_ylabel("Flux (Arbitrary Unit)")
    ax2.set_ylim(0., np.max(spmax) * 3.)
    ax2.set_xlim(np.min(lammin), np.max(lammax))
    ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    for j in range(aplength):
        ax1.plot(lams_sn[j], snr_val[j], color=twocolors_2[j % 2], lw=3.)
        ax1.scatter(lams_sn[j], snr_val[j], color=twocolors_2[j % 2], s=30)

    plt.savefig(outputpng)
    plt.close()


def plot_all_frames_flux_all_orders(splist, outputf, apnum, fnum):
    fig, ax1 = plt.subplots(figsize=(9, 6))

    divnum = fnum / maxfnum

    pp = PdfPages(outputf)

    leftfnum = fnum
    repeatnum = 0

    while leftfnum > 0:

        for i in range(len(apnum)):
            spxlist, spylist, frameno, spymed = [], [], [], []
            for k in range(min(leftfnum, maxfnum)):
                spx, spy, _, _, _ = openspecfits(splist[repeatnum * maxfnum + k][i])
                spxlist.append(spx)
                spylist.append(spy)
                spymed.append(np.median(spy))
                frameno.append("frame_NO%d" % (k + 1 + repeatnum * maxfnum))

            spectrum_plot(ax1, spxlist, spylist, fsrs[apnum[i]], [0, max(spymed) * 1.5], colors=colornames,
                          linew=[0.5], labels=frameno, xticks_label="Wavelength ($\AA$)", yticks_label="Flux",
                          legend_flag=True)

        plt.savefig(pp, format="pdf")
        plt.clf()

        repeatnum += 1
        leftfnum -= maxfnum

    pp.close()
    plt.close()


def peak_count_fwhm(splist, outputf, apnum, fnum, fwhm="INDEF"):
    averageorder = [45, 142, 171]  # WIDE, HIRES-J, HIRES-Y
    averagecount = []
    frameno = [i + 1 for i in range(fnum)]

    for i in range(fnum):
        for j in range(len(apnum)):
            if apnum[j] in averageorder:
                spx, spy, _, _, _ = openspecfits(splist[i][j])
                averagecount.append(np.average(spy))
                selectedorder = apnum[j]

    fig, ax1 = plt.subplots(figsize=(6, 4))

    ax1.set_xlabel("Frame number")
    ax1.set_ylabel("Average count at m=%d" % selectedorder)
    ax1.set_xticks(frameno)
    ax1.set_xticklabels(["No.%d" % frameno[i] for i in range(fnum)])
    ax1.set_ylim(0., np.max(averagecount) * 1.2)
    ax1.set_xlim(min(frameno) - 0.5, max(frameno) + 0.5)
    ax1.scatter(frameno, averagecount, color="r", s=30, label="Count")

    if fwhm != "INDEF":
        ax2 = ax1.twinx()
        ax2.scatter(frameno, fwhm, color="b", s=30, label="FWHM")

        ax2.set_ylabel("FWHM (pix)")
        ax2.set_ylim(0., np.max(fwhm) * 1.2)
        ax2.set_xlim(np.min(frameno) - 0.5, np.max(frameno) + 0.5)

        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax1.legend(h1 + h2, l1 + l2, loc='upper right')
    else:
        ax1.legend(loc='upper right')

    plt.savefig(outputf)
    savefig_small(outputf)
    plt.clf()
    plt.close()


def aperture_plot(datfile, imgfile, lowlim, upplim, bgregion, flagskysub):
    rf = open(datfile, "r")
    rl = rf.readlines()
    rf.close()

    med_plt_x = np.array([float(rl[i].split()[0]) for i in range(len(rl))])
    med_plt_y = np.array([float(rl[i].split()[1]) for i in range(len(rl))])

    xshift_fix = (lowlim + upplim) / 2.

    plt.figure()
    plt.scatter(med_plt_x, med_plt_y)
    plt.plot([xshift_fix, xshift_fix], [-2, 2], "g", label="Center")
    plt.plot([lowlim, lowlim], [-2, 2], "g--", label="Aperture")
    plt.plot([upplim, upplim], [-2, 2], "g--")
    if flagskysub and bgregion != "INDEF":
        regions = bgregion.split(",")
        bg_pixs = []
        for j in range(len(regions)):
            plt.plot([float(regions[j].split(":")[0]), float(regions[j].split(":")[1])], [0.5, 0.5], "r")
            bg_pixs.append(float(regions[j].split(":")[0]))
            bg_pixs.append(float(regions[j].split(":")[1]))

        avg_bg = np.average(bg_pixs)
        plt.text(avg_bg, 0.55, "background", color="r", va="bottom", ha="center")

    plt.ylim(-2.0, 2.0)
    plt.xlabel("Distance from defined aperture center (pix)")
    plt.ylabel("Normalized count")
    plt.legend()
    plt.grid()
    plt.savefig(imgfile, format="png")
    savefig_small(imgfile)
    plt.clf()
    plt.close()


def cosmicRay2dImages(crMask, crfig, apfile, bpmaskflat, xlim1=-30, xlim2=30):
    apset = apertureSet(apfile)
    maskap = apset.apmaskArray(lowlim=xlim1, upplim=xlim2)
    slitcoord = apset.slitcoordArray(lowlim=xlim1, upplim=xlim2)

    apy = np.arange(apset.arrayLength) + 1.
    Xarray, Yarray = np.meshgrid(apy - 1, apy - 1)

    crf = fits.open(crMask + ".fits")
    crdata = crf[0].data
    crf.close()

    bpmaskf = fits.open(bpmaskflat)
    bpfdata = bpmaskf[0].data
    bpmaskf.close()

    reqnotbpf = bpfdata == 0
    reqbpf = bpfdata != 0
    reqcr = (crdata != 0) & reqnotbpf

    pp = PdfPages(crfig)

    plt.figure(figsize=(10, 10))
    plt.scatter(Xarray[reqcr], Yarray[reqcr], s=5., marker=".", color="r",
                label="CR mask ({} pix)".format(np.sum(reqcr)))
    plt.scatter(Xarray[reqbpf], Yarray[reqbpf], s=5., marker=".", color="b",
                label="bad pixel mask ({} pix)".format(np.sum(reqbpf)))
    plt.legend()
    plt.title(crMask)
    plt.xlim(0., apset.arrayLength)
    plt.ylim(0., apset.arrayLength)
    plt.savefig(pp, format="pdf")
    plt.clf()

    plt.figure()
    plt.hist(slitcoord[reqcr], bins=20, histtype="bar")
    plt.xlabel("Slit coordinate (pix)")
    plt.ylabel("Cosmic ray pixels")
    plt.xlim(xlim1, xlim2)
    plt.title(crMask)
    plt.savefig(pp, format="pdf")
    plt.clf()

    plt.hist(maskap[reqcr], bins=20, histtype="bar")
    plt.xlabel("Echelle order m")
    plt.ylabel("Cosmic ray pixels")
    plt.xlim(min(apset.echelleOrders) - 1, max(apset.echelleOrders) + 1)
    plt.title(crMask)
    plt.savefig(pp, format="pdf")
    plt.clf()

    pp.close()


if __name__ == "__main__":
    fi = sys.argv[1:]
    plot_2dimages_sv(fi[0], fi[1])
