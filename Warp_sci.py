#!/usr/bin/env python
# -*- coding:utf-8 -*-

__version__ = "3.8.10"

from pyraf import iraf
import sys, shutil, os, glob, time
import numpy as np
import pathlib
from astropy.io import fits
import argparse

sys.path.append(os.path.dirname(__file__))

from warp.config import config, constant_str_length
from warp.logger import warpLog
from warp.aperture import apertureSet
from warp.centersearch_fortrans import centersearch_fortrans, make_slit_profile
from warp.Spec2Dtools import flatfielding, header_key_read
from warp.apscatter import pyapscatter
from warp.cutransform import cutransform
from warp.Spec1Dtools import pyapall, resample2Dspec, truncate, dispcor_single, cut_1dspec, PyScombine, openspecfits, FSR_angstrom
from warp.ccwaveshift import waveshift_oneorder, PySpecshift, waveshiftClip
from warp.SNratio_estimate import snestimate
from warp.PyContinuum import PyContinuum
from warp.vac2air_spec import vac2air_spec
from warp.plotframes import plot_all_frames_norm, plot_all_frames_flux, plot_all_frames_flux_BG, plot_2dimages_mask, \
    plot_2dimages, snr_plots, plot_combined_norm, plot_2dimages_sv, peak_count_fwhm, aperture_plot, cosmicRay2dImages
from warp.badpixmask import pyfixpix, cosmicRayMask
import tex_source_maker


class shortFile(str):
    def __new__(cls, filename, obj, subs="OBJECT", ext="fits"):
        if len(obj) > 15:
            self = super().__new__(cls, filename.replace(obj, subs))
        else:
            self = super().__new__(cls, filename)
        self.obj = obj
        self.ext = self + ".{}".format(ext)
        self.long = filename
        self.longext = filename + ".{}".format(ext)
        return self


def alternativequestion(question, anss, defans):
    flagans = False
    while not flagans:
        ansinput = input(question)
        if ansinput in anss:
            flagans = True
        else:
            print("Answers: ", anss)

    if ansinput != "":
        return ansinput
    else:
        return defans


def valueQuestion(question, lowlim, upplim, defans):
    flagans = False
    while not flagans:
        ansinput = input(question)
        if ansinput == "":
            flagans = True
        else:
            try:
                ansinputfloat = float(ansinput)
                if ansinputfloat > lowlim and ansinputfloat < upplim:
                    flagans = True
                else:
                    print("Please input float number bewteen {} and {}. (Your answer: {})".format(lowlim, upplim,
                                                                                                  ansinput))
            except:
                print("Please input float number bewteen {} and {}. (Your answer: {})".format(lowlim, upplim, ansinput))

    if ansinput == "":
        return defans
    else:
        return ansinputfloat


def remove_or_move(fname, targetdir, trashdir, flagsave):
    if flagsave:
        shutil.move(fname, targetdir)
    else:
        shutil.move(fname, trashdir)


def remove_or_move_sf(shortfile: shortFile, targetdir, trashdir, flagsave):
    if flagsave:
        shutil.move(shortfile.ext, shortfile.longext)
        shutil.move(shortfile.longext, targetdir)
    else:
        shutil.move(shortfile.ext, trashdir)


def absPathStr(path):
    pathobj = pathlib.Path(path)
    return str(pathobj.resolve()) + "/"


def Warp_sci(listfile, rawdatapath, calibpath, destpath, viewerpath="INDEF", query=False, save=False, parameterfile=None,
             oldformat=False, fastMode=False, autoCalib=False, noreport=False):
    pipeline_ver = __version__
    startTimeSec = time.time()
    startTimeStr = time.ctime()

    conf = config("status.txt")
    fsr = FSR_angstrom()
    abspathDir = os.path.dirname(os.path.abspath(__file__))
    logo = abspathDir + "/winered_logo.pdf"

    constant_str_length("Make the working directory and copy necessary files.")
    # check the paths and make destination directory

    currentpath = absPathStr("./")

    if parameterfile is not None:
        parampath = pathlib.Path(parameterfile)
        paramfile = str(parampath.resolve())
        if query:
            print("\033[31m ERROR: The -q and -p options cannot be set at the same time. \033[0m")
            sys.exit()

    if not os.path.exists(rawdatapath):
        print("\033[31m ERROR: " + rawdatapath + " does not exist. \033[0m")
        sys.exit()
    if not os.path.exists(viewerpath):
        print("\033[31m ERROR: " + viewerpath + " does not exist. \033[0m")
        sys.exit()
    if not os.path.exists(calibpath):
        print("\033[31m ERROR: " + calibpath + " does not exist. \033[0m")
        sys.exit()

    rawdatapath = absPathStr(rawdatapath)
    viewerpath = absPathStr(viewerpath)
    calibpath = absPathStr(calibpath)

    if destpath != currentpath:
        if not os.path.exists(destpath):
            try:
                os.makedirs(destpath)
                print("Created a working directory and jumped in.")
            except:
                print("\033[31m ERROR: tried to make \"%s\" directory, but failed. \033[0m")
                sys.exit()
        else:
            print("\033[31m ERROR: the destination directory already exists. \033[0m")
            sys.exit()

    destpath = absPathStr(destpath)
    shutil.copy(listfile, destpath)
    listfile = os.path.basename(listfile)
    os.chdir(destpath)
    conf.writeStatus(pipeline_ver, startTimeStr)

    # open and read the input file list

    conf.inputDataList(listfile, oldFormat=oldformat)
    conf.writeInputDataList()
    for i in conf.imagelist:
        shutil.copy("{}{}.fits".format(rawdatapath, i), ".")
        iraf.hedit(i, "PIPELINE", pipeline_ver, add="yes", verify="no")
    conf.readInputDataHeader()
    longObjectName = False
    objectNameLimit = 15
    for i in range(conf.objnum):
        if len(conf.objname_obj[i]) > objectNameLimit:
            longObjectName = True
            print("\033[31m WARNING: The object name ({}) of {} ({} characters) is longer than the limit ({} characters). \033[0m".format(
                conf.objname_obj[i], conf.objectlist[i], len(conf.objname_obj[i]), objectNameLimit))

    if autoCalib:
        calibCandidate = glob.glob(calibpath + "*")
        calibPathList = []
        calibStatus = []
        for cc in calibCandidate:
            if os.path.exists(cc + "/input_files.txt") and os.path.exists(cc + "/database"):
                confCal = config()
                confCal.readInputCalib(cc + "/input_files.txt")
                confCal.inputDataList(listfile, oldFormat=oldformat)
                confCal .readInputDataHeader()
                calibPathList.append(cc)
                calibStatus.append(confCal.checkDataStatus(showDetail=False, ignore=True))
        if sum(calibStatus) == 0:
            print("\033[31m ERROR: Appropriate calib data could not be found. \033[0m")
            conf.writeError("ERROR: Appropriate calib data could not be found.")
            sys.exit()
        elif sum(calibStatus) == 1:
            for i in range(len(calibPathList)):
                if calibStatus[i]:
                    calibpath = calibPathList[i] + "/"
            print("Below calib data was successfully selected.")
            print(calibpath)
        else:
            print("\033[31m ERROR: Multiple calib data was found. Please select with -c option. \033[0m")
            conf.writeError("ERROR: Multiple calib data was found. Please select with -c option.")
            sys.exit()

    if os.path.exists(calibpath + "input_files.txt") and os.path.exists(calibpath + "database"):
        calibfiles = glob.glob(calibpath + "*")
        for i in calibfiles:
            if not os.path.exists(destpath + i.split("/")[-1]):
                if os.path.isfile(i):
                    shutil.copy(i, destpath)
                elif os.path.isdir(i):
                    shutil.copytree(i, destpath + i.split("/")[-1])
        print("Calib files were successfully copied.")
    else:
        print("\033[31m ERROR: Invalid calibration path (%s). \033[0m" % calibpath)
        conf.writeError("ERROR: Invalid calibration path (%s)." % calibpath)
        sys.exit()

    # read calibration data names and pipeline parameters and modes.

    constant_str_length("Read setting.")

    conf.readInputCalib("input_files.txt")
    dataStatus = conf.checkDataStatus(statusOutput=True)

    if fastMode:
        conf.setFastModeParam()
    if query:
        if fastMode:
            print("\033[31m WARNING: fast mode setting was dismissed. \033[0m")
        conf.readParamQuery()
    if parameterfile is not None:
        if fastMode:
            print("\033[31m WARNING: fast mode setting was dismissed. \033[0m")
        conf.readParamFile(paramfile)

    # read fits headers

    if conf.flag_svimage:
        for i in conf.svfr_str:
            if not os.path.exists(viewerpath + i):
                conf.flag_svimage = False
        for i in conf.svfr_end:
            if not os.path.exists(viewerpath + i):
                conf.flag_svimage = False

    if conf.flag_svimage:
        for i in range(conf.imnum):
            shutil.copy(viewerpath + conf.svfr_str[i], ".")
            shutil.copy(viewerpath + conf.svfr_end[i], ".")

    conf.writeParam()

    # ヘッダー情報の取得・エラーチェック
    # read the aperture information from ap_file.

    apset = apertureSet(conf.ap_file)
    aplength = len(apset.echelleOrders)
    aptranslist = [conf.aptrans_file + "_{}trans".format(m) for m in apset.echelleOrders]
    cutlength = len(conf.cutrange_list)

    # make the comparison file name for each echelle order.

    comp_file = conf.comp_file.rstrip("fits").rstrip(".")
    comp_file_id = []
    for m in apset.echelleOrders:
        comp_file_id.append(comp_file + "." + apset.apertures[m].make_apnumtext())

    # make the output file names of object frame for...
    # sky subtraction
    # scattered light subtraction (option)
    # flat fielding
    # bad pixel interpolation
    # transformation
    # bad pixel interpolation (option)
    #

    constant_str_length("Reduction start.")

    obj_s_list = [shortFile(conf.objname_obj[i] + "_NO{}_s".format(i + 1), conf.objname_obj[i]) for i in
                  range(conf.objnum)]  # FITS file after sky subtraction: "star_NO1_s.fits"
    obj_s_mask_list = [shortFile("mask_" + obj_s_list[i].long, conf.objname_obj[i]) for i in range(conf.objnum)]
    obj_s_maskfig_list = ["mask_" + obj_s_list[i].long + ".pdf" for i in range(conf.objnum)]
    obj_s_maskflat_list = [shortFile("maskflat_" + obj_s_list[i].long, conf.objname_obj[i]) for i in range(conf.objnum)]
    obj_s_mf1_list = [shortFile("mf1_" + obj_s_list[i].long, conf.objname_obj[i]) for i in range(conf.objnum)]
    obj_s_mf2_list = [shortFile("mf2_" + obj_s_list[i].long, conf.objname_obj[i]) for i in range(conf.objnum)]
    obj_s_noise_list = [shortFile("noise_" + obj_s_list[i].long, conf.objname_obj[i]) for i in range(conf.objnum)]
    bpthres_list = [0. for i in range(conf.objnum)]
    bpnum_list = [0 for i in range(conf.objnum)]
    if conf.flag_apscatter:
        obj_ssc_list = [shortFile(obj_s_list[i].long + "sc", conf.objname_obj[i]) for i in
                        range(conf.objnum)]  # FITS file after scattered light subtraction: "star_NO1_ssc.fits"
        obj_ssc_scatter_list = [shortFile("scatter_" + obj_s_list[i].long, conf.objname_obj[i]) for i in
                        range(conf.objnum)]  # FITS file after scattered light subtraction: "star_NO1_ssc.fits"
    else:
        obj_ssc_list = obj_s_list
    obj_sscf_list = [shortFile(obj_ssc_list[i].long + "f", conf.objname_obj[i]) for i in
                     range(conf.objnum)]  # FITS file after flat fielding: "star*_NO1_s(sc)f.fits"
    obj_sscfm_list = [shortFile(obj_sscf_list[i].long + "m", conf.objname_obj[i]) for i in
                      range(conf.objnum)]  # FITS file after fixpix: "star_NO1_s(sc)fm.fits"
    obj_sscfm_trans_list = [[shortFile(obj_sscfm_list[i].long + "_m{}trans".format(m), conf.objname_obj[i])
                            for m in apset.echelleOrders]
                            for i in range(conf.objnum)]  # FITS file after transform: "star_NO1_s(sc)fm_**trans.fits"
    obj_sscfm_cut_list = [[shortFile(obj_sscfm_list[i].long + "_m{}cut".format(m), conf.objname_obj[i])
                            for m in apset.echelleOrders]
                            for i in range(conf.objnum)]  # FITS file after transform: "star_NO1_s(sc)fm_**trans.fits"

    img_cs_list = [[obj_sscfm_list[i].long + "_m{}trans.png".format(m) for m in apset.echelleOrders] for i in
                   range(conf.objnum)]
    dat_cs_list = [[obj_sscfm_list[i].long + "_m{}trans.dat".format(m) for m in apset.echelleOrders] for i in
                   range(conf.objnum)]

    # make the output file names of sky frame for...
    # flat fielding
    # bad pixel interpolation
    # transformation
    #

    sky_f_list = [shortFile(conf.objname_obj[i] + "_skyNO{}".format(i + 1) + "_f", conf.objname_obj[i]) for i in
                  range(conf.objnum)]  # sky FITS file after flat fielding: "HD***_skyNO1_f.fits"
    sky_fm_list = [shortFile(sky_f_list[i].long + "m", conf.objname_obj[i]) for i in
                   range(conf.objnum)]  # sky FITS file after fixpix: "HD***_skyNO1_fm.fits"
    sky_fm_cut_list = [[shortFile(sky_fm_list[i].long + "_m{}cut".format(m), conf.objname_obj[i])
                         for m in apset.echelleOrders]
                         for i in range(conf.objnum)]  # sky FITS file after transform: "HD***_skyNO1_fm_**trans.fits"
    sky_fm_trans_list = [[shortFile(sky_fm_list[i].long + "_m{}trans".format(m), conf.objname_obj[i])
                         for m in apset.echelleOrders]
                         for i in range(conf.objnum)]  # sky FITS file after transform: "HD***_skyNO1_fm_**trans.fits"

    # reduction from sky subtraction to cutransform

    apscatter_log = "apscatter_log.txt"
    cutransform_log = "cutransform_log.txt"
    centersearch_txt, centersearch_npz = "centersearch_log.txt", "centersearch_log.npz"
    aperture_txt, aperture_npz = "aperture_log.txt", "aperture_log.npz"
    waveshift_txt, waveshift_npz = "waveshift_log.txt", "waveshift_log.npz"
    badpix_txt, badpix_npz = "cosmicray_log.txt", "cosmicray_log.npz"

    log = warpLog(apset.echelleOrders, conf.objnum)

    try:
        for i in range(conf.objnum):
            iraf.imarith(conf.objectlist[i], "-", conf.skylist[i], obj_s_list[i])  # sky subtraction
            if conf.flag_bpmask:
                constant_str_length("Cosmic ray detection for {}".format(obj_s_list[i]))
                if conf.nodpos_obj[i].find("O") == -1:
                    bpnum_list[i], bpthres_list[i] = \
                        cosmicRayMask(obj_s_list[i], conf.objectlist[i], conf.skylist[i], obj_s_mask_list[i],
                                      obj_s_mf1_list[i], obj_s_mf2_list[i], conf.ap_file, conf.mask_file, True,
                                      noisefits=obj_s_noise_list[i], threshold=conf.CRthreshold, varatio=conf.CRvaratio,
                                      slitposratio=conf.CRslitposratio, maxsigma=conf.CRmaxsigma, fixsigma=conf.CRfixsigma)
                else:
                    bpnum_list[i], bpthres_list[i] = \
                        cosmicRayMask(obj_s_list[i], conf.objectlist[i], conf.skylist[i], obj_s_mask_list[i],
                                      obj_s_mf1_list[i], obj_s_mf2_list[i], conf.ap_file, conf.mask_file, False,
                                      noisefits=obj_s_noise_list[i], threshold=conf.CRthreshold, varatio=conf.CRvaratio,
                                      slitposratio=conf.CRslitposratio, maxsigma=conf.CRmaxsigma, fixsigma=conf.CRfixsigma)
                cosmicRay2dImages(obj_s_mask_list[i], obj_s_maskfig_list[i], conf.ap_file, conf.mask_file)

                # bpnum_list[i], bpthres_list[i] = badpixmask(obj_s_list[i], obj_s_mask_list[i], obj_s_mf1_list[i],
                #                                             obj_s_mf2_list[i])
                iraf.imarith(conf.mask_file, "+", obj_s_mask_list[i], obj_s_maskflat_list[i])

            if conf.flag_apscatter:
                # scattered light subtraction (option)
                pyapscatter(obj_s_list[i], obj_ssc_list[i], conf.apsc_maskfile, apscatter_log)

            flatfielding(obj_ssc_list[i], obj_sscf_list[i], conf.flat_file)  # (obj) flat fielding
            if conf.flag_bpmask:
                pyfixpix(obj_sscf_list[i], obj_sscfm_list[i], obj_s_maskflat_list[i])
            else:
                pyfixpix(obj_sscf_list[i], obj_sscfm_list[i], conf.mask_file)  # (obj) bad pixel interpolation

            constant_str_length("Transform {}".format(obj_sscfm_list[i]))

            cutransform(obj_sscfm_list[i], conf.ap_file, cutransform_log, conf.dyinput,
                        conf.fluxinput)  # (obj) transformation

            if conf.flag_skyemission:
                constant_str_length("Reduce {} too.".format(sky_f_list[i]))
                flatfielding(conf.skylist[i], sky_f_list[i], conf.flat_file)  # (sky) flat fielding
                if conf.flag_bpmask:
                    pyfixpix(sky_f_list[i], sky_fm_list[i], obj_s_maskflat_list[i])  # (sky) bad pixel interpolation
                else:
                    pyfixpix(sky_f_list[i], sky_fm_list[i], conf.mask_file)  # (sky) bad pixel interpolation
                cutransform(sky_fm_list[i], conf.ap_file, cutransform_log, conf.dyinput,
                            conf.fluxinput)  # (sky) transformation
    except Exception as e:
        conf.writeError("ERROR in 2d image reduction")
        conf.writeError(e)
        print(e)
        sys.exit()

    if conf.flag_bpmask:
        log.cosmicRayLog(bpthres_list, bpnum_list)
        log.writeCosmicRayLogText(badpix_txt)
        log.writeCosmicRayLogNpz(badpix_npz)

    # center search and bad pix detection

    obj_sscfm_transm_list = obj_sscfm_trans_list

    try:
        if not conf.flag_manual_aperture:
            constant_str_length("Center search of the PSF.")
            xcenter = [[] for i in range(conf.objnum)]
            fwhm = [[] for i in range(conf.objnum)]
            for i in range(conf.objnum):
                for j in range(aplength):
                    if conf.nodpos_obj[i].find("O") == -1 and conf.nodpos_obj[i] != "NA":
                        tmpx, tmpg = centersearch_fortrans(obj_sscfm_trans_list[i][j], aptranslist[j],
                                                           dat_cs_list[i][j],
                                                           abbaflag=True)
                    else:
                        tmpx, tmpg = centersearch_fortrans(obj_sscfm_trans_list[i][j], aptranslist[j],
                                                           dat_cs_list[i][j], abbaflag=False)
                    xcenter[i].append(tmpx)
                    fwhm[i].append(tmpg)

            log.psfLog(xcenter, fwhm)
            log.writePsfLogText(centersearch_txt)
            log.writePsfLogNpz(centersearch_npz)

            # setting the aperture range as 2 sigma.
            lowtrans = [[-1. * np.median(fwhm[i]) + np.median(xcenter[i]) for j in range(aplength)] for i in
                        range(conf.objnum)]
            hightrans = [[np.median(fwhm[i]) + np.median(xcenter[i]) for j in range(aplength)] for i in range(conf.objnum)]
            for i in range(conf.objnum):
                for j in range(aplength):
                    aperture_plot(dat_cs_list[i][j], img_cs_list[i][j], lowtrans[i][j], hightrans[i][j],
                                  conf.skysub_region[i], conf.flag_skysub)

        else:
            # setting aperture as input from input file list
            lowtrans = [[conf.lowlim_input[i] for j in range(aplength)] for i in range(conf.objnum)]
            hightrans = [[conf.upplim_input[i] for j in range(aplength)] for i in range(conf.objnum)]

            for i in range(conf.objnum):
                for j in range(aplength):
                    make_slit_profile(obj_sscfm_trans_list[i][j], aptranslist[j], dat_cs_list[i][j])
                    aperture_plot(dat_cs_list[i][j], img_cs_list[i][j], lowtrans[i][j], hightrans[i][j],
                                  conf.skysub_region[i], conf.flag_skysub)
    except Exception as e:
        conf.writeError("ERROR in PSF cetern search")
        conf.writeError(e)
        print(e)
        sys.exit()


    # setting background subtraction option as input from calibration parameter file

    bgsample = [[conf.skysub_region[i]] for i in range(conf.objnum)]
    if not conf.flag_skysub:
        conf.skysub_mode = "none"

    # make the output file names of object frame for...
    # extracted 1d spectrum (OBJ and SKY)
    # -> cut the edge
    # -> shift (only for OBJ)
    # -> apply the dispersion solution
    # -> cut at FSR * factor
    # -> VAC or AIR, flux or normalized
    #

    constant_str_length("Extract the spectrum from the transformed data.")

    obj_sscfm_transm_1d = [[shortFile(conf.objname_obj[i] + "_NO{}_m{}".format(i + 1, m), conf.objname_obj[i])
                            for m in apset.echelleOrders] for i in range(conf.objnum)]
    obj_sscfm_transm_1d_none = [[shortFile(conf.objname_obj[i] + "_NO{}_m{}_none".format(i + 1, m), conf.objname_obj[i])
                                for m in apset.echelleOrders]
                                for i in range(conf.objnum)]
    obj_sscfm_transm_1d_bg = [[shortFile(conf.objname_obj[i] + "_NO{}_m{}_bg".format(i + 1, m), conf.objname_obj[i])
                               for m in apset.echelleOrders]
                              for i in range(conf.objnum)]
    obj_sscfm_transm_1d_bgcut = [[shortFile(conf.objname_obj[i] + "_NO{}_m{}_bgc".format(i + 1, m), conf.objname_obj[i])
                                  for m in apset.echelleOrders]
                                 for i in range(conf.objnum)]
    obj_sscfm_transm_1d_bgcutw = [[shortFile(conf.objname_obj[i] + "_NO{}_m{}_bgcw".format(i + 1, m), conf.objname_obj[i])
                                   for m in apset.echelleOrders]
                                  for i in range(conf.objnum)]
    obj_sscfm_transm_1dap = [[] for i in range(conf.objnum)]
    obj_sscfm_transm_1d_noneap = [[] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcut = [[shortFile(obj_sscfm_transm_1d[i][j].long + "c", conf.objname_obj[i])
                               for j in range(aplength)] for i in range(conf.objnum)]
    file_matrix_forwaveshift = [[shortFile(obj_sscfm_transm_1d[i][j].long + "c", conf.objname_obj[i])
                                 for i in range(conf.objnum)] for j in range(aplength)]

    obj_sscfm_transm_1dcuts = [[shortFile(obj_sscfm_transm_1d[i][j].long + "cs", conf.objname_obj[i])
                                for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw = [[shortFile(obj_sscfm_transm_1d[i][j].long + "csw", conf.objname_obj[i])
                                 for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac = [
        [[shortFile(obj_sscfm_transm_1d[i][j].long + "_fsr%.2f_VAC" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac_norm = [
        [[shortFile(obj_sscfm_transm_1d[i][j].long + "_fsr%.2f_VAC_norm" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac_cont = [
        [[shortFile(obj_sscfm_transm_1d[i][j].long + "_fsr%.2f_VAC_cont" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw_fsr_air = [
        [[shortFile(obj_sscfm_transm_1d[i][j].long + "_fsr%.2f_AIR" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw_fsr_air_norm = [
        [[shortFile(obj_sscfm_transm_1d[i][j].long + "_fsr%.2f_AIR_norm" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw_fsr_air_cont = [
        [[shortFile(obj_sscfm_transm_1d[i][j].long + "_fsr%.2f_AIR_cont" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]

    sky_fm_trans_1d = [[shortFile(sky_fm_trans_list[i][j].long + "1d", conf.objname_obj[i]) for j in range(aplength)]
                       for i in range(conf.objnum)]
    sky_fm_trans_1dap = [[] for i in range(conf.objnum)]
    sky_fm_trans_1dcut = [[shortFile(sky_fm_trans_list[i][j].long + "1dcut", conf.objname_obj[i])
                           for j in range(aplength)] for i in range(conf.objnum)]
    sky_fm_trans_1dcutw = [[shortFile(sky_fm_trans_list[i][j].long + "1dcutw", conf.objname_obj[i])
                            for j in range(aplength)] for i in range(conf.objnum)]
    sky_fm_trans_1dcutw_fsr_vac = [
        [[shortFile(sky_fm_trans_list[i][j].long + "1dcutw_fsr%.2f_VAC" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]
    sky_fm_trans_1dcutw_fsr_air = [
        [[shortFile(sky_fm_trans_list[i][j].long + "1dcutw_fsr%.2f_AIR" % conf.cutrange_list[k], conf.objname_obj[i])
          for k in range(cutlength)] for j in range(aplength)] for i in range(conf.objnum)]

    # extracted 2d spectrum (only for OBJ)
    # -> cut the edge
    # -> shift
    # -> apply the dispersion solution
    # -> VAC or AIR
    #

    obj_sscfm_transm_2d = [[shortFile(obj_sscfm_transm_list[i][j].long + "2d", conf.objname_obj[i])
                            for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2d_resample = [[shortFile(obj_sscfm_transm_list[i][j].long + "2d_resample", conf.objname_obj[i])
                                     for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dap = [[] for i in range(conf.objnum)]
    obj_sscfm_transm_2dap_resample = [[] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcut = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcut", conf.objname_obj[i])
                               for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcut_resample = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcut_resample", conf.objname_obj[i])
                                        for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcuts = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcuts", conf.objname_obj[i])
                                for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcuts_resample = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcuts_resample", conf.objname_obj[i])
                                         for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcutsw_vac = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcutsw_VAC", conf.objname_obj[i])
                                     for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcutsw_resample_vac = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcutsw_resample_VAC", conf.objname_obj[i])
                                              for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcutsw_air = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcutsw_AIR", conf.objname_obj[i])
                                     for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_2dcutsw_resample_air = [[shortFile(obj_sscfm_transm_list[i][j].long + "2dcutsw_resample_AIR", conf.objname_obj[i])
                                              for j in range(aplength)] for i in range(conf.objnum)]

    aplow_log = [[] for i in range(conf.objnum)]
    aphigh_log = [[] for i in range(conf.objnum)]

    try:
        for i in range(conf.objnum):
            for j in range(aplength):
                apsetTrans = apertureSet(aptranslist[j])
                apTrans = apsetTrans.apertures[apset.echelleOrders[j]]

                # change the aperture range if the set range exceed limits.
                lowlim = apTrans.apLow if apTrans.apLow > lowtrans[i][j] else lowtrans[i][j]
                highlim = apTrans.apHigh if apTrans.apHigh < hightrans[i][j] else hightrans[i][j]
                apsetTrans.apertures[apset.echelleOrders[j]].apLow = lowlim
                apsetTrans.apertures[apset.echelleOrders[j]].apHigh = highlim

                # set aperture info.
                apsetTrans.write_apdatabase(obj_sscfm_trans_list[i][j], bgsample[i])

                # save aperture range for logging
                aplow_log[i].append(lowlim)
                aphigh_log[i].append(highlim)
                apname_trans = apTrans.make_apnumtext()

                # extract 1d spectrum (OBJ)
                pyapall(obj_sscfm_transm_list[i][j], obj_sscfm_transm_1d[i][j], obj_sscfm_trans_list[i][j],
                        conf.skysub_mode, "onedspec")
                if conf.flag_skysub:
                    pyapall(obj_sscfm_transm_list[i][j], obj_sscfm_transm_1d_none[i][j],
                            obj_sscfm_trans_list[i][j], "none", "onedspec")
                    iraf.sarith(obj_sscfm_transm_1d_none[i][j] + "." + apname_trans, "-",
                                obj_sscfm_transm_1d[i][j] + "." + apname_trans, obj_sscfm_transm_1d_bg[i][j])
                    truncate(obj_sscfm_transm_1d_bg[i][j], obj_sscfm_transm_1d_bgcut[i][j])
                    obj_sscfm_transm_1d_noneap[i].append(obj_sscfm_transm_1d_none[i][j] + "." + apname_trans)

                truncate(obj_sscfm_transm_1d[i][j] + "." + apname_trans, obj_sscfm_transm_1dcut[i][j])
                obj_sscfm_transm_1dap[i].append(shortFile(obj_sscfm_transm_1d[i][j] + "." + apname_trans, conf.objname_obj[i]))

                # extract 2d spectrum (OBJ)
                if conf.flag_extract2d:
                    pyapall(obj_sscfm_transm_list[i][j], obj_sscfm_transm_2d[i][j], obj_sscfm_trans_list[i][j],
                            conf.skysub_mode, "strip")
                    resampleFlag = resample2Dspec(obj_sscfm_transm_list[i][j], obj_sscfm_transm_2d_resample[i][j] + "." +
                                                  apname_trans, obj_sscfm_transm_2d[i][j] + "." + apname_trans, obj_sscfm_trans_list[i][j])
                    truncate(obj_sscfm_transm_2d[i][j] + "." + apname_trans, obj_sscfm_transm_2dcut[i][j])
                    obj_sscfm_transm_2dap[i].append(shortFile(obj_sscfm_transm_2d[i][j].long + "." + apname_trans, conf.objname_obj[i]))
                    obj_sscfm_transm_2dap_resample[i].append(shortFile(obj_sscfm_transm_2d_resample[i][j].long + "." + apname_trans, conf.objname_obj[i]))
                    if resampleFlag:
                        truncate(obj_sscfm_transm_2d_resample[i][j] + "." + apname_trans, obj_sscfm_transm_2dcut_resample[i][j])

                # extract 1d spectrum (SKY)
                if conf.flag_skyemission:
                    pyapall(sky_fm_trans_list[i][j], sky_fm_trans_1d[i][j], obj_sscfm_trans_list[i][j], "none",
                            "onedspec")
                    truncate(sky_fm_trans_1d[i][j] + "." + apname_trans, sky_fm_trans_1dcut[i][j])
                    sky_fm_trans_1dap[i].append(shortFile(sky_fm_trans_1d[i][j] + "." + apname_trans, conf.objname_obj[i]))
    except Exception as e:
        conf.writeError("ERROR in extraction")
        conf.writeError(e)
        print(e)
        sys.exit()


    log.apertureLog(aplow_log, aphigh_log)
    log.writeApertureLogText(aperture_txt)
    log.writeApertureLogNpz(aperture_npz)

    constant_str_length("Shift correction and continuum fitting.")

    # measuring the shift from the frame that have the highest count.
    try:
        if conf.objnum > 1 and conf.flag_wsmeasure:
            counts = [0. for i in range(conf.objnum)]
            for i in range(conf.objnum):
                for j in range(aplength):
                    spx, spy, _, _, _ = openspecfits(obj_sscfm_transm_1dcut[i][j])
                    counts[i] += np.median(spy)
            refid = np.argmax(counts) if conf.objnum > 2 else 0

            wshift_matrix = []
            for j in range(aplength):
                wshift_matrix.append(waveshift_oneorder(file_matrix_forwaveshift[j], refid))

            shift_average, shift_calcnum, shift_stddev, wshift_flag = waveshiftClip(wshift_matrix, conf.objnum, aplength)
            if conf.flag_wsmanual:
                shift_cor = conf.waveshift_man
            else:
                shift_cor = [shift_average[i] if wshift_flag[i] == 0 else 0. for i in range(conf.objnum)]
        else:
            wshift_matrix = [[0. for i in range(conf.objnum)] for j in range(aplength)]
            shift_average = [0. for i in conf.waveshift_man]
            shift_stddev = [0. for i in conf.waveshift_man]
            shift_calcnum = [0 for i in conf.waveshift_man]
            if conf.flag_wsmanual:
                shift_cor = conf.waveshift_man
            else:
                shift_cor = [0. for i in conf.waveshift_man]

        log.waveshiftLog(wshift_matrix, shift_average, shift_stddev, shift_calcnum, shift_cor)
        log.writeWaveshiftLogText(waveshift_txt)
        log.writeWaveshiftLogNpz(waveshift_npz)

        for i in range(conf.objnum):
            for j in range(aplength):
                # apply the dispersion solution to 1d spectra (OBJ)
                if conf.objnum > 1 and conf.flag_wscorrect:
                    # correct the shift
                    PySpecshift(obj_sscfm_transm_1dcut[i][j], obj_sscfm_transm_1dcuts[i][j], shift_cor[i])
                else:
                    iraf.scopy(obj_sscfm_transm_1dcut[i][j], obj_sscfm_transm_1dcuts[i][j])
                dispcor_single(obj_sscfm_transm_1dcuts[i][j], obj_sscfm_transm_1dcutsw[i][j], comp_file_id[j])
                if conf.flag_skysub:
                    dispcor_single(obj_sscfm_transm_1d_bgcut[i][j], obj_sscfm_transm_1d_bgcutw[i][j], comp_file_id[j])
                    iraf.hedit(obj_sscfm_transm_1d_bgcutw[i][j], "AIRORVAC", "vac", add="yes", verify="no")

                iraf.hedit(obj_sscfm_transm_1dcutsw[i][j], "AIRORVAC", "vac", add="yes", verify="no")

                # cut, normalize, convert to air wavelength for 1d spectra (OBJ)
                for k in range(cutlength):
                    cut_1dspec(obj_sscfm_transm_1dcutsw[i][j], obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k],
                               conf.cutrange_list[k], apset.echelleOrders[j], fsr)
                    if apset.echelleOrders[j] < 100.:
                        PyContinuum(obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k],
                                    obj_sscfm_transm_1dcutsw_fsr_vac_norm[i][j][k], 2., 3., 15, "spline3", "*", 10,
                                    continuumspec=obj_sscfm_transm_1dcutsw_fsr_vac_cont[i][j][k])
                    else:
                        PyContinuum(obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k],
                                    obj_sscfm_transm_1dcutsw_fsr_vac_norm[i][j][k], 2., 3., 7, "spline3", "*", 2,
                                    continuumspec=obj_sscfm_transm_1dcutsw_fsr_vac_cont[i][j][k])
                    vac2air_spec(obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k], obj_sscfm_transm_1dcutsw_fsr_air[i][j][k])
                    vac2air_spec(obj_sscfm_transm_1dcutsw_fsr_vac_norm[i][j][k],
                                 obj_sscfm_transm_1dcutsw_fsr_air_norm[i][j][k])
                    vac2air_spec(obj_sscfm_transm_1dcutsw_fsr_vac_cont[i][j][k],
                                 obj_sscfm_transm_1dcutsw_fsr_air_cont[i][j][k])

                # shift, apply dispersion solution, convert to air wavelength for 2d spectra (OBJ)
                if conf.flag_extract2d:
                    if conf.objnum > 1 and conf.flag_wscorrect:
                        PySpecshift(obj_sscfm_transm_2dcut[i][j], obj_sscfm_transm_2dcuts[i][j], shift_cor[i])
                        if os.path.exists(obj_sscfm_transm_2dcut_resample[i][j] + ".fits"):
                            PySpecshift(obj_sscfm_transm_2dcut_resample[i][j], obj_sscfm_transm_2dcuts_resample[i][j], shift_cor[i])
                    else:
                        iraf.scopy(obj_sscfm_transm_2dcut[i][j], obj_sscfm_transm_2dcuts[i][j])
                        if os.path.exists(obj_sscfm_transm_2dcut_resample[i][j] + ".fits"):
                            iraf.scopy(obj_sscfm_transm_2dcut_resample[i][j], obj_sscfm_transm_2dcuts_resample[i][j])
                    dispcor_single(obj_sscfm_transm_2dcuts[i][j], obj_sscfm_transm_2dcutsw_vac[i][j], comp_file_id[j])
                    iraf.hedit(obj_sscfm_transm_2dcutsw_vac[i][j], "AIRORVAC", "vac", add="yes", verify="no")
                    vac2air_spec(obj_sscfm_transm_2dcutsw_vac[i][j], obj_sscfm_transm_2dcutsw_air[i][j])
                    if os.path.exists(obj_sscfm_transm_2dcuts_resample[i][j] + ".fits"):
                        dispcor_single(obj_sscfm_transm_2dcuts_resample[i][j], obj_sscfm_transm_2dcutsw_resample_vac[i][j], comp_file_id[j])
                        iraf.hedit(obj_sscfm_transm_2dcutsw_resample_vac[i][j], "AIRORVAC", "vac", add="yes", verify="no")
                        vac2air_spec(obj_sscfm_transm_2dcutsw_resample_vac[i][j], obj_sscfm_transm_2dcutsw_resample_air[i][j])

                if conf.flag_skyemission:
                    # apply dispersion solution, cut, convert to air wavelength for 1d spectra (SKY)
                    dispcor_single(sky_fm_trans_1dcut[i][j], sky_fm_trans_1dcutw[i][j], comp_file_id[j])
                    iraf.hedit(sky_fm_trans_1dcutw[i][j], "AIRORVAC", "vac", add="yes", verify="no")

                    for k in range(cutlength):
                        cut_1dspec(sky_fm_trans_1dcutw[i][j], sky_fm_trans_1dcutw_fsr_vac[i][j][k],
                                   conf.cutrange_list[k], apset.echelleOrders[j], fsr)
                        vac2air_spec(sky_fm_trans_1dcutw_fsr_vac[i][j][k], sky_fm_trans_1dcutw_fsr_air[i][j][k])
    except Exception as e:
        conf.writeError("ERROR in reducing the 1D spectra")
        conf.writeError(e)
        print(e)
        sys.exit()

    # measuring signal-to-noize ratio, combine, normalize, conversion to air wavelength

    SNmatrix_fsr = [[[shortFile(obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k], conf.objname_obj[i])
                      for i in range(conf.objnum)]
                     for j in range(aplength)]
                    for k in range(cutlength)]
    SNoutput_fsr = [
        ["SNratio_m%d_fsr%.2f.dat" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for k in
        range(cutlength)]
    SN_png = ["SNratio_fsr%.2f.png" % conf.cutrange_list[k] for k in range(cutlength)]

    combined_spec_fsr_vac = [
        [shortFile(conf.objnameRep + "_sum_m%d_fsr%.2f_VAC" % (apset.echelleOrders[j], conf.cutrange_list[k]),
                   conf.objnameRep) for j in range(aplength)] for k in range(cutlength)]
    combined_spec_fsr_vac_norm = [
        [shortFile(conf.objnameRep + "_sum_m%d_fsr%.2f_VAC_norm" % (apset.echelleOrders[j], conf.cutrange_list[k]),
                   conf.objnameRep) for j in range(aplength)] for k in range(cutlength)]
    combined_spec_fsr_vac_cont = [
        [shortFile(conf.objnameRep + "_sum_m%d_fsr%.2f_VAC_cont" % (apset.echelleOrders[j], conf.cutrange_list[k]),
         conf.objnameRep) for j in range(aplength)] for k in range(cutlength)]
    combined_spec_fsr_air = [
        [shortFile(conf.objnameRep + "_sum_m%d_fsr%.2f_AIR" % (apset.echelleOrders[j], conf.cutrange_list[k]),
         conf.objnameRep) for j in range(aplength)] for k in range(cutlength)]
    combined_spec_fsr_air_norm = [
        [shortFile(conf.objnameRep + "_sum_m%d_fsr%.2f_AIR_norm" % (apset.echelleOrders[j], conf.cutrange_list[k]),
         conf.objnameRep) for j in range(aplength)] for k in range(cutlength)]
    combined_spec_fsr_air_cont = [
        [shortFile(conf.objnameRep + "_sum_m%d_fsr%.2f_AIR_cont" % (apset.echelleOrders[j], conf.cutrange_list[k]),
         conf.objnameRep) for j in range(aplength)] for k in range(cutlength)]
    combined_spec_fsr_vac_norm_combined = [shortFile(conf.objnameRep + "_sum_fsr%.2f_VAC_norm_comb" % (conf.cutrange_list[k]),
                                                     conf.objnameRep) for k in range(cutlength)]
    combined_spec_fsr_air_norm_combined = [shortFile(conf.objnameRep + "_sum_fsr%.2f_AIR_norm_comb" % (conf.cutrange_list[k]),
                                                     conf.objnameRep) for k in range(cutlength)]

    try:
        if conf.objnum > 1:
            lams_sn = [[] for k in range(cutlength)]
            snr_val = [[] for k in range(cutlength)]
            for k in range(cutlength):
                for j in range(aplength):
                    tmplam, tmpsnr = snestimate(SNmatrix_fsr[k][j], SNoutput_fsr[k][j])
                    lams_sn[k].append(tmplam)
                    snr_val[k].append(tmpsnr)

                    PyScombine(SNmatrix_fsr[k][j], combined_spec_fsr_vac[k][j])
                    if apset.echelleOrders[j] < 100:
                        PyContinuum(combined_spec_fsr_vac[k][j], combined_spec_fsr_vac_norm[k][j], 2., 3., 15, "spline3",
                                    "*", 10,
                                    continuumspec=combined_spec_fsr_vac_cont[k][j])  # the same parameter as pipeline ver2
                    else:
                        PyContinuum(combined_spec_fsr_vac[k][j], combined_spec_fsr_vac_norm[k][j], 2., 3., 7, "spline3",
                                    "*", 2, continuumspec=combined_spec_fsr_vac_cont[k][j])
                    vac2air_spec(combined_spec_fsr_vac[k][j], combined_spec_fsr_air[k][j])
                    vac2air_spec(combined_spec_fsr_vac_norm[k][j], combined_spec_fsr_air_norm[k][j])
                    vac2air_spec(combined_spec_fsr_vac_cont[k][j], combined_spec_fsr_air_cont[k][j])

                PyScombine(combined_spec_fsr_vac_norm[k], combined_spec_fsr_vac_norm_combined[k])
                PyScombine(combined_spec_fsr_air_norm[k], combined_spec_fsr_air_norm_combined[k])

        else:
            for k in range(cutlength):
                for j in range(aplength):
                    iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_air[0][j][k], combined_spec_fsr_air[k][j])
                    iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_air_norm[0][j][k], combined_spec_fsr_air_norm[k][j])
                    iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_air_cont[0][j][k], combined_spec_fsr_air_cont[k][j])
                    iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_vac[0][j][k], combined_spec_fsr_vac[k][j])
                    iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_vac_norm[0][j][k], combined_spec_fsr_vac_norm[k][j])
                    iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_vac_cont[0][j][k], combined_spec_fsr_vac_cont[k][j])
    except Exception as e:
        conf.writeError("ERROR in combine")
        conf.writeError(e)
        print(e)
        sys.exit()


    # make plots

    obj_sscfm_transm_1dcutsw_fsr_vac_norm_105 = [
        [obj_sscfm_transm_1dcutsw_fsr_vac_norm[i][j][0] for j in range(aplength)] for i in range(conf.objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac_flux_105 = [[obj_sscfm_transm_1dcutsw_fsr_vac[i][j][0] for j in range(aplength)]
                                                 for i in range(conf.objnum)]
    obj_comb_norm_png = [conf.objnameRep + "_AIRnorm_combined_m%d.png" % apset.echelleOrders[j] for j in
                         range(aplength)]

    plot_all_frames_norm(obj_sscfm_transm_1dcutsw_fsr_vac_norm_105, conf.objnameRep + "_VACnorm.pdf",
                         apset.echelleOrders, conf.objnum, fsr)
    if conf.flag_skysub:
        plot_all_frames_flux_BG(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, obj_sscfm_transm_1d_bgcutw,
                                conf.objnameRep + "_VACflux.pdf", apset.echelleOrders, conf.objnum, fsr)
    else:
        plot_all_frames_flux(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, conf.objnameRep + "_VACflux.pdf",
                             apset.echelleOrders, conf.objnum, fsr)

    plot_combined_norm(combined_spec_fsr_air_norm[0], obj_comb_norm_png, apset.echelleOrders, fsr)

    if conf.flag_bpmask:
        for i in range(conf.objnum):
            plot_2dimages_mask(obj_s_mask_list[i].ext, obj_s_mask_list[i].long + ".png")

    for i in range(conf.imnum):
        plot_2dimages(conf.imagelist[i] + ".fits", conf.imagelist[i] + ".png")
    for i in range(conf.objnum):
        plot_2dimages(obj_sscfm_list[i].ext, obj_sscfm_list[i].long + ".png")

    if conf.objnum > 1:
        for k in range(cutlength):
            snr_plots(lams_sn[k], snr_val[k], combined_spec_fsr_vac[k], aplength, SN_png[k])

    constant_str_length("Plot data.")

    countfwhmpng = "count_and_fwhm.png"
    if not conf.flag_manual_aperture:
        fwhmMedian = np.median(log.psfWidth, axis=1)
        peak_count_fwhm(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, countfwhmpng, apset.echelleOrders, conf.objnum,
                        fwhm=fwhmMedian)
    else:
        peak_count_fwhm(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, countfwhmpng, apset.echelleOrders, conf.objnum)

    if conf.flag_svimage:
        for i in range(conf.imnum):
            if conf.imobjid[i] != "sky":
                plot_2dimages_sv(conf.svfr_str[i], conf.imagelist[i] + "_expstart.png",
                                 lowlim=np.average(lowtrans[conf.imobjid[i]]),
                                 upplim=np.average(hightrans[conf.imobjid[i]]))
                plot_2dimages_sv(conf.svfr_end[i], conf.imagelist[i] + "_expend.png",
                                 lowlim=np.average(lowtrans[conf.imobjid[i]]),
                                 upplim=np.average(hightrans[conf.imobjid[i]]))
            else:
                plot_2dimages_sv(conf.svfr_str[i], conf.imagelist[i] + "_expstart.png")
                plot_2dimages_sv(conf.svfr_end[i], conf.imagelist[i] + "_expend.png")

    # make directries, and move files to appropriate directories
    # trash directory

    trashdir = "Trash"
    os.makedirs(trashdir)
    constant_str_length("Clean up the data.")

    # 1d spectra of OBJ

    onedspec_dirnames = ["AIR_flux", "AIR_norm", "AIR_cont", "VAC_flux", "VAC_norm", "VAC_cont"]

    onedspec_frames_dirs = [[["%s_NO%d/onedspec/%s/fsr%.2f/" % (
        conf.objname_obj[i], (i + 1), onedspec_dirnames[n], conf.cutrange_list[k]) for k in range(cutlength)]
                            for n in range(6)]
                            for i in range(conf.objnum)]
    onedspec_sum_dirs = [
        ["%s_sum/%s/fsr%.2f/" % (conf.objnameRep, onedspec_dirnames[n], conf.cutrange_list[k]) for k in
         range(cutlength)] for n in range(6)]

    for n in range(6):
        for k in range(cutlength):
            for i in range(conf.objnum):
                os.makedirs(onedspec_frames_dirs[i][n][k])
            os.makedirs(onedspec_sum_dirs[n][k])

    for j in range(aplength):
        for k in range(cutlength):
            for i in range(conf.objnum):
                remove_or_move_sf(obj_sscfm_transm_1dcutsw_fsr_air[i][j][k], onedspec_frames_dirs[i][0][k],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1dcutsw_fsr_air_norm[i][j][k], onedspec_frames_dirs[i][1][k],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1dcutsw_fsr_air_cont[i][j][k], onedspec_frames_dirs[i][2][k],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k], onedspec_frames_dirs[i][3][k],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1dcutsw_fsr_vac_norm[i][j][k], onedspec_frames_dirs[i][4][k],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1dcutsw_fsr_vac_cont[i][j][k], onedspec_frames_dirs[i][5][k],
                               trashdir, 1)

            remove_or_move_sf(combined_spec_fsr_air[k][j], onedspec_sum_dirs[0][k], trashdir, 1)
            remove_or_move_sf(combined_spec_fsr_air_norm[k][j], onedspec_sum_dirs[1][k], trashdir, 1)
            remove_or_move_sf(combined_spec_fsr_air_cont[k][j], onedspec_sum_dirs[2][k], trashdir, 1)
            remove_or_move_sf(combined_spec_fsr_vac[k][j], onedspec_sum_dirs[3][k], trashdir, 1)
            remove_or_move_sf(combined_spec_fsr_vac_norm[k][j], onedspec_sum_dirs[4][k], trashdir, 1)
            remove_or_move_sf(combined_spec_fsr_vac_cont[k][j], onedspec_sum_dirs[5][k], trashdir, 1)

    if conf.objnum > 1:
        for k in range(cutlength):
            remove_or_move_sf(combined_spec_fsr_vac_norm_combined[k], onedspec_sum_dirs[1][k], trashdir, 1)
            remove_or_move_sf(combined_spec_fsr_air_norm_combined[k], onedspec_sum_dirs[4][k], trashdir, 1)

    # signal-to-noize ratio

    SNRdat_frames_dirs = ["SNR_dat/fsr%.2f" % (conf.cutrange_list[k]) for k in range(cutlength)]

    if conf.objnum > 1:
        for k in range(cutlength):
            os.makedirs(SNRdat_frames_dirs[k])
            for j in range(aplength):
                remove_or_move("SNratio_m%d_fsr%.2f.dat" % (apset.echelleOrders[j], conf.cutrange_list[k]),
                               SNRdat_frames_dirs[k],
                               trashdir, 1)
            remove_or_move(SN_png[k], SNRdat_frames_dirs[k], trashdir, 1)

    # images

    images_frames_dirs_sp = ["%s_NO%d/images/spatial_profile/" % (conf.objname_obj[i], (i + 1)) for i in
                             range(conf.objnum)]

    for i in range(conf.objnum):
        os.makedirs(images_frames_dirs_sp[i])
        for j in range(aplength):
            remove_or_move(img_cs_list[i][j], images_frames_dirs_sp[i], trashdir, 1)
            remove_or_move(dat_cs_list[i][j], images_frames_dirs_sp[i], trashdir, 1)

    images_frames_dirs_2d = ["%s_NO%d/images/2d_image/" % (conf.objname_obj[i], (i + 1)) for i in range(conf.objnum)]

    for i in range(conf.objnum):
        os.makedirs(images_frames_dirs_2d[i])
        remove_or_move(obj_sscfm_list[i].long + ".png", images_frames_dirs_2d[i], trashdir, 1)

    if conf.flag_bpmask:
        images_frames_dirs_bp = ["%s_NO%d/images/badpixmask/" % (conf.objname_obj[i], (i + 1)) for i in
                                 range(conf.objnum)]

        for i in range(conf.objnum):
            os.makedirs(images_frames_dirs_bp[i])
            remove_or_move(obj_s_mask_list[i].long + ".png", images_frames_dirs_bp[i], trashdir, 1)
            remove_or_move(obj_s_maskfig_list[i], images_frames_dirs_bp[i], trashdir, 1)

    # 2d spectra of OBJ
    if conf.flag_extract2d:
        twodspec_dirnames = ["AIR", "VAC"]

        twodspec_frames_dirs = [
            ["%s_NO%d/twodspec/%s/" % (conf.objname_obj[i], (i + 1), twodspec_dirnames[n]) for n in range(2)]
            for i in range(conf.objnum)]

        for i in range(conf.objnum):
            os.makedirs(twodspec_frames_dirs[i][0])
            os.makedirs(twodspec_frames_dirs[i][1])
            for j in range(aplength):
                remove_or_move_sf(obj_sscfm_transm_2dcutsw_air[i][j], twodspec_frames_dirs[i][0], trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_2dcutsw_vac[i][j], twodspec_frames_dirs[i][1], trashdir, 1)
                if os.path.exists(obj_sscfm_transm_2dcutsw_resample_air[i][j].ext):
                    remove_or_move_sf(obj_sscfm_transm_2dcutsw_resample_air[i][j], twodspec_frames_dirs[i][0], trashdir, 1)
                    remove_or_move_sf(obj_sscfm_transm_2dcutsw_resample_vac[i][j], twodspec_frames_dirs[i][1], trashdir, 1)

    # 1d spectra of SKY

    if conf.flag_skyemission:
        skyemission_dirnames = ["AIR", "VAC"]

        skyemission_frames_dirs = [[["%s_NO%d/sky_emission/%s/fsr%.2f/" % (
            conf.objname_obj[i], i + 1, skyemission_dirnames[n], conf.cutrange_list[k]) for k in range(cutlength)] for n
                                    in range(2)] for i in range(conf.objnum)]

        for i in range(conf.objnum):
            for k in range(cutlength):
                os.makedirs(skyemission_frames_dirs[i][0][k])
                os.makedirs(skyemission_frames_dirs[i][1][k])
                for j in range(aplength):
                    remove_or_move_sf(sky_fm_trans_1dcutw_fsr_air[i][j][k], skyemission_frames_dirs[i][0][k],
                                   trashdir, 1)
                    remove_or_move_sf(sky_fm_trans_1dcutw_fsr_vac[i][j][k], skyemission_frames_dirs[i][1][k],
                                   trashdir, 1)

    # intermediate files of OBJ 2d images

    intermediate_obj_dirnames = ["1-OBJ_sky_subs", "2-OBJ_scatter_subs", "3-OBJ_flat", "4-OBJ_mask", "5-OBJ_cut",
                                 "6-OBJ_transform", "7-OBJ_mask2"]

    intermediate_obj_frames_dirs = [
        ["%s_NO%d/intermediate_files/OBJ/%s" % (conf.objname_obj[i], i + 1, intermediate_obj_dirnames[n]) for n in
         range(7)]
        for i in range(conf.objnum)]

    for i in range(conf.objnum):
        for n in range(7):
            os.makedirs(intermediate_obj_frames_dirs[i][n])
        remove_or_move_sf(obj_s_list[i], intermediate_obj_frames_dirs[i][0], trashdir, save)
        if conf.flag_apscatter:
            remove_or_move_sf(obj_ssc_list[i], intermediate_obj_frames_dirs[i][1], trashdir, save)
            remove_or_move_sf(obj_ssc_scatter_list[i], intermediate_obj_frames_dirs[i][1], trashdir, save)
        remove_or_move_sf(obj_sscf_list[i], intermediate_obj_frames_dirs[i][2], trashdir, save)
        remove_or_move_sf(obj_sscfm_list[i], intermediate_obj_frames_dirs[i][3], trashdir, 1)
        if conf.flag_bpmask:
            remove_or_move_sf(obj_s_mask_list[i], intermediate_obj_frames_dirs[i][3], trashdir, 1)  #
            remove_or_move_sf(obj_s_maskflat_list[i], intermediate_obj_frames_dirs[i][3], trashdir, 1)
            remove_or_move_sf(obj_s_mf1_list[i], intermediate_obj_frames_dirs[i][3], trashdir, 1)  #
            remove_or_move_sf(obj_s_mf2_list[i], intermediate_obj_frames_dirs[i][3], trashdir, 1)  #
            remove_or_move_sf(obj_s_noise_list[i], intermediate_obj_frames_dirs[i][3], trashdir, 1)  #
        for j in range(aplength):
            remove_or_move_sf(obj_sscfm_cut_list[i][j], intermediate_obj_frames_dirs[i][4], trashdir, save)
            remove_or_move_sf(obj_sscfm_trans_list[i][j], intermediate_obj_frames_dirs[i][5], trashdir, save)

    # intermediate files of OBJ 1d spectra

    intermediate_obj_1dspec_dirnames = ["1-OBJ-1DSPEC_extract", "2-OBJ-1DSPEC_truncate", "3-OBJ-1DSPEC_shift",
                                        "4-OBJ-1DSPEC_dispcor"]

    intermediate_obj_1dspec_frames_dirs = [["%s_NO%d/intermediate_files/OBJ/8A-OBJ-1DSPEC/%s" % (
        conf.objname_obj[i], i + 1, intermediate_obj_1dspec_dirnames[n]) for n in range(4)] for i in range(conf.objnum)]

    for i in range(conf.objnum):
        for n in range(4):
            os.makedirs(intermediate_obj_1dspec_frames_dirs[i][n])
        for j in range(aplength):
            remove_or_move_sf(obj_sscfm_transm_1dap[i][j], intermediate_obj_1dspec_frames_dirs[i][0], trashdir, 1)
            remove_or_move_sf(obj_sscfm_transm_1dcut[i][j], intermediate_obj_1dspec_frames_dirs[i][1], trashdir, 1)
            remove_or_move_sf(obj_sscfm_transm_1dcuts[i][j], intermediate_obj_1dspec_frames_dirs[i][2], trashdir, 1)
            remove_or_move_sf(obj_sscfm_transm_1dcutsw[i][j], intermediate_obj_1dspec_frames_dirs[i][3], trashdir, 1)
            if conf.flag_skysub:
                remove_or_move_sf(obj_sscfm_transm_1d_noneap[i][j], intermediate_obj_1dspec_frames_dirs[i][0],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1d_bg[i][j], intermediate_obj_1dspec_frames_dirs[i][0],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1d_bgcut[i][j], intermediate_obj_1dspec_frames_dirs[i][1],
                               trashdir, 1)
                remove_or_move_sf(obj_sscfm_transm_1d_bgcutw[i][j], intermediate_obj_1dspec_frames_dirs[i][3],
                               trashdir, 1)

    # intermediate files of OBJ 2d spectra

    if conf.flag_extract2d:
        intermediate_obj_2dspec_dirnames = ["1-OBJ-2DSPEC_extract", "2-OBJ-2DSPEC_truncate", "3-OBJ-2DSPEC_shift"]

        intermediate_obj_2dspec_frames_dirs = [["%s_NO%d/intermediate_files/OBJ/8B-OBJ-2DSPEC/%s" % (
            conf.objname_obj[i], i + 1, intermediate_obj_2dspec_dirnames[n]) for n in range(3)] for i in range(conf.objnum)]

        for i in range(conf.objnum):
            for n in range(3):
                os.makedirs(intermediate_obj_2dspec_frames_dirs[i][n])
            for j in range(aplength):
                remove_or_move_sf(obj_sscfm_transm_2dap[i][j], intermediate_obj_2dspec_frames_dirs[i][0], trashdir,
                               save)
                remove_or_move_sf(obj_sscfm_transm_2dcut[i][j], intermediate_obj_2dspec_frames_dirs[i][1], trashdir,
                               save)
                remove_or_move_sf(obj_sscfm_transm_2dcuts[i][j], intermediate_obj_2dspec_frames_dirs[i][2], trashdir,
                               1)
                if os.path.exists(obj_sscfm_transm_2dap_resample[i][j].ext):
                    remove_or_move_sf(obj_sscfm_transm_2dap_resample[i][j], intermediate_obj_2dspec_frames_dirs[i][0], trashdir,
                                   save)
                    remove_or_move_sf(obj_sscfm_transm_2dcut_resample[i][j], intermediate_obj_2dspec_frames_dirs[i][1], trashdir,
                                   save)
                    remove_or_move_sf(obj_sscfm_transm_2dcuts_resample[i][j], intermediate_obj_2dspec_frames_dirs[i][2], trashdir,
                                   1)

    # intermediate files of SKY 1d spectrum & 2d images

    if conf.flag_skyemission:

        intermediate_sky_dirnames = ["1-SKY_flat", "2-SKY_mask", "3-SKY_cut", "4-SKY_transform", "5-SKY_extract",
                                     "6-SKY_truncate", "7-SKY_dispcor"]

        intermediate_sky_frames_dirs = [
            ["%s_NO%d/intermediate_files/SKY/%s" % (conf.objname_obj[i], i + 1, intermediate_sky_dirnames[n]) for n in
             range(7)] for i in range(conf.objnum)]

        for i in range(conf.objnum):
            for n in range(7):
                os.makedirs(intermediate_sky_frames_dirs[i][n])
            remove_or_move_sf(sky_f_list[i], intermediate_sky_frames_dirs[i][0], trashdir, save)
            remove_or_move_sf(sky_fm_list[i], intermediate_sky_frames_dirs[i][1], trashdir, save)
            for j in range(aplength):
                remove_or_move_sf(sky_fm_cut_list[i][j], intermediate_sky_frames_dirs[i][2], trashdir, save)
                remove_or_move_sf(sky_fm_trans_list[i][j], intermediate_sky_frames_dirs[i][3], trashdir, save)
                remove_or_move_sf(sky_fm_trans_1dap[i][j], intermediate_sky_frames_dirs[i][4], trashdir, save)
                remove_or_move_sf(sky_fm_trans_1dcut[i][j], intermediate_sky_frames_dirs[i][5], trashdir, save)
                remove_or_move_sf(sky_fm_trans_1dcutw[i][j], intermediate_sky_frames_dirs[i][6], trashdir, 1)

    # Logs

    os.makedirs("reduction_log")
    if conf.flag_apscatter:
        remove_or_move(apscatter_log, "reduction_log", trashdir, 1)
    remove_or_move(cutransform_log, "reduction_log", trashdir, 1)
    if not conf.flag_manual_aperture:
        remove_or_move(centersearch_txt, "reduction_log", trashdir, 1)
        remove_or_move(centersearch_npz, "reduction_log", trashdir, 1)
    remove_or_move(aperture_txt, "reduction_log", trashdir, 1)
    remove_or_move(aperture_npz, "reduction_log", trashdir, 1)
    if conf.objnum > 1:
        remove_or_move(waveshift_txt, "reduction_log", trashdir, 1)
        remove_or_move(waveshift_npz, "reduction_log", trashdir, 1)
    if conf.flag_bpmask:
        remove_or_move(badpix_txt, "reduction_log", trashdir, 1)
        remove_or_move(badpix_npz, "reduction_log", trashdir, 1)
    remove_or_move(countfwhmpng, "reduction_log", trashdir, 1)

    # raw data

    os.makedirs("rawdata_image")
    for i in range(len(conf.imagelist)):
        remove_or_move(conf.imagelist[i] + ".fits", "rawdata_image", trashdir, 1)
        remove_or_move(conf.imagelist[i] + ".png", "rawdata_image", trashdir, 1)

    # spectra images

    os.makedirs("spectra_image")
    remove_or_move(conf.objnameRep + "_VACnorm.pdf", "spectra_image", trashdir, 1)
    remove_or_move(conf.objnameRep + "_VACflux.pdf", "spectra_image", trashdir, 1)
    for j in range(aplength):
        remove_or_move(obj_comb_norm_png[j], "spectra_image", trashdir, 1)

    # calibration data

    os.makedirs("calibration_data")
    remove_or_move(conf.flat_file, "calibration_data", trashdir, save)
    remove_or_move(conf.comp_file, "calibration_data", trashdir, 1)
    remove_or_move(conf.mask_file, "calibration_data", trashdir, save)
    remove_or_move("database", "calibration_data", trashdir, 1)
    remove_or_move("input_files.txt", "calibration_data", trashdir, 1)
    if os.path.exists("calibration_parameters.txt"):
        remove_or_move("calibration_parameters.txt", "calibration_data", trashdir, 1)
    if os.path.exists("readme.txt"):
        remove_or_move("readme.txt", "calibration_data", trashdir, 1)

    # slit viewer images

    if conf.flag_svimage:
        os.makedirs("slit_viewer")
        svlist = []
        for i in range(conf.imnum):
            if not conf.svfr_str[i] in svlist:
                remove_or_move(conf.svfr_str[i], "slit_viewer", trashdir, 1)
                svlist.append(conf.svfr_str[i])
            if not conf.svfr_end[i] in svlist:
                remove_or_move(conf.svfr_end[i], "slit_viewer", trashdir, 1)
                svlist.append(conf.svfr_end[i])
            remove_or_move(conf.imagelist[i] + "_expstart.png", "slit_viewer", trashdir, 1)
            remove_or_move(conf.imagelist[i] + "_expend.png", "slit_viewer", trashdir, 1)

    if os.path.exists("null"):
        remove_or_move("null", "reduction_log", trashdir, save)
    if os.path.exists("logfile"):
        remove_or_move("logfile", "reduction_log", trashdir, save)
    if os.path.exists("winered_logo.pdf"):
        remove_or_move("winered_logo.pdf", "reduction_log", trashdir, 1)

    shutil.rmtree(trashdir)
    elapsedTime = time.time() - startTimeSec
    endTimeStr = time.ctime()
    conf.writeStatus("reduction_log/status.txt", pipeline_ver, startTimeStr, endTimeStr, elapsedTime, "Successfully finished.")

    if not noreport:
        tex_source_maker.tex_source_make(conf, fsr, logo)

    # remove trash directory
    constant_str_length("Finished.")


if __name__ == "__main__":
    # option setting

    parser = argparse.ArgumentParser()
    parser.add_argument("listfile", type=str, help="fits file list to be reduced")
    parser.add_argument("-r", "--rawdatapath", type=str, default="../", help="directory path of input raw data")
    parser.add_argument("-v", "--viewerpath", type=str, default="../", help="directory path of slit vierer data")
    parser.add_argument("-c", "--calibpath", type=str, default="./", help="directory path of calibration data")
    parser.add_argument("-d", "--destpath", type=str, default="./", help="destination directory path")
    parser.add_argument("-q", "--query", action="store_true", help="query mode for setting pipeline parameters")
    # parser.add_argument("-m", "--mode", type=str, default="normal", help="reduction modes",
    #                     choices=["obs", "re", "normal"])
    parser.add_argument("-s", "--save", action="store_true", help="save all data")
    parser.add_argument("-p", "--parameterfile", type=str, help="pipeline parameter file")
    parser.add_argument("-o", "--oldformat", action="store_true", help="old (-ver3.6) input list format")
    parser.add_argument("-f", "--fastMode", action="store_true", help="Run WARP with the fast mode. (CR detection, wavelengthi shift are skipped.)")
    parser.add_argument("-a", "--autoCalib", action="store_true", help="Choose the appropriate calib data automatically.")
    parser.add_argument("--noreport", action="store_true", help="Not generate a pdf report.")


    # args = parser.parse_args()
    kwards = vars(parser.parse_args())

    # Warp_sci(listfile, rawdatapath, viewerpath, calibpath, destpath, flagquery, flagtrash, parameterfile,
    #          flagoldformat)
    Warp_sci(**kwards)
