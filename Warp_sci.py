#!/usr/bin/env python
# -*- coding:utf-8 -*-
# 2013.04.28 kondo version 1 作成#
# 簡易解析だけどまだ最終版ではないので注意
# (マイナーチェンジや機能追加もありうる)
# ホットピクセルは現段階でな未処理
# コンティニュームフィットは適当であくまで参考程度
# 本プログラムは独立に動くプログラムをつなぎあわせている
# 各プログラムはshellから引数を渡して独立で動かすことができる
# 将来的にはすべてpythonから処理する予定だが、エンジニアリングの段階では
# 統合化の予定はない
# 将来はなくなる機能としてゴーストライン確認も入っている
#
# ------------------------------------------------------------------------------
#
# 2015.06.19 hamano version 2.0 revised
# Added changes in this version
#   - The script of CUTRANSFORM is revised from v2.1 to v3.0
#   - The version of the pipeline is written in the reduced output files.
#   - The flat, comparison, and mask files are input from a plain text file.
#
# 2015.06.24 kondo version 2.1 minor update
#   - Intermediate files and calibration files are to be removed from output files.
#       (just for saving the HDD capacity of the server PC (merlot))
#   - If you want to take a look of the intermediate files, please use ver2.0.
# 2015.06.25 kondo version 2.2 minor update
#   - "apall_1d_v2.0.py" was updated to "apall_1d_v2.2.py"
#   - If an object is located on the edge of the slit,
#     a center finding function of the program  (v2.0) cannot work
#     and the program stops involuntarily. A bug in the function was fixed.
#     
# 2015.06.25 kondo version 2.2 minor update
#   - raw data are to be removed from output files
#
# 2015.07.06 kondo version 2.3 minor update
#   - raw data are to be removed from output files
#   - apall_1d_v2.2.py was updated to apall_1d_v2.3.py
#   - new parameter (center_position) was introduced for apall_1d_v2.3.py
#
# 2015.08.18 kondo version 2.4 minor update
#   - the programe was changed to make one cut5 figure
#
# ------------------------------------------------------------------------------
#
# 2017.01.25 Hamano version 3.0 MAJOR UPDATE!!!
#   - running other scripts by subprocess procedure are removed. the functions are imported instead.
#   - the parameters of WIDE mode are removed. all the parameters are read from calibration data.
#   - New functions:
#       * SNR estimate
#       * wavelength shift
#       * conversion from vacuum wavelength to air wavelength
#       * Bad pixel detection
#       * Air glow of night sky is extracted from sky frame.
#   - Some functions are optionized.
#       * subtraction of scattered light
#       * bad pix detection
#       * Aperture setting
#       * cutrange setting
#       * background subtraction mode setting
#
# 2017.02.01 Hamano version 3.1 updated:
#   - minor update of the vac2air_spec script
#   - comments are added
#   - Logs are added
#
# 2017.02.22 Hamano version 3.2 updated:
#   - plotframes_v1_0 is added. spectrum plotter.
#   - Until ver3.1, the aperture width was unintentionally set as 4*FWHMs, which is too large. From this ver., the width is set as 2*FWHMs.
#
# 2017.04.14 Hamano version 3.3 updated:
#   - bag fix
#   - When you do manual aperture setting in this ver., the region should be input from input file list, not from parameter file.
#   - cutransform is updated to v3.2.
#       - flux and dy can be set from calibration_parameters.txt.
#   - xshift parameter is removed to fix the bgsample region as set.
#   - function "remove_or_move" is newly installed and used instead of shutil.move.
#   - The dat and png file of spatial profiles are made from this ver.
#   - The default parameter of transform is changed from flux="yes" to flux="no".
#   - cut_v1 is updated.
#
# 2017.11.06 Hamano version 3.4 updated:
#   - apscatter is updated from 3.0 to 3.1.
#   - VACflux.pdf is drawn with plot_all_frames_flux.
#   - The parameters of PyContinuum are changed between WIDE and HIRES mode.
#   - Extraction of sky emission spectra is changed to be optional.
#   - The background regions are input from list file, not from calibration_parameters.txt.
#   - apall is updated from 1.1 to 1.2.
#   - centersearch_fortrans is updated from 2.0 to 2.1.
#   - TeX log file is made.
#   - badpixmask is updated from 1.0 to 1.1
#   - The form of waveshift log is changed.
#   - ccwaveshift_fortrans is updated from 1.3 to 2.0.
#   - plotframes are updated from 1.1 to 1.2
#   - specimage/ directory is made to contain the spectra images.
#
# 2017.12.21 Hamano version 3.5 updated:
#   - Bugs in the wavelength conversion of sum spectra are fixed
#   - tex_source_maker is updated from 1.0 to 1.2. (1.1 is not used.)
#   - The way to set calibration parameters is improved. Default-mode, query-mode, file-mode.
#   - The image of slit viewer is made.
#
# 2018.08.20 Hamano version 3.6 updated:
#   - Pipeline version is added in the fits header of raw images.
#   - Grid and info. are added in the SV image.
#   - Zscale is applied to the SV image.
#   - Aperture ranges and slit rectangle were plotted in the SV image.
#   - A plot of FWHM and average count was added.
#   - For the reading of fits image, astropy.io.fits is used instead of pyfits.
#   - Python 2.6 --> Python 3.5
#   - Free spectral ranges of HIRES-modes are replaced with more precise values.
#   - Reduce the file size of images using PIL.
#
# 2019.09.22 Hamano version 3.7.0 updated:
#   - Refactoring.
#   - the version number is removed from the name of the script.
#   - Clipping in the waveshift procedure is exported to ccwaveshift.py as the function "waveshiftClip".
#
# 2019.12.17 Hamano version 3.7.1 updated:
#   - The difference of pathlib package behavior is considered.
#
# 2020.01.13 Hamano version 3.7.2 updated:
#   - bug fixes in angle_measure and auto_ecidentify
#
# 2020.02.25 Hamano version 3.7.3 updated:
#   - improvement in auto_ecidentify
#   - The reference frame in waveshift analysis is changed to the highest flux frame.
#   - The parameters of wsmeasure, wscorrect, and wsmanual were added.
#   - The oldformat option was added.
#   - The log of waveshift is changed.
#   - New format of list file is incorpolated.
#
# 2020.03.06 Hamano version 3.7.4 updated:
#   - fixed bugs in plotframes.aperture_plot.
#   - make continuum spectra for the combined spectra.
#
# 2020.03.06 Hamano version 3.7.6 updated:
#   - bug fix
#
# 2020.06.10 Hamano version 3.7.7 updated:
#   - bug fix in SV image.
#   - bug fix in read_waveshift_log
#
# 2020.07.20 Hamano version 3.7.8 updated:
#   - New functions in aptools.
#   - Half size option in SV image (plotframes).
#   - Adjusting the figure design.
#
# 2021.11.19 Hamano version 3.7.9 updated:
#   - The letter '#' in the object name is replaced with '_'.
#   - ArrayMonitor, which produces dead pixel maps and a hot pixel map from flat data, is newly added.
#
# 2022.09.30 Hamano version 3.8.0 updated:
#   - Bug fix (Warp_calib: The dpmask files are removed in the aperture-replace mode.)
#   - Refactoring
#   - New function to detect cosmic rays was added.
#   - Tex section to show the cosmic ray statistics was added.
#   - Magellan telescope parameters were added.
#   - Parameters in auto_ecidentify were tuned.
#
# 2022.10.11 Hamano version 3.8.1 updated:
#   - Bug fix (ECtoID.py: Error case when IRAF/ecidentify output integer parameters with float type values.
#   - Bug fix (plotframes.plot_2dimage_sv: float conversion of SVSEEING value was added.)
#
# 2023.01.20 Hamano version 3.8.2 updated:
#   - Bug fix (plotframes.plot_2dimage_sv: offset parameter was added )
#   - config class: WARP parameter class
#   - warpLog class: WARP result log class
#   - apertureSet class: aperture information and methods class
#

__version__ = "3.8.2"

from pyraf import iraf
import sys, shutil, os, glob, time
import numpy as np
import pathlib
from astropy.io import fits
import argparse

sys.path.append(os.path.dirname(__file__))

from config import config
from logger import warpLog
from aperture import apertureSet
from centersearch_fortrans import centersearch_fortrans, make_slit_profile
from Spec2Dtools import flatfielding, header_key_read
from apscatter import pyapscatter
from cutransform import cutransform
from Spec1Dtools import pyapall, truncate, dispcor_single, cut_1dspec, PyScombine, FSR_angstrom, openspecfits
from ccwaveshift import waveshift_multiorder, waveshift_oneorder, PySpecshift, waveshiftClip
from SNratio_estimate import snestimate
from PyContinuum import PyContinuum
from vac2air_spec import vac2air_spec
from plotframes import plot_all_frames_norm, plot_all_frames_flux, plot_all_frames_flux_BG, plot_2dimages_mask, \
    plot_2dimages, snr_plots, plot_combined_norm, plot_2dimages_sv, peak_count_fwhm, aperture_plot, cosmicRay2dImages
from badpixmask import badpixmask, pyfixpix, cosmicRayMask
import tex_source_maker


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


def absPathStr(path):
    pathobj = pathlib.Path(path)
    return str(pathobj.resolve()) + "/"


def Warp_sci(listfile, rawdatapath, viewerpath, calibpath, destpath, flagquery, reductionmode, flagsave, parameterfile,
             flagoldformat):
    pipeline_ver = __version__
    conf = config()

    # check the paths and make destination directory

    currentpath = absPathStr("./")

    if parameterfile is not None:
        parampath = pathlib.Path(parameterfile)
        paramfile = str(parampath.resolve())

    if not os.path.exists(rawdatapath):
        print("Error: " + rawdatapath + "does not exist.")
        sys.exit()
    if not os.path.exists(viewerpath):
        print("Error: " + viewerpath + "does not exist.")
        sys.exit()
    if not os.path.exists(calibpath):
        print("Error: " + calibpath + "does not exist.")
        sys.exit()

    rawdatapath = absPathStr(rawdatapath)
    viewerpath = absPathStr(viewerpath)
    calibpath = absPathStr(calibpath)

    if destpath != currentpath:
        if not os.path.exists(destpath):
            try:
                os.makedirs(destpath)
            except:
                print("Error: tried to make \"%s\" directory, but failed. ")
                sys.exit()
        else:
            print("Error: the destination directory already exists.")
            sys.exit()

    destpath = absPathStr(destpath)
    shutil.copy(listfile, destpath)
    listfile = os.path.basename(listfile)
    os.chdir(destpath)

    if os.path.exists(calibpath + "input_files.txt") and os.path.exists(calibpath + "database"):
        calibfiles = glob.glob(calibpath + "*")
        for i in calibfiles:
            if not os.path.exists(destpath + i.split("/")[-1]):
                if os.path.isfile(i):
                    shutil.copy(i, destpath)
                elif os.path.isdir(i):
                    shutil.copytree(i, destpath + i.split("/")[-1])
    else:
        print("Error: Invalid calibration path (%s)." % calibpath)
        sys.exit()

    # read calibration data names and pipeline parameters and modes.

    conf.readInputCalib("input_files.txt")

    if flagquery:
        conf.readParamQuery()
    elif parameterfile is not None:
        conf.readParamFile(paramfile)

    # convert_flag_yesno = {1: "yes", 0: "no"}
    tfDict = {True: "yes", False: "no"}

    startTimeSec = time.time()
    startTimeStr = time.ctime()

    # open input file list

    rfile = open(listfile, "r")
    rlines = rfile.readlines()
    rfile.close()

    # read the file names of object frame and corresponding sky frame.

    objectlist = []
    skylist = []
    lowlim_input = []
    upplim_input = []
    skysub_region = []
    waveshift_man = []

    if flagoldformat:
        for i in range(len(rlines)):
            rl1 = rlines[i].split()
            for j in rl1:
                if j.find("ws=") != -1:
                    waveshift_man.append(float(j.lstrip("ws=")))
                    rl1.remove(j)
                    break
                if j == rl1[-1]:
                    waveshift_man.append(0.)
            objectlist.append(rl1[0])
            skylist.append(rl1[1])
            if conf.flag_manual_aperture:
                lowlim_input.append(float(rl1[2]))
                upplim_input.append(float(rl1[3]))
            if conf.flag_skysub:
                skysub_region.append(rl1[-1])
            else:
                skysub_region.append("INDEF")
    else:
        for i in range(len(rlines)):
            rl1 = rlines[i].split()
            objectlist.append(rl1[0])
            skylist.append(rl1[1])
            flagwsinput = False
            flagbginput = False
            for j in rl1[2:]:
                if j.find("ap=") != -1:
                    aptmp = j.lstrip("ap=")
                    lowlim_input.append(float(aptmp.split(":")[0]))
                    upplim_input.append(float(aptmp.split(":")[1]))
                if j.find("bg=") != -1:
                    skysub_region.append(j.lstrip("bg="))
                    flagbginput = True
                if j.find("ws=") != -1:
                    waveshift_man.append(float(j.lstrip("ws=")))
                    flagwsinput = True
            if not flagwsinput:
                waveshift_man.append(0.)
            if not flagbginput:
                skysub_region.append("INDEF")

    if len(objectlist) != len(lowlim_input) and conf.flag_manual_aperture:
        print("Aperture range parameter is not written in input list.")
        sys.exit()
    if len(objectlist) != len(skysub_region) and conf.flag_skysub:
        print("Background region parameter is not written in input list.")
        sys.exit()

    # synthesize the object frame list and sky frame list without overlaps.

    imagelist = list(set(objectlist + skylist))
    imlist_sort = sorted(imagelist)
    objnum = len(objectlist)
    imnum = len(imlist_sort)
    imobjid = []
    for i in range(imnum):
        if imlist_sort[i] in objectlist:
            imobjid.append(objectlist.index(imlist_sort[i]))
        else:
            imobjid.append("sky")

    for i in range(imnum):
        shutil.copy("%s%s.fits" % (rawdatapath, imagelist[i]), ".")
        iraf.hedit("%s" % imagelist[i], "PIPELINE", pipeline_ver, add="yes", verify="no")

    # read fits headers

    objname = []
    nodpos = []
    satupix = []
    svfr_str = []
    svfr_end = []
    for i in range(imnum):
        hdulist_obj = fits.open(imlist_sort[i] + ".fits")

        prihdr_obj = hdulist_obj[0].header
        data_obj = hdulist_obj[0].data

        hdulist_obj.close()

        data_obj[data_obj <= conf.saturation_thres] = 0
        data_obj[data_obj > conf.saturation_thres] = 1
        satupix.append(np.sum(data_obj))
        objname.append(
            header_key_read(prihdr_obj, "object").replace(" ", "_").replace(" ", "_").replace("'", "_").replace("\"",
                                                                                                                "_").replace(
                '#', '_'))

        nodpos.append(header_key_read(prihdr_obj, "NODPOS"))
        svfr_str.append(header_key_read(prihdr_obj, "SVFR-STR") + ".fits")
        svfr_end.append(header_key_read(prihdr_obj, "SVFR-END") + ".fits")

    flag_svimage = True
    for i in range(imnum):
        if flag_svimage:
            if svfr_str[i].find("N/A") != -1 and svfr_end[i].find("N/A") != -1:
                flag_svimage = False
            if not os.path.exists(viewerpath + svfr_str[i]):
                flag_svimage = False
            if not os.path.exists(viewerpath + svfr_end[i]):
                flag_svimage = False

    if flag_svimage:
        for i in range(imnum):
            shutil.copy(viewerpath + svfr_str[i], ".")
            shutil.copy(viewerpath + svfr_end[i], ".")

    objname_obj = []
    nodpos_obj = []
    for i in range(objnum):
        for j in range(imnum):
            if objectlist[i] == imlist_sort[j]:
                objname_obj.append(objname[j])
                nodpos_obj.append(nodpos[j])

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

    obj_s_list = [objname_obj[i] + "_NO{}_s".format(i + 1) for i in
                  range(objnum)]  # FITS file after sky subtraction: "star_NO1_s.fits"
    obj_s_mask_list = ["mask_" + obj_s_list[i] for i in range(objnum)]
    obj_s_maskfig_list = ["mask_" + obj_s_list[i] + ".pdf" for i in range(objnum)]
    obj_s_maskflat_list = ["maskflat_" + obj_s_list[i] for i in range(objnum)]
    obj_s_mf1_list = ["mf1_" + obj_s_list[i] for i in range(objnum)]
    obj_s_mf2_list = ["mf2_" + obj_s_list[i] for i in range(objnum)]
    bpthres_list = [0. for i in range(objnum)]
    bpnum_list = [0 for i in range(objnum)]
    if conf.flag_apscatter:
        obj_ssc_list = [obj_s_list[i] + "sc" for i in
                        range(objnum)]  # FITS file after scattered light subtraction: "star_NO1_ssc.fits"
    else:
        obj_ssc_list = obj_s_list
    obj_sscf_list = [obj_ssc_list[i] + "f" for i in
                     range(objnum)]  # FITS file after flat fielding: "star*_NO1_s(sc)f.fits"
    obj_sscfm_list = [obj_sscf_list[i] + "m" for i in range(objnum)]  # FITS file after fixpix: "star_NO1_s(sc)fm.fits"
    obj_sscfm_trans_list = [[obj_sscfm_list[i] + "_m{}trans".format(m) for m in apset.echelleOrders] for i in
                            range(objnum)]  # FITS file after transform: "star_NO1_s(sc)fm_**trans.fits"

    img_cs_list = [[obj_sscfm_list[i] + "_m{}trans.png".format(m) for m in apset.echelleOrders] for i in range(objnum)]
    dat_cs_list = [[obj_sscfm_list[i] + "_m{}trans.dat".format(m) for m in apset.echelleOrders] for i in range(objnum)]

    # make the output file names of sky frame for...
    # flat fielding
    # bad pixel interpolation
    # transformation
    #

    sky_f_list = [objname_obj[i] + "_skyNO{}".format(i + 1) + "_f" for i in
                  range(objnum)]  # sky FITS file after flat fielding: "HD***_skyNO1_f.fits"
    sky_fm_list = [sky_f_list[i] + "m" for i in range(objnum)]  # sky FITS file after fixpix: "HD***_skyNO1_fm.fits"
    sky_fm_trans_list = [[sky_fm_list[i] + "_m{}trans".format(m) for m in apset.echelleOrders] for i in
                         range(objnum)]  # sky FITS file after transform: "HD***_skyNO1_fm_**trans.fits"

    # reduction from sky subtraction to cutransform

    apscatter_log = "apscatter_log.txt"
    cutransform_log = "cutransform_log.txt"
    centersearch_txt, centersearch_npz = "centersearch_log.txt", "centersearch_log.npz"
    aperture_txt, aperture_npz = "aperture_log.txt", "aperture_log.npz"
    waveshift_txt, waveshift_npz = "waveshift_log.txt", "waveshift_log.npz"
    badpix_txt, badpix_npz = "cosmicray_log.txt", "cosmicray_log.npz"

    log = warpLog(apset.echelleOrders, objnum)

    if True:  # meaningless. just for debuging.
        for i in range(objnum):
            iraf.imarith(objectlist[i], "-", skylist[i], obj_s_list[i])  # sky subtraction
            if conf.flag_bpmask:
                print("Cosmic ray detection for {}...".format(obj_s_list[i]))
                if nodpos_obj[i].find("O") == -1:
                    bpnum_list[i], bpthres_list[i] = \
                        cosmicRayMask(obj_s_list[i], objectlist[i], skylist[i], obj_s_mask_list[i], obj_s_mf1_list[i],
                                      obj_s_mf2_list[i], conf.ap_file, conf.mask_file, True, threshold=conf.CRthreshold,
                                      varatio=conf.CRvaratio,
                                      slitposratio=conf.CRslitposratio, maxsigma=conf.CRmaxsigma,
                                      fixsigma=conf.CRfixsigma)
                else:
                    bpnum_list[i], bpthres_list[i] = \
                        cosmicRayMask(obj_s_list[i], objectlist[i], skylist[i], obj_s_mask_list[i], obj_s_mf1_list[i],
                                      obj_s_mf2_list[i], conf.ap_file, conf.mask_file, False,
                                      threshold=conf.CRthreshold,
                                      varatio=conf.CRvaratio, slitposratio=conf.CRslitposratio,
                                      maxsigma=conf.CRmaxsigma,
                                      fixsigma=conf.CRfixsigma)
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

            cutransform(obj_sscfm_list[i], conf.ap_file, cutransform_log, conf.dyinput,
                        conf.fluxinput)  # (obj) transformation

            if conf.flag_skyemission:
                flatfielding(skylist[i], sky_f_list[i], conf.flat_file)  # (sky) flat fielding
                if conf.flag_bpmask:
                    pyfixpix(sky_f_list[i], sky_fm_list[i], obj_s_maskflat_list[i])  # (sky) bad pixel interpolation
                else:
                    pyfixpix(sky_f_list[i], sky_fm_list[i], conf.mask_file)  # (sky) bad pixel interpolation
                cutransform(sky_fm_list[i], conf.ap_file, cutransform_log, conf.dyinput,
                            conf.fluxinput)  # (sky) transformation

        if conf.flag_bpmask:
            log.cosmicRayLog(bpthres_list, bpnum_list)
            log.writeCosmicRayLogText(badpix_txt)
            log.writeCosmicRayLogNpz(badpix_npz)
            # cosmicraymask_log(badpix_log, bpthres_list, bpnum_list, objnum)

    # center search and bad pix detection

    obj_sscfm_transm_list = obj_sscfm_trans_list

    if True:  # meaningless. just for debuging.
        if not conf.flag_manual_aperture:
            xs = [[] for i in range(objnum)]
            gw = [[] for i in range(objnum)]
            for i in range(objnum):
                for j in range(aplength):
                    if nodpos_obj[i].find("O") == -1 and nodpos_obj[i] != "NA":
                        tmpx, tmpg = centersearch_fortrans(obj_sscfm_trans_list[i][j], aptranslist[j],
                                                           dat_cs_list[i][j],
                                                           abbaflag=True)
                    else:
                        tmpx, tmpg = centersearch_fortrans(obj_sscfm_trans_list[i][j], aptranslist[j],
                                                           dat_cs_list[i][j], abbaflag=False)
                    xs[i].append(tmpx)
                    gw[i].append(tmpg)

            log.psfLog(xs, gw)
            log.writePsfLogText(centersearch_txt)
            log.writePsfLogNpz(centersearch_npz)
            # make_centersearch_log(centersearch_log, apnum, objnum, xs, gw)

            # setting the aperture range as 2 sigma.
            lowtrans = [[-1. * np.median(gw[i]) + np.median(xs[i]) for j in range(aplength)] for i in
                        range(objnum)]
            hightrans = [[np.median(gw[i]) + np.median(xs[i]) for j in range(aplength)] for i in range(objnum)]
            for i in range(objnum):
                for j in range(aplength):
                    aperture_plot(dat_cs_list[i][j], img_cs_list[i][j], lowtrans[i][j], hightrans[i][j],
                                  skysub_region[i], conf.flag_skysub)

        else:
            # setting aperture as input from input file list
            lowtrans = [[lowlim_input[i] for j in range(aplength)] for i in range(objnum)]
            hightrans = [[upplim_input[i] for j in range(aplength)] for i in range(objnum)]

            for i in range(objnum):
                for j in range(aplength):
                    make_slit_profile(obj_sscfm_trans_list[i][j], aptranslist[j], dat_cs_list[i][j])
                    aperture_plot(dat_cs_list[i][j], img_cs_list[i][j], lowtrans[i][j], hightrans[i][j],
                                  skysub_region[i], conf.flag_skysub)

    # setting background subtraction option as input from calibration parameter file

    bgsample = [[skysub_region[i]] for i in range(objnum)]
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

    obj_sscfm_transm_1d = [[objname_obj[i] + "_NO{}_m{}".format(i + 1, m) for m in apset.echelleOrders] for i in
                           range(objnum)]
    obj_sscfm_transm_1d_none = [[objname_obj[i] + "_NO{}_m{}_none".format(i + 1, m) for m in apset.echelleOrders] for i
                                in range(objnum)]
    obj_sscfm_transm_1d_bg = [[objname_obj[i] + "_NO{}_m{}_bg".format(i + 1, m) for m in apset.echelleOrders] for i in
                              range(objnum)]
    obj_sscfm_transm_1d_bgcut = [[objname_obj[i] + "_NO{}_m{}_bgc".format(i + 1, m) for m in apset.echelleOrders] for i
                                 in range(objnum)]
    obj_sscfm_transm_1d_bgcutw = [[objname_obj[i] + "_NO{}_m{}_bgcw".format(i + 1, m) for m in apset.echelleOrders] for
                                  i
                                  in range(objnum)]
    obj_sscfm_transm_1dap = [[] for i in range(objnum)]
    obj_sscfm_transm_1d_noneap = [[] for i in range(objnum)]
    obj_sscfm_transm_1dcut = [[obj_sscfm_transm_1d[i][j] + "c" for j in range(aplength)] for i in range(objnum)]
    file_matrix_forwaveshift = [[obj_sscfm_transm_1d[i][j] + "c" for i in range(objnum)] for j in range(aplength)]

    obj_sscfm_transm_1dcuts = [[obj_sscfm_transm_1d[i][j] + "cs" for j in range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw = [[obj_sscfm_transm_1d[i][j] + "csw" for j in range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac = [
        [[obj_sscfm_transm_1d[i][j] + "_fsr%.2f_VAC" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac_norm = [
        [[obj_sscfm_transm_1d[i][j] + "_fsr%.2f_VAC_norm" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac_cont = [
        [[obj_sscfm_transm_1d[i][j] + "_fsr%.2f_VAC_cont" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw_fsr_air = [
        [[obj_sscfm_transm_1d[i][j] + "_fsr%.2f_AIR" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw_fsr_air_norm = [
        [[obj_sscfm_transm_1d[i][j] + "_fsr%.2f_AIR_norm" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw_fsr_air_cont = [
        [[obj_sscfm_transm_1d[i][j] + "_fsr%.2f_AIR_cont" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]

    sky_fm_trans_1d = [[sky_fm_trans_list[i][j] + "1d" for j in range(aplength)] for i in range(objnum)]
    sky_fm_trans_1dap = [[] for i in range(objnum)]
    sky_fm_trans_1dcut = [[sky_fm_trans_list[i][j] + "1dcut" for j in range(aplength)] for i in range(objnum)]
    sky_fm_trans_1dcutw = [[sky_fm_trans_list[i][j] + "1dcutw" for j in range(aplength)] for i in range(objnum)]
    sky_fm_trans_1dcutw_fsr_vac = [
        [[sky_fm_trans_list[i][j] + "1dcutw_fsr%.2f_VAC" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]
    sky_fm_trans_1dcutw_fsr_air = [
        [[sky_fm_trans_list[i][j] + "1dcutw_fsr%.2f_AIR" % conf.cutrange_list[k] for k in range(cutlength)] for j in
         range(aplength)] for i in range(objnum)]

    # extracted 2d spectrum (only for OBJ)
    # -> cut the edge
    # -> shift
    # -> apply the dispersion solution
    # -> VAC or AIR
    #

    obj_sscfm_transm_2d = [[obj_sscfm_transm_list[i][j] + "2d" for j in range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_2dap = [[] for i in range(objnum)]
    obj_sscfm_transm_2dcut = [[obj_sscfm_transm_list[i][j] + "2dcut" for j in range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_2dcuts = [[obj_sscfm_transm_list[i][j] + "2dcuts" for j in range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_2dcutsw_vac = [[obj_sscfm_transm_list[i][j] + "2dcutsw_VAC" for j in range(aplength)] for i in
                                    range(objnum)]
    obj_sscfm_transm_2dcutsw_air = [[obj_sscfm_transm_list[i][j] + "2dcutsw_AIR" for j in range(aplength)] for i in
                                    range(objnum)]

    if True:  # meaningless. just for debuging.

        aplow_log = [[] for i in range(objnum)]
        aphigh_log = [[] for i in range(objnum)]

        for i in range(objnum):
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
                obj_sscfm_transm_1dap[i].append(obj_sscfm_transm_1d[i][j] + "." + apname_trans)

                # extract 2d spectrum (OBJ)
                pyapall(obj_sscfm_transm_list[i][j], obj_sscfm_transm_2d[i][j], obj_sscfm_trans_list[i][j],
                        conf.skysub_mode, "strip")
                truncate(obj_sscfm_transm_2d[i][j] + "." + apname_trans, obj_sscfm_transm_2dcut[i][j])
                obj_sscfm_transm_2dap[i].append(obj_sscfm_transm_2d[i][j] + "." + apname_trans)

                # extract 1d spectrum (SKY)
                if conf.flag_skyemission:
                    pyapall(sky_fm_trans_list[i][j], sky_fm_trans_1d[i][j], obj_sscfm_trans_list[i][j], "none",
                            "onedspec")
                    truncate(sky_fm_trans_1d[i][j] + "." + apname_trans, sky_fm_trans_1dcut[i][j])
                    sky_fm_trans_1dap[i].append(sky_fm_trans_1d[i][j] + "." + apname_trans)

        log.apertureLog(aplow_log, aphigh_log)
        log.writeApertureLogText(aperture_txt)
        log.writeApertureLogNpz(aperture_npz)

        # measuring the shift from the first frame
        if objnum > 1 and conf.flag_wsmeasure:
            counts = [0. for i in range(objnum)]
            for i in range(objnum):
                for j in range(aplength):
                    spx, spy, _, _, _ = openspecfits(obj_sscfm_transm_1dcut[i][j])
                    counts[i] += np.median(spy)
            refid = np.argmax(counts) if objnum > 2 else 0

            wshift_matrix = []
            for j in range(aplength):
                wshift_matrix.append(waveshift_oneorder(file_matrix_forwaveshift[j], refid))

            shift_average, shift_calcnum, shift_stddev, wshift_flag = waveshiftClip(wshift_matrix, objnum, aplength)
            if conf.flag_wsmanual:
                shift_cor = waveshift_man
            else:
                shift_cor = [shift_average[i] if wshift_flag[i] == 0 else 0. for i in range(objnum)]
        else:
            wshift_matrix = [[0. for i in range(objnum)] for j in range(aplength)]
            shift_average = [0. for i in waveshift_man]
            shift_stddev = [0. for i in waveshift_man]
            shift_calcnum = [0 for i in waveshift_man]
            if conf.flag_wsmanual:
                shift_cor = waveshift_man
            else:
                shift_cor = [0. for i in waveshift_man]

        log.waveshiftLog(wshift_matrix, shift_average, shift_stddev, shift_calcnum, shift_cor)
        log.writeWaveshiftLogText(waveshift_txt)
        log.writeWaveshiftLogNpz(waveshift_npz)

        for i in range(objnum):
            for j in range(aplength):

                # apply the dispersion solution to 1d spectra (OBJ)
                if objnum > 1 and conf.flag_wscorrect:
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
                               conf.cutrange_list[k], apset.echelleOrders[j])
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
                if objnum > 1:
                    PySpecshift(obj_sscfm_transm_2dcut[i][j], obj_sscfm_transm_2dcuts[i][j], shift_average[i])
                else:
                    iraf.scopy(obj_sscfm_transm_2dcut[i][j], obj_sscfm_transm_2dcuts[i][j])
                dispcor_single(obj_sscfm_transm_2dcuts[i][j], obj_sscfm_transm_2dcutsw_vac[i][j], comp_file_id[j])
                iraf.hedit(obj_sscfm_transm_2dcutsw_vac[i][j], "AIRORVAC", "vac", add="yes", verify="no")
                vac2air_spec(obj_sscfm_transm_2dcutsw_vac[i][j], obj_sscfm_transm_2dcutsw_air[i][j])

                if conf.flag_skyemission:
                    # apply dispersion solution, cut, convert to air wavelength for 1d spectra (SKY)
                    dispcor_single(sky_fm_trans_1dcut[i][j], sky_fm_trans_1dcutw[i][j], comp_file_id[j])
                    iraf.hedit(sky_fm_trans_1dcutw[i][j], "AIRORVAC", "vac", add="yes", verify="no")

                    for k in range(cutlength):
                        cut_1dspec(sky_fm_trans_1dcutw[i][j], sky_fm_trans_1dcutw_fsr_vac[i][j][k],
                                   conf.cutrange_list[k], apset.echelleOrders[j])
                        vac2air_spec(sky_fm_trans_1dcutw_fsr_vac[i][j][k], sky_fm_trans_1dcutw_fsr_air[i][j][k])

    # measuring signal-to-noize ratio, combine, normalize, conversion to air wavelength

    SNmatrix_fsr = [[[obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k] for i in range(objnum)] for j in range(aplength)] for k
                    in range(cutlength)]
    SNoutput_fsr = [["SNratio_m%d_fsr%.2f.dat" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for k in
                    range(cutlength)]
    SN_png = ["SNratio_fsr%.2f.png" % conf.cutrange_list[k] for k in range(cutlength)]

    combined_spec_fsr_vac = [
        [objname_obj[0] + "_sum_m%d_fsr%.2f_VAC" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for k in
        range(cutlength)]
    combined_spec_fsr_vac_norm = [
        [objname_obj[0] + "_sum_m%d_fsr%.2f_VAC_norm" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for
        k in
        range(cutlength)]
    combined_spec_fsr_vac_cont = [
        [objname_obj[0] + "_sum_m%d_fsr%.2f_VAC_cont" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for
        k in
        range(cutlength)]
    combined_spec_fsr_air = [
        [objname_obj[0] + "_sum_m%d_fsr%.2f_AIR" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for k in
        range(cutlength)]
    combined_spec_fsr_air_norm = [
        [objname_obj[0] + "_sum_m%d_fsr%.2f_AIR_norm" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for
        k in
        range(cutlength)]
    combined_spec_fsr_air_cont = [
        [objname_obj[0] + "_sum_m%d_fsr%.2f_AIR_cont" % (apset.echelleOrders[j], conf.cutrange_list[k]) for j in range(aplength)] for
        k in
        range(cutlength)]

    if objnum > 1:
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

            PyScombine(combined_spec_fsr_vac_norm[k],
                       objname_obj[0] + "_sum_fsr%.2f_VAC_norm_combine" % (conf.cutrange_list[k]))
            PyScombine(combined_spec_fsr_air_norm[k],
                       objname_obj[0] + "_sum_fsr%.2f_AIR_norm_combine" % (conf.cutrange_list[k]))

    else:
        for k in range(cutlength):
            for j in range(aplength):
                iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_air[0][j][k], combined_spec_fsr_air[k][j])
                iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_air_norm[0][j][k], combined_spec_fsr_air_norm[k][j])
                iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_air_cont[0][j][k], combined_spec_fsr_air_cont[k][j])
                iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_vac[0][j][k], combined_spec_fsr_vac[k][j])
                iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_vac_norm[0][j][k], combined_spec_fsr_vac_norm[k][j])
                iraf.scopy(obj_sscfm_transm_1dcutsw_fsr_vac_cont[0][j][k], combined_spec_fsr_vac_cont[k][j])

    # make plots

    obj_sscfm_transm_1dcutsw_fsr_vac_norm_105 = [
        [obj_sscfm_transm_1dcutsw_fsr_vac_norm[i][j][0] for j in range(aplength)] for i in range(objnum)]
    obj_sscfm_transm_1dcutsw_fsr_vac_flux_105 = [[obj_sscfm_transm_1dcutsw_fsr_vac[i][j][0] for j in range(aplength)]
                                                 for i in range(objnum)]
    obj_comb_norm_png = [objname_obj[0] + "_AIRnorm_combined_m%d.png" % apset.echelleOrders[j] for j in range(aplength)]

    plot_all_frames_norm(obj_sscfm_transm_1dcutsw_fsr_vac_norm_105, objname_obj[0] + "_VACnorm.pdf", apset.echelleOrders, objnum)
    if conf.flag_skysub:
        plot_all_frames_flux_BG(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, obj_sscfm_transm_1d_bgcutw,
                                objname_obj[0] + "_VACflux.pdf", apset.echelleOrders, objnum)
    else:
        plot_all_frames_flux(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, objname_obj[0] + "_VACflux.pdf", apset.echelleOrders, objnum)

    plot_combined_norm(combined_spec_fsr_air_norm[0], obj_comb_norm_png, apset.echelleOrders)

    if conf.flag_bpmask:
        for i in range(objnum):
            plot_2dimages_mask(obj_s_mask_list[i] + ".fits", obj_s_mask_list[i] + ".png")

    for i in range(imnum):
        plot_2dimages(imagelist[i] + ".fits", imagelist[i] + ".png")
    for i in range(objnum):
        plot_2dimages(obj_sscfm_list[i] + ".fits", obj_sscfm_list[i] + ".png")

    if objnum > 1:
        for k in range(cutlength):
            snr_plots(lams_sn[k], snr_val[k], combined_spec_fsr_vac[k], aplength, SN_png[k])

    countfwhmpng = "count_and_fwhm.png"
    if not conf.flag_manual_aperture:
        fwhm = [np.median(log.psfWidth) for i in range(objnum)]
        peak_count_fwhm(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, countfwhmpng, apset.echelleOrders, objnum, fwhm=fwhm)
    else:
        peak_count_fwhm(obj_sscfm_transm_1dcutsw_fsr_vac_flux_105, countfwhmpng, apset.echelleOrders, objnum)

    if flag_svimage:
        for i in range(imnum):
            if imobjid[i] != "sky":
                plot_2dimages_sv(svfr_str[i], imlist_sort[i] + "_expstart.png",
                                 lowlim=np.average(lowtrans[imobjid[i]]),
                                 upplim=np.average(hightrans[imobjid[i]]))
                plot_2dimages_sv(svfr_end[i], imlist_sort[i] + "_expend.png",
                                 lowlim=np.average(lowtrans[imobjid[i]]),
                                 upplim=np.average(hightrans[imobjid[i]]))
            else:
                plot_2dimages_sv(svfr_str[i], imlist_sort[i] + "_expstart.png")
                plot_2dimages_sv(svfr_end[i], imlist_sort[i] + "_expend.png")

    # make directries, and move files to appropriate directories
    # trash directory

    trashdir = "Trash"
    os.makedirs(trashdir)

    # 1d spectra of OBJ

    onedspec_dirnames = ["AIR_flux", "AIR_norm", "AIR_cont", "VAC_flux", "VAC_norm", "VAC_cont"]

    onedspec_frames_dirs = [[["%s_NO%d/onedspec/%s/fsr%.2f/" % (
        objname_obj[i], (i + 1), onedspec_dirnames[n], conf.cutrange_list[k]) for k in range(cutlength)] for n in
                             range(6)]
                            for i
                            in range(objnum)]
    onedspec_sum_dirs = [
        ["%s_sum/%s/fsr%.2f/" % (objname_obj[0], onedspec_dirnames[n], conf.cutrange_list[k]) for k in range(cutlength)]
        for
        n in range(6)]

    for n in range(6):
        for k in range(cutlength):
            for i in range(objnum):
                os.makedirs(onedspec_frames_dirs[i][n][k])
            os.makedirs(onedspec_sum_dirs[n][k])

    for j in range(aplength):
        for k in range(cutlength):
            for i in range(objnum):
                remove_or_move(obj_sscfm_transm_1dcutsw_fsr_air[i][j][k] + ".fits", onedspec_frames_dirs[i][0][k],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1dcutsw_fsr_air_norm[i][j][k] + ".fits", onedspec_frames_dirs[i][1][k],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1dcutsw_fsr_air_cont[i][j][k] + ".fits", onedspec_frames_dirs[i][2][k],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1dcutsw_fsr_vac[i][j][k] + ".fits", onedspec_frames_dirs[i][3][k],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1dcutsw_fsr_vac_norm[i][j][k] + ".fits", onedspec_frames_dirs[i][4][k],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1dcutsw_fsr_vac_cont[i][j][k] + ".fits", onedspec_frames_dirs[i][5][k],
                               trashdir, 1)

            remove_or_move(combined_spec_fsr_air[k][j] + ".fits", onedspec_sum_dirs[0][k], trashdir, 1)
            remove_or_move(combined_spec_fsr_air_norm[k][j] + ".fits", onedspec_sum_dirs[1][k], trashdir, 1)
            remove_or_move(combined_spec_fsr_air_cont[k][j] + ".fits", onedspec_sum_dirs[2][k], trashdir, 1)
            remove_or_move(combined_spec_fsr_vac[k][j] + ".fits", onedspec_sum_dirs[3][k], trashdir, 1)
            remove_or_move(combined_spec_fsr_vac_norm[k][j] + ".fits", onedspec_sum_dirs[4][k], trashdir, 1)
            remove_or_move(combined_spec_fsr_vac_cont[k][j] + ".fits", onedspec_sum_dirs[5][k], trashdir, 1)

    if objnum > 1:
        for k in range(cutlength):
            remove_or_move(objname_obj[0] + "_sum_fsr%.2f_AIR_norm_combine.fits" % (conf.cutrange_list[k]),
                           onedspec_sum_dirs[1][k], trashdir, 1)
            remove_or_move(objname_obj[0] + "_sum_fsr%.2f_VAC_norm_combine.fits" % (conf.cutrange_list[k]),
                           onedspec_sum_dirs[3][k], trashdir, 1)

    # signal-to-noize ratio

    SNRdat_frames_dirs = ["SNR_dat/fsr%.2f" % (conf.cutrange_list[k]) for k in range(cutlength)]

    if objnum > 1:
        for k in range(cutlength):
            os.makedirs(SNRdat_frames_dirs[k])
            for j in range(aplength):
                remove_or_move("SNratio_m%d_fsr%.2f.dat" % (apset.echelleOrders[j], conf.cutrange_list[k]), SNRdat_frames_dirs[k],
                               trashdir, 1)
            remove_or_move(SN_png[k], SNRdat_frames_dirs[k], trashdir, 1)

    # images

    images_frames_dirs_sp = ["%s_NO%d/images/spatial_profile/" % (objname_obj[i], (i + 1)) for i in range(objnum)]

    for i in range(objnum):
        os.makedirs(images_frames_dirs_sp[i])
        for j in range(aplength):
            remove_or_move(img_cs_list[i][j], images_frames_dirs_sp[i], trashdir, 1)
            remove_or_move(dat_cs_list[i][j], images_frames_dirs_sp[i], trashdir, 1)

    images_frames_dirs_2d = ["%s_NO%d/images/2d_image/" % (objname_obj[i], (i + 1)) for i in range(objnum)]

    for i in range(objnum):
        os.makedirs(images_frames_dirs_2d[i])
        remove_or_move(obj_sscfm_list[i] + ".png", images_frames_dirs_2d[i], trashdir, 1)

    if conf.flag_bpmask:
        images_frames_dirs_bp = ["%s_NO%d/images/badpixmask/" % (objname_obj[i], (i + 1)) for i in range(objnum)]

        for i in range(objnum):
            os.makedirs(images_frames_dirs_bp[i])
            remove_or_move(obj_s_mask_list[i] + ".png", images_frames_dirs_bp[i], trashdir, 1)
            remove_or_move(obj_s_maskfig_list[i], images_frames_dirs_bp[i], trashdir, 1)

    # 2d spectra of OBJ

    twodspec_dirnames = ["AIR", "VAC"]

    twodspec_frames_dirs = [["%s_NO%d/twodspec/%s/" % (objname_obj[i], (i + 1), twodspec_dirnames[n]) for n in range(2)]
                            for i in range(objnum)]

    for i in range(objnum):
        os.makedirs(twodspec_frames_dirs[i][0])
        os.makedirs(twodspec_frames_dirs[i][1])
        for j in range(aplength):
            remove_or_move(obj_sscfm_transm_2dcutsw_air[i][j] + ".fits", twodspec_frames_dirs[i][0], trashdir, 1)
            remove_or_move(obj_sscfm_transm_2dcutsw_vac[i][j] + ".fits", twodspec_frames_dirs[i][1], trashdir, 1)

    # 1d spectra of SKY

    if conf.flag_skyemission:
        skyemission_dirnames = ["AIR", "VAC"]

        skyemission_frames_dirs = [[["%s_NO%d/sky_emission/%s/fsr%.2f/" % (
            objname_obj[i], i + 1, skyemission_dirnames[n], conf.cutrange_list[k]) for k in range(cutlength)] for n in
                                    range(2)]
                                   for i in range(objnum)]

        for i in range(objnum):
            for k in range(cutlength):
                os.makedirs(skyemission_frames_dirs[i][0][k])
                os.makedirs(skyemission_frames_dirs[i][1][k])
                for j in range(aplength):
                    remove_or_move(sky_fm_trans_1dcutw_fsr_air[i][j][k] + ".fits", skyemission_frames_dirs[i][0][k],
                                   trashdir, 1)
                    remove_or_move(sky_fm_trans_1dcutw_fsr_vac[i][j][k] + ".fits", skyemission_frames_dirs[i][1][k],
                                   trashdir, 1)

    # intermediate files of OBJ 2d images

    intermediate_obj_dirnames = ["1-OBJ_sky_subs", "2-OBJ_scatter_subs", "3-OBJ_flat", "4-OBJ_mask", "5-OBJ_cut",
                                 "6-OBJ_transform", "7-OBJ_mask2"]

    intermediate_obj_frames_dirs = [
        ["%s_NO%d/intermediate_files/OBJ/%s" % (objname_obj[i], i + 1, intermediate_obj_dirnames[n]) for n in range(7)]
        for i in range(objnum)]

    for i in range(objnum):
        for n in range(7):
            os.makedirs(intermediate_obj_frames_dirs[i][n])
        remove_or_move(obj_s_list[i] + ".fits", intermediate_obj_frames_dirs[i][0], trashdir, flagsave)
        if conf.flag_apscatter:
            remove_or_move(obj_ssc_list[i] + ".fits", intermediate_obj_frames_dirs[i][1], trashdir, flagsave)
            remove_or_move("scatter_%s" % (obj_s_list[i] + ".fits"), intermediate_obj_frames_dirs[i][1], trashdir,
                           flagsave)
        remove_or_move(obj_sscf_list[i] + ".fits", intermediate_obj_frames_dirs[i][2], trashdir, flagsave)
        remove_or_move(obj_sscfm_list[i] + ".fits", intermediate_obj_frames_dirs[i][3], trashdir, 1)
        if conf.flag_bpmask:
            remove_or_move(obj_s_mask_list[i] + ".fits", intermediate_obj_frames_dirs[i][3], trashdir, 1)  #
            remove_or_move(obj_s_maskflat_list[i] + ".fits", intermediate_obj_frames_dirs[i][3], trashdir, 1)
            remove_or_move(obj_s_mf1_list[i] + ".fits", intermediate_obj_frames_dirs[i][3], trashdir, 1)  #
            remove_or_move(obj_s_mf2_list[i] + ".fits", intermediate_obj_frames_dirs[i][3], trashdir, 1)  #
        for j in range(aplength):
            remove_or_move(obj_sscfm_trans_list[i][j].replace("trans", "cut") + ".fits",
                           intermediate_obj_frames_dirs[i][4], trashdir, flagsave)
            remove_or_move(obj_sscfm_trans_list[i][j] + ".fits", intermediate_obj_frames_dirs[i][5], trashdir,
                           flagsave)

    # intermediate files of OBJ 1d spectra

    intermediate_obj_1dspec_dirnames = ["1-OBJ-1DSPEC_extract", "2-OBJ-1DSPEC_truncate", "3-OBJ-1DSPEC_shift",
                                        "4-OBJ-1DSPEC_dispcor"]

    intermediate_obj_1dspec_frames_dirs = [["%s_NO%d/intermediate_files/OBJ/8A-OBJ-1DSPEC/%s" % (
        objname_obj[i], i + 1, intermediate_obj_1dspec_dirnames[n]) for n in range(4)] for i in range(objnum)]

    for i in range(objnum):
        for n in range(4):
            os.makedirs(intermediate_obj_1dspec_frames_dirs[i][n])
        for j in range(aplength):
            remove_or_move(obj_sscfm_transm_1dap[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][0], trashdir,
                           1)
            remove_or_move(obj_sscfm_transm_1dcut[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][1], trashdir,
                           1)
            remove_or_move(obj_sscfm_transm_1dcuts[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][2], trashdir,
                           1)
            remove_or_move(obj_sscfm_transm_1dcutsw[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][3],
                           trashdir, 1)
            if conf.flag_skysub:
                remove_or_move(obj_sscfm_transm_1d_noneap[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][0],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1d_bg[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][0],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1d_bgcut[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][1],
                               trashdir, 1)
                remove_or_move(obj_sscfm_transm_1d_bgcutw[i][j] + ".fits", intermediate_obj_1dspec_frames_dirs[i][3],
                               trashdir, 1)

    # intermediate files of OBJ 2d spectra

    intermediate_obj_2dspec_dirnames = ["1-OBJ-2DSPEC_extract", "2-OBJ-2DSPEC_truncate", "3-OBJ-2DSPEC_shift"]

    intermediate_obj_2dspec_frames_dirs = [["%s_NO%d/intermediate_files/OBJ/8B-OBJ-2DSPEC/%s" % (
        objname_obj[i], i + 1, intermediate_obj_2dspec_dirnames[n]) for n in range(3)] for i in range(objnum)]

    for i in range(objnum):
        for n in range(3):
            os.makedirs(intermediate_obj_2dspec_frames_dirs[i][n])
        for j in range(aplength):
            remove_or_move(obj_sscfm_transm_2dap[i][j] + ".fits", intermediate_obj_2dspec_frames_dirs[i][0], trashdir,
                           flagsave)
            remove_or_move(obj_sscfm_transm_2dcut[i][j] + ".fits", intermediate_obj_2dspec_frames_dirs[i][1], trashdir,
                           flagsave)
            remove_or_move(obj_sscfm_transm_2dcuts[i][j] + ".fits", intermediate_obj_2dspec_frames_dirs[i][2], trashdir,
                           1)

    # intermediate files of SKY 1d spectrum & 2d images

    if conf.flag_skyemission:

        intermediate_sky_dirnames = ["1-SKY_flat", "2-SKY_mask", "3-SKY_cut", "4-SKY_transform", "5-SKY_extract",
                                     "6-SKY_truncate", "7-SKY_dispcor"]

        intermediate_sky_frames_dirs = [
            ["%s_NO%d/intermediate_files/SKY/%s" % (objname_obj[i], i + 1, intermediate_sky_dirnames[n]) for n in
             range(7)] for i in range(objnum)]

        for i in range(objnum):
            for n in range(7):
                os.makedirs(intermediate_sky_frames_dirs[i][n])
            remove_or_move(sky_f_list[i] + ".fits", intermediate_sky_frames_dirs[i][0], trashdir, flagsave)
            remove_or_move(sky_fm_list[i] + ".fits", intermediate_sky_frames_dirs[i][1], trashdir, flagsave)
            for j in range(aplength):
                remove_or_move(sky_fm_trans_list[i][j].replace("trans", "cut") + ".fits",
                               intermediate_sky_frames_dirs[i][2], trashdir, flagsave)
                remove_or_move(sky_fm_trans_list[i][j] + ".fits", intermediate_sky_frames_dirs[i][3], trashdir,
                               flagsave)
                remove_or_move(sky_fm_trans_1dap[i][j] + ".fits", intermediate_sky_frames_dirs[i][4], trashdir,
                               flagsave)
                remove_or_move(sky_fm_trans_1dcut[i][j] + ".fits", intermediate_sky_frames_dirs[i][5], trashdir,
                               flagsave)
                remove_or_move(sky_fm_trans_1dcutw[i][j] + ".fits", intermediate_sky_frames_dirs[i][6], trashdir, 1)

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
    if objnum > 1:
        remove_or_move(waveshift_txt, "reduction_log", trashdir, 1)
        remove_or_move(waveshift_npz, "reduction_log", trashdir, 1)
    if conf.flag_bpmask:
        remove_or_move(badpix_txt, "reduction_log", trashdir, 1)
        remove_or_move(badpix_npz, "reduction_log", trashdir, 1)
    remove_or_move(countfwhmpng, "reduction_log", trashdir, 1)

    # raw data

    os.makedirs("rawdata_image")
    for i in range(len(imagelist)):
        remove_or_move(imagelist[i] + ".fits", "rawdata_image", trashdir, 1)
        remove_or_move(imagelist[i] + ".png", "rawdata_image", trashdir, 1)

    # spectra images

    os.makedirs("spectra_image")
    remove_or_move(objname_obj[0] + "_VACnorm.pdf", "spectra_image", trashdir, 1)
    remove_or_move(objname_obj[0] + "_VACflux.pdf", "spectra_image", trashdir, 1)
    for j in range(aplength):
        remove_or_move(obj_comb_norm_png[j], "spectra_image", trashdir, 1)

    # calibration data

    os.makedirs("calibration_data")
    remove_or_move(conf.flat_file, "calibration_data", trashdir, flagsave)
    remove_or_move(conf.comp_file, "calibration_data", trashdir, 1)
    remove_or_move(conf.mask_file, "calibration_data", trashdir, flagsave)
    remove_or_move("database", "calibration_data", trashdir, 1)
    remove_or_move("input_files.txt", "calibration_data", trashdir, 1)
    if os.path.exists("calibration_parameters.txt"):
        remove_or_move("calibration_parameters.txt", "calibration_data", trashdir, 1)
    if os.path.exists("readme.txt"):
        remove_or_move("readme.txt", "calibration_data", trashdir, 1)

    # slit viewer images

    if flag_svimage:
        os.makedirs("slit_viewer")
        svlist = []
        for i in range(imnum):
            if not svfr_str[i] in svlist:
                remove_or_move(svfr_str[i], "slit_viewer", trashdir, 1)
                svlist.append(svfr_str[i])
            if not svfr_end[i] in svlist:
                remove_or_move(svfr_end[i], "slit_viewer", trashdir, 1)
                svlist.append(svfr_end[i])
            remove_or_move(imlist_sort[i] + "_expstart.png", "slit_viewer", trashdir, 1)
            remove_or_move(imlist_sort[i] + "_expend.png", "slit_viewer", trashdir, 1)

    if os.path.exists("null"):
        remove_or_move("null", "reduction_log", trashdir, flagsave)
    if os.path.exists("logfile"):
        remove_or_move("logfile", "reduction_log", trashdir, flagsave)
    if os.path.exists("winered_logo.eps"):
        remove_or_move("winered_logo.eps", "reduction_log", trashdir, 1)

    shutil.rmtree(trashdir)
    elapsedTime = time.time() - startTimeSec
    endTimeStr = time.ctime()
    conf.writeStatus("reduction_log/status.txt", pipeline_ver, startTimeStr, endTimeStr, elapsedTime)

    tex_source_maker.tex_source_make(imlist_sort, objectlist, skylist)

    # remove trash directory



if __name__ == "__main__":
    # option setting

    parser = argparse.ArgumentParser()
    parser.add_argument("listfile", type=str, help="fits file list to be reduced")
    parser.add_argument("-r", "--rawdatapath", type=str, default="../", help="directory path of input raw data")
    parser.add_argument("-v", "--viewerpath", type=str, default="../", help="directory path of slit vierer data")
    parser.add_argument("-c", "--calibpath", type=str, default="./", help="directory path of calibration data")
    parser.add_argument("-d", "--destpath", type=str, default="./", help="destination directory path")
    parser.add_argument("-q", "--query", action="store_true", help="query mode for setting pipeline parameters")
    parser.add_argument("-m", "--mode", type=str, default="normal", help="reduction modes",
                        choices=["obs", "re", "normal"])
    parser.add_argument("-s", "--save", action="store_true", help="save all data")
    parser.add_argument("-p", "--parameter", type=str, help="pipeline parameter file")
    parser.add_argument("-o", "--oldformat", action="store_true", help="old (-ver3.6) input list format")

    args = parser.parse_args()

    listfile = args.listfile
    rawdatapath = args.rawdatapath
    viewerpath = args.viewerpath
    calibpath = args.calibpath
    destpath = args.destpath
    flagquery = args.query
    reductionmode = args.mode
    flagtrash = args.save
    parameterfile = args.parameter
    flagoldformat = args.oldformat

    Warp_sci(listfile, rawdatapath, viewerpath, calibpath, destpath, flagquery, reductionmode, flagtrash, parameterfile,
             flagoldformat)
