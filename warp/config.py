#!/usr/bin/env python
# -*- coding:utf-8 -*-
import os.path

from astropy.io import fits
import sys
import numpy as np
from warp.Spec2Dtools import header_key_read
import traceback

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

def constant_str_length(comment):
    constlength = 72
    print("\033[31m\n=== %s %s\n\033[0m" % (comment, ("=" * (constlength - len(comment) - 4))))


class config:
    def __init__(self, statusFile='INDEF'):
        self.saturation_thres = 35000.
        self.flag_apscatter = True
        self.flag_manual_aperture = False
        self.flag_skysub = False
        self.flag_bpmask = True
        self.flag_skyemission = False
        self.flag_wsmeasure = True
        self.flag_wscorrect = True
        self.flag_wsmanual = False
        self.flag_extract2d = False
        self.skysub_mode = "none"
        self.skysubModeList = ["none", "average", "median", "minimum", "fit"]
        self.cutrange_list = [1.05, 1.30]
        self.fluxinput = "no"
        self.CRthreshold = 10.
        self.CRvaratio = 2.
        self.CRslitposratio = 1.5
        self.CRmaxsigma = 20.
        self.CRfixsigma = False
        self.status = statusFile
        self.frameNumberLimit = 28
        self.reduceFullData = True
        self.selectedOrders = []

    def inputDataList(self, listfile, oldFormat=False):
        # open input file list

        rfile = open(listfile, "r")
        rlines = rfile.readlines()
        rfile.close()

        # read the file names of object frame and corresponding sky frame.

        self.objectlist = []
        self.skylist = []
        self.lowlim_input = []
        self.upplim_input = []
        self.skysub_region = []
        self.waveshift_man = []

        if oldFormat:
            for i in range(len(rlines)):
                rl1 = rlines[i].split()
                for j in rl1:
                    if j.find("ws=") != -1:
                        self.waveshift_man.append(float(j.lstrip("ws=")))
                        rl1.remove(j)
                        break
                    if j == rl1[-1]:
                        self.waveshift_man.append(0.)
                self.objectlist.append(rl1[0])
                self.skylist.append(rl1[1])
                if self.flag_manual_aperture:
                    self.lowlim_input.append(float(rl1[2]))
                    self.upplim_input.append(float(rl1[3]))
                if self.flag_skysub:
                    self.skysub_region.append(rl1[-1])
                else:
                    self.skysub_region.append("INDEF")
        else:
            for i in range(len(rlines)):
                rl1 = rlines[i].split()
                self.objectlist.append(rl1[0])
                self.skylist.append(rl1[1])
                flagwsinput = False
                flagbginput = False
                for j in rl1[2:]:
                    if j.find("ap=") != -1:
                        aptmp = j.lstrip("ap=")
                        self.lowlim_input.append(float(aptmp.split(":")[0]))
                        self.upplim_input.append(float(aptmp.split(":")[1]))
                    if j.find("bg=") != -1:
                        self.skysub_region.append(j.lstrip("bg="))
                        flagbginput = True
                    if j.find("ws=") != -1:
                        self.waveshift_man.append(float(j.lstrip("ws=")))
                        flagwsinput = True
                if not flagwsinput:
                    self.waveshift_man.append(0.)
                if not flagbginput:
                    self.skysub_region.append("INDEF")

        if len(self.objectlist) != len(self.lowlim_input) and self.flag_manual_aperture:
            print("Aperture range parameter is not written in input list.")
            sys.exit()
        if len(self.objectlist) != len(self.skysub_region) and self.flag_skysub:
            print("Background region parameter is not written in input list.")
            sys.exit()

        # synthesize the object frame list and sky frame list without overlaps.

        self.imagelist = list(set(self.objectlist + self.skylist))
        self.imagelist.sort()
        self.imagelist = np.array(self.imagelist)
        self.objnum = len(self.objectlist)
        self.imnum = len(self.imagelist)
        self.imobjid = []
        for i in range(self.imnum):
            if self.imagelist[i] in self.objectlist:
                self.imobjid.append(self.objectlist.index(self.imagelist[i]))
            else:
                self.imobjid.append("sky")

    def readInputDataHeader(self):
        self.objname = np.array([])
        self.nodpos = np.array([])
        self.satupix = np.array([])
        self.svfr_str = np.array([])
        self.svfr_end = np.array([])
        self.period = np.array([])
        self.setting = np.array([])
        self.slit = np.array([])
        self.instmode = np.array([])
        for i in range(self.imnum):
            hdulist_obj = fits.open(self.imagelist[i] + ".fits")
            prihdr_obj = hdulist_obj[0].header
            data_obj = hdulist_obj[0].data
            hdulist_obj.close()
            data_obj[data_obj <= self.saturation_thres] = 0
            data_obj[data_obj > self.saturation_thres] = 1
            self.satupix = np.append(self.satupix, np.sum(data_obj))
            self.objname = np.append(self.objname,
                header_key_read(prihdr_obj, "OBJECT").replace(" ", "_").replace(" ", "_").replace("'", "_").replace(
                    "\"", "_").replace('#', '_').replace('/', '_'))
            self.nodpos = np.append(self.nodpos, header_key_read(prihdr_obj, "NODPOS"))
            self.svfr_str = np.append(self.svfr_str, header_key_read(prihdr_obj, "SVFR-STR") + ".fits")
            self.svfr_end = np.append(self.svfr_end, header_key_read(prihdr_obj, "SVFR-END") + ".fits")
            self.period = np.append(self.period, header_key_read(prihdr_obj, "PERIOD"))
            self.setting = np.append(self.setting, header_key_read(prihdr_obj, "SETTING"))
            self.slit = np.append(self.slit, header_key_read(prihdr_obj, "SLIT"))
            self.instmode = np.append(self.instmode, header_key_read(prihdr_obj, "INSTMODE"))
        self.objnameRep = self.objname[self.imagelist == self.objectlist[0]][0]

        self.flag_svimage = True
        for i in range(self.imnum):
            if self.svfr_str[i].find("N/A") != -1 and self.svfr_end[i].find("N/A") != -1:
                self.flag_svimage = False

        self.objname_obj = [self.objname[self.imagelist == i][0] for i in self.objectlist]
        self.nodpos_obj = [self.nodpos[self.imagelist == i][0] for i in self.objectlist]

    def readInputCalib(self, inputlist):
        para = []
        calpath = os.path.dirname(inputlist)
        if calpath != "":
            calpath += "/"
        for line in open(inputlist, "r"):
            items = line.split()
            para.append(str(items[0]))
        self.flat_file = calpath + para[1]
        self.mask_file = calpath + para[2]
        self.comp_file = calpath + para[3]
        self.ap_file = calpath + para[4]
        self.aptrans_file = calpath + para[5]
        self.apsc_maskfile = calpath + para[6]

        if os.path.exists(self.flat_file):
            flatf = fits.open(self.flat_file)
            prihdr_flat = flatf[0].header
            flatf.close()
            self.flatSetting = header_key_read(prihdr_flat, "SETTING")
            self.flatPeriod = header_key_read(prihdr_flat, "PERIOD")
            self.flatSlit = header_key_read(prihdr_flat, "SLIT")
            self.flatMode = header_key_read(prihdr_flat, "INSTMODE")

        if os.path.exists(self.comp_file):
            compf = fits.open(self.comp_file)
            prihdr_comp = compf[0].header
            compf.close()
            self.dyinput = prihdr_comp["CDELT1"]
            self.compSetting = header_key_read(prihdr_comp, "SETTING")
            self.compPeriod = header_key_read(prihdr_comp, "PERIOD")
            self.compSlit = header_key_read(prihdr_comp, "SLIT")
            self.compMode = header_key_read(prihdr_comp, "INSTMODE")

        return None

    def readParamQuery(self):
        print("======================================================")
        print("===    Please answer to the following questions.   ===")
        print("===                       OR                       ===")
        print("=== Just press enter key (adopt default settings). ===")
        print("======================================================")

        ynDict = {"yes": True, "no": False}
        tfDict = {True: "yes", False: "no"}

        # self.flagapscatter = ynDict[
        #     alternativequestion("Subtract scattered light? (def:{}) :".format(tfDict[self.flag_apscatter]),
        #                         ["yes", "no", ""], tfDict[self.flag_apscatter])]
        self.flag_manual_aperture = ynDict[
            alternativequestion("Adopt aperture ranges read from the input file? (def:{}) :".format(tfDict[self.flag_manual_aperture]),
                                ["yes", "no", ""], tfDict[self.flag_manual_aperture])]
        # self.skysub_mode = alternativequestion(
        #     "Subtract background spectra from object spectra? (def:{}) :".format(self.skysub_mode),
        #     self.skysubModeList + [""], self.skysub_mode)
        self.flag_skysub = ynDict[
            alternativequestion("Subtract background spectra from object spectra? (def:{}) :".format(tfDict[self.flag_skysub]),
            ["yes", "no", ""], tfDict[self.flag_skysub])]
        self.skysub_mode = "average" if self.flag_skysub else "none"
        self.flag_bpmask = ynDict[
            alternativequestion("Detect and interpolate the cosmic rays? (def:{}) :".format(tfDict[self.flag_bpmask]),
                                ["yes", "no", ""], tfDict[self.flag_bpmask])]

        # ans_query_CRparams = alternativequestion(
        #     "Change any parameters in the cosmic ray detection algorithm? (def:no) :",
        #     ["yes", "no", "detail", ""], "no")
        # if ans_query_CRparams == "yes":
        #     self.CRthreshold = valueQuestion(
        #         "Threshold for the cosmic ray detection (def: {} sigma) :".format(self.CRthreshold), 3., 100.,
        #         self.CRthreshold)
        #     self.CRfixsigma = ynDict[
        #         alternativequestion("Fix the threshold sigma (def: {}) :".format(tfDict[self.CRfixsigma]),
        #                             ["yes", "no", ""], tfDict[self.CRfixsigma])]
        # elif ans_query_CRparams == "detail":
        #     self.CRthreshold = valueQuestion(
        #         "Threshold for the cosmic ray detection (def: {} sigma) :".format(self.CRthreshold), 3., 100.,
        #         self.CRthreshold)
        #     self.CRfixsigma = ynDict[
        #         alternativequestion("Fix the threshold sigma (def: {}) :".format(tfDict[self.CRfixsigma]),
        #                             ["yes", "no", ""], tfDict[self.CRfixsigma])]
        #     if not self.CRfixsigma:
        #         self.CRmaxsigma = valueQuestion(
        #             "Maximum threshold for the cosmic ray detection (def: {}) :".format(self.CRmaxsigma), 3., 100.,
        #             self.CRmaxsigma)
        #         self.CRvaratio = valueQuestion(
        #             "Threshold for the variance / average of the cosmic ray distribution (def: {}) :".format(
        #                 self.CRvaratio), 1., 100., self.CRvaratio)
        #         self.CRslitposratio = valueQuestion(
        #             "Threshold for the cosmic ray number ratio between the slit positions (def: {}) :".format(
        #                 self.CRslitposratio), 1., 100., self.CRslitposratio)

        # self.flag_skyemission = ynDict[
        #     alternativequestion("Extract the spectra from sky frame? (def: {}) :".format(tfDict[self.flag_skyemission]),
        #                         ["yes", "no", ""], tfDict[self.flag_skyemission])]
        self.flag_wsmeasure = ynDict[
            alternativequestion("Measure the spectra offsets among multiple frames? (def: {}) :".format(tfDict[self.flag_wsmeasure]),
                                ["yes", "no", ""], tfDict[self.flag_wsmeasure])]
        self.flag_wscorrect = ynDict[
            alternativequestion("Correct the spectra offsets among multiple frames? (def: {}) :".format(tfDict[self.flag_wscorrect]),
                                ["yes", "no", ""], tfDict[self.flag_wscorrect])]
        if self.flag_wsmeasure and self.flag_wscorrect:
            self.flag_wsmanual = ynDict[alternativequestion(
                "Use the spectra offsets values written in list file? (def: {}) :".format(tfDict[self.flag_wsmanual]),
                ["yes", "no", ""], tfDict[self.flag_wsmanual])]
        elif not self.flag_wsmeasure and self.flag_wscorrect:
            self.flag_wsmanual = True
        else:
            self.flag_wsmanual = False

        # self.fluxinput = alternativequestion(
        #     "Conserve the flux in the transformation? (def: {}) :".format(self.fluxinput),
        #     ["yes", "no", ""], self.fluxinput)
        # self.flag_skysub = True if self.skysub_mode != "none" else False
        self.flag_extract2d = ynDict[
            alternativequestion("Extract two-dimensional spec (wavelength and spatial dimensions)? (def: {}) :".format(tfDict[self.flag_extract2d]),
                                ["yes", "no", ""], tfDict[self.flag_extract2d])]

        self.showAllParams()

        return None

    def readParamFile(self, inputlist):
        # read calibration_parameters.txt, in which the option and parameters for the pipeline is stored.

        ynDict = {"yes": True, "no": False}
        tfDict = {True: "yes", False: "no"}

        for line in open(inputlist, "r"):
            if line.find("Apscatter") != -1:
                self.flag_apscatter = ynDict[line.split(":")[1].split()[0]]
            if line.find("Manual Aperture") != -1:
                self.flag_manual_aperture = ynDict[line.split(":")[1].split()[0]]
            if line.find("Cosmic Ray Correction") != -1:
                self.flag_bpmask = ynDict[line.split(":")[1].split()[0]]
            if line.find("Extract all orders") != -1:
                self.reduceFullData = ynDict[line.split(":")[1].split()[0]]
            if line.find("Selected orders") != -1:
                if line.split(":")[1].split()[0] != "no":
                    tmpline = line.split(":")[1].split(",")
                    self.selectedOrders = [int(tmpline[i].split()[0]) for i in range(len(tmpline))]
            if line.find("Background Subtraction") != -1:
                if line.split(":")[1].split()[0] != "none":
                    self.flag_skysub = True
                    self.skysub_mode = line.split(":")[1].split()[0]
            if line.find("Set cut range") != -1:
                if line.split(":")[1].split()[0] != "no":
                    tmpline = line.split(":")[1].split(",")
                    self.cutrange_list = [float(tmpline[i].split()[0]) for i in range(len(tmpline))]
            if line.find("CUTRANSFORM flux") != -1:
                self.fluxinput = line.split(":")[1].split()[0]
            if line.find("Sky Emission") != -1:
                self.flag_skyemission = ynDict[line.split(":")[1].split()[0]]
            if line.find("Measure Shift") != -1:
                self.flag_wsmeasure = ynDict[line.split(":")[1].split()[0]]
            if line.find("Correct Shift") != -1:
                self.flag_wscorrect = ynDict[line.split(":")[1].split()[0]]
            if line.find("Manual Shift") != -1:
                self.flag_wsmanual = ynDict[line.split(":")[1].split()[0]]
            if line.find("Cosmic ray threshold sigma") != -1:
                try:
                    self.CRthreshold = float(line.split(":")[1].split()[0])
                except:
                    print("\033[31m WARNING: The input value could not be converted to float. value={} was set. \033[0m".format(
                        self.CRthreshold))
            if line.find("Cosmic ray maximum sigma") != -1:
                try:
                    self.CRmaxsigma = float(line.split(":")[1].split()[0])
                except:
                    print("\033[31m WARNING: The input value could not be converted to float. value={} was set. \033[0m".format(self.CRmaxsigma))
            if line.find("Cosmic ray Var/Ave ratio") != -1:
                try:
                    self.CRvaratio = float(line.split(":")[1].split()[0])
                except:
                    print("\033[31m WARNING: The input value could not be converted to float. value={} was set. \033[0m".format(self.CRvaratio))
            if line.find("Cosmic ray ratio between slit positions") != -1:
                try:
                    self.CRslitposratio = float(line.split(":")[1].split()[0])
                except:
                    print("\033[31m WARNING: The input value could not be converted to float. value={} was set. \033[0m".format(
                        self.CRslitposratio))
            if line.find("Cosmic ray fix sigma") != -1:
                self.CRfixsigma = ynDict[line.split(":")[1].split()[0]]
            if line.find("Extract 2d spectrum") != -1:
                self.flag_extract2d = ynDict[line.split(":")[1].split()[0]]

        self.showAllParams()

        return None

    def setFastModeParam(self):
        self.flag_apscatter = True
        self.flag_manual_aperture = False
        self.flag_skysub = False
        self.flag_bpmask = False
        self.flag_skyemission = False
        self.flag_wsmeasure = False
        self.flag_wscorrect = False
        self.flag_wsmanual = False
        self.flag_extract2d = False
        self.skysub_mode = "none"
        self.cutrange_list = [1.05]
        self.fluxinput = "no"
        self.CRthreshold = 10.
        self.CRvaratio = 2.
        self.CRslitposratio = 1.5
        self.CRmaxsigma = 20.
        self.CRfixsigma = False
        self.reduceFullData = True
        self.showAllParams()

    def writeStatus(self, pipelineVer, startTimeStr):
        status_file = open(self.status, "a")
        status_file.write("WARP ver.%s\n" % (pipelineVer))
        status_file.write("Starting time: %s\n" % startTimeStr)
        status_file.write("\n")
        status_file.close()

    def writeElapsedTime(self, endTimeStr, elapsedTime, status):
        status_file = open(self.status, "a")
        status_file.write("Termination time: %s\n" % (endTimeStr))
        status_file.write("Elapsed time: %dm%.1fs\n\n" % (int(elapsedTime / 60), elapsedTime % 60))
        status_file.write("Pipeline status: {}".format(status))
        status_file.close()

    def writePipelineSettings(self):
        status_file = open(self.status, "a")
        tfDict = {True: "yes", False: "no"}

        status_file.write("\nSettings:\n")
        status_file.write("     Apscatter: %s\n" % tfDict[self.flag_apscatter])
        status_file.write("     Manual Aperture: %s\n" % tfDict[self.flag_manual_aperture])
        status_file.write("     Background Subtraction: %s\n" % self.skysub_mode)
        status_file.write("     Extract all orders: %s\n" % tfDict[self.reduceFullData])
        if self.reduceFullData:
            status_file.write("     Selected orders: ")
            for m in self.selectedOrders:
                status_file.write(str(m))
                if m != self.cutrange_list[-1]:
                    status_file.write(", ")
        status_file.write("     Cosmic Ray Correction: %s\n" % tfDict[self.flag_bpmask])
        status_file.write("     Cosmic ray threshold sigma: %.1f\n" % self.CRthreshold)
        status_file.write("     Cosmic ray maximum sigma: %.1f\n" % self.CRmaxsigma)
        status_file.write("     Cosmic ray Var/Ave ratio: %.1f\n" % self.CRvaratio)
        status_file.write("     Cosmic ray ratio between slit positions: %.1f\n" % self.CRslitposratio)
        status_file.write("     Cosmic ray fix sigma: %s\n" % tfDict[self.CRfixsigma])
        status_file.write("     Set cut range: ")
        for i in self.cutrange_list:
            status_file.write(str(i))
            if i != self.cutrange_list[-1]:
                status_file.write(", ")
        status_file.write("\n")
        status_file.write("     Sky Emission: %s\n" % tfDict[self.flag_skyemission])
        status_file.write("     Measure Shift: %s\n" % tfDict[self.flag_wsmeasure])
        status_file.write("     Correct Shift: %s\n" % tfDict[self.flag_wscorrect])
        status_file.write("     Manual Shift: %s\n" % tfDict[self.flag_wsmanual])
        status_file.write("     CUTRANSFORM flux: %s\n" % self.fluxinput)
        status_file.write("     Extract 2d spectrum: %s\n" % tfDict[self.flag_extract2d])
        status_file.write("\n")
        status_file.close()

    def writeInputDataList(self):
        status_file = open(self.status, "a")
        status_file.write("Input data:\n")
        for i in range(len(self.objectlist)):
            dataLine = "    {} {}".format(self.objectlist[i], self.skylist[i])
            if len(self.lowlim_input) == len(self.objectlist) and self.flag_manual_aperture:
                dataLine += " ap={}:{}".format(self.lowlim_input[i], self.upplim_input[i])
            if len(self.skysub_region) == len(self.objectlist) and self.flag_skysub:
                dataLine += " bg={}".format(self.skysub_region[i])
            if len(self.waveshift_man) == len(self.objectlist) and self.flag_wsmanual:
                dataLine += " ws={}".format(self.waveshift_man[i])
            dataLine += "\n"
            status_file.write(dataLine)
        status_file.write("\n")
        status_file.close()

    def writeError(self, error):
        with open(self.status, "a") as f:
            f.write('{}\n'.format(error))
            traceback.print_exc(file=f)
            print('\033[31m {} \033[0m\n'.format(error))
            print(traceback.format_exc())
        print("\033[31m !!!Please send the \'status.txt\' to the developer to report the error.!!! \033[0m")

    def writeParam(self):
        status_file = open(self.status, "a")
        status_file.write("WARP Settings and Parameters:\n")
        status_file.write("    flag_apscatter: {}\n".format(self.flag_apscatter))
        status_file.write("    flag_manual_aperture: {}\n".format(self.flag_manual_aperture))
        status_file.write("    flag_skysub: {}\n".format(self.flag_skysub))
        status_file.write("    flag_bpmask: {}\n".format(self.flag_bpmask))
        status_file.write("    flag_skyemission: {}\n".format(self.flag_skyemission))
        status_file.write("    flag_wsmeasure: {}\n".format(self.flag_wsmeasure))
        status_file.write("    flag_wscorrect: {}\n".format(self.flag_wscorrect))
        status_file.write("    flag_wsmanual: {}\n".format(self.flag_wsmanual))
        status_file.write("    flag_extract2d: {}\n".format(self.flag_extract2d))
        status_file.write("    skysub_mode: {}\n".format(self.skysub_mode))
        status_file.write("    skysubModeList: {}\n".format(self.skysubModeList))
        status_file.write("    cutrange_list: {}\n".format(self.cutrange_list))
        status_file.write("    fluxinput: {}\n".format(self.fluxinput))
        status_file.write("    CRthreshold: {}\n".format(self.CRthreshold))
        status_file.write("    CRvaratio: {}\n".format(self.CRvaratio))
        status_file.write("    CRslitposratio: {}\n".format(self.CRslitposratio))
        status_file.write("    CRmaxsigma: {}\n".format(self.CRmaxsigma))
        status_file.write("    CRfixsigma: {}\n".format(self.CRfixsigma))
        status_file.write("    reduceFullData: {}\n".format(self.reduceFullData))
        status_file.write("    selectedOrders: {}\n".format(self.selectedOrders))
        status_file.write("\n")
        status_file.close()


    def readObservationInfo(self, hdulist):
        self.acqtime = [header_key_read(i, "ACQTIME1").split("-")[-1] for i in hdulist]
        self.acqdate = [header_key_read(hdulist[i], "ACQTIME1").split()[0].rstrip(self.acqtime[i]).rstrip("-") for i in
                   range(self.imnum)]
        self.exptime = [header_key_read(i, "EXPTIME") for i in hdulist]
        self.inttime = [header_key_read(i, "INTTIME") for i in hdulist]
        self.ra_hours = [header_key_read(i, "RA") for i in hdulist]
        self.dec_degree = [header_key_read(i, "DEC") for i in hdulist]
        self.modes = [header_key_read(i, "INSTMODE") for i in hdulist]
        self.teles = [header_key_read(i, "TELESCOP") for i in hdulist]
        self.seeing = [header_key_read(i, "SEEING") for i in hdulist]
        self.period = [header_key_read(i, "PERIOD") for i in hdulist]
        self.setting = [header_key_read(i, "SETTING") for i in hdulist]
        self.airmass = [header_key_read(i, "AIRMASS") for i in hdulist]
        self.airmass_start = [header_key_read(i, "AIRM-STR") for i in hdulist]
        self.airmass_end = [header_key_read(i, "AIRM-END") for i in hdulist]
        self.ut_start = [header_key_read(i, "UT-STR") for i in hdulist]
        self.ut_end = [header_key_read(i, "UT-END") for i in hdulist]
        self.observatory = [header_key_read(i, "OBSERVAT") for i in hdulist]
        self.observer = [header_key_read(i, "OBSERVER") for i in hdulist]
        self.humidity = [header_key_read(i, "OUT-HUM") for i in hdulist]
        self.temperature = [header_key_read(i, "OUT-TMP") for i in hdulist]
        self.air_pressure = [header_key_read(i, "OUT-PRS") for i in hdulist]
        self.wind_speed = [header_key_read(i, "OUT-WND") for i in hdulist]
        self.wodbtheme = [header_key_read(i, "WODBTHEM").replace("_", " ") for i in hdulist]
        self.wodbobsid = [header_key_read(i, "WODBOBS") for i in hdulist]
        self.wodbstd = [header_key_read(i, "WODBSTD") for i in hdulist]
        self.wodbpi = [header_key_read(i, "WODBPI") for i in hdulist]
        self.slitwidth = [header_key_read(i, "SLIT") for i in hdulist]
        self.agstatus = [header_key_read(i, "AUTOGUID") for i in hdulist]
        self.slitpa = [header_key_read(i, "SLT-PA") for i in hdulist]
        self.nodpat = [header_key_read(i, "NODPAT") for i in hdulist]
        self.nodamp = [header_key_read(i, "NODAMP") for i in hdulist]

    def checkDataStatus(self, showDetail=True, ignore=False, statusOutput=False):
        settingSet = set([self.setting[i] for i in range(self.imnum)] + [self.compSetting, self.flatSetting])
        periodSet = set([self.period[i] for i in range(self.imnum)] + [self.compPeriod, self.flatPeriod])
        slitSet = set([self.slit[i] for i in range(self.imnum)] + [self.compSlit, self.flatSlit])
        modeSet = set([self.instmode[i] for i in range(self.imnum)] + [self.compMode, self.flatMode])
        if showDetail:
            print("\n# Instrument setting IDs from fits header")
            print("Inputs: ", self.setting)
            print("Calibs (comp,flat): {}, {}".format(self.compSetting, self.flatSetting))
            print("\n# Period IDs from fits header")
            print("Inputs: ", self.period)
            print("Calibs (comp,flat): {}, {}".format(self.compPeriod, self.flatPeriod))
            print("\n# Slit from fits header")
            print("Inputs: ", self.slit)
            print("Calibs (comp, flat): {}, {}".format(self.compSlit, self.flatSlit))
            print("\n# Mode from fits header")
            print("Inputs: ", self.instmode)
            print("Calibs (comp, flat): {}, {}".format(self.compMode, self.flatMode))
        if statusOutput:
            status_file = open(self.status, "a")
            status_file.write("Data status:\n")
            status_file.write("# Instrument setting IDs from fits header\n")
            status_file.write("    Inputs: {}\n".format(self.setting))
            status_file.write("    Calibs (comp,flat): {}, {}\n".format(self.compSetting, self.flatSetting))
            status_file.write("# Period IDs from fits header\n")
            status_file.write("    Inputs: {}\n".format(self.period))
            status_file.write("    Calibs (comp,flat): {}, {}\n".format(self.compPeriod, self.flatPeriod))
            status_file.write("# Slit from fits header\n")
            status_file.write("    Inputs: {}\n".format(self.slit))
            status_file.write("    Calibs (comp, flat): {}, {}\n".format(self.compSlit, self.flatSlit))
            status_file.write("# Mode from fits header\n")
            status_file.write("    Inputs: {}\n".format(self.instmode))
            status_file.write("    Calibs (comp, flat): {}, {}\n".format(self.compMode, self.flatMode))
            status_file.write("\n")

            status_file.close()

        if len(settingSet) == 1 and len(periodSet) == 1 and len(slitSet) == 1 and len(modeSet) == 1:
            return True
        else:
            if not ignore:
                print("\033[31m WARNING: Multiple datatypes are mixed in the input data. \033[0m")
            return False

    def showAllParams(self):
        print("## WARP Settings and Parameters")
        print("flag_apscatter: ", self.flag_apscatter)
        print("flag_manual_aperture: ", self.flag_manual_aperture)
        print("flag_skysub: ", self.flag_skysub)
        print("flag_bpmask: ", self.flag_bpmask)
        print("flag_skyemission: ", self.flag_skyemission)
        print("flag_wsmeasure: ", self.flag_wsmeasure)
        print("flag_wscorrect: ", self.flag_wscorrect)
        print("flag_wsmanual: ", self.flag_wsmanual)
        print("flag_extract2d: ", self.flag_extract2d)
        print("skysub_mode: ", self.skysub_mode)
        print("skysubModeList: ", self.skysubModeList)
        print("cutrange_list: ", self.cutrange_list)
        print("fluxinput: ", self.fluxinput)
        print("CRthreshold: ", self.CRthreshold)
        print("CRvaratio: ", self.CRvaratio)
        print("CRslitposratio: ", self.CRslitposratio)
        print("CRmaxsigma: ", self.CRmaxsigma)
        print("CRfixsigma: ", self.CRfixsigma)
        print("reduceFullData: ", self.reduceFullData)
        print("selectedOrders: ", self.selectedOrders)
        print()

    def showAllCalibs(self):
        print("## Calibration files")
        print("flat_file: ", self.flat_file)
        print("mask_file: ", self.mask_file)
        print("comp_file: ", self.comp_file)
        print("ap_file: ", self.ap_file)
        print("aptrans_file: ", self.aptrans_file)
        print("apsc_maskfile: ", self.apsc_maskfile)
        print("dyinput: ", self.dyinput)
        print()
