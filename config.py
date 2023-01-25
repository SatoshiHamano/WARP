#!/usr/bin/env python
# -*- coding:utf-8 -*-

from astropy.io import fits


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


class config:
    def __init__(self):
        self.saturation_thres = 35000.
        self.flag_apscatter = True
        self.flag_manual_aperture = False
        self.flag_skysub = False
        self.flag_bpmask = True
        self.flag_skyemission = False
        self.flag_wsmeasure = True
        self.flag_wscorrect = True
        self.flag_wsmanual = False
        self.skysub_mode = "none"
        self.skysubModeList = ["none", "average", "median", "minimum", "fit"]
        self.cutrange_list = [1.05, 1.30]
        self.fluxinput = "no"
        self.CRthreshold = 10.
        self.CRvaratio = 2.
        self.CRslitposratio = 1.5
        self.CRmaxsigma = 20.
        self.CRfixsigma = False

    def readInputCalib(self, inputlist):
        para = []
        for line in open(inputlist, "r"):
            items = line.split()
            para.append(str(items[0]))
        self.flat_file = para[1]
        self.mask_file = para[2]
        self.comp_file = para[3]
        self.ap_file = para[4]
        self.aptrans_file = para[5]
        self.apsc_maskfile = para[6]

        compf = fits.open(self.comp_file)
        prihdr_comp = compf[0].header
        compf.close()
        self.dyinput = prihdr_comp["CDELT1"]

        return None

    def readParamQuery(self):
        print("======================================================")
        print("===    Please answer to the following questions.   ===")
        print("===                       OR                       ===")
        print("=== Just press Enter key (adopt default settings). ===")
        print("======================================================")

        ynDict = {"yes": True, "no": False}
        tfDict = {True: "yes", False: "no"}

        self.flag_apscatter = ynDict[
            alternativequestion("Subtract scattered light? (def:{}) :".format(tfDict[self.flag_apscatter]),
                                ["yes", "no", ""], tfDict[self.flag_apscatter])]
        self.flag_manual_aperture = ynDict[
            alternativequestion("Manually set aperture range? (def:{}) :".format(tfDict[self.flag_manual_aperture]),
                                ["yes", "no", ""], tfDict[self.flag_manual_aperture])]
        self.skysub_mode = alternativequestion(
            "Subtract background spectra from object spectra? (def:{}) :".format(self.skysub_mode),
            self.skysubModeList + [""], self.skysub_mode)
        self.flag_bpmask = ynDict[
            alternativequestion("Detect and interpolate the cosmic rays? (def:{}) :".format(tfDict[self.flag_bpmask]),
                                ["yes", "no", ""], tfDict[self.flag_bpmask])]

        ans_query_CRparams = alternativequestion(
            "Change any parameters in the cosmic ray detection algorithm? (def:no) :",
            ["yes", "no", "detail", ""], "no")
        if ans_query_CRparams == "yes":
            self.CRthreshold = valueQuestion(
                "Threshold for the cosmic ray detection (def: {} sigma) :".format(self.CRthreshold), 3., 100.,
                self.CRthreshold)
            self.CRfixsigma = ynDict[
                alternativequestion("Fix the threshold sigma (def: {}) :".format(tfDict[self.CRfixsigma]),
                                    ["yes", "no", ""], tfDict[self.CRfixsigma])]
            # ans_CRvaratio = df_CRvaratio
            # ans_CRmaxsigma = df_CRmaxsigma
            # ans_CRslitposratio = df_CRslitposratio
        elif ans_query_CRparams == "detail":
            self.CRthreshold = valueQuestion(
                "Threshold for the cosmic ray detection (def: {} sigma) :".format(self.CRthreshold), 3., 100.,
                self.CRthreshold)
            self.CRfixsigma = ynDict[
                alternativequestion("Fix the threshold sigma (def: {}) :".format(tfDict[self.CRfixsigma]),
                                    ["yes", "no", ""], tfDict[self.CRfixsigma])]
            if not self.CRfixsigma:
                self.CRmaxsigma = valueQuestion(
                    "Maximum threshold for the cosmic ray detection (def: {}) :".format(self.CRmaxsigma), 3., 100.,
                    self.CRmaxsigma)
                self.CRvaratio = valueQuestion(
                    "Threshold for the variance / average of the cosmic ray distribution (def: {}) :".format(
                        self.CRvaratio), 1., 100., self.CRvaratio)
                self.CRslitposratio = valueQuestion(
                    "Threshold for the cosmic ray number ratio between the slit positions (def: {}) :".format(
                        self.CRslitposratio), 1., 100., self.CRslitposratio)

        self.flag_skyemission = ynDict[
            alternativequestion("Extract the spectra from sky frame? (def: {}) :".format(tfDict[self.flag_skyemission]),
                                ["yes", "no", ""], tfDict[self.flag_skyemission])]
        self.flag_wsmeasure = ynDict[
            alternativequestion("Measure the pixel shifts? (def: {}) :".format(tfDict[self.flag_wsmeasure]),
                                ["yes", "no", ""], self.flag_wsmeasure)]
        self.flag_wscorrect = ynDict[
            alternativequestion("Correct the pixel shifts? (def: {}) :".format(tfDict[self.flag_wscorrect]),
                                ["yes", "no", ""], tfDict[self.flag_wscorrect])]
        if self.flag_wsmeasure and self.flag_wscorrect:
            self.flag_wsmanual = ynDict[alternativequestion(
                "Use the pixel shifts values written in list file? (def: {}) :".format(tfDict[self.flag_wsmanual]),
                ["yes", "no", ""], tfDict[self.flag_wsmanual])]
        elif not self.flag_wsmeasure and self.flag_wscorrect:
            self.flag_wsmanual = True
        else:
            self.flag_wsmanual = False

        self.fluxinput = alternativequestion(
            "Conserve the flux in the transformation? (def: {}) :".format(self.fluxinput),
            ["yes", "no", ""], self.fluxinput)

        self.flag_skysub = True if self.skysub_mode != "none" else False

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
                    print("WARNING: The input value could not be converted to float. value={} was set.".format(
                        self.CRthreshold))
            if line.find("Cosmic ray maximum sigma") != -1:
                try:
                    self.CRmaxsigma = float(line.split(":")[1].split()[0])
                except:
                    print("The input value could not be converted to float. value={} was set.".format(self.CRmaxsigma))
            if line.find("Cosmic ray Var/Ave ratio") != -1:
                try:
                    self.CRvaratio = float(line.split(":")[1].split()[0])
                except:
                    print("The input value could not be converted to float. value={} was set.".format(self.CRvaratio))
            if line.find("Cosmic ray ratio between slit positions") != -1:
                try:
                    self.CRslitposratio = float(line.split(":")[1].split()[0])
                except:
                    print("The input value could not be converted to float. value={} was set.".format(self.CRslitposratio))
            if line.find("Cosmic ray fix sigma") != -1:
                self.CRfixsigma = ynDict[line.split(":")[1].split()[0]]

        return None

    def writeStatus(self, wfile, pipelineVer, startTimeStr, endTimeStr, elapsedTime):
        tfDict = {True: "yes", False: "no"}

        status_file = open(wfile, "a")
        status_file.write("\nPIPELINE ver.%s\n\n" % (pipelineVer))
        status_file.write("\nStarting time: %s\n\n" % startTimeStr)
        status_file.write("\nSettings:\n")
        status_file.write("     Apscatter: %s\n" % tfDict[self.flag_apscatter])
        status_file.write("     Manual Aperture: %s\n" % tfDict[self.flag_manual_aperture])
        status_file.write("     Background Subtraction: %s\n" % self.skysub_mode)
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
        status_file.write("     CUTRANSFORM flux: %s\n\n\n" % self.fluxinput)
        status_file.write("Termination time: %s\n" % (endTimeStr))
        status_file.write("Elapsed time: %dm%.1fs\n\n" % (int(elapsedTime / 60), elapsedTime % 60))
        status_file.write("Pipeline status: Finished.")
        status_file.close()

    def showAllParams(self):
        print("flag_apscatter: ", self.flag_apscatter)
        print("flag_manual_aperture: ", self.flag_manual_aperture)
        print("flag_skysub: ", self.flag_skysub)
        print("flag_bpmask: ", self.flag_bpmask)
        print("flag_skyemission: ", self.flag_skyemission)
        print("flag_wsmeasure: ", self.flag_wsmeasure)
        print("flag_wscorrect: ", self.flag_wscorrect)
        print("flag_wsmanual: ", self.flag_wsmanual)
        print("skysub_mode: ", self.skysub_mode)
        print("skysubModeList: ", self.skysubModeList)
        print("cutrange_list: ", self.cutrange_list)
        print("fluxinput: ", self.fluxinput)
        print("CRthreshold: ", self.CRthreshold)
        print("CRvaratio: ", self.CRvaratio)
        print("CRslitposratio: ", self.CRslitposratio)
        print("CRmaxsigma: ", self.CRmaxsigma)
        print("CRfixsigma: ", self.CRfixsigma)

    def showAllCalibs(self):
        print("flat_file: ", self.flat_file)
        print("mask_file: ", self.mask_file)
        print("comp_file: ", self.comp_file)
        print("ap_file: ", self.ap_file)
        print("aptrans_file: ", self.aptrans_file)
        print("apsc_maskfile: ", self.apsc_maskfile)
        print("dyinput: ", self.dyinput)

