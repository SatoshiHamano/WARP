#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np


class warpLog:
    def __init__(self, orders, frameNum):
        self.psfCenter = []
        self.psfWidth = []
        self.apertureLow = []
        self.apertureUpp = []
        self.waveShift = []
        self.waveShiftAve = []
        self.waveShiftStd = []
        self.waveShiftNum = []
        self.waveShiftAdopted = []
        self.cosmicRaySigma = []
        self.cosmicRayNum = []
        self.signalToNoise = []
        self.echelleOrderArray = np.array([[m for m in orders] for _ in range(frameNum)])
        self.frameNumberArray = np.array([[n for _ in orders] for n in range(frameNum)])
        self.echelleOrderVector = orders.copy()
        self.frameNum = frameNum
        self.echelleOrderNum = len(orders)

    def log_title(self, wf, title, description, constlength=72):
        wf.write("%s\n" % ("#" * constlength))
        wf.write("### %s %s\n" % (title, ("#" * (constlength - len(title) - 5))))
        wf.write("%s\n\n" % ("#" * constlength))
        descsplit = description.split("\n")
        for i in descsplit:
            wf.write("# %s\n" % i)
        wf.write("\n")

    def log_file_header(self, wf):
        wf.write("Files\t")
        for i in range(self.frameNum):
            wf.write("No.%d\t" % (i + 1))
        wf.write("\n\n")

    def log_maketable(self, wf, values):
        for j in range(self.echelleOrderNum):
            wf.write("m=%d\t" % self.echelleOrderVector[j])
            for i in range(self.frameNum):
                wf.write(values[i][j] + "\t")
            wf.write("\n")

    def psfLog(self, center, width):
        self.psfCenter = np.array(center.copy())
        self.psfWidth = np.array(width.copy())

    def writePsfLogText(self, logfile):
        wf = open(logfile, "w")

        self.log_title(wf, "Log of Center Search",
                       "The profile is assumed to be Gaussian.\n" +
                       "Peak position (pix) from the center position in the aperture definition.")
        self.log_file_header(wf)
        self.log_maketable(wf, [["%.4f\t" % self.psfCenter[i][j] for j in range(self.echelleOrderNum)] for i in
                                range(self.frameNum)])

        wf.write("\nMedian:\t")
        for i in range(self.frameNum):
            wf.write("%.4f\t" % np.nanmedian(self.psfCenter[i]))

        wf.write("\n\n# FWHM (pix) of fitted Gaussian profiles. \n\nFiles\t")

        self.log_file_header(wf)
        self.log_maketable(wf, [["%.4f\t" % self.psfWidth[i][j] for j in range(self.echelleOrderNum)] for i in
                                range(self.frameNum)])

        wf.write("\nMedian:\t")
        for i in range(self.frameNum):
            wf.write("%.4f\t" % np.nanmedian(self.psfWidth[i]))

        wf.close()

    def readPsfLogText(self, centersearch_log):
        rf = open(centersearch_log, "r")
        rl = rf.readlines()
        rf.close()

        readlist = []

        counter = 0
        for i in range(len(rl)):
            if rl[i].find("Files") != -1:
                readlist.append([[] for n in range(self.frameNum)])
                for j in range(i + 2, i + self.echelleOrderNum + 2):
                    rl1comp = rl[j].split()
                    for k in range(self.frameNum):
                        readlist[counter][k].append(float(rl1comp[k + 1]))
                counter += 1

        [xs, gw] = readlist
        self.psfCenter = np.array(xs.copy())
        self.psfWidth = np.array(gw.copy())

    def writePsfLogNpz(self, npzfile):
        np.savez(npzfile, psfCenter=self.psfCenter, psfWidth=self.psfWidth, echelleOrder=self.echelleOrderVector,
                 frameNum=self.frameNum)

    def readPsfLogNpz(self, npzfile):
        p = np.load(npzfile)
        if (self.echelleOrderVector == p["echelleOrder"]).all():
            if (self.frameNum == p["frameNum"]).all():
                self.psfCenter = p["psfCenter"]
                self.psfWidth = p["psfWidth"]
            else:
                print("WARNING: The log was not read, because the frame number was different.")
        else:
            print("WARNING: The log was not read, because the echelle order was different.")

    def apertureLog(self, aplow, apupp):
        self.apertureLow = np.array(aplow.copy())
        self.apertureUpp = np.array(apupp.copy())

    def writeApertureLogText(self, logfile):
        wf = open(logfile, "w")

        self.log_title(wf, "Log of Apertures", "Lower and uppper limit (pix) of the apertures.")
        self.log_file_header(wf)
        self.log_maketable(wf,
                           [["%.2f:%.2f\t" % (self.apertureLow[i][j], self.apertureUpp[i][j]) for j in
                             range(self.echelleOrderNum)] for i in range(self.frameNum)])

        wf.close()

    def readApertureLogText(self, logfile):
        rf = open(logfile, "r")
        rl = rf.readlines()
        rf.close()

        aplow_log = [[] for i in range(self.frameNum)]
        aphigh_log = [[] for i in range(self.frameNum)]

        for i in range(len(rl)):
            if rl[i].find("Files") != -1:
                for j in range(i + 2, i + self.echelleOrderNum + 2):
                    rl1comp = rl[j].split()
                    for k in range(self.frameNum):
                        aplow_log[k].append(float(rl1comp[k + 1].split(":")[0]))
                        aphigh_log[k].append(float(rl1comp[k + 1].split(":")[1]))
                break

        self.apertureLow = np.array(aplow_log)
        self.apertureUpp = np.array(aphigh_log)

    def writeApertureLogNpz(self, npzfile):
        np.savez(npzfile, apertureLow=self.apertureLow, apertureUpp=self.apertureUpp,
                 echelleOrder=self.echelleOrderVector,
                 frameNum=self.frameNum)

    def readApertureLogNpz(self, npzfile):
        p = np.load(npzfile)
        if (self.echelleOrderVector == p["echelleOrder"]).all():
            if (self.frameNum == p["frameNum"]).all():
                self.apertureLow = p["apertureLow"]
                self.apertureUpp = p["apertureUpp"]
            else:
                print("WARNING: The log was not read, because the frame number was different.")
        else:
            print("WARNING: The log was not read, because the echelle order was different.")

    def waveshiftLog(self, wsmatrix, wsave, wsstd, wsnum, wsadopted):
        self.waveShift = np.array(wsmatrix.copy())
        self.waveShiftAve = np.array(wsave.copy())
        self.waveShiftStd = np.array(wsstd.copy())
        self.waveShiftNum = np.array(wsnum.copy())
        self.waveShiftAdopted = np.array(wsadopted.copy())

    def writeWaveshiftLogText(self, logfile):
        wf = open(logfile, "w")

        self.log_title(wf, "Log of Waveshift", "Pixel shift relative to 1st frame.")
        self.log_file_header(wf)
        self.log_maketable(wf,
                           [["%.4f\t" % self.waveShift[j][i] for j in
                             range(self.echelleOrderNum)] for i in range(self.frameNum)])

        wf.write("\nAverage:\t")
        for i in range(self.frameNum):
            wf.write("%.4f\t" % self.waveShiftAve[i])
        wf.write("\nStd dev:\t")
        for i in range(self.frameNum):
            wf.write("%.4f\t" % self.waveShiftStd[i])
        wf.write("\nNumber of orders:\t")
        for i in range(self.frameNum):
            wf.write("%d\t" % self.waveShiftNum[i])
        wf.write("\n\nAdopted shift:\t")
        for i in range(self.frameNum):
            wf.write("%.4f\t" % self.waveShiftAdopted[i])

        wf.close()

    def readWaveshiftLog(self, logfile):
        if self.frameNum == 1:
            self.waveShift = np.array([[0. for j in range(self.echelleOrderNum)] for i in range(self.frameNum)])
            self.waveShiftAve = np.array([0. for i in range(self.frameNum)])
            self.waveShiftStd = np.array([0. for i in range(self.frameNum)])
            self.waveShiftNum = np.array([0 for i in range(self.frameNum)])
            self.waveShiftAdopted = np.array([0. for i in range(self.frameNum)])

        rf = open(logfile, "r")
        rl = rf.readlines()
        rf.close()

        wsmatrix = [[] for i in range(self.frameNum)]

        wsave = []
        wsstd = []
        wsnum = []
        wsadopted = []

        for i in range(len(rl)):
            if rl[i].find("Files") != -1:
                for j in range(i + 2, i + self.echelleOrderNum + 2):
                    rl1comp = rl[j].split()
                    for k in range(self.frameNum):
                        wsmatrix[k].append(float(rl1comp[k + 1]))
            if rl[i].find("Average:") != -1:
                rl1comp = rl[i].split(":")[1].split()
                for k in range(self.frameNum):
                    wsave.append(float(rl1comp[k]))
            if rl[i].find("Std dev:") != -1:
                rl1comp = rl[i].split(":")[1].split()
                for k in range(self.frameNum):
                    wsstd.append(float(rl1comp[k]))
            if rl[i].find("Number of orders:") != -1:
                rl1comp = rl[i].split(":")[1].split()
                for k in range(self.frameNum):
                    wsnum.append(int(rl1comp[k]))
            if rl[i].find("Adopted shift:") != -1:
                rl1comp = rl[i].split(":")[1].split()
                for k in range(self.frameNum):
                    wsadopted.append(float(rl1comp[k]))

        self.waveShift = np.array(wsmatrix)
        self.waveShiftAve = np.array(wsave)
        self.waveShiftStd = np.array(wsstd)
        self.waveShiftNum = np.array(wsnum)
        self.waveShiftAdopted = np.array(wsadopted)

    def writeWaveshiftLogNpz(self, npzfile):
        np.savez(npzfile, waveShift=self.waveShift, waveShiftAve=self.waveShiftAve, waveShiftStd=self.waveShiftStd,
                 waveShiftNum=self.waveShiftNum, waveShiftAdopted=self.waveShiftAdopted,
                 echelleOrder=self.echelleOrderVector, frameNum=self.frameNum)

    def readWaveshiftLogNpz(self, npzfile):
        p = np.load(npzfile)
        if (self.echelleOrderVector == p["echelleOrder"]).all():
            if (self.frameNum == p["frameNum"]).all():
                self.waveShift = p["waveShift"]
                self.waveShiftAve = p["waveShiftAve"]
                self.waveShiftStd = p["waveShiftStd"]
                self.waveShiftNum = p["waveShiftNum"]
                self.waveShiftAdopted = p["waveShiftAdopted"]
            else:
                print("WARNING: The log was not read, because the frame number was different.")
        else:
            print("WARNING: The log was not read, because the echelle order was different.")

    def cosmicRayLog(self, crsigma, crnum):
        self.cosmicRaySigma = np.array(crsigma.copy())
        self.cosmicRayNum = np.array(crnum.copy())

    def writeCosmicRayLogText(self, logfile):
        wf = open(logfile, "w")

        self.log_title(wf, "Log of cosmic ray detection", "Number of cosmic rays and threshold.")
        for i in range(self.frameNum):
            wf.write("No.%d: %d\t%.1f\n" % (i + 1, self.cosmicRayNum[i], self.cosmicRaySigma[i]))

        wf.close()

    def readCosmicRayLogText(self, logfile):
        rf = open(logfile, "r")
        rl = rf.readlines()
        rf.close()

        crthres_list = []
        crnum_list = []

        for i in range(len(rl)):
            for j in range(self.frameNum):
                if rl[i].find("No.{:d}:".format(j + 1)) != -1:
                    rl1comp = rl[i].split(":")[1].split()
                    crnum_list.append(int(rl1comp[0]))
                    crthres_list.append(float(rl1comp[1]))

        self.cosmicRaySigma = np.array(crthres_list)
        self.cosmicRayNum = np.array(crnum_list)

    def writeCosmicRayLogNpz(self, npzfile):
        np.savez(npzfile, cosmicRaySigma=self.cosmicRaySigma, cosmicRayNum=self.cosmicRayNum,
                 echelleOrder=self.echelleOrderVector, frameNum=self.frameNum)

    def readCosmicRayLogNpz(self, npzfile):
        p = np.load(npzfile)
        if (self.echelleOrderVector == p["echelleOrder"]).all():
            if (self.frameNum == p["frameNum"]).all():
                self.cosmicRaySigma = p["cosmicRaySigma"]
                self.cosmicRayNum = p["cosmicRayNum"]
            else:
                print("WARNING: The log was not read, because the frame number was different.")
        else:
            print("WARNING: The log was not read, because the echelle order was different.")


    def signalToNoiseLog(self, snratio):
        self.signalToNoise = snratio.copy()
