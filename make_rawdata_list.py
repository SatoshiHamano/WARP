#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys, os, datetime, shutil
import glob
import numpy as np
from warp.Spec2Dtools import header_key_read
import argparse
import astropy.io.fits as fits

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("rawdatapath", type=str, help="directory path of input raw data")
    parser.add_argument("-t", "--uppTimeLimit", type=float, default=7200.,
                        help="Upper limit of time difference in a single list.")
    kwards = vars(parser.parse_args())
    rawdatapath = kwards["rawdatapath"]
    print("Raw data are in {}.".format(rawdatapath))

    uppTimeLimit = kwards["uppTimeLimit"]
    print("Upper limit of the time difference in a single list: {} sec".format(uppTimeLimit))
    rawfiles = glob.glob(rawdatapath + 'WINA00*fits')
    rawfiles.sort()
    hdrlist = []
    for f in rawfiles:
        with fits.open(f) as file:
            hdrlist.append(file[0].header)

    fname = np.array([i.split('/')[-1].rstrip("fits").rstrip(".") for i in rawfiles])
    object = np.array([header_key_read(i, "OBJECT").replace(" ", "_") for i in hdrlist])
    exptime = np.array([header_key_read(i, "EXPTIME") for i in hdrlist])
    nodpos = np.array([header_key_read(i, "NODPOS") for i in hdrlist])
    slit = np.array([header_key_read(i, "SLIT") for i in hdrlist])
    setting = np.array([header_key_read(i, "SETTING") for i in hdrlist])
    mode = np.array([header_key_read(i, "MODE") for i in hdrlist])
    period = np.array([header_key_read(i, "PERIOD") for i in hdrlist])
    acqtime = np.array([header_key_read(i, "ACQTIME1").split("-")[-1] for i in hdrlist])
    acqdate = np.array([header_key_read(hdrlist[i], "ACQTIME1").split()[0].rstrip(acqtime[i]).rstrip('-') for i in
                        range(len(hdrlist))])
    datatype = np.array([header_key_read(i, "DATA-TYP") for i in hdrlist])
    fitstime = np.array([])
    for i in range(len(acqtime)):
        [hour, minute, second] = acqtime[i].split(':')
        [obsyear, obsmonth, obsdate] = acqdate[i].split('-')
        fitstime = np.append(fitstime,
                             datetime.datetime(int(obsyear), int(obsmonth), int(obsdate), int(hour), int(minute),
                                               int(float(second))))

    object_list = set(object[datatype != 'TEST'])
    instwords = ['NA', 'COMPARISON', 'INSTFLAT', 'test', '']
    for w in instwords:
        object_list.discard(w)

    fileopen = False
    listlist = []
    for obj in object_list:
        counter = 1
        fnamelist = fname[object == obj]
        explist = exptime[object == obj]
        nodlist = nodpos[object == obj]
        slitlist = slit[object == obj]
        settinglist = setting[object == obj]
        modelist = mode[object == obj]
        periodlist = period[object == obj]
        acqlist = fitstime[object == obj]
        print(obj)

        settingSet = set(settinglist)
        modeSet = set(modelist)
        periodSet = set(periodlist)
        slitSet = set(slitlist)
        settingSet.discard('N/A')
        modeSet.discard('N/A')
        periodSet.discard('N/A')
        slitSet.discard('N/A')

        for s in settingSet:
            for m in modeSet:
                for p in periodSet:
                    for l in slitSet:
                        req = (settinglist == s) & (modelist == m) & (periodlist == p) & (slitlist == l)
                        acqobj = acqlist[req]
                        fnameobj = fnamelist[req]
                        expobj = explist[req]
                        nodobj = nodlist[req]
                        difacq = np.array([(np.roll(acqobj, -1)[i] - acqobj[i]).seconds for i in range(acqobj.size)])
                        for i in range(acqobj.size):
                            if not fileopen:
                                if not os.path.exists("{}_{}_{}{}_list.txt".format(obj, counter, m, l)):
                                    wfile = open("{}_{}_{}{}_list.txt".format(obj, counter, m, l), "w")
                                    listlist.append("{}_{}_{}{}_list.txt".format(obj, counter, m, l))
                                    fileopen = True
                                else:
                                    print("{}_{}_{}{}_list.txt already exists.".format(obj, counter, m, l))
                            if fileopen:
                                if nodobj[i].find("S") == -1:
                                    deltime = []
                                    kindex = []
                                    for k in range(len(fnameobj)):
                                        if k != i and nodobj[i][0] != nodobj[k][0] and expobj[i] == expobj[k]:
                                            if acqobj[i] > acqobj[k]:
                                                diftime = acqobj[i] - acqobj[k]
                                            else:
                                                diftime = acqobj[k] - acqobj[i]
                                            deltime.append(diftime.seconds)
                                            kindex.append(k)
                                    if len(kindex) == 0:
                                        print("Sky data was not found for {}. ({}, exp={:.2f}sec)".format(fnameobj[i],
                                                                                                         nodobj[i],
                                                                                                         expobj[i]))
                                    else:
                                        minindex = np.argmin(np.array(deltime))
                                        wfile.write("{} {}\n".format(fnameobj[i], fnameobj[kindex[minindex]]))
                                        print(
                                            "{} {}: position={}-{}, exptime={}, acqtime={}, delta_t={:.1f} sec\n".format(
                                                fnameobj[i], fnameobj[kindex[minindex]], nodobj[i],
                                                nodobj[kindex[minindex]], expobj[i], acqobj[i], np.amin(deltime)))

                            if difacq[i] > uppTimeLimit or i == acqobj.size - 1:
                                wfile.close()
                                fileopen = False

    print("\n\nFollowing files were created:\n")
    for l in listlist:
        print(l)