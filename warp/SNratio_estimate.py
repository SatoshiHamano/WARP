# -*- coding:utf-8 -*-

__version__ = "1.6"

import sys,os,datetime,glob
import numpy as np
import math,shutil
import matplotlib.pyplot as plt
from pyraf import iraf
import scipy.optimize
from warp.Spec1Dtools import openspecfits, binning_spec

iraf.noao()
iraf.onedspec()

#Description:
#
#   Using the difference spectra of the same order and same object,
#   the SN ratio of the combined spectrum can be estimated. Please
#   prepare the list file in which the names of fits files are
#   written.
#
#Usage:
#
#   $ python SNratio_estimate_v1_0.py <input> <output>
#
#   input -- The list of fits files of the same echelle order.
#   output -- The name of ascii file in which the results of SN estimate will be stored
#
#
#Updates:
#
#   ver1.0 made by S.Hamano in 2016/12/26
#   ver1.1-1.3 minor updates by S.Hamano
#   ver1.4
#       The measurement method of std. is changed from Gaussian fitting to the histgram to the simple statistical calculation.
#
#   ver1.5 updated by S.Hamano in 2018/04/19
#       - Import updated ver. of spectools (ver1.1).



def snestimate(inputlist, output):
    
    #set parameters for SN estimates

    clipsigma = 5   #sigma for clipping of outlier data in the difference spectrum.
    thres_low = 0.5 #the spectrum data lower than this threshold is not included in the SN calculation.
    thres_high = 1.3 #the spectrum data higher than this threshold is not included in the SN calculation.
    dividerange = 3 #the number of sections
    deleteflag = 1  #if 1, the intermediate files are deleted. If else, the intermediate files are left.
    binning_pix = 1 #the width of binning. If 1, there is no binning.
    binnum = 20 #the number of bins in the noise histgram.

    #setting the name of directory and the root of file name.

    direc = "tmpfits_%s-%d%d%d%d/" % (datetime.date.today(), datetime.datetime.now().hour, datetime.datetime.now().minute, datetime.datetime.now().second, datetime.datetime.now().microsecond)
    
    #creating the directory in which the intermediate files are stored.
    
    if not os.path.exists(direc):
        os.makedirs(direc)
    
    #open spectrum files


    spf = inputlist
    spnormf = []
    spcorf = []
    spmedian = []
    spoutput = []
    flist = ""
    for i in range(len(spf)):
        if i < 10:
            flist += spf[i]+","
        spnormf.append("%snorm%s" % (direc, spf[i].split("/")[-1]))
        tmplamx,tmpdata,rcrval1,rcdelt1,rcrpix1 = openspecfits(spf[i])
        spmedian.append(np.median(tmpdata))

    flist = flist.rstrip(",")


    #normalize the input spectrum files by the median of the spectrum itself

    for i in range(len(spf)):
        iraf.sarith(spf[i], "/", spmedian[i], spnormf[i])

    #combining and normalize the input spectrum files.

    iraf.scombine(flist, "%saverage_spec.fits" % direc, combine="average", group="all")
    iraf.continuum("%saverage_spec.fits" % direc, "%saverage_spec_norm.fits" % direc, order="10", high_rej="3", low_rej="2", function="spline3", interactive="no")


    #the created files are read and the binning is done.

    avelamx_tmp, avedata_tmp, avecrval1, avecdelt1, avecrpix1 = openspecfits("%saverage_spec_norm.fits" % (direc))
    avelamx, avedata = binning_spec(avelamx_tmp, avedata_tmp, binning_pix)

    #setting the ranges of the divided sections.

    minpix = 0
    maxpix = len(avelamx)-1
    widthpix = int((maxpix - minpix)/dividerange)

    minpix_ranges = [minpix + widthpix * i for i in range(dividerange)]
    maxpix_ranges = [minpix + widthpix * (i+1) for i in range(dividerange)]
    lamcen = [(avelamx[maxpix_ranges[i]] + avelamx[minpix_ranges[i]])/2. for i in range(dividerange)]


    #calculating the SN.

    difnormf = []
    difnormf_norm = []
    lamdif = []
    sn_vec = [[] for k in range(dividerange)]
    pairs = []
    stddev = [[] for i in range(len(spf))]

    for i in range(len(spf)):
        if i == len(spf) - 1:
            j = 0
        else:
            j = i+1

        difnormf.append(spnormf[i].rstrip("fits").rstrip(".") + "%d_%d.fits" % (i,j))
        difnormf_norm.append(spnormf[i].rstrip("fits").rstrip(".") + "%d_%dn.fits" % (i,j))
        iraf.sarith(spnormf[i], "-", spnormf[j], difnormf[i])
        iraf.continuum(difnormf[i], difnormf_norm[i], order="10", high_rej="3", low_rej="3", function="spline3", type="difference", interactive="no")

        lamx_tmp,normdata_tmp,crval1,cdelt1,crpix1 = openspecfits(difnormf_norm[i])
        lamx, normdata = binning_spec(lamx_tmp, normdata_tmp, binning_pix)
        normlength = len(lamx)
        pixx = np.array([l for l in range(normlength)])
        lamdif.append(lamx)

        pairs.append("%d--%d:\t\t" % ((i+1),(j+1)))

        for k in range(dividerange):
            lamdif_divide = lamdif[i][minpix_ranges[k]:min(maxpix_ranges[k],normlength)]
            normdata_divide = normdata[minpix_ranges[k]:min(maxpix_ranges[k],normlength)]
            divide_length = len(lamdif_divide)
            avedata_divide = avedata[minpix_ranges[k]:min(maxpix_ranges[k],normlength)]
        
            meddif = np.median(normdata_divide)
            stddif = np.std(normdata_divide)

            normdata_clipped = []
            lamdif_clipped = []
    
            for n in range(divide_length):
                if math.fabs(normdata_divide[n]-meddif) < stddif * clipsigma and thres_high > avedata_divide[n] > thres_low:
                    normdata_clipped.append(normdata_divide[n])
                    lamdif_clipped.append(lamdif_divide[n])

            std_cldata = np.std(normdata_clipped)

            print("%d -- %d (%d): sigma=%.1e" % (i,j,k,std_cldata))
            stddev[i].append(std_cldata)
            sn_vec[k].append(std_cldata**2)


    sn_sum = [sum(sn_vec[k])/2. for k in range(dividerange)]
    sntotal = [len(spf)/math.sqrt(sn_sum[k]) for k in range(dividerange)]
    
    #storing results into the output data

    wfile = open(output,"w")
    
    wfile.write("Time: %s\n\nInput: \t\n" % datetime.datetime.today())
    for i in range(len(inputlist)):
        wfile.write("\t\t%s\n" % inputlist[i])
    wfile.write("\nSigma for clipping: \t%d\n" % clipsigma)
    wfile.write("Threshold: \t%.2f-%.2f\n" % (thres_low, thres_high))
    wfile.write("# of division: \t%d\n" % dividerange)
    wfile.write("Binning size in pix: \t%d\n" % binning_pix)
    wfile.write("# of bins: \t%d\n" % binnum)
    
    wfile.write("\n\n")
    
    wfile.write("Frame pairs\t")
    for i in range(dividerange):
        wfile.write("%.1f\t" % lamcen[i])
    wfile.write("\n\n")
    
    for i in range(len(spf)):
        wfile.write(pairs[i])
        for k in range(dividerange):
            wfile.write("%.1e\t" % stddev[i][k])
        wfile.write("\n")

    wfile.write("\n\nS/N (total): \t")
    print("--------------------")
    for k in range(dividerange):
        print("S/N (total) at %.1f = %.1f" % (lamcen[k], sntotal[k]))
        wfile.write("%.1f\t" % sntotal[k])

    wfile.write("\n\n")
    wfile.close()

    #deleting the intermediate files

    if deleteflag == 1:
        flist = glob.glob("%s*" % direc)
        for i in range(len(flist)):
            os.remove(flist[i])
        os.rmdir(direc)

    #return results.

    return lamcen, sntotal



if __name__ == "__main__":

    filename = sys.argv[1:]

    spec1list = open(filename[0],"r")
    spec1lines = spec1list.readlines()
    spec1list.close()
    
    inputlist = []
    for i in range(len(spec1lines)):
        inputlist.append(spec1lines[i].split()[0])

    lamcen, sntotal = snestimate(inputlist, filename[1])

