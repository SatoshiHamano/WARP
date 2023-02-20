import sys
from pyraf import iraf
from warp.aperture import *
from warp.config import constant_str_length

__version__ = "3.5"

iraf.twodspec()
iraf.longslit()
iraf.images()
iraf.imutil()
iraf.imgeom()

# made by hamano (2014-02-20)
# revised by hamano (2014-2-24): ver 2
#   add a "hedit" command to define the dispersion axis before "transform".
# revised by hamano (2015-6-12): ver 3.0
#   the part of calculation of the transformation matrix is separated
#   from this script to APtoFC_conversion.py.
#   - the calculation time become short.
#   - the usage is changed (see below).
#   - output file ##transshift.fits is removed. You will use ##trans.fits
# revised by hamano (2017-01): ver 3.1
#   - Log is added
#   - functionized
#   - some parts are moved to aptools.py script. some functions are called from aptools in this script.
#   - flux option is set to "yes"
#
# revised by Hamano (2017-04-13): ver3.2
#   - set flux option and dy values changeable
#
# revised by Hamano (2017-10-10): ver3.3
#   - In command-line mode, the dy parameter can also be input from commmand line from this ver.
#     In the previous ver(3.2), it was not possible to call this script in the command line because
#     of the lack of dy input.
#   - In command-line mode, the flux parameter is always set as "no".
#
# revised by Hamano (2018-05-24): ver3.4
#   - Usable with Python 3 series.
#
# revised by Hamano (2019-08-31): ver3.5
#   - The version number was removed from the name of the script.
#   - cutransform script for calibration mode was integrated. ("calibmode")
#

# usage:
#   python cutransform_v3_2.py <file1> <file2> <log>

#   <file1>
#       FITS file to be cut and transformed.

#   <file2>
#       FITS file with which aperture trace functions are obtained with "aptrace".
#       This script will read "database/ap<file2>".

#   <log>
#       ASCII file in which the parameters of transform will be stored.
#

#   example:
#       If you want to cut and transform "star.fits" using aperture function
#       defined with "galaxy.fits" conserving flux, you should run this script with:
#
#       python CUTRANSFORM.py star galaxy log.txt
#
#       Then, in your working directory, below files are generated:
#           1) star_<o>cut.fits        :fits images for order <o> that is just cut,
#           2) star_<o>trans.fits      :fits images transformed for order <o>,
#       In general use, 1) is not needed.
#       You will use only 2).


def cutransform_log(logf, input_file, referencename, xmargin, interpfunc, ypixsize, fluxopt):
    logfile = open(logf, "a")
    logfile.write("##########################\n")
    logfile.write("### Log of cutransform ###\n")
    logfile.write("##########################\n\n")

    logfile.write("Files: \n")
    logfile.write("\tinput:\t %s\n" % input_file)
    logfile.write("\treference:\t %s\n\n" % referencename)

    logfile.write("Parameters: \n")
    logfile.write("\txmargin=%d\n" % xmargin)
    logfile.write("\tinterptype=%s\n" % interpfunc)
    logfile.write("\tdy=%s\n" % ypixsize)
    logfile.write("\tflux=%s\n" % fluxopt)

    logfile.write("\n\n\n")

    logfile.close()


def cutransform(inputfile, referencename, logf, dyinput, fluxinput, calibflag=False):
    arraylength = 2048
    xmargin = 100
    interpfunc = "spline3"
    ypixsize = dyinput
    fluxopt = fluxinput

    if inputfile.split(".")[-1].find("fits") != -1:
        inputfile = inputfile.rstrip("fits").rstrip(".")
    if referencename.split(".")[-1].find("fits") != -1:
        referencename = referencename.rstrip("fits").rstrip(".")

    cutransform_log(logf, inputfile, referencename, xmargin, interpfunc, ypixsize, fluxopt)

    # read parameters defined with aptrace

    apset = apertureSet(referencename)

    xmin = []
    xmax = []
    for m in apset.echelleOrders:
        ap = apset.apertures[m]
        if int(ap.tracex[0] - xmargin) < 1:
            xmin.append(1)
        else:
            xmin.append(int(ap.tracex[0] - xmargin))
        if int(ap.tracex[-1] + xmargin) > arraylength:
            xmax.append(arraylength)
        else:
            xmax.append(int(ap.tracex[-1] + xmargin))

    iraf.hedit("%s" % (inputfile), "DISPAXIS", "2", verify="no", add="yes")

    for i in range(len(apset.echelleOrders)):
        m = apset.echelleOrders[i]
        constant_str_length("aperture %d is being reduced." % m)
        if calibflag:
            cutname, transname = "%s_%dcut" % (inputfile, m), "%s_%dtrans" % (inputfile, m)
        else:
            cutname, transname = "%s_m%dcut" % (inputfile, m), "%s_m%dtrans" % (inputfile, m)
        iraf.imcopy("%s[%d:%d,*]" % (inputfile, xmin[i], xmax[i]), cutname)
        print("\nimcopy: %s[%d:%d,] --> %s_%dcut\n" % (inputfile, xmin[i], xmax[i], inputfile, m))
        iraf.transform(cutname, transname, "%s_%d" % (referencename, m), flux="%s" % fluxopt, dy=ypixsize,
                       interptype=interpfunc)
        iraf.hedit(transname, "DISPAXIS", "2", verify="no")


if __name__ == "__main__":
    inputfile = str(sys.argv[1])
    referencename = str(sys.argv[2])
    logf = str(sys.argv[3])
    dyinput = float(sys.argv[4])
    fluxinput = "no"
    cutransform(inputfile, referencename, logf, dyinput, fluxinput)
