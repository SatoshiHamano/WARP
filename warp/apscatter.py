# -*- coding:utf-8 -*-

__version__ = "3.2"

# usage:
#   python apscatter_v3_0.py <file1> <file2> <file3>

#   <file1>
#       FITS file to be "apscatter"
#
#   <file2>
#       FITS file after the "apscatter".
#
#   <file3>
#       FITS file, with which the mask information was defined.
#

# 2013.04.27 kondo version 1 作成#
# 2015.06.16 hamano version 2.0 revised
#   - base is changed to be input from command line.
# 2017.01.18 hamano version 3.0 revised
#   - Whole script is changed to function.
#   - apextract.apscatter.apscat2.order=15 is added.
# 2017.07.29 hamano version 3.1 revised
#   - Fitting function is changed from spline3 to legendre.
# 2019.08.31 hamano version 3.2 revised
#   - The name of the script is changed.
#

import sys
from pyraf import iraf
from iraf import twodspec
from iraf import apextract


def apscat_log(logf, input_file, output_file, base, ap1sample, ap1order, ap1low, ap1high, ap1nite, ap1function,
               ap2order, ap2sample, ap2nite, ap2function):
    logfile = open(logf, "a")
    logfile.write("########################\n")
    logfile.write("### Log of apscatter ###\n")
    logfile.write("########################\n\n")

    logfile.write("Files: \n")
    logfile.write("\tinput:\t %s\n" % input_file)
    logfile.write("\toutput:\t %s\n" % output_file)
    logfile.write("\treference:\t %s\n\n" % base)

    logfile.write("Parameters: \n")
    logfile.write("\tapextract.apscatter.apscat1.order=%d\n" % ap1order)
    logfile.write("\tapextract.apscatter.apscat1.sample=%s\n" % ap1sample)
    logfile.write("\tapextract.apscatter.apscat1.low_reject=%.1f\n" % ap1low)
    logfile.write("\tapextract.apscatter.apscat1.high_reject=%.1f\n" % ap1high)
    logfile.write("\tapextract.apscatter.apscat1.niterate=%d\n" % ap1nite)
    logfile.write("\tapextract.apscatter.apscat1.function=%s\n" % ap1function)
    logfile.write("\tapextract.apscatter.apscat2.order=%d\n" % ap2order)
    logfile.write("\tapextract.apscatter.apscat2.sample=%s\n" % ap2sample)
    logfile.write("\tapextract.apscatter.apscat2.niterate=%d\n" % ap2nite)
    logfile.write("\tapextract.apscatter.apscat2.function=%s\n" % ap2function)

    logfile.write("\n\n\n")

    logfile.close()


def pyapscatter(input_file, output_file, base, logf):
    ap1sample = "10:2000"
    ap1order = 3
    ap1low = 3.0
    ap1high = 2.0
    ap1nite = 100
    ap1function = "legendre"
    ap2order = 5
    ap2sample = "10:2000"
    ap2nite = 20
    ap2function = "legendre"

    apscat_log(logf, input_file, output_file, base, ap1sample, ap1order, ap1low, ap1high, ap1nite, ap1function,
               ap2order, ap2sample, ap2nite, ap2function)

    ##### base para #####
    apextract.apscatter.apscat1.order = ap1order
    apextract.apscatter.apscat1.sample = ap1sample
    apextract.apscatter.apscat1.low_reject = ap1low
    apextract.apscatter.apscat1.high_reject = ap1high
    apextract.apscatter.apscat1.niterate = ap1nite
    apextract.apscatter.apscat1.function = ap1function
    apextract.apscatter.apscat2.order = ap2order
    apextract.apscatter.apscat2.sample = ap2sample
    apextract.apscatter.apscat2.niterate = ap2nite
    apextract.apscatter.apscat2.function = ap2function

    ##### apscatter#####
    apextract.apscatter.refe = base
    apextract.apscatter.inter = "no"
    apextract.apscatter.find = "no"
    apextract.apscatter.rece = "no"
    apextract.apscatter.resize = "no"
    apextract.apscatter.edit = "no"
    apextract.apscatter.trace = "no"
    apextract.apscatter.fittrac = "no"
    apextract.apscatter.scatter = "scatter_" + input_file
    apextract.apscatter(input_file, output=output_file)

    iraf.hedit(output_file, "APSCMASK", base, add="yes", verify="no")


if __name__ == "__main__":
    input_file = str(sys.argv[1])
    output_file = str(sys.argv[2])
    base = str(sys.argv[3])
    logf = str(sys.argv[4])

    pyapscatter(input_file, output_file, base, logf)
