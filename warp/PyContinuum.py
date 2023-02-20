import sys
from pyraf import iraf

__version__ = "1.3"

# Description:
#   This "PyContinuum.py" script was made by Satoshi Hamano in 2016/04/28.
#
#   This script enables you to continuum fit with arbitraly parameters.
#
# Usage:
#
#   $ python PyContinuum.py <input> <output> <continuum parameter>
#
#   input -- input fits file or the list of fits files
#   output -- output fits file or the list fo fits files
#   continuum parameter -- the list of parameters of continuum fitting with the following format
#       ===starting format===
#       echelle_order    low_rej    high_rej   fitting_order     fitting_fuction   not_sampled_region
#       .
#       .
#       .
#       ===ending format===
#
#       echelle_order: the echelle order of spectrum. from m=42 to 61 for WINERED WIDE mode.
#       low_rej: low_rej parameter for IRAF/continuum.
#       high_rej: high_rej parameter for IRAF/continuum.
#       fitting_order: order parameter for IRAF/continuum.
#       fitting_fuction: fitting function for IRAF/continuum.
#       not_sampled_region: specify the region (in angstrom) not included in the continuum fitting.
#                           (e.g.,) 10010:10200,10440:10560
#                           Then, two regions, 10010-10200A and 10440-10560A, are NOT included in fitting.
#                           If you want to include whole range of the spectrum, please type "no" instead.
#
#
# Update:
#
# 2017.09.15 Hamano ver1.1 updated:
#
#   - The variable of "niterate" is added to PyContinuum function.
#
# 2018.05.24 Hamano ver1.2 updated:
#   - usable in python 3
#
# 2019.08.31 hamano ver1.3 updated:
#   - The name of the script was changed.
#   - Continuum spectrum can be created optionally.
#

def PyContinuum(input, output, continuum_lowrej, continuum_highrej, continuum_order, continuum_func, continuum_sample,
                continuum_nite, continuumspec="INDEF"):
    iraf.continuum(input, output, interactive="no", low_rej=continuum_lowrej, high_rej=continuum_highrej, grow=1.,
                   order=continuum_order, sample=continuum_sample, function=continuum_func, niterate=continuum_nite)
    if continuumspec != "INDEF":
        iraf.sarith(input, "/", output, continuumspec)



def read_continuum_parameters(paramfile, minlam, maxlam):
    contparam = open(paramfile, "r")
    contlines = contparam.readlines()
    contparam.close()
    continuum_lowrej = []
    continuum_highrej = []
    continuum_order = []
    continuum_func = []
    continuum_nosample = []
    for i in range(len(contlines)):
        continuum_lowrej.append(int(contlines[i].split()[1]))
        continuum_highrej.append(int(contlines[i].split()[2]))
        continuum_order.append(int(contlines[i].split()[3]))
        continuum_func.append(contlines[i].split()[4])
        continuum_nosample.append(contlines[i].split()[5])

    continuum_sample = []
    for i in range(len(continuum_nosample)):

        if continuum_nosample[i] == "no":
            continuum_sample.append("*")

        else:
            tmpconline = continuum_nosample[i].split(",")
            tmpconstart = []
            tmpconend = []
            tmpsample = ""
            for k in range(len(tmpconline)):
                tmpconstart.append(float(tmpconline[k].split(":")[0]))
                tmpconend.append(float(tmpconline[k].split(":")[1]))
            if tmpconstart[0] > minlam[i] and tmpconstart[0] < maxlam[i]:
                tmpsample += "%d:%d," % (minlam[i], tmpconstart[0])
            for k in range(len(tmpconline)):
                if k == len(tmpconline) - 1:
                    if tmpconend[k] < maxlam[i] and minlam[i] < tmpconend[k]:
                        tmpsample += "%d:%d," % (tmpconend[k], maxlam[i])

                else:
                    if tmpconstart[k + 1] < maxlam[i] and minlam[i] < tmpconend[k]:
                        tmpsample += "%d:%d," % (tmpconend[k], tmpconstart[k + 1])
                    elif tmpconend[k] > minlam[i] and tmpconend[k] < maxlam[i] and tmpconstart[k + 1] > maxlam[i]:
                        tmpsample += "%d:%d," % (tmpconend[k], minlam[i])
                    elif tmpconend[k] > minlam[i] and tmpconstart[k + 1] > minlam[i] and tmpconstart[k + 1] < maxlam[i]:
                        tmpsample += "%d:%d," % (minlam[i], tmpconstart[k + 1])

                if tmpsample == "":
                    continuum_sample.append("*")
                else:
                    continuum_sample.append(tmpsample.rstrip(","))

    return continuum_lowrej, continuum_highrej, continuum_order, continuum_func, continuum_sample


if __name__ == "__main__":

    filename = sys.argv[1:]
    input = filename[0]
    output = filename[1]

    continuum_lowrej, continuum_highrej, continuum_order, continuum_func, continuum_sample = read_continuum_parameters(
        filename[2])

    if input.find("fits"):
        if len(continuum_lowrej) > 1:
            print("Too much continuum parameters for input list.")
            sys.exit()
        PyContinuum(input, output, continuum_lowrej[0], continuum_highrej[0], continuum_order[0], continuum_func[0],
                    continuum_sample[0])

    else:
        inputf = open(input, "r")
        outputf = open(output, "r")
        inputl = inputf.readlines()
        outputl = outputl.readlines()
        if len(inputl) != len(outputl):
            print("Error: the lengths of input and output list is different.")
            sys.exit()
        if len(inputl) != len(continuum_lowrej):
            print("Error: the lengths of input and continuum parameter list is different.")

    for i in range(len(inputl)):
        PyContinuum(inputl[i].split()[0], outputl[i].split()[0], continuum_lowrej[i], continuum_highrej[i],
                    continuum_order[i], continuum_func[i], continuum_sample[i])
