# -*- coding:utf-8 -*-

__version__ = "1.0"

import sys
import astropy.io.fits as fits
from pyraf import iraf
from iraf import onedspec
from iraf import twodspec
from iraf import apextract
import os
import shutil
import glob

sys.path.append(os.path.dirname(__file__))

from aperture_auto_trace import auto_aptrace
from aperture import *
from badpixmask import badpixmask_flaton, badpixmask_flatoff, pyfixpix, deadpixmap_inter
from Spec2Dtools import header_key_read



def constant_str_length(comment):
    constlength = 72
    print("\033[31m\n=== %s %s\n\033[0m" % (comment, ("=" * (constlength - len(comment) - 4))))


def remove_if_exist(filename):
    if os.path.exists(filename):
        os.remove(filename)


def main():
    filename = sys.argv[1:]
    root_path = os.path.dirname(os.path.abspath(__file__))

    ###################################################################################################
    # Read input files. ###############################################################################
    ###################################################################################################

    inputlist = open(filename[0], "r")
    inputlines = inputlist.readlines()
    inputlist.close()

    constant_str_length("Reading %s" % filename[0])

    inputfiles = {}
    inputkeywords = ["pinhole", "flaton", "flatoff"]
    inputflag = [False for i in range(len(inputkeywords))]
    for i in inputlines:
        strline = i.split(":")
        if strline[0] in inputkeywords:
            inputfiles[strline[0]] = strline[1].split()[0].rstrip("fits").rstrip(".")
            inputflag[inputkeywords.index(strline[0])] = True

    if False in inputflag:
        for i in range(len(inputkeywords)):
            if not inputflag[inputkeywords[i]]:
                print("ERROR: \"%s\" file is not listed in %s." % (inputkeywords[i], filename[0]))
        sys.exit()

    [pinholefile, flatonfile, flatofffile] = [inputfiles[inputkeywords[i]] for i in range(3)]

    for i in range(3):
        print("%s: %s" % (inputkeywords[i], inputfiles[inputkeywords[i]]))

    transformonly = False
    aperturereplace = False
    if len(filename) > 1:
        for i in range(1, len(filename)):
            if filename[i] == "transform":
                transformonly = True
                constant_str_length("Transform mode")
        for i in range(1, len(filename)):
            if filename[i] == "aperture-replace":
                aperturereplace = True
                constant_str_length("Aperture was replaced")

    ###################################################################################################
    # Read reference files. ###########################################################################
    ###################################################################################################

    winered_mode = ["WIDE", "HIRES-Y", "HIRES-J"]
    reference_files = {winered_mode[0]: "wide_reference_files.txt",
                       winered_mode[1]: "hiresy_reference_files.txt",
                       winered_mode[2]: "hiresj_reference_files.txt"}

    if os.path.exists(flatonfile + ".fits"):
        flatonfits = fits.open(flatonfile + ".fits")
        flatonhdr = flatonfits[0].header
        flatonfits.close()
        inputmode = header_key_read(flatonhdr, "INSTMODE")
        inputslit = header_key_read(flatonhdr, "SLIT")
        flatdate = header_key_read(flatonhdr, "DATE-OBS").replace("-", "")
    else:
        sys.exit("ERROR: \"%s\" does not exist." % (flatonfile + ".fits"))

    refpath = root_path + "/reference/" + inputmode + "/"

    if not inputmode in winered_mode:
        print("Caution: WINERED has %s, %s and %s modes." % (winered_mode[0], winered_mode[1], winered_mode[2]))
        sys.exit("ERROR: \"%s\" is not in the list of WINERED observational modes" % inputmode)

    constant_str_length("Reading %s" % reference_files[inputmode])

    if os.path.exists(refpath + reference_files[inputmode]):
        reffile = open(refpath + reference_files[inputmode], "r")
        reflines = reffile.readlines()
        reffile.close()
    else:
        sys.exit("ERROR: \"%s\" does not exist." % (reference_files[inputmode]))

    refkeywords = ["order", "aperture(pinhole)", "aperture(flat)", "linelist(comp)", "aperture_position",
                   "aperture_position(transform)", "angle_function", "transform_dy"]
    refflag = [False for i in range(len(refkeywords))]
    reffiles = {}
    for i in reflines:
        strline = i.split(":")
        if strline[0] in refkeywords:
            reffiles[strline[0]] = strline[1].split()[0]
            refflag[refkeywords.index(strline[0])] = True

    if False in refflag:
        for i in range(len(refkeywords)):
            if not inputflag[refkeywords[i]]:
                print("ERROR: \"%s\" file is not listed in %s." % (refkeywords[i], reference_files[inputmode]))
        sys.exit()

    [order_word, apfile_pinhole, apfile_flat, linelist, ap_npz, ap_trans_npz_core, anglefile_npz, transdy] = [
        reffiles[refkeywords[i]] for i in range(8)]

    linelist = refpath + linelist
    ap_npz = refpath + ap_npz
    ap_trans_npz_core = refpath + ap_trans_npz_core
    anglefile_npz = refpath + anglefile_npz

    for i in range(8):
        print("%s: %s" % (refkeywords[i], reffiles[refkeywords[i]]))

    transdy = float(transdy)
    apnum = range(int(order_word.split("-")[0]), int(order_word.split("-")[1]) + 1)

    if not os.path.exists("./database"):
        os.makedirs("./database")
    databasefile = glob.glob(refpath + "database/*")
    for i in databasefile:
        shutil.copy(i, "./database")

    ###################################################################################################
    # Output file name ################################################################################
    ###################################################################################################

    if os.path.exists(flatonfile + ".fits"):
        flatfits = fits.open(flatonfile + ".fits")
        flathdr = flatfits[0].header
        flatfits.close()
    else:
        sys.exit("ERROR: \"%s\" does not exist." % (flatonfile + ".fits"))

    print("\nMODE: %s" % inputmode)
    print("SLIT: %s" % inputslit)
    print("DATE (flat): %s" % flatdate)

    flat_op = "flat_%s%s_%s" % (inputmode, inputslit, flatdate)

    mask_flatoff = "mask_flatoff_%s%s_%s" % (inputmode, inputslit, flatdate)
    mask_flaton = "mask_flaton_%s%s_%s" % (inputmode, inputslit, flatdate)
    mask_flat = "mask_flat_%s%s_%s" % (inputmode, inputslit, flatdate)
    deadpixmap_interorder = "dpmask_interorder_%s%s_%s" % (inputmode, inputslit, flatdate)
    deadpixmap_intraorder = "dpmask_intraorder_%s%s_%s" % (inputmode, inputslit, flatdate)

    onefitsimage = root_path + "/one.fits"

    aperturemask = "aperture_%s%s_%s.pl" % (inputmode, inputslit, flatdate)
    wideapmask = "wideap_%s%s_%s.pl" % (inputmode, inputslit, flatdate)
    wideapmaskfits = "wideap_%s%s_%s.fits" % (inputmode, inputslit, flatdate)

    ###################################################################################################
    # Aperture trace ##################################################################################
    ###################################################################################################

    aplength = len(apnum)

    if not os.path.exists("database/ap" + pinholefile):
        constant_str_length("Tracing %s" % pinholefile)
        apsetFlat = auto_aptrace(pinholefile, apfile_flat, ap_npz, 30)
        aperturereplace = True
    else:
        constant_str_length("SKIP Tracing %s" % pinholefile)
        constant_str_length("Reading %s" % pinholefile)
        apsetFlat = apertureSet(pinholefile)

    if aperturereplace:
        remove_if_exist("database/ap%s" % flat_op)
        remove_if_exist("database/ap%s" % flatonfile)

    if not os.path.exists("database/ap" + flat_op):
        bgregion_flat = ["-10:-5,5:10" for _ in apsetFlat.echelleOrders]

        for i in apsetFlat.echelleOrders:
            if i != apsetFlat.echelleOrders[-1]:
                apsetFlat.apertures[i].apLow = -32
                apsetFlat.apertures[i].apHigh = 32
            else:
                apsetFlat.apertures[i].apLow = -32
                apsetFlat.apertures[i].apHigh = 500

        apsetFlat.write_apdatabase(flat_op, bgregion_flat)
        apsetFlat.write_apdatabase(flatofffile, bgregion_flat)
        apsetFlat.renewOrders(apnum)
        bgregion = ["-10:-5,5:10" for _ in apsetFlat.echelleOrders]
        for i in apsetFlat.echelleOrders:
            apsetFlat.apertures[i].apLow = -27
            apsetFlat.apertures[i].apHigh = 26
        apsetFlat.write_apdatabase(flatonfile, bgregion)

        for i in apsetFlat.echelleOrders:
            apsetFlat.apertures[i].apLow = -28
            apsetFlat.apertures[i].apHigh = 28
        apsetFlat.write_apdatabase(pinholefile, bgregion)

    ###################################################################################################
    # Make aperture mask ##############################################################################
    ###################################################################################################

    constant_str_length("Making aperture mask")
    iraf.apmask(pinholefile, aperturemask, interactive="no", find="no", recenter="no", resize="no", edit="no",
                trace="no", mask="yes")
    iraf.apmask(flatofffile, wideapmask, interactive="no", find="no", recenter="no", resize="no", edit="no",
                trace="no", mask="yes")
    iraf.imarith(onefitsimage, "*", wideapmask, wideapmaskfits)

    ###################################################################################################
    # Make bad pix mask and flat ######################################################################
    ###################################################################################################

    constant_str_length("Make flat from %s/%s" % (flatonfile, flatofffile))

    badpixmask_flatoff(flatofffile, mask_flatoff)
    badpixmask_flaton(flatonfile, mask_flaton, deadpixmap_intraorder)

    deadpixmap_inter(flatonfile, wideapmaskfits, deadpixmap_interorder)

    iraf.imarith(mask_flaton, "+", mask_flatoff, mask_flat)


if __name__ == "__main__":
    main()
