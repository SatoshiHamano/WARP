# -*- coding:utf-8 -*-

__version__ = "1.5"

import sys
import astropy.io.fits as fits
from pyraf import iraf
from iraf import onedspec
import os
import shutil
import subprocess
import glob
import numpy as np
import argparse

sys.path.append(os.path.dirname(__file__))

from aperture_auto_trace import auto_aptrace
from aperture import apertureSet
from angle_measure import comparison_extract, auto_angle_measurement
from AP_FC_conversion import AP_to_FC
from cutransform import cutransform
from badpixmask import badpixmask_flaton, badpixmask_flatoff, pyfixpix, deadpixmap_inter
from apscatter import pyapscatter
from auto_ecidentify import auto_ecidentify
from Resolution import resolution_measure
from spectrum_map import spectrum_mapping
from Spec2Dtools import header_key_read
from Spec1Dtools import pyapall, truncate, dispcor_single

# from ECtoID_v1_3 import ECtoID


def constant_str_length(comment):
    constlength = 72
    print("\033[31m\n=== %s %s\n\033[0m" % (comment, ("=" * (constlength - len(comment) - 4))))


def remove_if_exist(filename):
    if os.path.exists(filename):
        os.remove(filename)


def main(inputfile, aperturereplace, transformonly):
    root_path = os.path.dirname(os.path.abspath(__file__))

    ###################################################################################################
    # Read input files. ###############################################################################
    ###################################################################################################

    inputlist = open(inputfile, "r")
    inputlines = inputlist.readlines()
    inputlist.close()

    constant_str_length("Reading %s" % inputfile)

    inputfiles = {}
    inputkeywords = ["pinhole", "flaton", "flatoff", "comp"]
    inputflag = [False for i in range(len(inputkeywords))]
    for i in inputlines:
        strline = i.split(":")
        if strline[0] in inputkeywords:
            inputfiles[strline[0]] = strline[1].split()[0].rstrip("fits").rstrip(".")
            inputflag[inputkeywords.index(strline[0])] = True

    if False in inputflag:
        for i in range(len(inputkeywords)):
            if not inputflag[inputkeywords[i]]:
                print("ERROR: \"%s\" file is not listed in %s." % (inputkeywords[i], inputfile))
        sys.exit()

    [pinholefile, flatonfile, flatofffile, compfile] = [inputfiles[inputkeywords[i]] for i in range(4)]

    for i in range(4):
        print("%s: %s" % (inputkeywords[i], inputfiles[inputkeywords[i]]))

    if transformonly:
        constant_str_length("Transform mode")
    if aperturereplace:
        constant_str_length("Aperture was replaced")

    ###################################################################################################
    # Read reference files. ###########################################################################
    ###################################################################################################

    winered_mode = ["WIDE", "HIRES-Y", "HIRES-J"]
    reference_files = {winered_mode[0]: "wide_reference_files.txt",
                       winered_mode[1]: "hiresy_reference_files.txt",
                       winered_mode[2]: "hiresj_reference_files.txt"}

    if os.path.exists(compfile + ".fits"):
        compfits = fits.open(compfile + ".fits")
        comphdr = compfits[0].header
        compfits.close()
        inputmode = header_key_read(comphdr, "INSTMODE")
        inputslit = header_key_read(comphdr, "SLIT")
        compdate = header_key_read(comphdr, "DATE-OBS").replace("-", "")
    else:
        sys.exit("ERROR: \"%s\" does not exist." % (compfile + ".fits"))

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

    flatdate = header_key_read(flathdr, "DATE-OBS").replace("-", "")
    flatmode = header_key_read(flathdr, "INSTMODE")
    flatslit = header_key_read(flathdr, "SLIT")

    print("\nMODE: %s" % inputmode)
    print("SLIT (comp): %s" % inputslit)
    print("SLIT (flat): %s" % flatslit)
    print("DATE (comp): %s" % compdate)
    print("DATE (flat): %s" % flatdate)

    if inputslit != flatslit:
        print("Warning: The slit widths of flat and comp data do not match!")

    flat_op = "flat_%s%s_%s" % (flatmode, flatslit, flatdate)

    mask_flatoff = "mask_flatoff_%s%s_%s" % (flatmode, flatslit, flatdate)
    mask_flaton = "mask_flaton_%s%s_%s" % (flatmode, flatslit, flatdate)
    mask_flat = "mask_flat_%s%s_%s" % (flatmode, flatslit, flatdate)
    deadpixmap_interorder = "dpmask_interorder_%s%s_%s" % (flatmode, flatslit, flatdate)
    deadpixmap_intraorder = "dpmask_intraorder_%s%s_%s" % (flatmode, flatslit, flatdate)

    onefitsimage = root_path + "/one.fits"

    aperturemask = "aperture_%s%s_%s.pl" % (flatmode, flatslit, flatdate)
    wideapmask = "wideap_%s%s_%s.pl" % (flatmode, flatslit, flatdate)
    wideapmaskfits = "wideap_%s%s_%s.fits" % (flatmode, flatslit, flatdate)

    flatoff_m = flatofffile + "_m"
    flaton_m = flatonfile + "_m"

    flat_op_m = flat_op + "_m"
    flat_op_msc = flat_op + "_msc"
    flat_op_mscm = flat_op + "_mscm"
    flat_op_mscmn = flat_op + "_mscmn"

    comp_op = "comp_%s%s_%s" % (inputmode, inputslit, compdate)
    comp_f = comp_op + "_f"
    comp_fm = comp_op + "_fm"

    comp_eclist = "comp_%s%s_%s_fm_ecall.list" % (inputmode, inputslit, compdate)
    comp_ecfits = "comp_%s%s_%s_fm_ecall.fits" % (inputmode, inputslit, compdate)
    comp_idfits = ["comp_%s%s_%s_fm_ecall_m%d.fits" % (inputmode, inputslit, compdate, apnum[i]) for i in
                   range(len(apnum))]
    comp_idfits_w = ["comp_%s%s_%s_fm_ecall_m%dw.fits" % (inputmode, inputslit, compdate, apnum[i]) for i in
                     range(len(apnum))]
    # comp_idref = ["comp_%s%s_%s_fm_ecall.%s" % (inputmode, inputslit, compdate, make_apnumtext(apnum[i])) for i in
    #               range(len(apnum))]

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

    if not transformonly:
        if aperturereplace:
            remove_if_exist(aperturemask)
            remove_if_exist(wideapmask)
            remove_if_exist(wideapmaskfits)

        if not os.path.exists(aperturemask):
            constant_str_length("Making aperture mask")
            iraf.apmask(pinholefile, aperturemask, interactive="no", find="no", recenter="no", resize="no", edit="no",
                        trace="no", mask="yes")
            iraf.apmask(flatofffile, wideapmask, interactive="no", find="no", recenter="no", resize="no", edit="no",
                        trace="no", mask="yes")
            iraf.imarith(onefitsimage, "*", wideapmask, wideapmaskfits)

        else:
            constant_str_length("SKIP Making aperture mask")

    ###################################################################################################
    # Make bad pix mask and flat ######################################################################
    ###################################################################################################

    if not transformonly:
        if aperturereplace:
            remove_if_exist(mask_flatoff + ".fits")
            remove_if_exist(mask_flaton + ".fits")
            remove_if_exist(mask_flat + ".fits")
            remove_if_exist(flatoff_m + ".fits")
            remove_if_exist(flaton_m + ".fits")
            remove_if_exist(flatofffile + "_medfilter.fits")
            remove_if_exist(flatonfile + "_norm.fits")
            remove_if_exist(flat_op_m + ".fits")
            remove_if_exist(flat_op_msc + ".fits")
            remove_if_exist("scatter_" + flat_op_m + ".fits")
            remove_if_exist(flat_op_mscm + ".fits")
            remove_if_exist(flat_op_mscmn + ".fits")
            remove_if_exist(deadpixmap_interorder)
            remove_if_exist(deadpixmap_intraorder)

        if not os.path.exists(flat_op_mscmn + ".fits"):
            constant_str_length("Make flat from %s/%s" % (flatonfile, flatofffile))

            badpixmask_flatoff(flatofffile, mask_flatoff)
            badpixmask_flaton(flatonfile, mask_flaton, deadpixmap_intraorder)

            deadpixmap_inter(flatonfile, wideapmaskfits, deadpixmap_interorder)

            iraf.imarith(mask_flaton, "+", mask_flatoff, mask_flat)

            pyfixpix(flatofffile, flatoff_m, mask_flatoff)
            pyfixpix(flatonfile, flaton_m, mask_flaton)

            iraf.imarith(flaton_m, "-", flatoff_m, flat_op_m)

            pyapscatter(flat_op_m, flat_op_msc, flat_op, "apscatter_log.txt")

            iraf.imarith(flat_op_msc, "*", aperturemask, flat_op_mscm)
            iraf.imarith(flat_op_mscm, "/", 10000, flat_op_mscmn)
        else:
            constant_str_length("SKIP Make flat from %s/%s" % (flatonfile, flatofffile))

    ###################################################################################################
    # Transform function and transform ################################################################
    ###################################################################################################

    ns = 15
    shift = [-21 + 3 * i for i in range(ns)]
    apfs = ["shift%d.fits" % shift[i] for i in range(int((ns + 1) / 2))] + \
           ["shift+%d.fits" % shift[i] for i in range(int((ns + 1) / 2), ns)]
    apfile_trans = [apfile_pinhole % apnum[i] for i in range(aplength)]
    pinhole_cut = [pinholefile + "_%dcut" % apnum[i] for i in range(aplength)]
    pinhole_trans = [pinholefile + "_%dtrans" % apnum[i] for i in range(aplength)]
    pinhole_fc = [pinholefile + "_%d" % apnum[i] for i in range(aplength)]

    if aperturereplace:
        for i in range(aplength):
            remove_if_exist("database/fc" + pinhole_fc[i])

    if not os.path.exists("database/fc" + pinhole_fc[0]):
        for i in range(len(apfs)):
            remove_if_exist(apfs[i])
        constant_str_length("Extract %s spectra" % compfile)
        comparison_extract(compfile, pinholefile, shift, apfs)

        constant_str_length("Measure the angle of slit images")

        [t0, t1, t2], [ap0, ap1, ap2] = auto_angle_measurement(compfile, shift, apfs, apnum, anglefile_npz)

        print("[t0, t1, t2], [ap0, ap1, ap2] = [%.3e, %.3e, %.3e], [%.3e, %.3e, %.3e]" % (t0, t1, t2, ap0, ap1, ap2))
        constant_str_length("Calculate transform function")

        AP_to_FC(pinholefile, [t0, t1, t2], [ap0, ap1, ap2])
    else:
        constant_str_length("SKIP Calculate transform function")
        params = np.load(anglefile_npz)
        [t0, t1, t2] = params["p0_yn"]
        [ap0, ap1, ap2] = params["p0_m"]

    if aperturereplace:
        for i in range(aplength):
            remove_if_exist(pinhole_cut[i] + ".fits")
            remove_if_exist(pinhole_trans[i] + ".fits")

    if not os.path.exists(pinhole_trans[0] + ".fits"):
        constant_str_length("Transform %s" % pinholefile)
        cutransform(pinholefile, pinholefile, "cutransform_log.txt", transdy, "no", calibflag=True)
    else:
        constant_str_length("SKIP Transform %s" % pinholefile)

    if transformonly:
        constant_str_length("Thank you :) (Transform mode)")
        sys.exit()

    ###################################################################################################
    # Aperture trace for transformed data #############################################################
    ###################################################################################################

    ap_trans_npz = [ap_trans_npz_core % apnum[i] for i in range(aplength)]

    if aperturereplace:
        for i in range(aplength):
            remove_if_exist("database/ap" + pinhole_trans[i])

    if not os.path.exists("database/ap" + pinhole_trans[-1]):
        constant_str_length("Trace transformed %s" % pinholefile)
        for i in range(aplength):
            apsetTrans = auto_aptrace(pinhole_trans[i], apfile_trans[i], ap_trans_npz[i], 5)
            for j in apsetTrans.echelleOrders:
                apsetTrans.apertures[j].apLow = -28
                apsetTrans.apertures[j].apHigh = 28
            bgregion_trans = ["-10:-5,5:10" for j in apsetTrans.echelleOrders]
            apsetTrans.write_apdatabase(pinhole_trans[i], bgregion_trans)
    else:
        constant_str_length("SKIP Trace transformed %s" % pinholefile)

    ###################################################################################################
    # Transform comparison image ######################################################################
    ###################################################################################################

    compfile_cut = [comp_fm + "_%dcut" % apnum[i] for i in range(aplength)]
    compfile_trans = [comp_fm + "_%dtrans" % apnum[i] for i in range(aplength)]
    compfile_trans_ec = [compfile_trans[i] + ".ec" for i in range(aplength)]
    compfile_trans_eccut = [compfile_trans[i] + ".eccut" for i in range(aplength)]

    if aperturereplace:
        remove_if_exist(comp_f + ".fits")
        remove_if_exist(comp_fm + ".fits")
        for i in range(aplength):
            remove_if_exist(compfile_cut[i] + ".fits")
            remove_if_exist(compfile_trans[i] + ".fits")
            remove_if_exist(compfile_trans_ec[i] + ".fits")
            remove_if_exist(compfile_trans_eccut[i] + ".fits")

    if not os.path.exists(comp_fm + ".fits"):
        constant_str_length("Transform %s" % compfile)

        iraf.imarith(compfile, "/", flat_op_mscmn, comp_f)
        pyfixpix(comp_f, comp_fm, mask_flat)
        cutransform(comp_fm, pinholefile, "cutransform_log.txt", transdy, "no", calibflag=True)
    else:
        constant_str_length("SKIP Transform %s" % compfile)

    ###################################################################################################
    # Extract comparison spectra ######################################################################
    ###################################################################################################

    if aperturereplace:
        remove_if_exist(comp_ecfits)
        for i in range(aplength):
            remove_if_exist(compfile_trans_ec[i] + ".fits")
            remove_if_exist(compfile_trans_eccut[i] + ".fits")
            remove_if_exist(comp_idfits[i])
            remove_if_exist(comp_idfits_w[i])

    if not os.path.exists(comp_ecfits):
        constant_str_length("Extract comparison spectra %s" % comp_ecfits)
        wf = open(comp_eclist, "w")
        for i in range(aplength):
            print(compfile_trans[i], compfile_trans_ec[i], pinhole_trans[i])
            pyapall(compfile_trans[i], compfile_trans_ec[i], pinhole_trans[i], "none", "echelle")
            truncate(compfile_trans_ec[i], compfile_trans_eccut[i])
            wf.write(compfile_trans_eccut[i] + ".fits\n")
        wf.close()
        iraf.scopy("@" + comp_eclist, comp_ecfits, renumber="yes", offset=min(apnum) - 1)
    else:
        constant_str_length("SKIP Extract comparison spectra %s" % comp_ecfits)

    ###################################################################################################
    # Ecidentify ######################################################################################
    ###################################################################################################

    if aperturereplace:
        remove_if_exist("database/ec%s" % comp_ecfits.rstrip("fits").rstrip("."))

    if not os.path.exists("database/ec%s" % comp_ecfits.rstrip("fits").rstrip(".")):
        constant_str_length("Identify emission lines in %s" % comp_ecfits)
        auto_ecidentify(comp_ecfits, inputslit, linelist)
        constant_str_length("Next steps:")
        print("\033[34m1) Copy the following commands to your terminal and run it.\n\n\033[0m")
        print("\033[34mcd calib_data_for_pipeline_%s \npython %s/ECtoID.py %s  \n \033[0m" % (compdate, root_path, comp_ecfits))
        print("\033[34m2) Centering the identified features using IRAF/ecidentify task.\033[0m")
        print("\033[34m3) Fit dispersion solution & add/delete features if necessary.\033[0m")
        # proc = subprocess.run(["python %s/ECtoID.py %s" % (root_path, comp_ecfits)],stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # resolution_list = open("resolution.list", "w")
        # for i in range(aplength):
        #     iraf.scopy(comp_ecfits + "[*,%d]" % (i + 1), comp_idfits[i])
        #     dispcor_single(comp_idfits[i], comp_idfits_w[i], comp_idref[i])
        #     resolution_list.write(comp_idfits_w[i] + "\n")
        # resolution_list.close()
        #
        # resolution_measure("resolution.list", linelist, "resolution", inputmode, inputslit)

        # ECtoID(comp_ecfits)
    else:
        constant_str_length("SKIP Identify emission lines in %s" % comp_ecfits)

    ###################################################################################################
    # File copy #######################################################################################
    ###################################################################################################

    outputdir = "calib_data_for_pipeline_" + compdate

    if os.path.exists(outputdir):
        shutil.rmtree(outputdir)

    os.makedirs(outputdir)

    inputf = open(outputdir + "/input_files.txt", "w")
    inputf.write("WINERED_pipeline\n")

    inputf.write("%s\t#flat file\n" % (flat_op_mscmn + ".fits"))
    shutil.copy(flat_op_mscmn + ".fits", outputdir)

    inputf.write("%s\t#mask file\n" % (mask_flat + ".fits"))
    shutil.copy(mask_flat + ".fits", outputdir)

    inputf.write("%s\t#comp file\n" % comp_ecfits)
    shutil.copy(comp_ecfits, outputdir)

    inputf.write("%s\t#ap file\n" % pinholefile)
    inputf.write("%s\t#aptrans file\n" % pinholefile)
    inputf.write("%s\t#ap file for apscatter" % flat_op)
    if os.path.exists(outputdir + "/database"):
        shutil.rmtree(outputdir + "/database")
    shutil.copytree("./database", "./%s/database" % outputdir)

    inputf.close()

    # spectrum_mapping(outputdir + ".pdf", comp_ecfits, pinholefile, transdy, [[t0, t1, t2], [ap0, ap1, ap2]])

    constant_str_length("Calib data is made.")
    print("\033[34m\t Calib data is in ./%s.\n\033[0m" % outputdir)
    constant_str_length("Thank you :)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("listfile", type=str, help="fits file list to be reduced")
    parser.add_argument("-a", "--aperturereplace", action="store_true", help="Aperture replace mode")
    parser.add_argument("-t", "--transformonly", action="store_true", help="Transform only mode")

    args = parser.parse_args()

    listfile = args.listfile
    aperturereplace = args.aperturereplace
    transformonly = args.transformonly

    main(listfile, aperturereplace, transformonly)
