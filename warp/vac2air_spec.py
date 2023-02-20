#!/usr/bin/env python
# -*- coding:utf-8 -*-
# 2015 9/3 made by S. Kondo
# ver1.1 updated by Hamano in 2017.01.25
#   - functionized
# ver1.2 updated by Hamano in 2018.05.24
#   - pyfits --> astropy.io.fits
#   - usable with Python 3
# ver1.3 updated by Hamano in 2019.09.08
#   - air2vac_spec is newly developed.
#   - refractive_index is newly developed.
#   - main function is newly written.
#   - the version number is removed from the name of the script.
#
# The conversion of the spectrum wavelength from vacuum wavelength to wavelength in air.
# The conversion formula is the same as APOGEE.
#

import sys
import math
import numpy as np
from astropy.io import fits
import scipy.constants
from pyraf import iraf
from iraf import onedspec
from warp.Spec1Dtools import openspecfits

__version__ = "1.3"


def refractive_index(lamvac):
    s = pow(10, 4) / lamvac
    n = 1.0 + 5.792105e-2 / (238.0185 - s ** 2) + 1.67917e-3 / (57.362 - s ** 2)

    return n


def vac2air_spec(input, output):
    if input.find(".fits") == -1:
        input += ".fits"

    spx, spy, _, _, _ = openspecfits(input)

    vac = np.average(spx)
    n = refractive_index(vac)

    vshift = (1. / n - 1) * (scipy.constants.c / 1000.)

    onedspec.dopcor(input, output, redshift=-vshift, isveloc="yes")
    iraf.images.imutil.hedit(output, "vac2air", vshift, add="yes", verify="no")  #
    iraf.images.imutil.hedit(output, "AIRORVAC", "air", add="yes", verify="no")  #


def air2vac_spec(input, output):
    if input.find(".fits") == -1:
        input += ".fits"

    spx, spy, _, _, _ = openspecfits(input)

    vac = np.average(spx)
    n = refractive_index(vac)

    vshift = (n - 1.) * (scipy.constants.c / 1000.)

    onedspec.dopcor(input, output, redshift=-vshift, isveloc="yes")
    iraf.images.imutil.hedit(output, "air2vac", vshift, add="yes", verify="no")  #
    iraf.images.imutil.hedit(output, "AIRORVAC", "vac", add="yes", verify="no")  #

# VALD3 stores vacuum wavelengths of all the transitions. This allows uniform handling of extraction across the whole spectral range. On the other hand, the selection tools in VALD3 include options for returning the air wavelengths. Furthermore, some of the original line lists include wavelengths measured in the air. This calls for conversion tools. Such tools must be uniformly accurate and reversible across the whole spectral range.

# For the vacuum to air conversion the formula from Donald Morton (2000, ApJ. Suppl., 130, 403) is used for the refraction index, which is also the IAU standard:

# n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2), where s = 10^4 / λvac and λvac is in Ångströms.

# The conversion is then: λair = λvac / n.

# This formula comes from Birch and Downs (1994, Metrologia, 31, 315) and applies to dry air at 1 atm pressure and 15ºC with 0.045% CO2 by volume. The corrections to Edlén (1953, J. Opt. Soc. Am., 43, 339) are less than 0.0001 Å at 2000 Å and less than 0.001 Å at 30000 Å.

if __name__ == "__main__":
    filename = sys.argv[1:]
    vac2air_spec(filename[0], filename[1])