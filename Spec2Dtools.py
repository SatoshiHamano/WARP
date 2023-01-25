# -*- coding:utf-8 -*-

__version__ = "1.2"

import sys
from pyraf import iraf
import astropy.io.fits as fits

# 2019.08.31 S.Hamano created.
#   imarith was integrated.
#   savefitsimage was added from badpixmask.py

def flatfielding(input, output, flatfile):

    iraf.imarith(input, "/", flatfile, output)
    iraf.hedit(output, "FLAT", flatfile, add="yes", verify="no")

def header_key_read(hdulist, keyword):
    hdr_value = "N/A"

    try:
        hdr_value = hdulist[keyword]
    except:
        print(("No header value \"%s\" is found." % keyword))

    return hdr_value


def savefitsimage(inputdata, outputfits):

    hdu = fits.PrimaryHDU(inputdata)
    img = fits.HDUList([hdu])
    img.writeto(outputfits+".fits")
    img.close()

