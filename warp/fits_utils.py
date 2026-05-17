# -*- coding:utf-8 -*-


def header_key_read(hdulist, keyword):
    hdr_value = "N/A"

    try:
        hdr_value = hdulist[keyword]
    except:
        print(("No header value \"%s\" is found." % keyword))

    return hdr_value
