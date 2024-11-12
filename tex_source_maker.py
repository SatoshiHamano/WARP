# -*- coding:utf-8 -*-

# ver1.1 updated 2017.12
#   - The symbol of underscore("_") can be included in the file name.

# ver1.2 updated by Hamano 2017.12.25
#   - tex_images is renamed to tex_spprofile.
#   - tex_slitviewer is added.

# ver1.3 updated by Hamano 2018.01.15?
#   - Exception handling was added to function "table_1line_float_range".

# ver1.4 updated by Hamano 2018.04
#   - Function "slitlength_scale_identifier" was added.
#   - Table style was slightly adjusted.
#   - \clearpage was added to prevent the collaption of figure positions.
#   - Some new parameters were added to the table in the cover page.

# ver2.0 updated by Hamano 2019.12.01
#   - Function "tex_source_make" is newly added.
#
# ver2.1 updated by Hamano 2020.02.24
#   - Addition of pipeline parameters (flag_wsmeasure, flag_wscorrect, flag_wsmanual) were incorpolated.
#   - Refactoring
#

import time, datetime
import numpy as np
import shutil
import sys
import os
import subprocess
from  astropy.io import fits

from warp.aperture import *
from warp.Spec2Dtools import header_key_read
import Warp_sci
from warp.logger import warpLog
from warp.config import config

__version__ = "2.2"


def shortening_list(inputlist):
    shortlist = []

    for i in range(len(inputlist)):
        flag = 0
        for j in range(len(shortlist)):
            if shortlist[j] == inputlist[i]:
                flag += 1
        if flag == 0:
            shortlist.append(inputlist[i])

    return shortlist


############################################

def shortening_float_list(inputlist, strnum):
    shortlist = []

    for i in range(len(inputlist)):
        if inputlist[i] == "N/A":
            shortlist.append(inputlist[i])
        elif len(str(inputlist[i])) < strnum:
            shortlist.append(str(inputlist[i]))
        else:
            shortlist.append(str(inputlist[i]).split()[0][0:strnum])

    return shortlist


############################################
def table_1line(wf, key, value_list):
    wf.write("%s & " % key)
    for i in range(len(value_list)):
        wf.write("%s" % str(value_list[i]).replace("_", " "))
        if i != len(value_list) - 1:
            wf.write(",")
    wf.write("\\\\\n")


############################################

def table_1line_float_range(wf, key, value_list):
    value_float = []
    for i in range(len(value_list)):
        if value_list[i] != "N/A":
            try:
                value_float.append(float(value_list[i]))
            except:
                print("Failed to convert %s to float!" % key)

    wf.write("%s & " % key)
    if value_float == []:
        wf.write("N/A \\\\\n")
    else:
        wf.write("%.3f -- %.3f \\\\\n" % (min(value_float), max(value_float)))


############################################

def table_1line_time_range(wf, key, value_list):
    values = []
    for i in range(len(value_list)):
        if value_list[i] != "N/A":
            values.append(value_list[i])

    wf.write("%s & " % key)
    if values == []:
        wf.write("N/A \\\\\n")
    else:
        wf.write("%s to %s \\\\\n" % (values[0], values[len(values) - 1]))


############################################

def table_1line_sum(wf, key, value_list):
    value_sum = 0
    for i in range(len(value_list)):
        if value_list[i] != "N/A":
            value_sum += float(value_list[i])

    wf.write("%s & " % key)
    if value_sum == 0:
        wf.write("N/A \\\\\n")
    else:
        wf.write("%.0f \\\\\n" % value_sum)


############################################

def pix_scale_identifier(obsdate):
    pix_scale = "N/A"

    thres_date1 = datetime.date(2015, 1, 1)
    thres_date2 = datetime.date(2017, 1, 1)
    thres_date3 = datetime.date(2021, 1, 1)
    pix_scale_Araki1 = "0.82"
    pix_scale_Araki2 = "0.75"
    pix_scale_NTT = "0.27"
    pix_scale_Clay = "0.144"
    if len(obsdate.split("-")) == 3:
        [year_obs, month_obs, day_obs] = obsdate.split("-")
        obsdate_date = datetime.date(int(year_obs), int(month_obs), int(day_obs))
        if obsdate_date < thres_date1:
            pix_scale = pix_scale_Araki1
        elif thres_date1 < obsdate_date < thres_date2:
            pix_scale = pix_scale_Araki2
        elif thres_date2 < obsdate_date < thres_date3:
            pix_scale = pix_scale_NTT
        elif thres_date3 < obsdate_date:
            pix_scale = pix_scale_Clay

    return pix_scale


############################################

def slitlength_scale_identifier(obsdate):
    slitlength = "N/A"

    thres_date1 = datetime.date(2015, 1, 1)
    thres_date2 = datetime.date(2017, 1, 1)
    thres_date3 = datetime.date(2021, 1, 1)
    slitlength_Araki1 = "49.63"
    slitlength_Araki2 = "45.39"
    slitlength_NTT = "16.34"
    slitlength_Clay = "8.7"
    if len(obsdate.split("-")) == 3:
        [year_obs, month_obs, day_obs] = obsdate.split("-")
        obsdate_date = datetime.date(int(year_obs), int(month_obs), int(day_obs))
        if obsdate_date < thres_date1:
            slitlength = slitlength_Araki1
        elif thres_date1 < obsdate_date < thres_date2:
            slitlength = slitlength_Araki2
        elif thres_date2 < obsdate_date < thres_date3:
            slitlength = slitlength_NTT
        elif thres_date3 < obsdate_date:
            slitlength = slitlength_Clay

    return slitlength


############################################

def exptime_identifier(obsdate, exptime, inttime):
    exptime_selected = exptime

    thres_date1 = datetime.date(2015, 1, 1)
    if len(obsdate.split("-")) == 3:
        [year_obs, month_obs, day_obs] = obsdate.split("-")
        obsdate_date = datetime.date(int(year_obs), int(month_obs), int(day_obs))
        if obsdate_date < thres_date1:
            exptime_selected = inttime
        elif thres_date1 < obsdate_date:
            exptime_selected = exptime

    return exptime_selected


############################################

def slitwidth_converter(slitwidthlist, obsdate):
    pix_scale = pix_scale_identifier(obsdate)
    slitwidth_arcsec = ["N/A" for i in range(len(slitwidthlist))]
    mic_pix_ratio = 50.

    if pix_scale != "N/A":
        pix_scale_value = float(pix_scale.split()[0])
        for i in range(len(slitwidthlist)):
            if slitwidthlist[i] != "N/A":
                slitwidth_arcsec[i] = "%.2f" % (pix_scale_value * float(slitwidthlist[i]) / mic_pix_ratio)

    return slitwidth_arcsec


############################################


def tex_cover_and_obsinfo(texfile, conf: config):
    wf = open(texfile, "w")

    wf.write(r"""\documentclass{article}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{fancyhdr}
\usepackage{longtable}
\usepackage{grffile}
\usepackage{color} % Required for color customization

\lhead{\begin{minipage}{10mm}
\vspace*{-10mm}
\includegraphics[width=30mm]{winered_logo.pdf}
\end{minipage}}
\rhead{\leftmark}

\addtolength{\oddsidemargin}{-70pt}
\addtolength{\marginparwidth}{-50pt}
\addtolength{\textwidth}{140pt}

\pagestyle{fancy}


\title{WINERED pipeline reduction log}
\author{This file was created automatically by the pipeline.}
        """)

    wf.write("\date{Report created on %s.}\n" % time.ctime())

    wf.write(r"""
        
\begin{document}
        
\maketitle

\section*{Observation information}

\begin{table}[htb]
\centering
\begin{tabular}{p{15em}p{15em}} \hline
Key & Value \\ \hline
""")

    objname_short = shortening_list(conf.objname)
    modes_short = shortening_list(conf.modes)
    teles_short = shortening_list(conf.teles)
    period_short = shortening_list(conf.period)
    setting_short = shortening_list(conf.setting)
    acqdate_short = shortening_list(conf.acqdate)
    observer_short = shortening_list(conf.observer)
    observatory_short = shortening_list(conf.observatory)
    wodbobsid_short = shortening_list(conf.wodbobsid)
    wodbstd_short = shortening_list(conf.wodbstd)
    wodbpi_short = shortening_list(conf.wodbpi)
    wodbtheme_short = shortening_list(conf.wodbtheme)
    slitwidth_short = shortening_list(conf.slitwidth)
    agstatus_short = shortening_list(conf.agstatus)
    nodpat_short = shortening_list(conf.nodpat)
    nodamp_short = shortening_list(conf.nodamp)

    airmass_start_short = shortening_float_list(conf.airmass_start, 5)
    airmass_end_short = shortening_float_list(conf.airmass_end, 5)
    seeing_short = shortening_float_list(conf.seeing, 4)
    air_pressure_short = shortening_float_list(conf.air_pressure, 5)
    temperature_short = shortening_float_list(conf.temperature, 5)
    wind_speed_short = shortening_float_list(conf.wind_speed, 3)
    slitpa_short = shortening_float_list(conf.slitpa, 5)

    exptime = exptime_identifier(acqdate_short[0], conf.exptime, conf.inttime)

    nodpos_obj, nodpos_sky = [], []
    exptime_obj, exptime_sky = [], []
    slitpa_obj, slitpa_sky = [], []
    satupix_obj = []
    for i in range(len(conf.objectlist)):
        for j in range(len(conf.imagelist)):
            if conf.objectlist[i] == conf.imagelist[j]:
                exptime_obj.append(exptime[j])
                nodpos_obj.append(conf.nodpos[j])
                slitpa_obj.append(slitpa_short[j])
                satupix_obj.append(conf.satupix[j])
            elif conf.skylist[i] == conf.imagelist[j]:
                exptime_sky.append(exptime[j])
                nodpos_sky.append(conf.nodpos[j])
                slitpa_sky.append(slitpa_short[j])

    table_1line(wf, "Object", objname_short)
    wf.write("R.A. & %s\\\\\n" % conf.ra_hours[0])
    wf.write("Dec. & %s\\\\\n" % conf.dec_degree[0])
    table_1line(wf, "WODB OBS ID", wodbobsid_short)
    table_1line(wf, "WODB STD ID", wodbstd_short)
    table_1line(wf, "WODB PI", wodbpi_short)
    table_1line(wf, "WODB Theme", wodbtheme_short)
    wf.write(" & \\\\\n")

    table_1line(wf, "Mode", modes_short)
    table_1line(wf, "Observatory", observatory_short)
    table_1line(wf, "Telescope", teles_short)
    table_1line(wf, "Period", period_short)
    table_1line(wf, "Setting", setting_short)
    table_1line(wf, "Pixel scale (arcsec/pix)", [pix_scale_identifier(acqdate_short[0])])
    table_1line(wf, "Slit width (arcsec)", slitwidth_converter(slitwidth_short, acqdate_short[0]))
    table_1line(wf, "Slit length (arcsec)", [slitlength_scale_identifier(acqdate_short[0])])
    wf.write(" & \\\\\n")

    table_1line(wf, "Date (UT)", acqdate_short)
    table_1line_time_range(wf, "Time (UT)", conf.acqtime)
    table_1line(wf, "Observer", observer_short)
    table_1line_float_range(wf, "Airmass", conf.airmass)
    table_1line_float_range(wf, "Seeing (arcsec)", conf.seeing)
    table_1line_sum(wf, "Exptime (sec)", exptime_obj)
    table_1line(wf, "Auto guider status", agstatus_short)
    table_1line(wf, "Nodding pattern", nodpat_short)
    table_1line(wf, "Nodding amplitude (arcsec)", nodamp_short)
    table_1line_float_range(wf, "Slit PA (deg.)", conf.slitpa)

    wf.write(r"""\hline
\end{tabular}
\end{table}

\newpage

\section*{Reduced data}

Detailed information of each frame.

\scriptsize
\begin{longtable}[c]{cccccccccccccc}
\hline
Fits & Slit & I.T. & UT &            & R.A. & Dec. & Airmass &       & Seeing & RH$^a$ & $P_\text{air}$ & $T_\text{amb}$        & $V_\text{wind}$ \\
     &         & sec  & (start) & (end) &   &         & (start) & (end) & arcsec & \%      & mbar         & $^\circ\mathrm{C}$ & m/sec      \\
\endfirsthead
\multicolumn{8}{l}{\small\slshape continued from previous page} \\
\hline
Fits & Slit & I.T. & UT &            & R.A. & Dec. & Airmass &       & Seeing & RH$^a$ & $P_\text{air}$ &   $T_\text{amb}$       & $V_\text{wind}$ \\
     &         & sec  & (start) & (end) &   &         & (start) & (end) & arcsec & \%      & mbar         & $^\circ\mathrm{C}$ & m/sec      \\
\hline
\endhead
\hline
\multicolumn{8}{r}{\small\sl table continued on next page} \\
\endfoot
\hline
\endlastfoot
\hline
""")
    for i in range(len(conf.imagelist)):
        wf.write(" %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n" % (
            conf.imagelist[i].replace("_", "\_"), conf.nodpos[i], exptime[i], conf.ut_start[i].split('.')[0],
            conf.ut_end[i].split('.')[0],
            conf.ra_hours[i][0:-2], conf.dec_degree[i][0:-2], airmass_start_short[i], airmass_end_short[i], seeing_short[i],
            conf.humidity[i], air_pressure_short[i], temperature_short[i], wind_speed_short[i]))

    wf.write(r"""\hline
\end{longtable}
\vspace{-10pt}
$^a$Relative humidity.

\normalsize

\vspace{40pt}

\noindent
Frame pairs.

\begin{longtable}[c]{ccccccccc}
\hline
& OBJECT & & & & & SKY & & \\
Frame & Fits & I.T. & Position & Slit PA & $>$35000 (pix)$^a$ & Fits &  I.T. & Position \\
\endfirsthead
\multicolumn{8}{l}{\small\slshape continued from previous page} \\
\hline
& OBJECT & & & & & SKY & & \\
Frame & Fits & I.T. & Position & Slit PA & $>$35000 (pix)$^a$ & Fits &  I.T. & Position \\
\hline
\endhead
\hline
\multicolumn{8}{r}{\small\sl table continued on next page} \\
\endfoot
\hline
\endlastfoot
\hline
""")
    for i in range(len(conf.objectlist)):
        wf.write("No.%d & %s & %d sec & %s & %s & %d & %s & %d sec & %s \\\\\n" % (
            i + 1, conf.objectlist[i].replace("_", "\_"), exptime_obj[i], nodpos_obj[i], slitpa_obj[i], satupix_obj[i],
            conf.skylist[i].replace("_", "\_"), exptime_sky[i], nodpos_sky[i]))

    wf.write(r"""\hline
\end{longtable}
\vspace{-10pt}
$^a$Number of pixels with $>$35000 counts, which are possibly saturated.

%------------------------------------------------------------------------------------------------------------------------------------------------
\newpage
""")

    wf.close()


############################################

def tex_pipeline_para(texfile, pipeline_ver, conf):
    wf = open(texfile, "a")
    wf.write(r"""\section*{Pipeline information}

\begin{table}[htb]
\centering
\begin{tabular}{p{18em}p{7em}} \hline
Parameter & Value \\ \hline
""")
    wf.write("ver. & %s\\\\\n" % pipeline_ver)
    yes_or_no = {True: "yes", False: "no"}
    wf.write("Scattered light subtraction & %s \\\\\n" % yes_or_no[conf.flag_apscatter])
    wf.write("Manual aperture range setting & %s \\\\\n" % yes_or_no[conf.flag_manual_aperture])
    wf.write("Background subtraction mode & %s \\\\\n" % conf.skysub_mode)
    # wf.write("(background region) & %s \\\\\n" % skysubs_region)
    wf.write("Cosmic ray correction & %s \\\\\n" % yes_or_no[conf.flag_bpmask])
    wf.write("Cosmic ray threshold sigma & %.1f\\\\\n" % conf.CRthreshold)
    wf.write("Cosmic ray maximum sigma & %.1f\\\\\n" % conf.CRmaxsigma)
    wf.write("Cosmic ray Var/Ave ratio & %.1f\\\\\n" % conf.CRvaratio)
    wf.write("Cosmic ray ratio between slit positions & %.1f\\\\\n" % conf.CRslitposratio)
    wf.write("Cosmic ray fix sigma & %s\\\\\n" % yes_or_no[conf.CRfixsigma])

    wf.write("Sky emission spectra & %s \\\\\n" % yes_or_no[conf.flag_skyemission])
    wf.write("Measure pixel shifts & %s \\\\\n" % yes_or_no[conf.flag_wsmeasure])
    wf.write("Correct pixel shifts & %s \\\\\n" % yes_or_no[conf.flag_wscorrect])
    wf.write("Use the pixel shifts in list & %s \\\\\n" % yes_or_no[conf.flag_wsmanual])
    wf.write("Transform dy & %s \\\\\n" % conf.dyinput)
    wf.write("Transform flux & %s \\\\\n" % conf.fluxinput)
    wf.write("Cut ranges (x FSRs) & ")
    for i in range(len(conf.cutrange_list)):
        wf.write("%.2f" % conf.cutrange_list[i])
        if i != len(conf.cutrange_list) - 1:
            wf.write(",")
    wf.write("\\\\\n")
    wf.write(r"""\hline
\end{tabular}
\end{table}

%------------------------------------------------------------------------------------------------------------------------------------------------
""")

    # \section*{Free spectral ranges}
    #
    # \begin{table}[htb]
    # \centering
    # \begin{tabular}{cp{7em}p{7em}} \hline
    # Echelle order & $\lambda _\text{start}$ (\r{A}) & $\lambda _\text{end}$ (\r{A}) \\ \hline
    # """)

    #    for i in range(len(apnum)):
    #        wf.write("%d & %.1f & %.1f \\\\\n" % (apnum[i], fsr_angstrom[apnum[i]][0], fsr_angstrom[apnum[i]][1]))
    #    wf.write(r"""\hline
    # \end{tabular}
    # \end{table}

    wf.close()


def tex_calibration_data(texfile, conf):
    wf = open(texfile, "a")

    wf.write(r"""\section*{Calibration data}
        
\begin{table}[htb]
\centering
\begin{tabular}{p{14em}p{21em}} \hline
Calibration data & Name \\ \hline
        """)
    wf.write("Flat fielding image& %s\\\\\n" % conf.flat_file.replace("_", "\_"))
    yes_or_no = ["no", "yes"]
    wf.write("Bad pixel mask & %s \\\\\n" % conf.mask_file.replace("_", "\_"))
    wf.write("Comparison image & %s \\\\\n" % conf.comp_file.replace("_", "\_"))
    wf.write("Aperture trace image & %s \\\\\n" % conf.ap_file.replace("_", "\_"))
    wf.write("Aperture mask for apscatter & %s \\\\\n" % conf.apsc_maskfile.replace("_", "\_"))
    wf.write(r"""\hline
\end{tabular}
\end{table}

%------------------------------------------------------------------------------------------------------------------------------------------------
\clearpage
\newpage
""")

    wf.close()


############################################

def tex_spprofile(texfile, images_frames_dirs_sp, obj_sscfm_list, img_cs_list):
    wf = open(texfile, "a")

    wf.write(r"""
\section*{Spatial profiles}

The spatial profiles of the object in the slit. These profiles were obtained by integrating the 2D spectra along the dispersion direction.
The images and plotted data are stored in "object\_NO1/images/spatial\_profile/".

""")

    for i in range(len(images_frames_dirs_sp)):
        wf.write("\subsection*{No.%d}\n" % (i + 1))
        wf.write(r"""\begin{figure}[!h]""")
        wf.write(r"""\includegraphics[width=9.3cm]""")
        wf.write("{%s%s}\n" % (images_frames_dirs_sp[i], img_cs_list[i][int(len(img_cs_list[i])/2)]))
        wf.write(r"""\end{figure}

             
""")
        if i % 2 == 1:
            wf.write(r"""\clearpage
\newpage
            
""")


############################################

def tex_spec_images(texfile, spectra_image, obj_comb_norm_png):
    wf = open(texfile, "a")

    wf.write(r"""
\section*{Spectra}

Continuum normalized spectra. The plotted data is stored in "object\_sum/AIR\_norm/".

""")

    for i in range(len(obj_comb_norm_png)):
        wf.write(r"""\begin{figure}[!h]
\includegraphics[width=15.5cm]""")
        wf.write("{%s/%s}\n" % (spectra_image, obj_comb_norm_png[len(obj_comb_norm_png) - i - 1]))
        wf.write(r"""\end{figure}
""")

        if (i + 1) % 4 == 0:
            wf.write(r"""\newpage \
""")

    wf.write(r"""%-----------------------------------------------------------------------------------------------------------------------------------------------
""")


############################################

def tex_cosmicrays(texfile, warpLog):
    wf = open(texfile, "a")

    wf.write(r"""
    \section*{Cosmic ray detection}

    \begin{table}[htb]
    \centering
    \begin{tabular}{p{6em}|p{6em}p{6em}p{6em}} \hline
    Frame & Cosmic ray (pix) & Threshold (sigma) \\ \hline
    """)

    for i in range(warpLog.frameNum):
        wf.write("No.%d & %.d & %.1f  \\\\\n" % ((i + 1), warpLog.cosmicRayNum[i], warpLog.cosmicRaySigma[i]))

    wf.write(r"""\hline
    \end{tabular}
    \end{table}
    """)

    wf.write(r"""%-----------------------------------------------------------------------------------------------------------------------------------------------
    """)

    wf.close()


############################################

def tex_apertures(texfile, warpLog):
    wf = open(texfile, "a")

    wf.write(r"""
\section*{Aperture range and Linear shift of spectra}

\begin{table}[htb]
\centering
\begin{tabular}{p{6em}|p{6em}p{6em}p{6em}|p{11em}p{11em}} \hline
Frame & From (pix) & To (pix) & Width (pix) & Measured shift (pix) & Corrected shift (pix) \\ \hline
""")

    if warpLog.frameNum > 1:
        for i in range(warpLog.frameNum):
            wf.write("No.%d & %.2f & %.2f & %.2f  & %.2f$\pm$%.2f ($N=$%d) & %.2f \\\\\n" % (
                (i + 1), warpLog.apertureLow[i][0], warpLog.apertureUpp[i][0],
                warpLog.apertureUpp[i][0] - warpLog.apertureLow[i][0], warpLog.waveShiftAve[i],
                warpLog.waveShiftStd[i], warpLog.waveShiftNum[i], warpLog.waveShiftAdopted[i]))
    else:
        wf.write("No.%d & %.2f & %.2f & %.2f  & %.2f$\pm$%.2f ($N=$%d) & %.2f \\\\\n" % (
            (1), warpLog.apertureLow[0][0], warpLog.apertureUpp[0][0],
            warpLog.apertureUpp[0][0] - warpLog.apertureLow[0][0], 0.,
            0., 0., 0.))

    wf.write(r"""\hline
\end{tabular}
\end{table}
""")

    wf.write(r"""%-----------------------------------------------------------------------------------------------------------------------------------------------
""")

    wf.close()


############################################


def tex_psf_center_width(texfile, warpLog, fsr_angstrom):
    wf = open(texfile, "a")

    maxfrms = 5

    leftfnum = warpLog.frameNum
    repeatnum = 0

    wf.write(r"""\section*{Central positions and full widths at half maximum}

The central positions and FWHMs of the object in the slit. 
The values are written in the slit coordinate, which is defined from the calibration file. 
The free spectral ranges are also shown.

""")

    while leftfnum > 0:
        wf.write(r"""\begin{table}[htb]
\centering
\scriptsize
\begin{tabular}{""")

        wf.write("ccc")
        for i in range(min(leftfnum, maxfrms)):
            wf.write("cc")
        wf.write("} \\hline \n")

        wf.write(r""" & Free spectral range &  """)
        for i in range(min(leftfnum, maxfrms)):
            wf.write("& No.%d & " % (i + 1 + repeatnum * maxfrms))
        wf.write(" \\\\ \n")
        wf.write(r"""$m$ & $\lambda _\text{start}$ (\r{A}) & $\lambda _\text{end}$ (\r{A}) """)
        for i in range(min(leftfnum, maxfrms)):
            wf.write("& $x_c$ & FWHM ")
        wf.write(" \\\\ \\hline \n")
        for i in range(warpLog.echelleOrderNum):
            wf.write("%d & %.1f & %.1f " % (warpLog.echelleOrderVector[i], fsr_angstrom[warpLog.echelleOrderVector[i]][0], fsr_angstrom[warpLog.echelleOrderVector[i]][1]))
            for j in range(min(leftfnum, maxfrms)):
                wf.write("& %.2f & %.2f " % (warpLog.psfCenter[j + repeatnum * maxfrms][i], warpLog.psfWidth[j + repeatnum * maxfrms][i]))
            wf.write("\\\\ \n")

        wf.write(r"""\hline
Median: """)
        wf.write(" & & ")
        for i in range(min(leftfnum, maxfrms)):
            wf.write("& %.2f & %.2f " % (
                np.nanmedian(warpLog.psfCenter[i + repeatnum * maxfrms]), np.nanmedian(warpLog.psfWidth[i + repeatnum * maxfrms])))

        wf.write(r""" \\ \hline
\end{tabular}
\end{table}
""")

        repeatnum += 1
        leftfnum -= maxfrms

        # \subsubsection*{Full width at half maximum}
    #
    # \begin{table}[htb]
    # \centering
    # \scriptsize
    # \begin{tabular}{c""")
    #
    #    for i in range(objnum):
    #        wf.write("c")
    #    wf.write("} \\hline \n")
    #
    #    wf.write("Order ")
    #    for i in range(objnum):
    #        wf.write("& No.%d " % (i+1))
    #    wf.write(" \\\\ \\hline \n")
    #    for i in range(len(apnum)):
    #        wf.write("%d " % apnum[i])
    #        for j in range(objnum):
    #            wf.write("& %.2f " % gw[j][i])
    #        wf.write("\\\\ \n")
    #
    #    wf.write(r"""\hline
    # Median: """)
    #    for i in range(objnum):
    #        wf.write("& %.2f" % numpy.median(gw[i]))
    #
    #    wf.write(r""" \\ \hline
    # \end{tabular}
    # \end{table}
    # """)
    wf.write(r"""%-----------------------------------------------------------------------------------------------------------------------------------------------
\clearpage
\newpage
""")

    wf.close()


############################################


def tex_snr(texfile, SNRdat_frames_dirs, SN_png):
    wf = open(texfile, "a")

    wf.write(r"""
\section*{Signal to noise ratio}


\textbf{S/N} indicated by the blue line graph is
based on measurements of the consistency between
two adjasent spectra; differences between two normalized spectra
are measured and those from all the pairs are combined to
calculate the S/N values (per pixel) around
three sections of each order.
The red curve indicates \textbf{the combined spectra from the all orders 
\underline{before the continuum normalization}}. The flux scale is uncalibrated.
The plotted data is stored in "object\_sum/AIR\_flux/".

""")

    wf.write(r"""\begin{figure}[!h]
\centering
\includegraphics[width=15cm]""")
    wf.write("{%s/%s}\n" % (SNRdat_frames_dirs, SN_png))
    wf.write(r"""\end{figure}

\newpage


""")

    wf.close()


def tex_fwhm_count(texfile, logdir, FWHMpng):
    wf = open(texfile, "a")

    wf.write(r"""
\section*{Counts and FWHMs}

The average counts of the spectrum (red circles) and FWHMs of the object in the slit profile (blue circles). 
The echelle order used in the calculation of average counts is indicated in the left axis. 

""")

    wf.write(r"""\begin{figure}[!h]
\centering
\includegraphics[width=15cm]""")
    wf.write("{%s/%s}\n" % (logdir, FWHMpng))
    wf.write(r"""\end{figure}

\clearpage
\newpage
""")

    wf.close()


def tex_slitviewer(texfile, svdir, imlist_sort):
    wf = open(texfile, "a")
    wf.write(r"""
\section*{Slit viewer images}

The left and right images show the slit viewer images taken at the starting and ending time, respectively, of each exposure. 
The png and fits images are stored in "slit\_viewer/".

\vspace{20pt}

""")
    for i in range(len(imlist_sort)):
        if i % 2 == 0 and i != 0:
            wf.write(r"""
\newpage
""")
        wf.write("\\noindent %s\n\n" % imlist_sort[i].replace("_", "\_"))

        wf.write(r"""\begin{figure}[!h]
        \includegraphics[width=8.5cm]""")
        wf.write("{%s/%s}\n" % (svdir, imlist_sort[i] + "_expstart.png"))
        wf.write("""\includegraphics[width=8.5cm]""")
        wf.write("{%s/%s}\n" % (svdir, imlist_sort[i] + "_expend.png"))
        wf.write(r"""\end{figure}

""")

    wf.close()


def tex_closing(texfile):
    wf = open(texfile, "a")
    wf.write("\end{document}")
    wf.close()


def tex_source_make(conf: config, fsr, logo, redpath="./"):
    # read fits headers
    hdulist = []
    for i in range(conf.imnum):
        hdulist_obj = fits.open(redpath + "rawdata_image/" + conf.imagelist[i] + ".fits")
        hdulist.append(hdulist_obj[0].header)
        hdulist_obj.close()

    conf.readObservationInfo(hdulist)
    texfile = "%s_%s.tex" % (conf.objnameRep, conf.acqdate[0])
    tex_cover_and_obsinfo(texfile, conf)
    pipeline_ver = Warp_sci.__version__

    conf.readParamFile("reduction_log/status.txt")
    os.chdir("calibration_data/")
    conf.readInputCalib("input_files.txt")
    apset = apertureSet(conf.ap_file, reduceFullData=conf.reduceFullData, selectedOrders=conf.selectedOrders)
    cutlength = len(conf.cutrange_list)
    aplength = len(apset.echelleOrders)
    os.chdir("../")
    log = warpLog(apset.echelleOrders, conf.objnum)
    tex_pipeline_para(texfile, pipeline_ver, conf)
    tex_calibration_data(texfile, conf)
    if conf.flag_bpmask:
        log.readCosmicRayLogNpz("reduction_log/cosmicray_log.npz")
        tex_cosmicrays(texfile, log)
    log.readApertureLogNpz("reduction_log/aperture_log.npz")
    if conf.objnum > 1:
        log.readWaveshiftLogNpz("reduction_log/waveshift_log.npz")
    tex_apertures(texfile, log)

    if not conf.flag_manual_aperture:
        log.readPsfLogNpz("reduction_log/centersearch_log.npz")
        tex_psf_center_width(texfile, log, fsr)

    obj_comb_norm_png = [conf.objnameRep + "_AIRnorm_combined_m%d.png" % apset.echelleOrders[j] for j in range(aplength)]
    SNRdat_frames_dirs = ["SNR_dat/fsr%.2f" % (conf.cutrange_list[k]) for k in range(cutlength)]
    SN_png = ["SNratio_fsr%.2f.png" % conf.cutrange_list[k] for k in range(cutlength)]

    tex_spec_images(texfile, "spectra_image", obj_comb_norm_png)
    if 1 < conf.objnum <= conf.frameNumberLimit:
        tex_snr(texfile, SNRdat_frames_dirs[0], SN_png[0])
    tex_fwhm_count(texfile, "reduction_log", "count_and_fwhm.png")

    images_frames_dirs_sp = ["%s_NO%d/images/spatial_profile/" % (conf.objname_obj[i], (i + 1)) for i in range(conf.objnum)]

    obj_s_list = [conf.objname_obj[i] + "_NO%d_s" % (i + 1) for i in
                  range(conf.objnum)]  # FITS file after sky subtruction: "star_NO1_s.fits"
    if conf.flag_apscatter:
        obj_ssc_list = [obj_s_list[i] + "sc" for i in
                        range(conf.objnum)]  # FITS file after scattered light subtruction: "star_NO1_ssc.fits"
    else:
        obj_ssc_list = obj_s_list
    obj_sscf_list = [obj_ssc_list[i] + "f" for i in
                     range(conf.objnum)]  # FITS file after flat fielding: "star*_NO1_s(sc)f.fits"
    obj_sscfm_list = [obj_sscf_list[i] + "m" for i in range(conf.objnum)]  # FITS file after fixpix: "star_NO1_s(sc)fm.fits"

    img_cs_list = [[obj_sscfm_list[i] + "_m%dtrans.png" % apset.echelleOrders[j] for j in range(aplength)] for i in range(conf.objnum)]

    tex_spprofile(texfile, images_frames_dirs_sp, obj_sscfm_list, img_cs_list)
    if conf.flag_svimage:
        tex_slitviewer(texfile, "slit_viewer", conf.imagelist)
    tex_closing(texfile)

    if not os.path.exists("./winered_logo.pdf"):
        shutil.copy(logo, ".")

    try:
        for _ in range(2):
            subprocess.run("pdflatex %s" % texfile, shell=True, check=True, timeout=10)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        print(e)
        print("WARNING: Reduction report could not be generated by a Tex compile error. "
              "Reduction processes were successfully finished.")

    shutil.move(texfile, "reduction_log/" + texfile)
    shutil.move(texfile.rstrip("tex") + "aux", "reduction_log/" + texfile.rstrip("tex") + "aux")
    shutil.move(texfile.rstrip("tex") + "log", "reduction_log/" + texfile.rstrip("tex") + "log")
    shutil.move("winered_logo.pdf", "reduction_log")


if __name__ == "__main__":
    listfile = sys.argv[1]

    rfile = open(listfile, "r")
    rlines = rfile.readlines()
    rfile.close()

    objectlist = [i.split()[0] for i in rlines]
    skylist = [i.split()[1] for i in rlines]

    imagelist = list(set(objectlist + skylist))
    imlist_sort = sorted(imagelist)

    tex_source_make(imlist_sort, objectlist, skylist)
