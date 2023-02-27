#!/bin/sh

python Warp_sci.py ./TEST/WIDE/4_Ari_list.txt -r ./TEST/WIDE/ -c ./TEST/WIDE/WINERED_calibration_LCO22b_WIDE100_20220914_v2/ -v ./TEST/WIDE/ -d ./TEST/4_Ari_test
python Warp_sci.py ./TEST/HIRES-J/21_Peg_list.txt -r ./TEST/HIRES-J/ -c ./TEST/HIRES-J/WINERED_calibration_NTTrun_HIRESJ_20170728/ -v ./TEST/HIRES-J/ -d ./TEST/21_Peg_test
python Warp_sci.py ./TEST/HIRES-Y/HD_163336_list.txt -r ./TEST/HIRES-Y/ -c ./TEST/HIRES-Y/WINERED_calibration_NTTrun_HIRESY_20170727/ -v ./TEST/HIRES-Y/ -d ./TEST/HD_163336_test
