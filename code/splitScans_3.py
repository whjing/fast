#!/usr/bin/env python

################################################################################
#  This is the code to calibrate polarization data observed with the FAST      #
#  telescope                                                                   # 
#      The code:                                                               #
#         - measures calibration parameters from the input calibration scans   #
#         - calibrate ALL rpf files in the current directory                   #
#         - bin the data in df-MHz wide bands                                  #
#         - write sdfits output file.                                          #
#                                                                              #
#      It assumes the following scans are available:                           #
#         - an unpolarised flux calibrator (1934, 0407, and HydA only for now) #
#         - a polarised calibrator (0043-424 and 3C138 only for now).          #
################################################################################


from astropy.io import fits as pf
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
import numpy as np
from openpyxl import load_workbook
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, ITRS, EarthLocation, AltAz
from astropy.stats import sigma_clip
from scipy import interpolate
import csv
import sys
import os
import glob
import re as RegEx
from astropy import wcs
from scipy.interpolate import griddata
from scipy import ndimage
import gc


def splitData(infList, scanList, oufDir, oufTemplate):
    """
    input file list infList, scan list scanList, oufDir, oufTemplate
    
    split scans based on the MJD in scanList
    """

    inf = infList[0]
    print(f'reading file {inf}')
    hdul = pf.open(inf)
    nrows = hdul[1].data.shape[0]

    Nrows_total = nrows 

    for i in range(1, len(infList)):
        inf_1 = infList[i]
        print(f'reading file {inf_1}')
        hdul_1 = pf.open(inf_1)
        nrows_1 = hdul_1[1].data.shape[0]
        hdu = pf.BinTableHDU.from_columns(hdul[1].columns, nrows=Nrows_total+nrows_1)
        for colname in hdul[1].columns.names:
            hdu.data[colname][Nrows_total:] = hdul_1[1].data[colname]
        hdul[1] = hdu
        Nrows_total += nrows_1

    tmpfits = f'{oufDir}{inf.split("/")[-1][:-12]}.fits'#.
    hdul.writeto(tmpfits, overwrite=True)
    print(f'{tmpfits} has been saved!')
    hdul.close()

    print("We have read all files. Now is time to generate new files! ")

    for i in range(len(scanList)):

        hdul = pf.open(tmpfits) 
        hdu = hdul[1]

        mask1 = hdu.data['UTOBS'] >= scanList[i][0]
        mask2 = hdu.data['UTOBS'] <= scanList[i][1]
        mask = mask1 & mask2
        hdu.data = hdu.data[mask]

        if hdu.data.shape[0] > 0:
            hdul_n = pf.HDUList([hdul[0], hdu])
            outFile = oufDir + oufTemplate.replace('XXX',str(i+1))
            hdul_n.writeto(outFile, overwrite=True)
            print(f'The file {outFile} has been saved')
        
        else:
            print(f'hdu.data.shape[0]{hdu.data.shape[0]} <= 0')

        hdul.close()



def main():

    dataDir = "/share/home/whjing/fast/orig/G182_3/"
    dateList = ['20221205']
    beamList = np.arange(1,20)

    scanList = np.arange(1,6)
    print(f'the scan list is {scanList}')
    outputDir = '/share/home/whjing/fast/file_redo-day-3/'
    scan_pts_file = outputDir + 'G182_3_scan_points.dat'
    f = open(scan_pts_file,'w')
    for scan in scanList:
        print(f'Now we are working on scan {scan}')
        f1 = open(outputDir +'../file_redo-day-3/G182_3_2022_12_05_00_58_R%d.csv'%scan)
        print(f'Now we are reading {f1}')
        allLines = f1.readlines()
        f1.close()
        s1 = allLines[1].split(',')[0]
        s2 = allLines[-1].split(',')[0]
        f.write(s1+'  '+s2+'\n')
    f.close()


    
    f = open(scan_pts_file, "r")
    allLines = f.readlines()
    f.close()

    scanList = []
    for term in allLines:
        x = term.split()
        scanList.append([float(x[0]), float(x[1])])

    oufDir = outputDir + 'file_redo-day-3/'
    oufTemplate = 'G182_RXXX-MYYY_0001.fits'

    for date in dateList:
        infDir = dataDir +date+'/' 
        if not os.path.exists(oufDir): os.mkdir(oufDir)

        for beam in beamList:
            infList = glob.glob(infDir + '*M%02d*.fits' % beam)
            infList.sort()

            splitData(infList, scanList, oufDir, oufTemplate.replace('YYY','%02d'% beam))
            
if __name__ == "__main__":
    main()

