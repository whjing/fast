#!/bin/bash
#SBATCH --job-name=20241204_3C138_multibeamcalibration-M07
#SBATCH --output=logs/20241204_3C138_multibeamcalibration-M07_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/3C138/20241204/3C138_multibeamcalibration-M07_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241204_3C138_multibeamcalibration-M07.log