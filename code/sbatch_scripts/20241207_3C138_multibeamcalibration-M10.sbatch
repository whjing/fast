#!/bin/bash
#SBATCH --job-name=20241207_3C138_multibeamcalibration-M10
#SBATCH --output=logs/20241207_3C138_multibeamcalibration-M10_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/3C138/20241207/3C138_multibeamcalibration-M10_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241207_3C138_multibeamcalibration-M10.log