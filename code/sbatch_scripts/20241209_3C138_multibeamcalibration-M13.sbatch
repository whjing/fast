#!/bin/bash
#SBATCH --job-name=20241209_3C138_multibeamcalibration-M13
#SBATCH --output=logs/20241209_3C138_multibeamcalibration-M13_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/3C138/20241209/3C138_multibeamcalibration-M13_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241209_3C138_multibeamcalibration-M13.log