#!/bin/bash
#SBATCH --job-name=20241201_3C138_multibeamcalibration-M05
#SBATCH --output=logs/20241201_3C138_multibeamcalibration-M05_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/3C138/20241201/3C138_multibeamcalibration-M05_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241201_3C138_multibeamcalibration-M05.log