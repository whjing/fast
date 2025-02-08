#!/bin/bash
#SBATCH --job-name=20241210_3C138_multibeamcalibration-M03
#SBATCH --output=logs/20241210_3C138_multibeamcalibration-M03_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/3C138/20241210/3C138_multibeamcalibration-M03_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241210_3C138_multibeamcalibration-M03.log