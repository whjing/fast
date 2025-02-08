#!/bin/bash
#SBATCH --job-name=20241202_3C138_multibeamcalibration-M08
#SBATCH --output=logs/20241202_3C138_multibeamcalibration-M08_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/3C138/20241202/3C138_multibeamcalibration-M08_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241202_3C138_multibeamcalibration-M08.log