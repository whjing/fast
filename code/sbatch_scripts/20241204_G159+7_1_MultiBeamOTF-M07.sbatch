#!/bin/bash
#SBATCH --job-name=20241204_G159+7_1_MultiBeamOTF-M07
#SBATCH --output=logs/20241204_G159+7_1_MultiBeamOTF-M07_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/G159+7_1/20241204/G159+7_1_MultiBeamOTF-M07_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241204_G159+7_1_MultiBeamOTF-M07.log