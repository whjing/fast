#!/bin/bash
#SBATCH --job-name=20241209_G159+7_5_MultiBeamOTF-M07
#SBATCH --output=logs/20241209_G159+7_5_MultiBeamOTF-M07_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/G159+7_5/20241209/G159+7_5_MultiBeamOTF-M07_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241209_G159+7_5_MultiBeamOTF-M07.log