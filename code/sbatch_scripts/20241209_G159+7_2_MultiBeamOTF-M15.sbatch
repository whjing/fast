#!/bin/bash
#SBATCH --job-name=20241209_G159+7_2_MultiBeamOTF-M15
#SBATCH --output=logs/20241209_G159+7_2_MultiBeamOTF-M15_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/G159+7_2/20241209/G159+7_2_MultiBeamOTF-M15_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241209_G159+7_2_MultiBeamOTF-M15.log