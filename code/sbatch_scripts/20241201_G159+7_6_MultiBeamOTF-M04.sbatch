#!/bin/bash
#SBATCH --job-name=20241201_G159+7_6_MultiBeamOTF-M04
#SBATCH --output=logs/20241201_G159+7_6_MultiBeamOTF-M04_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/G159+7_6/20241201/G159+7_6_MultiBeamOTF-M04_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241201_G159+7_6_MultiBeamOTF-M04.log