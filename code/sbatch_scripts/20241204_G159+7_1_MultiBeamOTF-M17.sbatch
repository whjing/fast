#!/bin/bash
#SBATCH --job-name=20241204_G159+7_1_MultiBeamOTF-M17
#SBATCH --output=logs/20241204_G159+7_1_MultiBeamOTF-M17_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/G159+7_1/20241204/G159+7_1_MultiBeamOTF-M17_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241204_G159+7_1_MultiBeamOTF-M17.log