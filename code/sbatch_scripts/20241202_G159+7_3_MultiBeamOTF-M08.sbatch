#!/bin/bash
#SBATCH --job-name=20241202_G159+7_3_MultiBeamOTF-M08
#SBATCH --output=logs/20241202_G159+7_3_MultiBeamOTF-M08_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/G159+7_3/20241202/G159+7_3_MultiBeamOTF-M08_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241202_G159+7_3_MultiBeamOTF-M08.log