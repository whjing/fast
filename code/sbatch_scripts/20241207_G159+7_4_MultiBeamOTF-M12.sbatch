#!/bin/bash
#SBATCH --job-name=20241207_G159+7_4_MultiBeamOTF-M12
#SBATCH --output=logs/20241207_G159+7_4_MultiBeamOTF-M12_%j.out

python -c "from loadFits import load_fits; load_fits('../data_orig/G159+7_4/20241207/G159+7_4_MultiBeamOTF-M12_W_0001.fits', '../data_ave-64', avechans=64)" > logs/20241207_G159+7_4_MultiBeamOTF-M12.log