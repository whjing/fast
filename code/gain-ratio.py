#!/usr/bin/env python
#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from scipy.signal import medfilt


# inf1 = '3c48_20210813.cal.dat'
# inf2 = '3c48_20210815.cal.dat'
# inf3 = '3c48_20210818.cal.dat'
# outname = 'gain_ratio_R.txt'
dir = '/share/home/whjing/fast/3c138_cal_0725'
inf1 = '/share/home/whjing/fast/3c138_cal_0725/3c138_20221201_all.cal.dat'
outname = '/share/home/whjing/fast/3c138_cal_0725/gain_ratio_R.txt'

Ncol = 13
data1 = np.loadtxt(inf1)#, sep='\n')
print(f'orignal data shape {data1.shape}')
# data1 = data1.reshape((data1.size//Ncol, Ncol))
# print(f'reshaped data shape {data1.shape}')
#data2 = np.fromfile(inf2, sep='\n')
#data2 = data2.reshape((data2.size//Ncol, Ncol))
#data3 = np.fromfile(inf3, sep='\n')
#data3 = data3.reshape((data3.size//Ncol, Ncol))
#data4 = np.fromfile(inf4, sep='\n')
#data4 = data4.reshape((data4.size//Ncol, Ncol))


freq = data1[:,-1][::19]
print(freq[1])

print(f'freq data shape {freq.shape}')
gain1 = 1./data1[:,2]
print(f'gain1 data shape {gain1.shape}')
#gain2 = 1./data2[:,2]
#gain3 = 1./data3[:,2]
#gain4 = 1./data4[:,2]

N_beams = 19

ouf = open(outname, 'w')
ouf.write('1, 0\n')

# y0 = (gain1[0::19] + gain2[0::19] + gain3[0::19]) / 3.
y0 = gain1[0::19]
print(f'y0 shape {y0.shape}')

fig = plt.figure(figsize = (10,20))


for beam in range(1, N_beams):
    beam_name = beam+1
    y1 = gain1[beam::19]
    # y = y0 / y1 
    y = y1 / y0 
    

    p = models.Polynomial1D(1)
    pfit = fitting.LinearLSQFitter()

    
    model = pfit(p, freq, y)
    yf = model(freq)
    ss = '%.2f, %.2f' % (model.c0.value, model.c1.value)
    ouf.write(ss+'\n')

    ymedfit = medfilt(y, kernel_size= 99)


    # pos = f'0306{beam}'
    ax = fig.add_subplot(9,2,beam)
    ax.plot(freq, y, '.', color = 'lightblue')
    ax.plot(freq, yf, label = f'{ss}')
    ax.plot(freq, ymedfit, color = 'magenta', alpha = 0.5, label = 'medfit')
    ax.set_title(f'Beam {beam_name}')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel(r'$G_{B}/G_{B1}$')
    plt.tight_layout()
    plt.legend()
    plt.savefig(f'{dir}/gain-ratio_rev.png', dpi=100)
    

# plt.close()


# %%
