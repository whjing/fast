#!/usr/bin/env python
#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from scipy.signal import medfilt


def obtain_gainratio(file_in, file_out, plot=True):


    Ncol = 13
    data1 = np.loadtxt(file_in)#, sep='\n')
    print(f'orignal data shape {data1.shape}')

    freq = data1[:,-1][::19]
    print(freq[1])

    print(f'freq data shape {freq.shape}')
    gain1 = 1./data1[:,2]
    print(f'gain1 data shape {gain1.shape}')
    N_beams = 19

    ouf = open(file_out, 'w')
    ouf.write('1, 0\n')

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

        ax = fig.add_subplot(9,2,beam)
        ax.plot(freq, y, '.', color = 'lightblue')
        ax.plot(freq, yf, label = f'{ss}')
        ax.plot(freq, ymedfit, color = 'magenta', alpha = 0.5, label = 'medfit')
        ax.set_title(f'Beam {beam_name}')
        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel(r'$G_{B}/G_{B1}$')
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"{file_out}.png", dpi=100)


dir = '../data_ave-64/3C138_20241207/'
file_in = f'{dir}/gain_ratio1.cal.dat'
file_out = f'{dir}/gain_ratio.dat'
obtain_gainratio(file_in, file_out, plot=True)
# plt.close()


# %%
