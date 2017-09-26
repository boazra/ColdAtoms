# -*- coding: utf-8 -*-
"""
Created on Sun May 28 16:18:37 2017

@author: Boaz
"""

import numpy as np
from matplotlib.pylab import plot,show,figure,legend,rc,xlabel,ylabel,title

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 26}
rc('font', **font)

path = 'C:/Users/Boaz/Documents/Nir Lab/battery exp/0.3Hz/scope_7.csv'
data = np.loadtxt(path, delimiter=",")
data[:,0] = data[:,0] - min(data[:,0])
plot(data[:,0],data[:,1]*10,data[:,0],data[:,2]*-100)
title('Battery current VS Time')
xlabel('Time (s)')
ylabel('Current(A)')
show()
figure()
NoiseSpecturm = np.fft.rfft(data[:,2])

dt = (data[-1,0] - data[0,0] )/len(data)
df = 1/dt
Freq = range(len(NoiseSpecturm))*df / 1.0e6
plot(Freq,abs(NoiseSpecturm))
title('Nosie Frequency spectrum')
xlabel('Frequency (MHz)')
ylabel('Amplitude (a.u)')
show()