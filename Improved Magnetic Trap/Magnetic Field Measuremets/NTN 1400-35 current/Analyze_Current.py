# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:12:07 2017

@author: Boaz
"""
import os
import matplotlib.pylab as plt

path = r'C:\Users\Boaz\Documents\Nir Lab\Improved Magnetic Trap\battery exp'
files = os.listdir(path)

font = {'family' : 'DejaVu Sans',        
        'size'   : 26}
plt.rc('font', **font)

for file in files:
    if not '.csv' in file:
        continue
    plt.close('all')
    plt.figure('Current Response',(20,10))
    data = np.genfromtxt(os.path.join(path,file), delimiter=',')[2:,:]
    ratio = data[:,2].ptp()/data[:,1].ptp()
    plt.plot(data[:,0], data[:,1],data[:,0], data[:,2]/ratio)
    plt.xlabel('Time(s)')
    plt.legend(['Current Output', 'Signal Voltage'], loc=2)
    plt.savefig(os.path.join(path,'{}.jpg'.format(file[:-4])))
    print('finished file : {}'.format(file))    