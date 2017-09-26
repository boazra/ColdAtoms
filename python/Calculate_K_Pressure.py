# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 13:28:46 2017

@author: Boaz
"""

from numpy import pi,log

lambda0 = 766.7e-9
sigma0 = 3*lambda0**2/2.0/pi
Isat = 1.75
Ibeam = 6**2 * Isat
sigma = sigma0 * 1/(1+(Ibeam/Isat))
L = 13e-2
V = L
I0 = 200.0
I = 137.0
OD = -log(I/I0)
n = OD/V/sigma
k = 1.3806e-23
T = 300
P = n*k*T
P_torr = P * 0.00750062
print(P_torr)