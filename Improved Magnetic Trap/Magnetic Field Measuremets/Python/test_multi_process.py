# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 13:29:49 2017

@author: nirlab
"""

from multiprocessing import Pool

def func(a):
    return a**2

p = Pool(5)
x = p.map(func,range(100))