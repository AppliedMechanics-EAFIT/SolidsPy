# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 18:43:29 2016

@author: casierraa
"""
import numpy as np
import matplotlib.pyplot as plt
#
f=np.cos
fx = lambda x: 1.0/(2.0*np.cos(x)*np.sin(x));
xx = np.linspace(10, 80, 100)
zz = np.zeros((100))
for i in range(100):
    xxx=xx[i]*np.pi/180.0
    zz[i]=fx(xxx)
plt.plot(xx,zz)
plt.ylim([0.0,2.0])
