############################################################################
# Computes radial loads for the ring problem                               #   
#                                                                          #
#  (Last updated December 2015)                                            #
############################################################################
"""
Created on Tue Dec 22 18:16:19 2015

@author: eafit
"""
import math as mt
import numpy as np
#
nodes    = np.loadtxt('nodloads.txt')
nn =len(nodes[:,0])
F=1.0
Fl=np.zeros([nn,2],float)

COORD  =np.zeros([nn,2],dtype=np.float)
COORD[:,0]=nodes[:,1]
COORD[:,1]=nodes[:,2]
for i in range(0,nn):
    x=COORD[i,0]
    y=COORD[i,1]
    theta=mt.atan2(y,x)
    Fl[i,0]=F*mt.cos(theta)
    Fl[i,1]=F*mt.sin(theta)
np.savetxt("KNODES.txt", Fl,    fmt='%5.2f', delimiter=' ')