# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:36:56 2015

@author: eafit
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import mis_funciones as mf
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from sympy import *
from sympy import init_printing
init_printing()
#
#Computes the solution for a symmetric wedge subject to constant surface shear
# see Principles of solid mechanics. Rowland Richards, Jr. CRC Press.
#
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['image.cmap'] = "jet"
#
nodes    = np.loadtxt('nodes.txt')
nn =len(nodes[:,0])
coords=np.zeros([nn,2])
UX=np.zeros([nn])
UY=np.zeros([nn])
coords[:,0]=nodes[:,1]
coords[:,1]=nodes[:,2]
#
l=np.sqrt(8)
phi=np.pi/4
E=1.0
nu=0.30
S=0.45
#
for i in range(0,nn):
    x=coords[i,0]
    y=coords[i,1]
    ux, uy =mf.cunia(x,y,phi,l,nu,E,S)
    UX[i]=ux
    UY[i]=uy
#
grid_x, grid_y = np.mgrid[0:4:100j, -2:2:100j]
grid_z0 = griddata(coords, UX, (grid_x, grid_y), method='linear')
grid_z1 = griddata(coords, UY, (grid_x, grid_y), method='linear')
#
plt.imshow(grid_z0.T,aspect='equal', extent=(0,4,-2,2), origin='lower')
plt.title('Horizontal field')
plt.colorbar(orientation='vertical')
plt.grid()
plt.savefig('cuniahorizo.pdf')
plt.show()
#
plt.imshow(grid_z1.T,aspect='equal', extent=(0,4,-2,2), origin='lower')
plt.title('Vertical field')
plt.colorbar(orientation='vertical')
plt.grid()
plt.savefig('cuniavertic.pdf')
plt.show()
#