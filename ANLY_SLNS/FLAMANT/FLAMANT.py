# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:36:56 2015

@author: eafit
"""
from mpl_toolkits.mplot3d import Axes3D 

import mis_funciones as mf
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

plt.close("all")

##############################################################################
# Flamant solution corresponding to a load and moment in a half-space. In this
# particular case we are applying two line loads representing the footings
# of a building.
#
# Input data

ly=10.#...................base of the computational domain (dist from -yi to yf)
lx=10.#...................Height of the computational domain (dist from 0  to xf)
vf=[1.,1] #...............Surface loads vector
vbeta=[90.,90.] #.........Vector  that stores the angle in which the forces are applied
vm=[1.,1.] #............. Surface moment vector
vcy=[0,1]
##############################################################################
nx=100 #................................Number of data points along the x-axis
ny=nx #.................................Number of data points along the y-axis

##############################################################################

ejex= np.linspace(0., lx, nx)#............ X-axis 
ejey= np.linspace(-ly/2., ly/2., ny)#..... Y-axis

x, y = np.meshgrid(ejex, ejey) #.......... x and y coordinates matrix
sp=np.zeros((2,2))
mquiz=np.zeros((ny,nx))
msx=np.zeros((ny,nx))
msy=np.zeros((ny,nx))
mtxy=np.zeros((ny,nx))
msmin=np.zeros((ny,nx))
msmax=np.zeros((ny,nx))
mtmax=np.zeros((ny,nx))

for k in range (len(vf)): 
    for i in range (len(ejey)):
        for j in range (len(ejex)):
            xp=x[i,j]
            yp=y[i,j]-vcy[k]
            rp=(xp**2.+yp**2.)**0.5
            tetap=mf.grados(np.arcsin(yp/rp)) 
            sc=mf.tensor_cart_m(rp,tetap,vf[k],vm[k],vbeta[k])
            msx[i,j]=sc[0,0]+msx[i,j]
            msy[i,j]=sc[1,1]+msy[i,j]
            mtxy[i,j]=sc[0,1]+mtxy[i,j]
# superimposes the maximum values
for i in range (len(ejey)):
   for j in range (len(ejex)):
        # Cartesian tensor superimposed        
        sp[0,0]=msx[i,j] 
        sp[1,1]=msy[i,j]
        sp[0,1]=mtxy[i,j]
        sp[1,1]=mtxy[i,j]
        # find the maximum values
        la,ve = eigh(sp)
        msmin[i,j]=la[0]
        msmax[i,j]=la[1]
        mtmax[i,j]=abs(la[0]-la[1])/2.
        
        
levels=np.linspace(-5, 5, 20)
############################################################################3
plt.contourf(y,-x,mtmax ,levels, alpha=.75, cmap='jet')
plt.colorbar()
plt.grid()
#


## plots the building
px=np.zeros(2)
py=np.zeros(2)
for i in range (len(vf)):
    px[0]=vcy[i]
    px[1]=vcy[i]    
    py[0]=0
    py[1]=0.3
    plt.plot(px,py)
#
#
###levels=np.linspace(-5, 0, 5)
##CS4 = plt.contour(y, -x, msx, levels, colors=('k',),linewidths=(3,))
##plt.clabel(CS4, fmt='%2.1f', colors='w', fontsize=14)
###
##plt.title('Con todo')
#
##
##figura_3d = plt.figure()
##ax = figura_3d.gca(projection = '3d')
##ax.plot_surface(x, y, mtmax, rstride=1, cstride=1, cmap='coolwarm',linewidth=0, antialiased=False)
###ax.set_zlim(-1, 1)
##
#























































