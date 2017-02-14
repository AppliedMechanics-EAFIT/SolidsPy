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

ly=20.#...................Base of the computational domain (dist from -yi to yf)
lx=10.#...................Height of the computational domain (dist from 0  to xf)
vf=[1.0,1.0] #............Surface loads vector 
vbeta=[90.,90.] #.........Vector con la inclinacion de las cargas respectoa la horizontal  
vm=[0.,0.] #............. Surface moment vector
vcy=[-5.,5.] #..............Vector con la coordenada "Y"  en que se aplica cada carga 

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
        sp[1,0]=mtxy[i,j]
        # find the maximum values
        la,ve = eigh(sp)
        msmin[i,j]=la[0]
        msmax[i,j]=la[1]
        mtmax[i,j]=abs(la[0]-la[1])/2.
        
        
levels=np.linspace(-0.5, 0.5, 10.0)
############################################################################
### ojo X va hacia abajo y Y hacia la derecha
#######
# Se calculó msx=(Matriz con los esfuerzoz axiales Sigma X de cada punto)
# Se calculó msy=(Matriz con los esfuerzoz axiales Sigma Y de cada punto)
# Se calculó mtxy=(Matriz con los esfuerzoz tangenciales TaoXY de cada punto)
# Se calculó mtxy=(Matriz con los esfuerzoz tangenciales TaoXY de cada punto)
# Se calculó msmin=(Matriz con los esfuerzoz axiales mínimos en cada punto)
# Se calculó msmax=(Matriz con los esfuerzoz axiales máximos en cada punto)
# Se calculó mtmax=(Matriz con los esfuerzoz tangenciales máximos en cada punto)


plt.figure(1)
plt.contourf(y,-x,msx ,levels, alpha=.75, cmap='jet')
plt.title('Esfuerzos axiales en X')
plt.colorbar()
plt.grid()

plt.figure(2)
plt.contourf(y,-x,msy ,levels, alpha=.75, cmap='jet')
plt.title('Esfuerzos axiales en Y')
plt.colorbar()
plt.grid()

plt.figure(3)
plt.contourf(y,-x,mtxy ,levels, alpha=.75, cmap='jet')
plt.title('Esfuerzos cortantes caras X y Y ')
plt.colorbar()
plt.grid()

plt.figure(4)
plt.contourf(y,-x,msmin ,levels, alpha=.75, cmap='jet')
plt.title('Esfuerzos axiales Minimos')
plt.colorbar()
plt.grid()

plt.figure(5)
plt.contourf(y,-x,msmax ,levels, alpha=.75, cmap='jet')
plt.title('Esfuerzos axiales Maximos')
plt.colorbar()
plt.grid()

plt.figure(6)
plt.contourf(y,-x,mtmax ,levels, alpha=.75, cmap='jet')
plt.title('Esfuerzos cortantes Maximos')
plt.colorbar()
plt.grid()







































