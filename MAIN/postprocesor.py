# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 10:02:49 2015

@author: eafit
"""
from __future__ import division
import numpy as np
import femutil as fe
import preprocesor as pre
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from sympy import *
from sympy import init_printing
init_printing()
#
def plotdis(IBC,UG,nodes,nn, xmin,xmax,ymin,ymax):
#
# Plots the nodal displacements solution using the Python
# function griddata()
#
    plt.close()
    points     = np.zeros([nn,2],dtype=np.float)
    for i in range(0,nn):
        points[i,0]=nodes[i,1]
        points[i,1]=nodes[i,2]
    grid_x, grid_y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    #
    UC     = np.zeros([nn,2],dtype=np.float)
    for i in range(0,nn):
        for j in range(0,2):
            kk=IBC[i,j]
            if kk==-1:
                UC[i,j]=0.0
            else:
                UC[i,j]=UG[kk]
    #
    grid_z0 = griddata(points, UC[:,0], (grid_x, grid_y), method='linear')
    grid_z1 = griddata(points, UC[:,1], (grid_x, grid_y), method='linear')
    
    plt.imshow(grid_z0.T,aspect='equal', extent=(xmin,xmax,ymin,ymax), origin='lower')
    plt.title('Horizontal field')
    plt.colorbar(orientation='vertical')
    plt.grid()
#    plt.savefig('numhorizo.pdf')
    plt.show()
    #
    plt.imshow(grid_z1.T,aspect='equal', extent=(xmin,xmax,ymin,ymax), origin='lower')
    plt.title('Vertical field')
    plt.colorbar(orientation='vertical')
    plt.grid()
#    plt.savefig('numvertic.pdf')
    plt.show()
#
def grafmat(k):
#
#   Plots matrix k
#
    plt.close()
    plt.figure()
    plt.spy(k), plt.title("Stiffness matrix")
    plt.ylabel(r"$i$ index", size=10)
    plt.xlabel(r"$j$ index", size=10)
    plt.show()
#    
    return xc
#
def scatter(DME,UG,ne,neq,elements):
#
# Scatters the nodal displacements vector UG
# over each element.
#
    iet=elements[0,1]
    ndof,nnodes, ngpts =fe.eletype(iet)
    UU=np.zeros([ne,ndof],dtype=np.float)
    for i in range(0,ne):
        for ii in range(0,ndof):
            kk=DME[i,ii]
            if kk<>-1:
                UU[i,ii]=UG[kk]
#
    return UU
#
def plotstrain(EG,XS,ne,xmin,xmax,ymin,ymax):
#
# Using griddata plots the strain solution over the full
# domain defined by the integration points. The integration
# points physical coordinates are stored in XS[] while the
# strain solution is stored in EG[].
#
    plt.close()
    grid_x, grid_y = np.mgrid[xmin:xmax:20j,ymin:ymax:20j]
    #
    grid_z0 = griddata(XS, EG[:,0], (grid_x, grid_y), method='linear')
    grid_z1 = griddata(XS, EG[:,1], (grid_x, grid_y), method='linear')
    grid_z2 = griddata(XS, EG[:,2], (grid_x, grid_y), method='linear')
    
    plt.imshow(grid_z0.T,aspect='equal', extent=(xmin,xmax,ymin,ymax), origin='lower')
    plt.title('Epsilon-xx')
    plt.colorbar(orientation='vertical')
    plt.grid()
#    plt.savefig('numepsixx.pdf')
    plt.show()
    #
    plt.imshow(grid_z1.T,aspect='equal', extent=(xmin,xmax,ymin,ymax), origin='lower')
    plt.title('Epsilon-yy')
    plt.colorbar(orientation='vertical')
    plt.grid()
#    plt.savefig('numepsiyy.pdf')
    plt.show()
#
    plt.imshow(grid_z2.T,aspect='equal', extent=(xmin,xmax,ymin,ymax), origin='lower')
    plt.title('Gamma-xy')
    plt.colorbar(orientation='vertical')
    plt.grid()
#    plt.savefig('numgamaxy.pdf')
    plt.show()
#
    return
#
def xstrain(IELCON,nodes,ne,hh):
#
# Computes the physical coordinates of the complete
# domain integration points.
#
    XS=np.zeros([4*ne,2],dtype=np.float)
    xl=np.zeros([4,2],dtype=np.float)
    for i in range(0,ne):
        idp=IELCON[i,0]
        xp=nodes[idp,1]
        yp=nodes[idp,2]
        xl[0,0]=xp+hh/2
        xl[1,0]=xp+3*hh/2
        xl[2,0]=xl[0,0]
        xl[3,0]=xl[1,0]
        xl[0,1]=yp+3*hh/2
        xl[1,1]=xl[0,1]
        xl[2,1]=yp+hh/2
        xl[3,1]=xl[2,1]
        for j in range(0,4):
            XS[4*i+j,0]=xl[j,0]
            XS[4*i+j,1]=xl[j,1]           
#
    return XS
#
def strainGLO(IELCON,UU,ne,COORD,elements):
#
# Computes the strain solution for all the elements
# in the domain and the physical coordinates of the complete
# domain integration points. It then assembles all the element strains
# into a global strains vector EG[].
#
    iet=elements[0,1]
    ndof,nnodes, ngpts=fe.eletype(iet)
#
    XS=np.zeros([ngpts*ne,2],dtype=np.float)
    elcoor=np.zeros([nnodes,2],dtype=np.float)
    EG=np.zeros([ngpts*ne,3],dtype=np.float)
    ul=np.zeros([ndof],dtype=np.float)
    for i in range(0,ne):
        for j in range(0,nnodes):
            elcoor[j,0]=COORD[IELCON[i,j],0]
            elcoor[j,1]=COORD[IELCON[i,j],1]
        for j in range(0,ndof):
            ul[j]=UU[i,j]
        if iet==1:
            epsG,xl =fe.str_el4(elcoor,ul)
        elif iet==2:
            epsG,xl =fe.str_el6(elcoor,ul)
        elif iet==3:
            epsG,xl =fe.str_el3(elcoor,ul)
        
        for j in range(ngpts):
            XS[ngpts*i+j,0]=xl[j,0]
            XS[ngpts*i+j,1]=xl[j,1]
            for k in range(0,3):
                EG[ngpts*i+j,k]=epsG[j,k]
            
#
    return EG,XS
#
def axisscale(COORD,nn):
#
#   Using the nodal coordinates it retrives
#   minimum and maximum values along each
#   direction for plotting purposes.
#
    xs=np.zeros([nn,2],dtype=np.float)
    ys=np.zeros([nn,2],dtype=np.float)
    for i in range(0,nn):
        xs[i]=COORD[i,0]
        ys[i]=COORD[i,1]
    xmin=np.amin(xs)
    xmax=np.amax(xs)
    ymin=np.amin(ys)
    ymax=np.amax(ys)
    
    return xmin,xmax,ymin,ymax
#
def plotstraincontours(EG,XS,ne,elements):
#
# Using Python function contourf() it plots the
# strain solution.
#
    plt.close()
    iet=elements[0,1]
    ndof,nnodes, ngpts=pre.eletype(iet)
    epsx=np.zeros([ne*ngpts,ne*ngpts],dtype=np.float)
    epsy=np.zeros([ne*ngpts,ne*ngpts],dtype=np.float)
    gamx=np.zeros([ne*ngpts,ne*ngpts],dtype=np.float)
    for i in range(0,ne*ngpts):
        for j in range(0,ne*ngpts):
            epsx[i,j]=EG[j,0]
            epsy[i,j]=EG[j,1]
            gamx[i,j]=EG[j,2]
    
    X,Y=np.meshgrid(XS[:,0],XS[:,1])
#    levels=np.linspace(-0.60,-20.0, 100)
    plt.contourf(X,Y, epsx.T,alpha=0.75, cmap='jet')
    plt.title('Epsilon-xx')
    plt.colorbar(orientation='horizontal')
    plt.grid()
    plt.show()
#
    plt.contourf(X,Y, epsy.T,alpha=0.75, cmap='jet')
    plt.title('Epsilon-yy')
    plt.colorbar(orientation='horizontal')
    plt.grid()
    plt.show()
#
    plt.contourf(X,Y, gamx.T,alpha=0.75, cmap='jet')
    plt.title('Gamma-xy')
    plt.colorbar(orientation='horizontal')
    plt.grid()
    plt.show()
#
def locstrain4nQ(ul,coord,enu,Emod):
#
# Plot the strain and the stress field for a linear 4-noded quad
# element in the natural space.
#
    r=np.zeros([4,2])
    eG=np.zeros([3,4])
    sG=np.zeros([3,4])
    eps=np.zeros([3,4])
    e=np.zeros([4])
    s=np.zeros([4])
    C=np.zeros([3,3])
#
    r[0,0]=-1.0
    r[0,1]=-1.0
    r[1,0]= 1.0
    r[1,1]=-1.0
    r[2,0]= 1.0
    r[2,1]= 1.0
    r[3,0]=-1.0
    r[3,1]= 1.0
    C=fe.umat(enu,Emod)
    for i in range(0,4):
        ri=r[i,0]
        si=r[i,1]
        ddet, B = fe.stdm4NQ(ri,si,coord)
        eps=B*ul
        sig=C*eps
        eG[:,i]=eps[:]
        sG[:,i]=sig[:]
    grid_x, grid_y = np.mgrid[-1:1:100j, -1:1:100j]
#
    print 'strain field'
    for j in range(0,3):
        for i in range(0,4):
            e[i]=eG[j,i]
        grid_z0 = griddata(r, e, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T,aspect='equal', extent=(-1.0,1.0,-1.0,1.0), origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()
        plt.show()
#
    print 'stress field'
    for j in range(0,3):
        for i in range(0,4):
            s[i]=sG[j,i]
        grid_z0 = griddata(r, s, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T,aspect='equal', extent=(-1.0,1.0,-1.0,1.0), origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()
        plt.show()        
#
    return
#
def locstrain3nT(ul,coord,enu,Emod):
#
# Plot the strain and the stress field for a linear 4-noded quad
# element in the natural space.
#
    r=np.zeros([3,2])
    eG=np.zeros([3,3])
    sG=np.zeros([3,3])
    eps=np.zeros([3,3])
    e=np.zeros([3])
    s=np.zeros([3])
    C=np.zeros([3,3])
#
    r[0,0]= 0.0
    r[0,1]= 0.0
    r[1,0]= 1.0
    r[1,1]= 0.0
    r[2,0]= 0.0
    r[2,1]= 1.0
    C=fe.umat(enu,Emod)
    for i in range(0,3):
        ri=r[i,0]
        si=r[i,1]
        ddet, B = fe.stdm3NT(ri,si,coord)
        eps=B*ul
        sig=C*eps
        eG[:,i]=eps[:]
        sG[:,i]=sig[:]
    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:100j]
#
    print 'strain field'
    for j in range(0,3):
        for i in range(0,3):
            e[i]=eG[j,i]
        grid_z0 = griddata(r, e, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T,aspect='equal', extent=(0.0,1.0,0.0,1.0), origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()
        plt.show()
#
    print 'stress field'
    for j in range(0,3):
        for i in range(0,3):
            s[i]=sG[j,i]
        grid_z0 = griddata(r, s, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T,aspect='equal', extent=(0.0,1.0,0.0,1.0), origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()
        plt.show()        
#
    return
#
def gmeshpost(IBC,nn,UG):
#
# Stores the nodal displacements solution vector
# into the file salida.txt required to produce
# Gmesh readable files.
#
    UR=np.zeros([nn,2])
    for i in range(0,nn):
        for j in range(0,2):
            k=IBC[i,j]
            if k==-1:
                UR[i,j]=0.0
            else:
                UR[i,j]=UG[k]
    nomfile1='../MESHUTIL/salida.txt'    
    np.savetxt(nomfile1,UR, fmt='%.18e', delimiter=' ')
    
    return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    