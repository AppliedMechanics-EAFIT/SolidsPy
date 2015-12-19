# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 10:02:49 2015

@author: eafit
"""
from __future__ import division
from sympy import *
import gaussutil as gau
import numpy as np
from sympy import init_printing
init_printing()
#
def eletype(iet):
#
# According to iet assigns number of degrees of freedom,
# number of nodes and minimum required number of integration
# points.
#

    if iet==1:
        ndof=8
        nnodes=4
        ngpts=4
    if iet==2:
        ndof=12
        nnodes=6
        ngpts=7
    if iet==3:
        ndof=6
        nnodes=3
        ngpts=3
#    
    return ndof,nnodes,ngpts
#
def sha4(x,y):
#
# Shape functions for a 4-noded quad element.
#
    H=zeros(4)
    N = zeros(2,8)
    H = S(1)/4*Matrix([(1 - x)*(1 - y),
         (1 + x)*(1 - y),
         (1 + x)*(1 + y),
         (1 - x)*(1 + y)])
#
    for i in range(4):
        N[0, 2*i] = H[i]
        N[1, 2*i + 1] = H[i]
#
    return N
#
def sha6(x,y):
#
# Shape functions for a 6-noded triang element.
#
    H=zeros(6)
    N = zeros(2,12)
#
    H=Matrix([(1-x-y)-2*x*(1-x - y)-2*y*(1-x-y),
         x-2*x*(1-x-y)-2*x*y,
         y-2*x*y-2*y*(1-x-y),
         4*x*(1 -x-y),
         4*x*y,
         4*y*(1-x-y)])         
#
    for i in range(6):
        N[0, 2*i] = H[i]
        N[1, 2*i + 1] = H[i]
#
    return N
#
def sha3(x,y):
#
# Shape functions for a 3-noded triang element.
#
    H=zeros(3)
    N = zeros(2,6)
#
    H=Matrix([(1-x-y),
         x,
         y])         
#
    for i in range(3):
        N[0, 2*i] = H[i]
        N[1, 2*i + 1] = H[i]
#
    return N
#
def stdm4NQ(r,s,coord):
#
# Strain-displacement interpolator B[] for
# a 4-noded quad element.
#
    rr , ss = symbols('rr ss')
    nn=4
    N=zeros(nn)
    B = zeros(3,2*nn)
    dhdx=zeros(2,nn)
    DNR=zeros(2,nn)
    N = S(1)/4*Matrix([(1 - rr)*(1 - ss),
         (1 + rr)*(1 - ss),
         (1 + rr)*(1 + ss),
         (1 - rr)*(1 + ss)])    
    for i in range(nn):
        dhdx[0,i]=diff(N[i],rr)
        dhdx[1,i]=diff(N[i],ss)
    DNR=dhdx.subs([(rr,r),(ss, s)])
#
    xj=jacoper(DNR,coord,nn)
    ddet=np.linalg.det(xj)
    xi=np.linalg.inv(xj)
    aux1=xi*DNR
    for i in range(nn):
        B[0, 2*i  ] = aux1[0, i]
        B[1, 2*i+1] = aux1[1, i]
        B[2, 2*i  ] = aux1[1, i]
        B[2, 2*i+1] = aux1[0, i]
    return ddet, B
#
def stdm6NT(r,s,coord):
#
# Strain-displacement interpolator B[] for
# a 6-noded triang element.
#
    rr , ss = symbols('rr ss')
    nn=6
    N=zeros(nn)
    B = zeros(3,2*nn)
    dhdx=zeros(2,nn)
    DNR=zeros(2,nn)
    N=Matrix([(1 - rr-ss)-2*rr*(1-rr - ss)-2*ss*(1-rr-ss),
         rr-2*rr*(1-rr-ss)-2*rr*ss,
         ss-2*rr*ss-2*ss*(1-rr-ss),
         4*rr*(1 -rr-ss),
         4*rr*ss,
         4*ss*(1-rr-ss)])
    for i in range(nn):
        dhdx[0,i]=diff(N[i],rr)
        dhdx[1,i]=diff(N[i],ss)
    DNR=dhdx.subs([(rr,r),(ss, s)])
#
    xj=jacoper(DNR,coord,nn)
    ddet=np.linalg.det(xj)
    xi=np.linalg.inv(xj)
    aux1=xi*DNR
    for i in range(nn):
        B[0, 2*i  ] = aux1[0, i]
        B[1, 2*i+1] = aux1[1, i]
        B[2, 2*i  ] = aux1[1, i]
        B[2, 2*i+1] = aux1[0, i]
    return ddet, B    
#
def stdm3NT(r,s,coord):
#
# Strain-displacement interpolator B[] for
# a 3-noded triang element.
#
    rr , ss = symbols('rr ss')
    nn=3
    N=zeros(nn)
    B = zeros(3,2*nn)
    dhdx=zeros(2,nn)
    DNR=zeros(2,nn)
    N=Matrix([(1 - rr-ss),
         rr,
         ss])
    for i in range(nn):
        dhdx[0,i]=diff(N[i],rr)
        dhdx[1,i]=diff(N[i],ss)
    DNR=dhdx.subs([(rr,r),(ss, s)])
#
    xj=jacoper(DNR,coord,nn)
    xi=np.linalg.inv(xj)
    ddet=np.linalg.det(xj)    
    aux1=xi*DNR
    for i in range(nn):
        B[0, 2*i  ] = aux1[0, i]
        B[1, 2*i+1] = aux1[1, i]
        B[2, 2*i  ] = aux1[1, i]
        B[2, 2*i+1] = aux1[0, i]
    return ddet, B      
#
def jacoper(dhdx,coord,nn):
    dum=0.0
    xja=np.zeros([2,2],dtype=np.float)
    for k in range(0,2):
        for j in range(0,2):
            for i in range(0,nn):
                dum=dum+dhdx[k,i]*coord[i,j]
            xja[k,j]=dum
            dum=0.0
    
    return xja
#
def umat(nu,E):
#
# 2D Elasticity consitutive matrix.
# For plane strain use effective properties.
#
    C=zeros(3,3)
    enu=E/(1-nu**2)
    mnu=(1-nu)/2.0
    C[0,0]=enu
    C[0,1]=nu*enu
    C[1,0]=C[0,1]
    C[1,1]=enu
    C[2,2]=enu*mnu
#
    return C       
#
def str_el4(coord,ul):
#
# Computes the strains at each element integration point.
# 4-noded quad.
#
    epsl=np.zeros([3],dtype=np.float)
    epsG=np.zeros([3,4],dtype=np.float)
    epsGT=np.zeros([4,3],dtype=np.float)
    xl=np.zeros([4,2],dtype=np.float)
    x, y = symbols('x, y')
#    
    XW, XP=gau.gpoints2x2()
    for i in range(0,4):
        ri=XP[i,0]
        si=XP[i,1]
        ddet, B = stdm4NQ(ri,si,coord)
        epsl=B*ul
        epsG[:,i]=epsl[:]
        N=sha4(ri,si)
        NN=N.subs([(x, ri), (y, si)])
        xl[i,0]=NN[0,0]*coord[0,0]+NN[0,2]*coord[1,0]+NN[0,4]*coord[2,0]+NN[0,6]*coord[3,0]
        xl[i,1]=NN[0,0]*coord[0,1]+NN[0,2]*coord[1,1]+NN[0,4]*coord[2,1]+NN[0,6]*coord[3,1]
    epsGT=epsG.T
    return epsGT,xl
#
def str_el6(coord,ul):
#
# 6-noded triang.
#
    epsl=np.zeros([3],dtype=np.float)
    epsG=np.zeros([3,7],dtype=np.float)
    epsGT=np.zeros([7,3],dtype=np.float)
    xl=np.zeros([7,2],dtype=np.float)
    x, y = symbols('x, y')
#    
    XW, XP=gau.gpoints7()
    for i in range(0,7):
        ri=XP[i,0]
        si=XP[i,1]
        ddet, B = stdm6NT(ri,si,coord)
        epsl=B*ul
        epsG[:,i]=epsl[:]
        N=sha6(ri,si)
        NN=N.subs([(x, ri), (y, si)])
        xl[i,0]=NN[0,0]*coord[0,0]+NN[0,2]*coord[1,0]+NN[0,4]*coord[2,0]+NN[0,6]*coord[3,0]+NN[0,8]*coord[4,0]+NN[0,10]*coord[5,0]
        xl[i,1]=NN[0,0]*coord[0,1]+NN[0,2]*coord[1,1]+NN[0,4]*coord[2,1]+NN[0,6]*coord[3,1]+NN[0,8]*coord[4,1]+NN[0,10]*coord[5,1]
    epsGT=epsG.T
    return epsGT,xl
#
#
def str_el3(coord,ul):
#
# 3-noded triang.
#
    epsl=np.zeros([3],dtype=np.float)
    epsG=np.zeros([3,3],dtype=np.float)
    epsGT=np.zeros([3,3],dtype=np.float)
    xl=np.zeros([3,2],dtype=np.float)
    x, y = symbols('x, y')
#    
    XW, XP=gau.gpoints3()
    for i in range(0,3):
        ri=XP[i,0]
        si=XP[i,1]
        ddet, B = stdm3NT(ri,si,coord)
        epsl=B*ul
        epsG[:,i]=epsl[:]
        N=sha3(ri,si)
        NN=N.subs([(x, ri), (y, si)])
        xl[i,0]=NN[0,0]*coord[0,0]+NN[0,2]*coord[1,0]+NN[0,4]*coord[2,0]
        xl[i,1]=NN[0,0]*coord[0,1]+NN[0,2]*coord[1,1]+NN[0,4]*coord[2,1]
    epsGT=epsG.T
    return epsGT,xl     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        