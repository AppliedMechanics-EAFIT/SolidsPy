# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 10:02:49 2015

@author: eafit
"""
from __future__ import division
import numpy as np
from sympy import *
import uelutil as ue
import femutil as fem
import preprocesor as pre
from sympy import init_printing
init_printing()
#
def eqcounter(nn,nodes):
#
# Counts active equations and creates
# BCs array IBC[]
#
    IBC    = np.zeros([nn,2],dtype=np.integer)
    neq=0
    for i in range(0,nn):
        for j in range(0,2):
            IBC[i,j]=int(nodes[i,j+3])
            if IBC[i,j]==0:
                IBC[i,j]=neq
                neq=neq+1
#    
    return neq, IBC  
#
def DME(IBC,ne, elements):    
#
# Creates the assembly operator DME[] and
# processes the element connectivity array
# IELCON[]
#
    IELCON = np.zeros([ne,9],dtype=np.integer)
    DME    = np.zeros([ne,18],dtype=np.integer)
#
    for i in range(0,ne):
        iet=elements[i,1]
        ndof,nnodes, ngpts =fem.eletype(iet)  
        for j in range(0,nnodes):
            IELCON[i,j]=elements[i,j+3]
            kk=IELCON[i,j]
            for l in range(0,2):
                DME[i,2*j+l]=IBC[kk,l]
#
    return DME, IELCON
#
def matassem(IBC,mats,elements,nn,ne,neq,COORD,DME,IELCON):
#
# Assembles the global stiffness matrix KG
#
    KG =np.zeros([neq,neq],dtype=np.float)
    for i in range(0,ne):
        iet=elements[i,1]
        ndof,nnodes, ngpts =fem.eletype(iet)
        elcoor=np.zeros([nnodes,2],dtype=np.float)
        kloc=np.zeros([ndof,ndof],dtype=np.float)
        im=elements[i,2]                 
        emod=mats[im,0]
        enu =mats[im,1]
        dme=np.zeros([ndof],dtype=np.integer)
        for j in range(0,nnodes):
            elcoor[j,0]=COORD[IELCON[i,j],0]
            elcoor[j,1]=COORD[IELCON[i,j],1]
        if iet==1:
            kloc=ue.uel4nquad(elcoor,enu,emod)
        elif iet==2:
            kloc=ue.uel6ntrian(elcoor,enu,emod)
        elif iet==3:
            kloc=ue.uel3ntrian(elcoor,enu,emod)
        for ii in range(0,ndof):
            dme[ii]=DME[i,ii]
        for ii in range(0,ndof):
            kk=dme[ii]
            if kk<>-1:
                for jj in range(0,ndof):
                    ll=dme[jj]
                    if ll<>-1:
                        KG[kk,ll]=KG[kk,ll]+kloc[ii,jj]
#    
    return KG 
#
def loadasem(loads,IBC,neq,nl):
#
# Assembles the global Right Hand Side Vector RHSG
#
    RHSG=np.zeros([neq],dtype=np.float)                   
    for i in range(0,nl):
        il=int(loads[i,0])
        ilx=IBC[il,0]
        ily=IBC[il,1]
        RHSG[ilx]=1*loads[i,1]
        RHSG[ily]=1*loads[i,2]
#
    return RHSG     
#
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    