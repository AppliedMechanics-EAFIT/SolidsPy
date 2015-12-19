# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:49:18 2015

@author: eafit
"""
#
from __future__ import division
import numpy as np
import femutil as fem
import gaussutil as gau
from sympy import *
from sympy import init_printing
init_printing()
#
############################################################################
#-------------------------ELEMENT SUBROUTINES----------------------------- #
# Each UEL subroutine computes the local stiffness matrix for a given      #
# finite element.                                                          #
# New elements can be added by including additional subroutines            #
#                                                                          #
############################################################################
def uel4nquad(coord,enu,Emod):
    kl=np.zeros([8,8],dtype=np.float)
    C=fem.umat(enu,Emod)
    XW,XP =gau.gpoints2x2()
    ngpts=4
    for i in range(0,ngpts):
        ri=XP[i,0]
        si=XP[i,1]
        alf=XW[i]
        ddet, B = fem.stdm4NQ(ri,si,coord)
        kl=kl+B.T*C*B*alf*ddet
    return kl
#
def uel6ntrian(coord,enu,Emod):
    kl=np.zeros([12,12],dtype=np.float)
    C=fem.umat(enu,Emod)
    XW,XP =gau.gpoints7()
    ngpts=7
    for i in range(0,ngpts):
        ri=XP[i,0]
        si=XP[i,1]
        alf=XW[i]
        ddet, B = fem.stdm6NT(ri,si,coord)
        kl=kl+0.5*B.T*C*B*alf*ddet
    return kl    
#
def uel3ntrian(coord,enu,Emod):
    kl=np.zeros([6,6],dtype=np.float)
    C=fem.umat(enu,Emod)
    XW,XP =gau.gpoints3()
    ngpts=3
    for i in range(0,ngpts):
        ri=XP[i,0]
        si=XP[i,1]
        alf=XW[i]
        ddet, B = fem.stdm3NT(ri,si,coord)
        kl=kl+0.5*B.T*C*B*alf*ddet
#
    return kl
#