# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 18:52:09 2015

@author: eafit
"""

from __future__ import division
import femutil as fe
import uelutil as ue
import postprocesor as pos
import numpy as np
from sympy import *
from sympy import init_printing
init_printing()
#
xxx=np.zeros([4,2],dtype=np.float)
vl=np.zeros([8],dtype=np.float)
xxx[0,0]=0.0
xxx[0,1]=0.0
xxx[1,0]=2.0
xxx[1,1]=0.0
xxx[2,0]=2.0
xxx[2,1]=2.0
xxx[3,0]=0.0
xxx[3,1]=2.0
print xxx
#
vl[0]=0.0
vl[1]=0.0
vl[2]=0.0
vl[2]=-6.55275061e-04
vl[3]=-9.01011415e-04
vl[4]=-6.46540321e-05
vl[5]=-8.09443450e-04
vl[6]=0.0
vl[7]=0.0
#
xdet,BB=fe.stdm4NQ(-1,1,xxx)
pos.locstrain4nQ(vl,xxx,0.30,2700000)