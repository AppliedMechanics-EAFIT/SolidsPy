# -*- coding: utf-8 -*-
"""
Created on Nov 26 16:45:00 2015

@author: 

Ruttine for mesh by gmesh 
"""
from __future__ import division
from sympy import init_printing
init_printing()
# ***** Import Package ****
import os
# ***** end *****
os.system ('/Applications/Gmsh.app/Contents/MacOS/gmsh dam.geo -2 -order 1')
#
print('End Proccess')