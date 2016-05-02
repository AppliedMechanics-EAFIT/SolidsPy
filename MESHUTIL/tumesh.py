# -*- coding: utf-8 -*-
"""
Created on Nov 26 16:45:00 2015

@author: 

Ruttine for mesh by gmesh 
"""
from __future__ import division
import os

os.system ('/Applications/Gmsh.app/Contents/MacOS/gmsh tap2.geo -2 -order 2')
print('End Proccess')
