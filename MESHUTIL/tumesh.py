# -*- coding: utf-8 -*-
"""
Created on Nov 26 16:45:00 2015

@author: 

Ruttine for mesh by gmesh 
"""
from __future__ import division
import os

os.system ('/Applications/Gmsh.app/Contents/MacOS/gmsh tunel.geo -2 -order 1')
print('End Proccess')
