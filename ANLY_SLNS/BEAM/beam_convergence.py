# -*- coding: utf-8 -*-
"""

"""
from __future__ import division
from os import sys
sys.path.append("../../MAIN/")
import numpy as np
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
import matplotlib.pyplot as plt
from beam import beam_sln

nx = 24
ny = 8
P = -50
E = 1000
nu = 0.3
L = 24
h = 8
I = 42.67
x = np.linspace(0, L, nx)
y = np.linspace(-h/2, h/2, ny)
x, y = np.mgrid[0:L:nx*1j, -h/2:h/2:ny*1j]


nodes = np.zeros((nx*ny, 5))
mats = np.array([1, 0.3])
loads = np.zeros((ny, 3))
elements = np.zeros(((nx - 1)*(ny - 1), 6))
u, v, _, _, _  = beam_sln(x, y, nu, P, E, I, L, h)

#ne, nn, nm, nl, COORD = pre.proini(nodes, mats, elements, loads)
#neq, IBC = ass.eqcounter(nn, nodes)
#DME, IELCON = ass.DME(IBC, ne, elements)
#KG = ass.matassem(IBC, mats, elements, nn, ne, neq, COORD, DME, IELCON)
#
#RHSG = ass.loadasem(loads, IBC, neq, nl)