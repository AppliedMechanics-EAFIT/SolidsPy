# -*- coding: utf-8 -*-
"""
Elasticity solution for a cantilever loaded at the tip (See Timoshenko
and Young).

The script uses the inout files nodes.txt and mater.txt.

The paramters for the mater.txt input file are Poisson's ratio, tip
load, Young's modulus moment of inertia of the cross section, length
and heigth of the beam.
"""

from __future__ import division
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
<<<<<<< HEAD
from sympy import *
from sympy import init_printing
init_printing()
=======
>>>>>>> 75aa1daec17d6a185a063a28cd03e0706b1ba33f
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
<<<<<<< HEAD
#rcParams['image.cmap'] = "YlGnBu_r"
rcParams['image.cmap'] = "jet"
#
# Elasticity solution for a cantilever loaded at the tip (See Timoshenko and Young).
# The script uses the inout files nodes.txt and mater.txt.
# The paramters for the mater.txt input file are Poisson's ratio, tip load, Young's modulus
# moment of inertia of the cross section, length and heigth of the beam.
=======
rcParams['image.cmap'] = "YlGnBu_r"

>>>>>>> 75aa1daec17d6a185a063a28cd03e0706b1ba33f
plt.close("all")
nodes = np.loadtxt('nodes.txt')
mater = np.loadtxt('mater.txt')
nu = mater[0]
P = mater[1]
E = mater[2]
G = E/(2*(1+nu))
I = mater[3]
l = mater[4]
h = mater[5]
c = h/2
nn = len(nodes[:,0])
xx = nodes[:, 1:3]
#
# Assign symbols
#
##
C1 = -P/(2*E*I)
C2 = -(nu*P)/(6*E*I)
C3 = P/(6*I*G)
C4 = (P*l**2)/(2*E*I)
C5 = -(P*c**2)/(2*I*G)
C6 = C4 + C5
#
C7 = (nu*P)/(2*E*I)
C8 = P/(6*E*I)
C9 = -(P*l**2)/(2*E*I)
C10 = (P*l**3)/(3*E*I)
#
B1 = -P/(E*I)
B2 = (nu*P)/(E*I)
B3 = P/(2*I*G)
#
u = np.zeros([nn])
v = np.zeros([nn])
exx = np.zeros([nn])
eyy = np.zeros([nn])
gammaxy = np.zeros([nn])
grid_x, grid_y = np.mgrid[0:24:100j, -4:4:100j]
x = xx[:, 0]
y = xx[:, 1]
u = C1*y*x**2 + C2*y**3 + C3*y**3 + C6*y
v = C7*x*y**2 + C8*x**3 + C9*x + C10
exx = B1*x*y
eyy = B2*x*y
gammaxy = B3*(y**2 - c**2)

#
grid_z0 = griddata(xx, u, (grid_x, grid_y), method='linear')
grid_z1 = griddata(xx, v, (grid_x, grid_y), method='linear')
grid_z2 = griddata(xx, exx, (grid_x, grid_y), method='linear')
grid_z3 = griddata(xx, eyy, (grid_x, grid_y), method='linear')
grid_z4 = griddata(xx, gammaxy, (grid_x, grid_y), method='linear')
#
plt.imshow(grid_z0.T,aspect='equal', extent=(0,24,-4,4), origin='lower')
plt.title('Horizontal field')
plt.colorbar(orientation='horizontal')
plt.grid()
plt.savefig('anahorizo.pdf')
plt.show()
#
plt.imshow(grid_z1.T,aspect='equal', extent=(0,24,-4,4), origin='lower')
plt.title('Vertical field')
plt.colorbar(orientation='horizontal')
plt.grid()
plt.savefig('anavertic.pdf')
plt.show()
#
plt.imshow(grid_z2.T,aspect='equal', extent=(0,24,-4,4), origin='lower')
plt.title('Epsilon xx')
plt.colorbar(orientation='horizontal')
plt.grid()
plt.savefig('anaepsixx.pdf')
plt.show()
#
plt.imshow(grid_z3.T,aspect='equal', extent=(0,24,-4,4), origin='lower')
plt.title('Epsilon yy')
plt.colorbar(orientation='horizontal')
plt.grid()
plt.savefig('anaepsiyy.pdf')
plt.show()
#
plt.imshow(grid_z4.T,aspect='equal', extent=(0,24,-4,4), origin='lower')
plt.title('Gamma xy')
plt.colorbar(orientation='horizontal')
plt.grid()
plt.savefig('anagamaxy.pdf')
plt.show()