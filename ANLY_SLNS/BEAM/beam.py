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
<<<<<<< HEAD
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
=======
rcParams['font.size'] = 12
>>>>>>> 87dcacff7624933db8ea0217b03966b128e9680b
rcParams['image.cmap'] = "YlGnBu_r"

>>>>>>> 75aa1daec17d6a185a063a28cd03e0706b1ba33f
plt.close("all")
nodes = np.loadtxt('nodes.txt')
mater = np.loadtxt('mater.txt')
nu, P, E, I, L, h = mater
G = E/(2*(1 + nu))
c = h/2
xx = nodes[:, 1:3]
#
# Assign symbols
C1 = -P/(2*E*I)
C2 = -(nu*P)/(6*E*I)
C3 = P/(6*I*G)
C4 = (P*L**2)/(2*E*I)
C5 = -(P*c**2)/(2*I*G)
C6 = C4 + C5
#
C7 = (nu*P)/(2*E*I)
C8 = P/(6*E*I)
C9 = -(P*L**2)/(2*E*I)
C10 = (P*L**3)/(3*E*I)
#
B1 = -P/(E*I)
B2 = (nu*P)/(E*I)
B3 = P/(2*I*G)
#
grid_x, grid_y = np.mgrid[0:24:100j, -4:4:100j]
x = xx[:, 0]
y = xx[:, 1]
u = C1*y*x**2 + C2*y**3 + C3*y**3 + C6*y
v = C7*x*y**2 + C8*x**3 + C9*x + C10
exx = B1*x*y
eyy = B2*x*y
gammaxy = B3*(y**2 - c**2)

#%% Plot results
def plot_field(xx, field, grid_x, grid_y, fname, title):
    plt.figure()
    grid_z = griddata(xx, field, (grid_x, grid_y), method='linear')
    plt.imshow(grid_z.T, aspect='equal', extent=(0,24,-4,4),
               origin='lower')
    plt.title(title)
    plt.colorbar(orientation='horizontal')
    plt.grid()
    plt.savefig(fname)


plot_field(xx, u, grid_x, grid_y, 'anahorizo.pdf', 'Horizontal field')
plot_field(xx, v, grid_x, grid_y, 'anavertic.pdf', 'Vertical field')
plot_field(xx, exx, grid_x, grid_y, 'anaepsixx.pdf', 'Epsilon xx')
plot_field(xx, eyy, grid_x, grid_y, 'anaepsiyy.pdf', 'Epsilon yy')
plot_field(xx, gammaxy, grid_x, grid_y, 'anagamaxy.pdf', 'Gamma xy')
plt.show()