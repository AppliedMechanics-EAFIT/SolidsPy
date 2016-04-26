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
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 12
rcParams['image.cmap'] = "YlGnBu_r"


def beam_sln(x, y, nu, P, E, I, L, h):
    """Compute the solution for a cantilever beam

    Parameters
    ----------
    x : ndarray (float)
        Array with x coordinates.
    y : ndarray (float)
        Array with y coordinates.
    nu : float, (-1, 0.5)
        Poisson coefficient.
    P : float
        Applied force at the end of the beam.
    E : float, >0
        Young modulus.
    I : float, >0
        Moment of inertia.
    L : float, >0
        Length of the beam.
    h : float, >0
        Height of the beam.

    Returns
    -------
    u : ndarray (float)
        Horizontal displacement at the nodes.
    v : ndarray (float)
        Vertical displacement at the nodes.
    exx : ndarray (float)
        xx component of the strain tensor.
    eyy : ndarray (float)
        yy component of the strain tensor.
    gammaxy : ndarray (float)
        xy component of the strain tensor.
        
    References
    ----------
    .. [1] Timoshenko, S. & Goodier, J., 1970. Theory of Elasticity,
        McGraw-Hill, 3rd Ed.

    """
    G = E/(2*(1 + nu))
    c = h/2
    C1 = -P/(2*E*I)
    C2 = -(nu*P)/(6*E*I)
    C3 = P/(6*I*G)
    C4 = (P*L**2)/(2*E*I)
    C5 = -(P*c**2)/(2*I*G)
    C6 = C4 + C5
    C7 = (nu*P)/(2*E*I)
    C8 = P/(6*E*I)
    C9 = -(P*L**2)/(2*E*I)
    C10 = (P*L**3)/(3*E*I)
    B1 = -P/(E*I)
    B2 = (nu*P)/(E*I)
    B3 = P/(2*I*G)
    u = C1*y*x**2 + C2*y**3 + C3*y**3 + C6*y
    v = C7*x*y**2 + C8*x**3 + C9*x + C10
    exx = B1*x*y
    eyy = B2*x*y
    gammaxy = B3*(y**2 - c**2)

    return u, v, exx, eyy, gammaxy    


def plot_field(xx, field, grid_x, grid_y, fname, title):
    """Plot contours of the analytical solution"""    
    plt.figure()
    grid_z = griddata(xx, field, (grid_x, grid_y), method='cubic')
    plt.contourf(grid_x, grid_y, grid_z, 10)
    plt.axis("image")
    plt.title(title)
    plt.colorbar(orientation='horizontal')
    plt.grid()
    plt.savefig(fname)


if __name__ == "__main__":
    #%% Material definition
    nodes = np.loadtxt('nodes.txt')
    mater = np.loadtxt('mater.txt')
    nu, P, E, I, L, h = mater
    
    xx = nodes[:, 1:3]
    grid_x, grid_y = np.mgrid[0:24:100j, -4:4:100j]
    u, v, exx, eyy, gammaxy  = beam_sln(xx[:, 0], xx[:, 1], nu, P, E, I, L, h)
    
    #%% Plot results
    plt.close("all")
    plot_field(xx, u, grid_x, grid_y, 'img/anahorizo.pdf', 'Horizontal field')
    plot_field(xx, v, grid_x, grid_y, 'img/anavertic.pdf', 'Vertical field')
    plot_field(xx, exx, grid_x, grid_y, 'img/anaepsixx.pdf', 'Epsilon xx')
    plot_field(xx, eyy, grid_x, grid_y, 'img/anaepsiyy.pdf', 'Epsilon yy')
    plot_field(xx, gammaxy, grid_x, grid_y, 'img/anagamaxy.pdf', 'Gamma xy')
    plt.show()