# -*- coding: utf-8 -*-
"""
Postprocessot subroutines
-------------------------

@author: eafit
"""
from __future__ import division
import numpy as np
from sympy import *
import femutil as fe
import preprocesor as pre
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['image.cmap'] = "YlGnBu_r"


def plotdis(IBC, UG, nodes, nn, xmin, xmax, ymin, ymax, savefigs=False):
    """Plot the nodal displacement solution using `griddata()`

    Parameters
    ----------
    IBC : ndarray (int)
      IBC (Indicator of Boundary Conditions) indicates if the nodes
      has any type of boundary conditions applied to it.
    UG : ndarray (float)
      Array with the computed displacements.
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    nn : int
      Number of nodes.
    xmin : float
      Minimum x value for the grid.
    xmax : float
      Maximum x value for the grid.
    ymin : float
      Minimum y value for the grid.
    ymax : float
      Maximum y value for the grid.

    """
    points = nodes[:, 1:3]
    grid_x, grid_y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

    UC = np.zeros([nn, 2], dtype=np.float)
    for i in range(nn):
        for j in range(2):
            kk = IBC[i, j]
            if kk == -1:
                UC[i, j] = 0.0
            else:
                UC[i, j] = UG[kk]

    grid_z0 = griddata(points, UC[:, 0], (grid_x, grid_y), method='linear')
    grid_z1 = griddata(points, UC[:, 1], (grid_x, grid_y), method='linear')

    plt.figure("Solution: Horizontal displacement")
    plt.imshow(grid_z0.T, aspect='equal', extent=(xmin, xmax, ymin, ymax),
               origin='lower')
    plt.title(r'$u_x$')
    plt.colorbar(orientation='vertical')
    plt.grid()
    if savefigs:
        plt.savefig('numhorizo.pdf')

    plt.figure("Solution: Vertical displacement")
    plt.imshow(grid_z1.T, aspect='equal', extent=(xmin, xmax, ymin, ymax),
               origin='lower')
    plt.title(r'$u_y$')
    plt.colorbar(orientation='vertical')
    plt.grid()
    if savefigs:
        plt.savefig('numvertic.pdf')


def grafmat(k):
    """Plot stiffness matrix sparsity

    Parameters
    ----------
    k : ndarray (int)
      Stiffness matrix of the system.

    """
    plt.figure("Stiffness matrix")
    plt.spy(k)
    plt.title("Stiffness matrix")
    plt.ylabel(r"$i$ index", size=10)
    plt.xlabel(r"$j$ index", size=10)


def scatter(DME, UG, ne, neq, elements):
    """Scatter the nodal displacements vector `UG` over each element

    Parameters
    ----------
    DME : ndarray (int)
      Array that shows the connectivity of degrees of freedom.
    UG : ndarray (float)
      Array with the computed displacements.
    ne : int
      Number of elements.
    neq : int
      Number of equations (degrees of freedom).
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    Returns
    -------
    UU : ndarray (float)
      Array with the displacements. This one contains both, the
      computed and imposed values.

    """
    iet = elements[0, 1]
    ndof, nnodes, ngpts = fe.eletype(iet)
    UU = np.zeros([ne, ndof], dtype=np.float)
    for i in range(ne):
        for ii in range(ndof):
            kk = DME[i, ii]
            if kk != -1:
                UU[i, ii] = UG[kk]

    return UU


def plotstrain(EG, XS, xmin, xmax, ymin, ymax, savefigs=False):
    """Plot the strain solution over the full domain

    Using griddata plots the strain solution over the full
    domain defined by the integration points. The integration
    points physical coordinates are stored in XS[] while the
    strain solution is stored in EG[].

    Parameters
    ----------
    EG : ndarray (float)
      Array that contains the strain solution for each integration
      point in physical coordinates.
    XS : ndarray (float)
      Array with the coordinates of the integration points.
    xmin : float
      Minimum x value for the grid.
    xmax : float
      Maximum x value for the grid.
    ymin : float
      Minimum y value for the grid.
    ymax : float
      Maximum y value for the grid.

    """
    grid_x, grid_y = np.mgrid[xmin:xmax:20j, ymin:ymax:20j]
    grid_z0 = griddata(XS, EG[:, 0], (grid_x, grid_y), method='linear')
    grid_z1 = griddata(XS, EG[:, 1], (grid_x, grid_y), method='linear')
    grid_z2 = griddata(XS, EG[:, 2], (grid_x, grid_y), method='linear')

    plt.figure("Solution: epsilon-xx strain")
    plt.imshow(grid_z0.T, aspect='equal', extent=(xmin, xmax, ymin, ymax),
               origin='lower')
    plt.title(r'$\epsilon_{xx}$')
    plt.colorbar(orientation='vertical')
    plt.grid()
    if savefigs:
        plt.savefig('numepsixx.pdf')

    plt.figure("Solution: epsilon-yy strain")
    plt.imshow(grid_z1.T, aspect='equal', extent=(xmin, xmax, ymin, ymax),
               origin='lower')
    plt.title(r'$\epsilon_{yy}$')
    plt.colorbar(orientation='vertical')
    plt.grid()
    if savefigs:
        plt.savefig('numepsiyy.pdf')

    plt.figure("Solution: gamma-xy strain")
    plt.imshow(grid_z2.T, aspect='equal', extent=(xmin, xmax, ymin, ymax),
               origin='lower')
    plt.title(r'$\gamma_{xy}$')
    plt.colorbar(orientation='vertical')
    plt.grid()
    if savefigs:
        plt.savefig('numgamaxy.pdf')


def xstrain(IELCON, nodes, ne, hh):
    """Compute physical coordinates of the domain integration points

    Parameters
    ----------
    IELCON : ndarray (int)
      Array with the nodes numbers for each element.
    nodes : ndarray (float)
      Array with nodes coordinates.
    ne : int
      Number of elements.
    hh : float
      Description pending...

    """
    XS = np.zeros([4*ne, 2], dtype=np.float)
    xl = np.zeros([4, 2], dtype=np.float)
    for i in range(ne):
        idp = IELCON[i, 0]
        xp = nodes[idp, 1]
        yp = nodes[idp, 2]
        xl[0, 0] = xp + hh/2
        xl[1, 0] = xp + 3*hh/2
        xl[2, 0] = xl[0, 0]
        xl[3, 0] = xl[1, 0]
        xl[0, 1] = yp + 3*hh/2
        xl[1, 1] = xl[0, 1]
        xl[2, 1] = yp + hh/2
        xl[3, 1] = xl[2, 1]
        XS[4*i:4*i + 4, 0] = xl[:, 0]
        XS[4*i:4*i + 4, 1] = xl[:, 1]

    return XS


def strainGLO(IELCON, UU, ne, COORD, elements):
    """Compute the strain solution for all the elements

    Computes the strain solution for all the elements
    in the domain and the physical coordinates of the complete
    domain integration points. It then assembles all the element strains
    into a global strains vector EG[].

    Parameters
    ----------
    IELCON : ndarray (int)
      Array with the nodes numbers for each element.
    UU : ndarray (float)
      Array with the displacements. This one contains both, the
      computed and imposed values.
    ne : int
      Number of elements.
    COORD : ndarray (float).
      Array with nodes coordinates.
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    Returns
    -------
    EG : ndarray (float)
      Array that contains the strain solution for each integration
      point in physical coordinates.
    XS : ndarray (float)
      Array with the coordinates of the integration points.

    """
    iet = elements[0, 1]
    ndof, nnodes, ngpts = fe.eletype(iet)

    XS = np.zeros([ngpts*ne, 2], dtype=np.float)
    elcoor = np.zeros([nnodes, 2], dtype=np.float)
    EG = np.zeros([ngpts*ne, 3], dtype=np.float)
    ul = np.zeros([ndof], dtype=np.float)
    for i in range(ne):
        for j in range(nnodes):
            elcoor[j, 0] = COORD[IELCON[i, j], 0]
            elcoor[j, 1] = COORD[IELCON[i, j], 1]
        for j in range(ndof):
            ul[j] = UU[i, j]
        if iet == 1:
            epsG, xl = fe.str_el4(elcoor, ul)
        elif iet == 2:
            epsG, xl = fe.str_el6(elcoor, ul)
        elif iet == 3:
            epsG, xl = fe.str_el3(elcoor, ul)

        for j in range(ngpts):
            XS[ngpts*i + j, 0] = xl[j, 0]
            XS[ngpts*i + j, 1] = xl[j, 1]
            for k in range(3):
                EG[ngpts*i + j, k] = epsG[j, k]

    return EG, XS


def axisscale(COORD, nn):
    """Determine plotting range

    Using the nodal coordinates it retrives minimum and maximum values
    along each direction for plotting purposes.

    Parameters
    ----------
    COORD : ndarray (float).
      Array with nodes coordinates.
    nn : int
      Number of nodes.

    Returns
    -------
    xmin : float
      Minimum value for x.
    xmax : float
      Maximum value for x.
    ymin : float
      Minimum value for y.
    ymax : float
      Maximum value for y.

    """
    xmin, ymin = np.min(COORD, axis=0)
    xmax, ymax = np.max(COORD, axis=0)

    return xmin, xmax, ymin, ymax


def plotstraincontours(EG, XS, ne, elements):
    """Plot contours for the strain solution

    Using Python function contourf() it plots the strain solution.

    Parameters
    ----------
    EG : ndarray (float)
      Array that contains the strain solution for each integration
      point in physical coordinates.
    XS : ndarray (float)
      Array with the coordinates of the integration points.
    ne : int
      Number of elements.
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    """
    iet = elements[0, 1]
    ndof, nnodes, ngpts = pre.eletype(iet)
    epsx = EG[:, 0]
    epsy = EG[:, 1]
    gamx = EG[:, 2]
    epsx.shape = (ne*ngpts, ne*ngpts)
    epsy.shape = (ne*ngpts, ne*ngpts)
    gamx.shape = (ne*ngpts, ne*ngpts)

    X, Y = np.meshgrid(XS[:, 0], XS[:, 1])
    plt.contourf(X, Y, epsx.T, alpha=0.75)
    plt.title(r'$\epsilon_{xx}$')
    plt.colorbar(orientation='horizontal')
    plt.grid()

    plt.contourf(X, Y, epsy.T, alpha=0.75)
    plt.title(r'$\epsilon_{yy}$')
    plt.colorbar(orientation='horizontal')
    plt.grid()

    plt.contourf(X, Y, gamx.T, alpha=0.75)
    plt.title(r'$\gamma_{xy}$')
    plt.colorbar(orientation='horizontal')
    plt.grid()


def locstrain4nQ(ul, coord, enu, Emod):
    """
    Plot strain and stress for a bilinear 4-noded element in the
    natural space

    Parameters
    ----------
    ul : ndarray (float)
      Array with the displacements for the element.
    coord : ndarray (float)
      Coordinates for the nodes of the element.
    enu : float
      Poisson coefficient in the range (-1, 0.5).
    Emod : float
      Young modulus, should be greater than zero.

    """
    r = np.zeros([4, 2])
    eG = np.zeros([3, 4])
    sG = np.zeros([3, 4])
    eps = np.zeros([3, 4])
    e = np.zeros([4])
    s = np.zeros([4])
    C = np.zeros([3, 3])

    r[0, 0] = -1.0
    r[0, 1] = -1.0
    r[1, 0] = 1.0
    r[1, 1] = -1.0
    r[2, 0] = 1.0
    r[2, 1] = 1.0
    r[3, 0] = -1.0
    r[3, 1] = 1.0
    C = fe.umat(enu, Emod)
    for i in range(4):
        ri = r[i, 0]
        si = r[i, 1]
        ddet, B = fe.stdm4NQ(ri, si, coord)
        eps = B*ul
        sig = C*eps
        eG[:, i] = eps[:]
        sG[:, i] = sig[:]
    grid_x, grid_y = np.mgrid[-1:1:100j, -1:1:100j]

    print('Strain field')
    for j in range(3):
        for i in range(4):
            e[i] = eG[j, i]
        grid_z0 = griddata(r, e, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T, aspect='equal', extent=(-1.0, 1.0, -1.0, 1.0),
                   origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()
        plt.show()

    print('Stress field')
    for j in range(3):
        for i in range(4):
            s[i] = sG[j, i]
        grid_z0 = griddata(r, s, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T, aspect='equal', extent=(-1.0, 1.0, -1.0, 1.0),
                   origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()
        plt.show()


def locstrain3nT(ul, coord, enu, Emod):
    """
    Plot strain and stress for a linear 3-noded element in the
    natural space

    Parameters
    ----------
    ul : ndarray (float)
      Array with the displacements for the element.
    coord : ndarray (float)
      Coordinates for the nodes of the element.
    enu : float
      Poisson coefficient in the range (-1, 0.5).
    Emod : float
      Young modulus, should be greater than zero.

    """
    r = np.zeros([3, 2])
    eG = np.zeros([3, 3])
    sG = np.zeros([3, 3])
    eps = np.zeros([3, 3])
    e = np.zeros([3])
    s = np.zeros([3])
    C = np.zeros([3, 3])

    r[0, 0] = 0.0
    r[0, 1] = 0.0
    r[1, 0] = 1.0
    r[1, 1] = 0.0
    r[2, 0] = 0.0
    r[2, 1] = 1.0
    C = fe.umat(enu, Emod)
    for i in range(3):
        ri = r[i, 0]
        si = r[i, 1]
        ddet, B = fe.stdm3NT(ri, si, coord)
        eps = B*ul
        sig = C*eps
        eG[:, i] = eps[:]
        sG[:, i] = sig[:]
    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:100j]

    print('Strain field')
    for j in range(3):
        for i in range(3):
            e[i] = eG[j, i]
        plt.figure()
        grid_z0 = griddata(r, e, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T, aspect='equal', extent=(0.0, 1.0, 0.0, 1.0),
                   origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()

    print('Stress field')
    for j in range(3):
        for i in range(3):
            s[i] = sG[j, i]
        plt.figure()
        grid_z0 = griddata(r, s, (grid_x, grid_y), method='linear')
        plt.imshow(grid_z0.T, aspect='equal', extent=(0.0, 1.0, 0.0, 1.0),
                   origin='lower')
        plt.colorbar(orientation='vertical')
        plt.grid()


def gmeshpost(IBC, nn, UG):
    """Export the nodal displacement solution

    Stores the nodal displacements solution vector into the file
    `out.txt` required to produce Gmesh readable files.

    """
    UR = np.zeros([nn, 2])
    for i in range(nn):
        for j in range(2):
            k = IBC[i, j]
            if k == -1:
                UR[i, j] = 0.0
            else:
                UR[i, j] = UG[k]
    nomfile1 = '../MESHUTIL/out.txt'
    np.savetxt(nomfile1, UR, fmt='%.18e', delimiter=' ')
