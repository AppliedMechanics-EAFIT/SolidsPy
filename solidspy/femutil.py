 # -*- coding: utf-8 -*-
"""
femutil.py
----------

Functions to compute kinematics variables for the Finite
Element Analysis.

The elements included are:
    1. 4 node bilinear quadrilateral.
    2. 6 node quadratic triangle.
    3. 3 node linear triangle.

The notation used is similar to the one used by Bathe [1]_.


References
----------
.. [1] Bathe, Klaus-Jürgen. Finite element procedures. Prentice Hall,
   Pearson Education, 2006.

"""
from __future__ import division, print_function
import gaussutil as gau
import numpy as np


def eletype(iet):
    """Assigns number to degrees of freedom

    According to iet assigns number of degrees of freedom, number of
    nodes and minimum required number of integration points.

    Parameters
    ----------
    iet :  int
      Type of element. These are:
        1. 4 node bilinear quadrilateral.
        2. 6 node quadratic triangle.
        3. 3 node linear triangle.
        5. 2 node spring.
        6. 2 node truss element.
        7. 2 node beam (3 DOF per node).

    Returns
    -------
    ndof : int
      Number of degrees of freedom for the selected element.
    nnodes : int
      Number of nodes for the selected element.
    ngpts : int
      Number of Gauss points for the selected element.

    """
    if iet == 1:
        ndof = 8
        nnodes = 4
        ngpts = 4
    if iet == 2:
        ndof = 12
        nnodes = 6
        ngpts = 7
    if iet == 3:
        ndof = 6
        nnodes = 3
        ngpts = 3
    if iet == 5:
        ndof = 4
        nnodes = 2
        ngpts = 3
    if iet == 6:
        ndof = 4
        nnodes = 2
        ngpts = 3
    if iet == 7:
        ndof = 6
        nnodes = 2
        ngpts = 3

    return ndof, nnodes, ngpts


#%% Shape functions and derivatives
def sha4(x, y):
    """Shape functions for a 4-noded quad element

    Parameters
    ----------
    x : float
      x coordinate for a point within the element.
    y : float
      y coordinate for a point within the element.

    Returns
    -------
    N : Numpy array
      Array of interpolation functions.

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (1, 1). Thus

    >>> N = sha4(0, 0)
    >>> N_ex = np.array([
    ...    [1/4, 0, 1/4, 0, 1/4, 0, 1/4, 0],
    ...    [0, 1/4, 0, 1/4, 0, 1/4, 0, 1/4]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N = sha4(1, 1)
    >>> N_ex = np.array([
    ...    [0, 0, 0, 0, 1, 0, 0, 0],
    ...    [0, 0, 0, 0, 0, 1, 0, 0]])
    >>> np.allclose(N, N_ex)
    True

    """
    N = np.zeros((2, 8))
    H = 0.25*np.array(
        [(1 - x)*(1 - y),
         (1 + x)*(1 - y),
         (1 + x)*(1 + y),
         (1 - x)*(1 + y)])
    N[0, ::2] = H
    N[1, 1::2] = H

    return N


def sha6(x, y):
    """Shape functions for a 6-noded triangular element

    Parameters
    ----------
    x : float
      x coordinate for a point within the element.
    y : float
      y coordinate for a point within the element.

    Returns
    -------
    N : Numpy array
      Array of interpolation functions.

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (0.5, 0.5). Thus

    >>> N = sha6(0, 0)
    >>> N_ex = np.array([
    ...    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ...    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N = sha6(1/2, 1/2)
    >>> N_ex = np.array([
    ...     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    ...     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]])
    >>> np.allclose(N, N_ex)
    True

    """
    N = np.zeros((2, 12))
    H = np.array(
        [(1 - x - y) - 2*x*(1 - x - y) - 2*y*(1 - x - y),
         x - 2*x*(1 - x - y) - 2*x*y,
         y - 2*x*y - 2*y*(1-x-y),
         4*x*(1 - x - y),
         4*x*y,
         4*y*(1 - x - y)])
    N[0, ::2] = H
    N[1, 1::2] = H

    return N


def sha3(x, y):
    """Shape functions for a 3-noded triangular element

    Parameters
    ----------
    x : float
      x coordinate for a point within the element.
    y : float
      y coordinate for a point within the element.

    Returns
    -------
    N : Numpy array
      Array of interpolation functions.

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (0, 0.5). Thus

    >>> N = sha3(0, 0)
    >>> N_ex = np.array([
    ...    [1, 0, 0, 0, 0, 0],
    ...    [0, 1, 0, 0, 0, 0]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N = sha3(1/2, 1/2)
    >>> N_ex = np.array([
    ...    [0, 0, 1/2, 0, 1/2, 0],
    ...    [0, 0, 0, 1/2, 0, 1/2]])
    >>> np.allclose(N, N_ex)
    True

    """
    N = np.zeros((2, 6))
    H = np.array([
        (1 - x - y),
         x,
         y])
    N[0, ::2] = H
    N[1, 1::2] = H

    return N



def stdm4NQ(r, s, coord):
    """Strain-displacement interpolator B for a 4-noded quad element

    Parameters
    ----------
    r : float
      r component in the natural space.
    s : float
      s component in the natural space.
    coord : ndarray
      Coordinates of the nodes of the element (4, 2).

    Returns
    -------
    ddet : float
      Determinant evaluated at `(r, s)`.
    B : ndarray
      Strain-displacement interpolator evaluated at `(r, s)`.

    """
    nn = 4
    B = np.zeros((3, 2*nn))
    dhdx = 0.25*np.array([
            [s - 1, -s + 1, s + 1, -s - 1],
            [r - 1, -r - 1, r + 1, -r + 1]])
    det, jaco_inv = jacoper(dhdx, coord)
    dhdx = np.dot(jaco_inv, dhdx)
    B[0, ::2] = dhdx[0, :]
    B[1, 1::2] = dhdx[1, :]
    B[2, ::2] = dhdx[1, :]
    B[2, 1::2] = dhdx[0, :]
    return det, B


def stdm6NT(r, s, coord):
    """Strain-displacement interpolator B for a 6-noded triang element

    Parameters
    ----------
    r : float
      r component in the natural space.
    s : float
      s component in the natural space.
    coord : ndarray
      Coordinates of the nodes of the element (6, 2).

    Returns
    -------
    ddet : float
      Determinant evaluated at `(r, s)`.
    B : ndarray
      Strain-displacement interpolator evaluated at `(r, s)`.

    """
    nn = 6
    B = np.zeros((3, 2*nn))
    dhdx = np.array([
        [4*r + 4*s - 3, 4*r - 1, 0, -8*r - 4*s + 4, 4*s,  -4*s],
        [4*r + 4*s - 3, 0, 4*s - 1,  -4*r, 4*r, -4*r - 8*s + 4]])
    det, jaco_inv = jacoper(dhdx, coord)
    dhdx = np.dot(jaco_inv, dhdx)
    B[0, ::2] = dhdx[0, :]
    B[1, 1::2] = dhdx[1, :]
    B[2, ::2] = dhdx[1, :]
    B[2, 1::2] = dhdx[0, :]
    return det, B


def stdm3NT(r, s, coord):
    """Strain-displacement interpolator B for a 3-noded triang element

    Parameters
    ----------
    r : float
      r component in the natural space.
    s : float
      s component in the natural space.
    coord : ndarray
      Coordinates of the nodes of the element (3, 2).

    Returns
    -------
    det : float
      Determinant evaluated at `(r, s)`.
    B : ndarray
      Strain-displacement interpolator evaluated at `(r, s)`.

    """
    nn = 3
    B = np.zeros((3, 2*nn))
    dhdx = np.array([
            [-1, 1, 0],
            [-1, 0, 1]])
    det, jaco_inv = jacoper(dhdx, coord)
    dhdx = np.dot(jaco_inv, dhdx)
    B[0, ::2] = dhdx[0, :]
    B[1, 1::2] = dhdx[1, :]
    B[2, ::2] = dhdx[1, :]
    B[2, 1::2] = dhdx[0, :]
    return det, B


def jacoper(dhdx, coord):
    """

    Parameters
    ----------
    dhdx : ndarray
      Derivatives of the interpolation function with respect to the
      natural coordinates.
    coord : ndarray
      Coordinates of the nodes of the element (nn, 2).

    Returns
    -------
    xja : ndarray (2, 2)
      Jacobian of the transformation evaluated at `(r, s)`.

    """
    jaco = dhdx.dot(coord)
    det = np.linalg.det(jaco)
    jaco_inv = np.linalg.inv(jaco)
    return det, jaco_inv

#%% Material routines
def umat(nu, E):
    """2D Elasticity consitutive matrix in plane stress

    For plane strain use effective properties.

    Parameters
    ----------
    nu : float
      Poisson coefficient (-1, 0.5).
    E : float
      Young modulus (>0).

    Returns
    -------
    C : ndarray
      Constitutive tensor in Voigt notation.

    Examples
    --------

    >>> C = umat(1/3, 8/3)
    >>> C_ex = np.array([
    ...    [3, 1, 0],
    ...    [1, 3, 0],
    ...    [0, 0, 1]])
    >>> np.allclose(C, C_ex)
    True

    """
    C = np.zeros((3, 3))
    enu = E/(1 - nu**2)
    mnu = (1 - nu)/2
    C[0, 0] = enu
    C[0, 1] = nu*enu
    C[1, 0] = C[0, 1]
    C[1, 1] = enu
    C[2, 2] = enu*mnu

    return C

#%% Elemental strains
def str_el4(coord, ul):
    """Compute the strains at each element integration point

    This one is used for 4-noded quadrilateral elements.

    Parameters
    ----------
    coord : ndarray
      Coordinates of the nodes of the element (4, 2).
    ul : ndarray
      Array with displacements for the element.

    Returns
    -------
    epsGT : ndarray
      Strain components for the Gauss points.
    xl : ndarray
      Configuration of the Gauss points after deformation.

    """
    epsl = np.zeros([3])
    epsG = np.zeros([3, 4])
    xl = np.zeros([4, 2])
    XW, XP = gau.gpoints2x2()
    for i in range(4):
        ri = XP[i, 0]
        si = XP[i, 1]
        ddet, B = stdm4NQ(ri, si, coord)
        epsl = np.dot(B, ul)
        epsG[:, i] = epsl[:]
        N = sha4(ri, si)
        xl[i, 0] = sum(N[0, 2*i]*coord[i, 0] for i in range(4))
        xl[i, 1] = sum(N[0, 2*i]*coord[i, 1] for i in range(4))
    return epsG.T, xl


def str_el6(coord, ul):
    """Compute the strains at each element integration point

    This one is used for 6-noded triangular elements.

    Parameters
    ----------
    coord : ndarray
      Coordinates of the nodes of the element (6, 2).
    ul : ndarray
      Array with displacements for the element.

    Returns
    -------
    epsGT : ndarray
      Strain components for the Gauss points.
    xl : ndarray
      Configuration of the Gauss points after deformation.

    """
    epsl = np.zeros([3])
    epsG = np.zeros([3, 7])
    xl = np.zeros([7, 2])
    XW, XP = gau.gpoints7()
    for i in range(7):
        ri = XP[i, 0]
        si = XP[i, 1]
        ddet, B = stdm6NT(ri, si, coord)
        epsl = np.dot(B, ul)
        epsG[:, i] = epsl[:]
        N = sha6(ri, si)
        xl[i, 0] = sum(N[0, 2*i]*coord[i, 0] for i in range(6))
        xl[i, 1] = sum(N[0, 2*i]*coord[i, 1] for i in range(6))
    return epsG.T, xl


def str_el3(coord, ul):
    """Compute the strains at each element integration point

    This one is used for 3-noded triangular elements.

    Parameters
    ----------
    coord : ndarray
      Coordinates of the nodes of the element (nn, 2).
    ul : ndarray
      Array with displacements for the element.

    Returns
    -------
    epsGT : ndarray
      Strain components for the Gauss points.
    xl : ndarray
      Configuration of the Gauss points after deformation.

    """
    epsl = np.zeros([3])
    epsG = np.zeros([3, 3])
    xl = np.zeros([3, 2])
    XW, XP = gau.gpoints3()
    for i in range(3):
        ri = XP[i, 0]
        si = XP[i, 1]
        ddet, B = stdm3NT(ri, si, coord)
        epsl = np.dot(B, ul)
        epsG[:, i] = epsl
        N = sha3(ri, si)
        xl[i, 0] = sum(N[0, 2*i]*coord[i, 0] for i in range(3))
        xl[i, 1] = sum(N[0, 2*i]*coord[i, 1] for i in range(3))
    return epsG.T, xl


if __name__ == "__main__":
    import doctest
    doctest.testmod()
