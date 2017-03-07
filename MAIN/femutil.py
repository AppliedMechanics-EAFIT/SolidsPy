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
import sympy as sym


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

    return ndof, nnodes, ngpts


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
    N : Matrix
      Matrix of interpolation functions.

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (1, 1). Thus

    >>> N = sha4(0, 0)
    >>> print(N[0, :])
    Matrix([[1/4, 0, 1/4, 0, 1/4, 0, 1/4, 0]])
    >>> print(N[1, :])
    Matrix([[0, 1/4, 0, 1/4, 0, 1/4, 0, 1/4]])

    and

    >>> N = sha4(1, 1)
    >>> print(N[0, :])
    Matrix([[0, 0, 0, 0, 1, 0, 0, 0]])
    >>> print(N[1, :])
    Matrix([[0, 0, 0, 0, 0, 1, 0, 0]])

    """
    N = sym.zeros(2, 8)
    H = sym.S(1)/4*sym.Matrix(
        [(1 - x)*(1 - y),
         (1 + x)*(1 - y),
         (1 + x)*(1 + y),
         (1 - x)*(1 + y)])

    for i in range(4):
        N[0, 2*i] = H[i]
        N[1, 2*i + 1] = H[i]

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
    N : Matrix
      Matrix of interpolation functions.

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (0, 0.5). Thus

    >>> N = sha6(0, 0)
    >>> print(N[0, :])
    Matrix([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> print(N[1, :])
    Matrix([[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    and

    >>> N = sha6(sym.S(1)/2, sym.S(1)/2)
    >>> print(N[0, :])
    Matrix([[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]])
    >>> print(N[1, :])
    Matrix([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]])

    """
    N = sym.zeros(2, 12)
    H = sym.Matrix(
        [(1 - x - y) - 2*x*(1 - x - y) - 2*y*(1 - x - y),
         x - 2*x*(1 - x - y) - 2*x*y,
         y - 2*x*y - 2*y*(1-x-y),
         4*x*(1 - x - y),
         4*x*y,
         4*y*(1 - x - y)])

    for i in range(6):
        N[0, 2*i] = H[i]
        N[1, 2*i + 1] = H[i]

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
    N : Matrix
      Matrix of interpolation functions.

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (0, 0.5). Thus

    >>> N = sha3(0, 0)
    >>> print(N[0, :])
    Matrix([[1, 0, 0, 0, 0, 0]])
    >>> print(N[1, :])
    Matrix([[0, 1, 0, 0, 0, 0]])

    and

    >>> N = sha3(sym.S(1)/2, sym.S(1)/2)
    >>> print(N[0, :])
    Matrix([[0, 0, 1/2, 0, 1/2, 0]])
    >>> print(N[1, :])
    Matrix([[0, 0, 0, 1/2, 0, 1/2]])

    """
    N = sym.zeros(2, 6)
    H = sym.Matrix(
        [(1 - x - y),
         x,
         y])

    for i in range(3):
        N[0, 2*i] = H[i]
        N[1, 2*i + 1] = H[i]

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
    rr, ss = sym.symbols('rr ss')
    nn = 4
    B = sym.zeros(3, 2*nn)
    dhdx = sym.zeros(2, nn)
    DNR = sym.zeros(2, nn)
    N = sym.S(1)/4*sym.Matrix(
        [(1 - rr)*(1 - ss),
         (1 + rr)*(1 - ss),
         (1 + rr)*(1 + ss),
         (1 - rr)*(1 + ss)])
    for i in range(nn):
        dhdx[0, i] = sym.diff(N[i], rr)
        dhdx[1, i] = sym.diff(N[i], ss)
    DNR = dhdx.subs([(rr, r), (ss, s)])

    xj = jacoper(DNR, coord, nn)
    ddet = np.linalg.det(xj)
    xi = np.linalg.inv(xj)
    aux1 = xi*DNR
    for i in range(nn):
        B[0, 2*i] = aux1[0, i]
        B[1, 2*i+1] = aux1[1, i]
        B[2, 2*i] = aux1[1, i]
        B[2, 2*i+1] = aux1[0, i]
    return ddet, B


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
    rr, ss = sym.symbols('rr ss')
    nn = 6
    N = sym.zeros(nn)
    B = sym.zeros(3, 2*nn)
    dhdx = sym.zeros(2, nn)
    DNR = sym.zeros(2, nn)
    N = sym.Matrix(
        [(1 - rr-ss) - 2*rr*(1 - rr - ss) - 2*ss*(1 - rr - ss),
         rr - 2*rr*(1 - rr - ss) - 2*rr*ss,
         ss - 2*rr*ss - 2*ss*(1 - rr - ss),
         4*rr*(1 - rr - ss),
         4*rr*ss,
         4*ss*(1 - rr - ss)])
    for i in range(nn):
        dhdx[0, i] = sym.diff(N[i], rr)
        dhdx[1, i] = sym.diff(N[i], ss)
    DNR = dhdx.subs([(rr, r), (ss, s)])

    xj = jacoper(DNR, coord, nn)
    ddet = np.linalg.det(xj)
    xi = np.linalg.inv(xj)
    aux1 = xi*DNR
    for i in range(nn):
        B[0, 2*i] = aux1[0, i]
        B[1, 2*i+1] = aux1[1, i]
        B[2, 2*i] = aux1[1, i]
        B[2, 2*i+1] = aux1[0, i]
    return ddet, B


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
    ddet : float
      Determinant evaluated at `(r, s)`.
    B : ndarray
      Strain-displacement interpolator evaluated at `(r, s)`.

    """
    rr, ss = sym.symbols('rr ss')
    nn = 3
    N = sym.zeros(nn)
    B = sym.zeros(3, 2*nn)
    dhdx = sym.zeros(2, nn)
    DNR = sym.zeros(2, nn)
    N = sym.Matrix(
        [(1 - rr - ss),
         rr,
         ss])
    for i in range(nn):
        dhdx[0, i] = sym.diff(N[i], rr)
        dhdx[1, i] = sym.diff(N[i], ss)
    DNR = dhdx.subs([(rr, r), (ss, s)])

    xj = jacoper(DNR, coord, nn)
    xi = np.linalg.inv(xj)
    ddet = np.linalg.det(xj)
    aux1 = xi*DNR
    for i in range(nn):
        B[0, 2*i] = aux1[0, i]
        B[1, 2*i+1] = aux1[1, i]
        B[2, 2*i] = aux1[1, i]
        B[2, 2*i+1] = aux1[0, i]
    return ddet, B


def jacoper(dhdx, coord, nn):
    """

    Parameters
    ----------
    dhdx : ndarray
      Derivatives of the interpolation function with respect to the
      natural coordinates.
    coord : ndarray
      Coordinates of the nodes of the element (nn, 2).
    nn : int
      Number of nodes in the element

    Returns
    -------
    xja : ndarray (2, 2)
      Jacobian of the transformation evaluated at `(r, s)`.

    """
    dum = 0.0
    xja = np.zeros([2, 2])
    for k in range(2):
        for j in range(2):
            for i in range(nn):
                dum = dum + dhdx[k, i]*coord[i, j]
            xja[k, j] = dum
            dum = 0.0

    return xja


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

    >>> C = umat(sym.S(1)/3, sym.S(8)/3)
    >>> print(C)
    Matrix([[3, 1, 0], [1, 3, 0], [0, 0, 1]])

    """
    C = sym.zeros(3, 3)
    enu = E/(1 - nu**2)
    mnu = (1 - nu)/2
    C[0, 0] = enu
    C[0, 1] = nu*enu
    C[1, 0] = C[0, 1]
    C[1, 1] = enu
    C[2, 2] = enu*mnu

    return C


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
    epsGT = np.zeros([4, 3])
    xl = np.zeros([4, 2])
    x, y = sym.symbols('x y')

    XW, XP = gau.gpoints2x2()
    for i in range(4):
        ri = XP[i, 0]
        si = XP[i, 1]
        ddet, B = stdm4NQ(ri, si, coord)
        epsl = B*ul
        epsG[:, i] = epsl[:]
        N = sha4(ri, si)
        NN = N.subs([(x, ri), (y, si)])
        xl[i, 0] = sum(NN[0, 2*i]*coord[i, 0] for i in range(4))
        xl[i, 1] = sum(NN[0, 2*i]*coord[i, 1] for i in range(4))
    epsGT = epsG.T
    return epsGT, xl


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
    epsGT = np.zeros([7, 3])
    xl = np.zeros([7, 2])
    x, y = sym.symbols('x y')

    XW, XP = gau.gpoints7()
    for i in range(7):
        ri = XP[i, 0]
        si = XP[i, 1]
        ddet, B = stdm6NT(ri, si, coord)
        epsl = B*ul
        epsG[:, i] = epsl[:]
        N = sha6(ri, si)
        NN = N.subs([(x, ri), (y, si)])
        xl[i, 0] = sum(NN[0, 2*i]*coord[i, 0] for i in range(6))
        xl[i, 1] = sum(NN[0, 2*i]*coord[i, 1] for i in range(6))
    epsGT = epsG.T
    return epsGT, xl


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
    epsGT = np.zeros([3, 3])
    xl = np.zeros([3, 2])
    x, y = sym.symbols('x y')
    XW, XP = gau.gpoints3()
    for i in range(3):
        ri = XP[i, 0]
        si = XP[i, 1]
        ddet, B = stdm3NT(ri, si, coord)
        epsl = B*ul
        epsG[:, i] = epsl[:]
        N = sha3(ri, si)
        NN = N.subs([(x, ri), (y, si)])
        xl[i, 0] = sum(NN[0, 2*i]*coord[i, 0] for i in range(3))
        xl[i, 1] = sum(NN[0, 2*i]*coord[i, 1] for i in range(3))
    epsGT = epsG.T
    return epsGT, xl
