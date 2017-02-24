# -*- coding: utf-8 -*-
"""
Element subroutines
-------------------
Each UEL subroutine computes the local stiffness matrix for a given
finite element.

New elements can be added by including additional subroutines.

"""
from __future__ import division
import numpy as np
import femutil as fem
import gaussutil as gau
from sympy import Matrix, S

def uel4nquad(coord, enu, Emod):
    """Quadrilateral element with 4 nodes

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (4, 2).
    enu : float
      Poisson coefficient (-1, 0.5).
    Emod : float
      Young modulus (>0).

    Returns
    -------
    kl : ndarray
      Local stiffness matrix for the element (8, 8).

    Examples
    --------
    
    >>> coord = Matrix([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    >>> stiff = uel4nquad(coord, S(1)/3, S(8)/3)
    >>> stiff_ex = 1/6 * Matrix([
    ...             [ 8,  3, -5,  0, -4, -3,  1,  0],
    ...             [ 3,  8,  0,  1, -3, -4,  0, -5],
    ...             [-5,  0,  8, -3,  1,  0, -4,  3],
    ...             [ 0,  1, -3,  8,  0, -5,  3, -4],
    ...             [-4, -3,  1,  0,  8,  3, -5,  0],
    ...             [-3, -4,  0, -5,  3,  8,  0,  1],
    ...             [ 1,  0, -4,  3, -5,  0,  8, -3],
    ...             [ 0, -5,  3, -4,  0,  1, -3,  8]])
    >>> (stiff - stiff_ex).norm()/stiff_ex.norm() < 1e-6
    True

    """
    kl = np.zeros([8, 8])
    C = fem.umat(enu, Emod)
    XW, XP = gau.gpoints2x2()
    ngpts = 4
    for i in range(0, ngpts):
        ri = XP[i, 0]
        si = XP[i, 1]
        alf = XW[i]
        ddet, B = fem.stdm4NQ(ri, si, coord)
        kl = kl + B.T*C*B*alf*ddet
    return kl


def uel6ntrian(coord, enu, Emod):
    """Triangular element with 6 nodes

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (6, 2).
    enu : float
      Poisson coefficient (-1, 0.5).
    Emod : float
      Young modulus (>0).

    Returns
    -------
    kl : ndarray
      Local stiffness matrix for the element (12, 12).
 
    Examples
    --------
    
    >>> coord = Matrix([
    ...         [0, 0],
    ...         [1, 0],
    ...         [0, 1],
    ...         [0.5, 0],
    ...         [0.5, 0.5],
    ...         [0, 0.5]])
    >>> stiff = uel6ntrian(coord, S(1)/3, S(8)/3)
    >>> stiff_ex = 1/6 * Matrix([
    ...            [12, 6, 3, 1, 1, 1, -12, -4, 0, 0, -4, -4],
    ...            [6, 12, 1, 1, 1, 3, -4, -4, 0, 0, -4, -12],
    ...            [3, 1, 9, 0, 0, -1, -12, -4, 0, 4, 0, 0],
    ...            [1, 1, 0, 3, -1, 0, -4, -4, 4, 0, 0, 0],
    ...            [1, 1, 0, -1, 3, 0, 0, 0, 0, 4, -4, -4],
    ...            [1, 3, -1, 0, 0, 9, 0, 0, 4, 0, -4, -12],
    ...            [-12, -4, -12, -4, 0, 0, 32, 8, -8, -8, 0, 8],
    ...            [-4, -4, -4, -4, 0, 0, 8, 32, -8, -24, 8, 0],
    ...            [0, 0, 0, 4, 0, 4, -8, -8, 32, 8, -24, -8],
    ...            [0, 0, 4, 0, 4, 0, -8, -24, 8, 32, -8, -8],
    ...            [-4, -4, 0, 0, -4, -4, 0, 8, -24, -8, 32, 8],
    ...            [-4, -12, 0, 0, -4, -12, 8, 0, -8, -8, 8, 32]])
    >>> (stiff - stiff_ex).norm()/stiff_ex.norm() < 1e-6
    True
    
    """
    kl = np.zeros([12, 12])
    C = fem.umat(enu, Emod)
    XW, XP = gau.gpoints7()
    ngpts = 7
    for i in range(ngpts):
        ri = XP[i, 0]
        si = XP[i, 1]
        alf = XW[i]
        ddet, B = fem.stdm6NT(ri, si, coord)
        kl = kl + 0.5*B.T*C*B*alf*ddet
    return kl


def uel3ntrian(coord, enu, Emod):
    """Triangular element with 3 nodes

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (3, 2).
    enu : float
      Poisson coefficient (-1, 0.5).
    Emod : float
      Young modulus (>0).

    Returns
    -------
    kl : ndarray
      Local stiffness matrix for the element (6, 6).

    Examples
    --------
    
    >>> coord = Matrix([
    ...         [0, 0],
    ...         [1, 0],
    ...         [0, 1]])
    >>> stiff = uel3ntrian(coord, S(1)/3, S(8)/3)
    >>> stiff_ex = 1/2 * Matrix([
    ...            [4, 2, -3, -1, -1, -1],
    ...            [2, 4, -1, -1, -1, -3],
    ...            [-3, -1, 3, 0, 0, 1],
    ...            [-1, -1, 0, 1, 1, 0],
    ...            [-1, -1, 0, 1, 1, 0],
    ...            [-1, -3, 1, 0, 0, 3]])
    >>> (stiff - stiff_ex).norm()/stiff_ex.norm() < 1e-6
    True

    """
    kl = np.zeros([6, 6])
    C = fem.umat(enu, Emod)
    XW, XP = gau.gpoints3()
    ngpts = 3
    for i in range(ngpts):
        ri = XP[i, 0]
        si = XP[i, 1]
        alf = XW[i]
        ddet, B = fem.stdm3NT(ri, si, coord)
        kl = kl + 0.5*B.T*C*B*alf*ddet
    return kl


def uelspring(coord, enu, Emod):
    """1D-2-noded Spring element

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (3, 2).
    enu : float
      Fictitious parameter.
    Emod : float
      Stiffness coefficient (>0).

    Returns
    -------
    kl : ndarray
      Local stiffness matrix for the element (2, 2).


    """
    kl = np.zeros([4, 4])
    kl[0,0]= Emod
    kl[0,2]=-Emod
    kl[2,0]=-Emod
    kl[2,2]= Emod
    return kl




if __name__ == "__main__":
    import doctest
    doctest.testmod()