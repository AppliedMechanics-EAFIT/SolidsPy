# -*- coding: utf-8 -*-
"""
Element subroutines
-------------------
Each UEL subroutine computes the local stiffness matrix for a given
finite element.

New elements can be added by including additional subroutines.

"""
from __future__ import division, print_function
import numpy as np
import femutil as fem
import gaussutil as gau


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

    >>> coord = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    >>> stiff = uel4nquad(coord, 1/3, 8/3)
    >>> stiff_ex = 1/6 * np.array([
    ...             [ 8,  3, -5,  0, -4, -3,  1,  0],
    ...             [ 3,  8,  0,  1, -3, -4,  0, -5],
    ...             [-5,  0,  8, -3,  1,  0, -4,  3],
    ...             [ 0,  1, -3,  8,  0, -5,  3, -4],
    ...             [-4, -3,  1,  0,  8,  3, -5,  0],
    ...             [-3, -4,  0, -5,  3,  8,  0,  1],
    ...             [ 1,  0, -4,  3, -5,  0,  8, -3],
    ...             [ 0, -5,  3, -4,  0,  1, -3,  8]])
    >>> np.allclose(stiff, stiff_ex)
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
        kl = kl + np.dot(np.dot(B.T,C), B)*alf*ddet
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

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0],
    ...         [0, 1],
    ...         [0.5, 0],
    ...         [0.5, 0.5],
    ...         [0, 0.5]])
    >>> stiff = uel6ntrian(coord,1/3, 8/3)
    >>> stiff_ex = 1/6 * np.array([
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
    >>> np.allclose(stiff, stiff_ex)
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
        kl = kl + 0.5*np.dot(np.dot(B.T,C), B)*alf*ddet
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

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0],
    ...         [0, 1]])
    >>> stiff = uel3ntrian(coord, 1/3, 8/3)
    >>> stiff_ex = 1/2 * np.array([
    ...            [4, 2, -3, -1, -1, -1],
    ...            [2, 4, -1, -1, -1, -3],
    ...            [-3, -1, 3, 0, 0, 1],
    ...            [-1, -1, 0, 1, 1, 0],
    ...            [-1, -1, 0, 1, 1, 0],
    ...            [-1, -3, 1, 0, 0, 3]])
    >>> np.allclose(stiff, stiff_ex)
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
        kl = kl + 0.5*np.dot(np.dot(B.T,C), B)*alf*ddet
    return kl


def uelspring(coord, enu, Emod):
    """1D-2-noded Spring element

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    enu : float
      Fictitious parameter.
    Emod : float
      Stiffness coefficient (>0).

    Returns
    -------
    kl : ndarray
      Local stiffness matrix for the element (4, 4).

    Examples
    --------

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0]])
    >>> stiff = uelspring(coord, 1/3, 8/3)
    >>> stiff_ex = 8/3 * np.array([
    ...    [-1, 0, 1, 0],
    ...    [0, 0, 0, 0],
    ...    [1, 0, -1, 0],
    ...    [0, 0, 0, 0]])
    >>> np.allclose(stiff, stiff_ex)
    True

    """
    vec = coord[1, :] - coord[0, :]
    nx = vec[0]/np.linalg.norm(vec)
    ny = vec[1]/np.linalg.norm(vec)
    Q = np.array([
        [nx, ny , 0 , 0],
        [0,  0, nx , ny]])
    kl = Emod * np.array([
        [1, -1],
        [-1, 1]])
    kG = np.dot(np.dot(Q.T, kl), Q)
    return kG


def ueltruss2D(coord, A, Emod):
    """2D-2-noded truss element

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    A : float
      Cross section area.
    Emod : float
      Young modulus (>0).

    Returns
    -------
    kl : ndarray
      Local stiffness matrix for the element (4, 4).

    Examples
    --------

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0]])
    >>> stiff = ueltruss2D(coord, 1.0 , 1000.0)
    >>> stiff_ex = 8/3 * np.array([
    ...    [-1, 0, 1, 0],
    ...    [0, 0, 0, 0],
    ...    [1, 0, -1, 0],
    ...    [0, 0, 0, 0]])
    >>> np.allclose(stiff, stiff_ex)
    True

    """
    vec = coord[1, :] - coord[0, :]
    nx = vec[0]/np.linalg.norm(vec)
    ny = vec[1]/np.linalg.norm(vec)
    length = np.linalg.norm(vec) 
    Q = np.array([
        [nx, ny , 0 , 0],
        [0,  0, nx , ny]])
    kl =(A*Emod/length) * np.array([
        [1, -1],
        [-1, 1]])
    kG = np.dot(np.dot(Q.T, kl), Q)
    return kG

def uelbeam2DU(coord, I , Emod):
    """2D-2-noded bam element
       without axial dformation

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    A : float
      Cross section area.
    Emod : float
      Young modulus (>0).

    Returns
    -------
    kl : ndarray
      Local stiffness matrix for the element (4, 4).

    Examples
    --------

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0]])
    >>> stiff = ueltruss2D(coord, 1.0 , 1000.0)
    >>> stiff_ex = 8/3 * np.array([
    ...    [-1, 0, 1, 0],
    ...    [0, 0, 0, 0],
    ...    [1, 0, -1, 0],
    ...    [0, 0, 0, 0]])
    >>> np.allclose(stiff, stiff_ex)
    True

    """
    vec = coord[1, :] - coord[0, :]
    nx = vec[0]/np.linalg.norm(vec)
    ny = vec[1]/np.linalg.norm(vec)
    L = np.linalg.norm(vec) 
    Q = np.array([
        [-ny , nx ,   0 ,  0 ,  0 , 0 ],
        [  0 ,  0 , 1.0 ,  0 ,  0 , 0 ],
        [  0 ,  0 ,   0 ,-ny , nx , 0 ],
        [  0 ,  0  ,  0 ,  0 ,  0 , 1.0 ]])
    kl =(I*Emod/(L*L*L)) * np.array([
        [12.0, 6*L , -12.0 , 6*L],
        [6*L,  4*L*L , -6*L , 2*L*L],
        [-12.0,  -6*L , 12.0 , -6*L],
        [6*L,  2*L*L , -6*L , 4*L*L]])
    kG = np.dot(np.dot(Q.T, kl), Q)
    return kG

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    