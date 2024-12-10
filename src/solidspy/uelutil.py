# -*- coding: utf-8 -*-
"""
Element subroutines
-------------------
Each UEL subroutine computes the local stiffness matrix for a given
finite element.

New elements can be added by including additional subroutines.

"""
import numpy as np
from . import femutil as fem
from . import gaussutil as gau
from numpy import ndarray
from typing import Tuple, Sequence


#%% Continuum elements

# Triangles
def elast_tri3(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """Triangular element with 3 nodes

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (3, 2).
    params : Sequence[float]
      Material parameters in the following order:

          young : float
            Young modulus (>0).
          poisson : float
            Poisson coefficient (-1, 0.5).
          dens : float, optional
            Density (>0).

    Returns
    -------
    stiff_mat : ndarray
      Local stiffness matrix for the element (6, 6).
    mass_mat : ndarray
      Local mass matrix for the element (6, 6).

    Examples
    --------

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0],
    ...         [0, 1]])
    >>> params = [8/3, 1/3]
    >>> stiff, mass = elast_tri3(coord, params)
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
    stiff_mat = np.zeros([6, 6])
    mass_mat = np.zeros([6, 6])
    C = fem.umat(params[:2])
    if len(params) == 2:
        dens = 1
    else:
        dens = params[-1]
    gpts, gwts = gau.gauss_tri(order=2)
    for cont in range(gpts.shape[0]):
        r, s = gpts[cont, :]
        H, B, det = fem.elast_diff_2d(r, s, coord, fem.shape_tri3)
        factor = det * gwts[cont]
        stiff_mat += 0.5 * factor * (B.T @ C @ B)
        mass_mat += 0.5 * dens * factor * (H.T @ H)
    return stiff_mat, mass_mat


def elast_tri6(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """Triangular element with 6 nodes

    Parameters
    ----------
    coord : ndarray
        Coordinates for the nodes of the element (6, 2).
    params : Sequence[float]
        Material parameters in the following order:

            young : float
                Young modulus (>0).
            poisson : float
                Poisson coefficient (-1, 0.5).
            dens : float, optional
                Density (>0).

    Returns
    -------
    stiff_mat : ndarray
        Local stiffness matrix for the element (12, 12).
    mass_mat : ndarray
        Local mass matrix for the element (12, 12).

    Examples
    --------

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0],
    ...         [0, 1],
    ...         [0.5, 0],
    ...         [0.5, 0.5],
    ...         [0, 0.5]])
    >>> params = [8/3, 1/3]
    >>> stiff, mass = elast_tri6(coord, params)
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
    stiff_mat = np.zeros([12, 12])
    mass_mat = np.zeros([12, 12])
    C = fem.umat(params[:2])
    if len(params) == 2:
        dens = 1
    else:
        dens = params[-1]
    gpts, gwts = gau.gauss_tri(order=3)
    for cont in range(gpts.shape[0]):
        r, s = gpts[cont, :]
        H, B, det = fem.elast_diff_2d(r, s, coord, fem.shape_tri6)
        factor = gwts[cont] * det
        stiff_mat += 0.5 * factor * (B.T @ C @ B)
        mass_mat += 0.5 * dens * factor * (H.T @ H)
    return stiff_mat, mass_mat


# Quadrilaterals
def elast_quad4(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """Quadrilateral element with 4 nodes

    Parameters
    ----------
    coord : ndarray
        Coordinates for the nodes of the element (4, 2).
    params : Sequence[float]
        Material parameters in the following order:

            young : float
                Young modulus (>0).
            poisson : float
                Poisson coefficient (-1, 0.5).
            dens : float, optional
                Density (>0).

    Returns
    -------
    stiff_mat : ndarray
        Local stiffness matrix for the element (8, 8).
    mass_mat : ndarray
        Local mass matrix for the element (8, 8).

    Examples
    --------

    >>> coord = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    >>> params = [8/3, 1/3, 1]
    >>> stiff, mass = elast_quad4(coord, params)
    >>> stiff_ex = 1/6 * np.array([
    ...             [ 8,  3, -5,  0, -4, -3,  1,  0],
    ...             [ 3,  8,  0,  1, -3, -4,  0, -5],
    ...             [-5,  0,  8, -3,  1,  0, -4,  3],
    ...             [ 0,  1, -3,  8,  0, -5,  3, -4],
    ...             [-4, -3,  1,  0,  8,  3, -5,  0],
    ...             [-3, -4,  0, -5,  3,  8,  0,  1],
    ...             [ 1,  0, -4,  3, -5,  0,  8, -3],
    ...             [ 0, -5,  3, -4,  0,  1, -3,  8]])
    >>> mass_ex = 1/9 * np.array([
    ...             [4, 0, 2, 0, 1, 0, 2, 0],
    ...             [0, 4, 0, 2, 0, 1, 0, 2],
    ...             [2, 0, 4, 0, 2, 0, 1, 0],
    ...             [0, 2, 0, 4, 0, 2, 0, 1],
    ...             [1, 0, 2, 0, 4, 0, 2, 0],
    ...             [0, 1, 0, 2, 0, 4, 0, 2],
    ...             [2, 0, 1, 0, 2, 0, 4, 0],
    ...             [0, 2, 0, 1, 0, 2, 0, 4]])
    >>> np.allclose(stiff, stiff_ex)
    True
    >>> np.allclose(mass, mass_ex)
    True

    """
    stiff_mat = np.zeros([8, 8])
    mass_mat = np.zeros([8, 8])
    C = fem.umat(params[:2])
    if len(params) == 2:
        dens = 1
    else:
        dens = params[-1]
    gpts, gwts = gau.gauss_nd(2)
    for cont in range(gpts.shape[0]): # pylint: disable=E1136  # pylint/issues/3139
        r, s = gpts[cont, :]
        H, B, det = fem.elast_diff_2d(r, s, coord, fem.shape_quad4)
        factor = det * gwts[cont]
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += dens * factor * (H.T @ H)
    return stiff_mat, mass_mat


def elast_quad9(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """
    Quadrilateral element with 9 nodes for classic elasticity
    under plane-strain

    Parameters
    ----------
    coord : ndarray
        Coordinates of the element.
    params : Sequence[float]
        List with material parameters in the following order:
        [Young modulus, Poisson coefficient, density].

    Returns
    -------
    stiff_mat : ndarray (float)
        Local stiffness matrix.
    mass_mat : ndarray (float)
        Local mass matrix.
    """
    stiff_mat = np.zeros((18, 18))
    mass_mat = np.zeros((18, 18))
    C = fem.umat(params[:2])
    if len(params) == 2:
        dens = 1
    else:
        dens = params[-1]
    gpts, gwts = gau.gauss_nd(3, ndim=2)
    for cont in range(gpts.shape[0]): # pylint: disable=E1136  # pylint/issues/3139
        r, s = gpts[cont, :]
        H, B, det = fem.elast_diff_2d(r, s, coord, fem.shape_quad9)
        factor = det * gwts[cont]
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += dens * factor * (H.T @ H)
    return stiff_mat, mass_mat


def elast_quad8(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """
    Quadrilateral element with 8 nodes for classic elasticity
    under plane-strain

    Parameters
    ----------
    coord : ndarray
        Coordinates of the element.
    params : Sequence[float]
        List with material parameters in the following order:
        [Young modulus, Poisson coefficient, density].

    Returns
    -------
    stiff_mat : ndarray (float)
        Local stiffness matrix.
    mass_mat : ndarray (float)
        Local mass matrix.
    """
    
    E, nu, rho = params
    C = fem.umat((E, nu))
    stiff_mat = np.zeros((16, 16))
    mass_mat = np.zeros((16, 16))
    gpts, gwts = gau.gauss_nd(3, ndim=2)
    for cont in range(gpts.shape[0]):
        r = gpts[cont, 0]
        s = gpts[cont, 1]
        H, B, det = fem.elast_diff_2d(r, s, coord, fem.shape_quad8)
        factor = det * gwts[cont]
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += rho * factor * (H.T @ H)
    return stiff_mat, mass_mat


## Axisymmetric
def elast_axi_quad9(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """
    Quadrilateral element with 9 nodes for classic elasticity
    under axisymmetric conditions

    Parameters
    ----------
    coord : ndarray
        Coordinates of the element.
    params : Sequence[float]
        List with material parameters in the following order:
        [Young modulus, Poisson coefficient, density].

    Returns
    -------
    stiff_mat : ndarray (float)
        Local stiffness matrix.
    mass_mat : ndarray (float)
        Local mass matrix.
    """
    
    E, nu, rho = params
    C = fem.elast_mat_axi((E, nu))
    stiff_mat = np.zeros((18, 18))
    mass_mat = np.zeros((18, 18))
    gpts, gwts = gau.gauss_nd(3, ndim=2)
    for cont in range(gpts.shape[0]):
        r = gpts[cont, 0]
        s = gpts[cont, 1]
        H, B, det = fem.elast_diff_axi(r, s, coord, fem.shape_quad9)
        factor = det * gwts[cont]
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += rho * factor * (H.T @ H)
    return stiff_mat, mass_mat


## 3D elements
def elast_tet4(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """Tetraedral element with 4 nodes for classic elasticity

    Parameters
    ----------
    coord : ndarray
        Coordinates of the element.
    params : Sequence[float]
        List with material parameters in the following order:
        [Young modulus, Poisson coefficient, density].

    Returns
    -------
    stiff_mat : ndarray (float)
        Local stiffness matrix.
    mass_mat : ndarray (float)
        Local mass matrix.
    """
    E, nu, rho = params
    C = fem.elast_mat((E, nu))
    stiff_mat = np.zeros((12, 12))
    mass_mat = np.zeros((12, 12))
    gpts, gwts = gau.gauss_tet(3)
    for cont in range(gpts.shape[0]):
        r = gpts[cont, 0]
        s = gpts[cont, 1]
        t = gpts[cont, 2]
        H, B, det = fem.elast_diff_3d(r, s, t, coord, fem.shape_tet4)
        factor = det * gwts[cont] / 6
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += rho * factor * (H.T @ H)
    return stiff_mat, mass_mat


def elast_hex8(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """Hexaedral element with 8 nodes for classic elasticity

    Parameters
    ----------
    coord : ndarray
        Coordinates of the element.
    params : Sequence[float]
        List with material parameters in the following order:
        [Young modulus, Poisson coefficient, density].

    Returns
    -------
    stiff_mat : ndarray (float)
        Local stiffness matrix.
    mass_mat : ndarray (float)
        Local mass matrix.
    """
    E, nu, rho = params
    C = fem.elast_mat((E, nu))
    stiff_mat = np.zeros((24, 24))
    mass_mat = np.zeros((24, 24))
    gpts, gwts = gau.gauss_nd(2, ndim=3)
    for cont in range(gpts.shape[0]):
        r = gpts[cont, 0]
        s = gpts[cont, 1]
        t = gpts[cont, 2]
        H, B, det = fem.elast_diff_3d(r, s, t, coord)
        factor = det * gwts[cont]
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += rho * factor * (H.T @ H)
    return stiff_mat, mass_mat


#%% Structural elements
def spring(
    coord: ndarray, 
    stiff: float
) -> Tuple[ndarray, ndarray]:
    """1D 2-noded Spring element

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    stiff : float
      Spring stiffness (>0).

    Returns
    -------
    stiff_mat : ndarray
      Local stiffness matrix for the element (4, 4).
    mass_mat : ndarray
      Local mass matrix for the element (4, 4). For now it
      returns a zero matrix and it is used for compatibility
      reasons.

    Examples
    --------

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0]])
    >>> stiff = 8/3
    >>> stiff, mass = spring(coord, stiff)
    >>> stiff_ex = 8/3 * np.array([
    ...    [1, 0, -1, 0],
    ...    [0, 0, 0, 0],
    ...    [-1, 0, 1, 0],
    ...    [0, 0, 0, 0]])
    >>> np.allclose(stiff, stiff_ex)
    True

    """
    vec = coord[1, :] - coord[0, :]
    nx = vec[0]/np.linalg.norm(vec)
    ny = vec[1]/np.linalg.norm(vec)
    Q = np.array([
        [nx, ny, 0, 0],
        [0, 0, nx, ny]])
    stiff_mat = stiff * np.array([
        [1, -1],
        [-1, 1]])
    stiff_mat = Q.T @ stiff_mat @ Q
    return stiff_mat, np.zeros(4)


def truss2D(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """2D 2-noded truss element

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    params : Sequence[float]
      Element parameters in the following order:

          young : float
            Young modulus (>0).
          area : float
            Cross-sectional area (>0).
          dens : float, optional
            Density (>0).

    Returns
    -------
    stiff_mat : ndarray
      Local stiffness matrix for the element (4, 4).
    mass_mat : ndarray
      Local mass matrix for the element (4, 4).

    Examples
    --------

    >>> coord = np.array([
    ...         [0, 0],
    ...         [1, 0]])
    >>> params = [1.0 , 1.0]
    >>> stiff, mass = truss2D(coord, params)
    >>> stiff_ex =  np.array([
    ...    [1, 0, -1, 0],
    ...    [0, 0, 0, 0],
    ...    [-1, 0, 1, 0],
    ...    [0, 0, 0, 0]])
    >>> np.allclose(stiff, stiff_ex)
    True

    """
    vec = coord[1, :] - coord[0, :]
    length = np.linalg.norm(vec)
    nx = vec[0]/length
    ny = vec[1]/length
    Q = np.array([
        [nx, ny, 0, 0],
        [0, 0, nx, ny]])
    young, area  = params[:2]
    stiff = area * young / length
    stiff_mat = stiff * np.array([
        [1, -1],
        [-1, 1]])
    stiff_mat = Q.T @ stiff_mat @ Q
    if len(params) == 2:
        dens = 1.0
    else:
        dens = params[-1]
    mass = area * length * dens
    mass_mat = mass / 6 * np.array([
            [2, 0, 1, 0],
            [0, 2, 0, 1],
            [1, 0, 2, 0],
            [0, 1, 0, 2]])
    return stiff_mat, mass_mat


def beam2DU(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """2D 2-noded beam element without axial deformation

    Parameters
    ----------
    coord : ndarray
        Coordinates for the nodes of the element (2, 2).
    params : Sequence[float]
        Element parameters in the following order:
    
            young : float
              Young modulus (>0).
            area_moment : float
              Second moment of area (>0).
            dens : float
              Density (>0).
            area : float, optional
              Cross-sectional area (>0).

    Returns
    -------
    stiff_mat : ndarray
      Local stiffness matrix for the element (6, 6).
    mass_mat : ndarray
      Local mass matrix for the element (6, 6).
    """
    vec = coord[1, :] - coord[0, :]
    nx = vec[0]/np.linalg.norm(vec)
    ny = vec[1]/np.linalg.norm(vec)
    L = np.linalg.norm(vec)
    Q = np.array([
        [ny, nx, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, ny, nx, 0],
        [0, 0, 0, 0, 0, 1]])
    young, area_moment = params[:2]
    bending_stiff = area_moment * young / L**3
    stiff_mat = bending_stiff * np.array([
        [12, 6*L, -12, 6*L],
        [6*L, 4*L**2, -6*L, 2*L**2],
        [-12, -6*L, 12, -6*L],
        [6*L, 2*L**2, -6*L, 4*L**2]])
    if len(params) == 2:
        dens = 1.0
        area = 1.0
    else:
        dens, area = params[2:]
    mass = area * L * dens
    mass_mat = mass / 420 * np.array([
            [156, 22*L, 54, -13*L],
            [22*L, 4*L**2, 13*L, -3*L**2],
            [54, 13*L, 156, -22*L],
            [-13*L, -3*L**2, -22*L, 4*L**2]])
    stiff_mat = Q.T @ stiff_mat @ Q
    mass_mat = Q.T @ mass_mat @ Q
    return stiff_mat, mass_mat


def beam2D(
    coord: ndarray, 
    params: Sequence[float]
) -> Tuple[ndarray, ndarray]:
    """2D 2-noded beam element with axial deformation

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    params : Sequence[float]
        Element parameters in the following order:
    
            young : float
              Young modulus (>0).
            area_moment : float
              Second moment of area (>0).
            area : float, optional
              Cross-sectional area (>0).
            dens : float
              Density (>0).

    Returns
    -------
    stiff_mat : ndarray
      Local stiffness matrix for the element (6, 6).
    mass_mat : ndarray
      Local mass matrix for the element (6, 6).
    """
    vec = coord[1, :] - coord[0, :]
    nx = vec[0]/np.linalg.norm(vec)
    ny = vec[1]/np.linalg.norm(vec)
    L = np.linalg.norm(vec)
    Q = np.array([
        [nx, -ny, 0, 0, 0, 0],
        [ny, nx, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, nx, -ny, 0],
        [0, 0, 0, ny, nx, 0],
        [0, 0, 0, 0, 0, 1]])
    young, area_moment, area = params[:3]
    bending_stiff = area_moment * young
    ratio = area / area_moment
    stiff_mat = bending_stiff / L**3 * np.array([
        [ratio*L**2, 0, 0, -ratio*L**2, 0, 0],
        [0, 12, 6*L, 0, -12, 6*L],
        [0, 6*L, 4*L**2, 0, -6*L, 2*L**2],
        [-ratio*L**2, 0, 0, ratio*L**2, 0, 0],
        [0, -12, -6*L, 0, 12, -6*L],
        [0, 6*L, 2*L**2, 0, -6*L, 4*L**2]])
    
    if len(params) == 3:
        dens = 1.0
    else:
        dens = params[3:]
    mass = area * L * dens
    mass_mat = mass / 420 * np.array([
            [140, 0, 0, 70, 0, 0],
            [0, 156, 22*L, 0, 54, -13*L],
            [0, 22*L, 4*L**2, 0, 13*L, -3*L**2],
            [70, 0, 0, 140, 0, 0],
            [0, 54, 13*L, 0, 156, -22*L],
            [0, -13*L, -3*L**2, 0, -22*L, 4*L**2]])
    stiff_mat = Q.T @ stiff_mat @ Q
    mass_mat = Q.T @ mass_mat @ Q
    return stiff_mat, mass_mat


if __name__ == "__main__":
    import doctest
    doctest.testmod()
