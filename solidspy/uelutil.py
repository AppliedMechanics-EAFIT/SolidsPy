# -*- coding: utf-8 -*-
"""
Element subroutines
-------------------
Each UEL subroutine computes the local stiffness matrix for a given
finite element.

New elements can be added by including additional subroutines.

"""
import numpy as np
import solidspy.femutil as fem
import solidspy.gaussutil as gau


#%% Continuum elements
def uel4nquad(coord, params):
    """Quadrilateral element with 4 nodes

    Parameters
    ----------
    coord : ndarray
        Coordinates for the nodes of the element (4, 2).
    params : tuple
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
    >>> stiff, mass = uel4nquad(coord, params)
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
    for i in range(gpts.shape[0]):
        ri = gpts[i, 0]
        si = gpts[i, 1]
        ddet, B = fem.stdm4NQ(ri, si, coord)
        N = fem.sha4(ri, si)
        factor = ddet * gwts[i]
        stiff_mat += factor * (B.T @ C @ B)
        mass_mat += dens*factor* (N.T @ N)
    return stiff_mat, mass_mat


def uel6ntrian(coord, params):
    """Triangular element with 6 nodes

    Parameters
    ----------
    coord : ndarray
        Coordinates for the nodes of the element (6, 2).
    params : tuple
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
    >>> stiff, mass = uel6ntrian(coord, params)
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
    for i in range(gpts.shape[0]):
        ri = gpts[i, 0]
        si = gpts[i, 1]
        ddet, B = fem.stdm6NT(ri, si, coord)
        N = fem.sha6(ri, si)
        factor = gwts[i] * ddet
        stiff_mat += 0.5*factor*(B.T @ C @ B)
        mass_mat += dens*factor* (N.T @ N)
    return stiff_mat, mass_mat


def uel3ntrian(coord, params):
    """Triangular element with 3 nodes

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (3, 2).
    params : tuple
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
    >>> stiff, mass = uel3ntrian(coord, params)
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
    gpts, gwts = gau.gauss_tri(order=1)
    for i in range(gpts.shape[0]):
        ri = gpts[i, 0]
        si = gpts[i, 1]
        ddet, B = fem.stdm3NT(ri, si, coord)
        N = fem.sha3(ri, si)
        factor = ddet * gwts[i]
        stiff_mat += 0.5*factor*(B.T @ C @ B)
        mass_mat += dens*factor* (N.T @ N)
    return stiff_mat, mass_mat


#%% Structural elements
def uelspring(coord, stiff):
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
    >>> stiff, mass = uelspring(coord, 8/3)
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


def ueltruss2D(coord, params):
    """2D 2-noded truss element

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    params : tuple
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
    >>> stiff, mass = ueltruss2D(coord, params)
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
    stiff = area*young/length
    stiff_mat = stiff * np.array([
        [1, -1],
        [-1, 1]])
    stiff_mat = Q.T @ stiff_mat @ Q
    if len(params) == 2:
        dens = 1.0
    else:
        dens = params[-1]
    mass = area * length * dens
    mass_mat = mass/6*np.array([
            [2, 0, 1, 0],
            [0, 2, 0, 1],
            [1, 0, 2, 0],
            [0, 1, 0, 2]])
    return stiff_mat, mass_mat


def uelbeam2DU(coord, params):
    """2D 2-noded beam element without axial deformation

    Parameters
    ----------
    coord : ndarray
        Coordinates for the nodes of the element (2, 2).
    params : tuple
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
    bending_stiff = area_moment*young/L**3
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
    mass_mat = mass/420*np.array([
            [156, 22*L, 54, -13*L],
            [22*L, 4*L**2, 13*L, -3*L**2],
            [54, 13*L, 156, -22*L],
            [-13*L, -3*L**2, -22*L, 4*L**2]])
    stiff_mat = Q.T @ stiff_mat @ Q
    mass_mat = Q.T @ mass_mat @ Q
    return stiff_mat, mass_mat


def uelbeam2D(coord, params):
    """2D 2-noded beam element with axial deformation

    Parameters
    ----------
    coord : ndarray
      Coordinates for the nodes of the element (2, 2).
    params : tuple
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
    ratio = area/area_moment
    stiff_mat = bending_stiff/L**3 * np.array([
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
    mass_mat = mass/420*np.array([
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
