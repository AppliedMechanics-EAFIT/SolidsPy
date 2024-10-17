# -*- coding: utf-8 -*-
"""
FEM routines
------------
    
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
import numpy as np
import solidspy.gaussutil as gau
from typing import Tuple, Callable
from numpy.typing import NDArray


def eletype(eletype: int) -> Tuple[int, int, int]:
    """Assigns number to degrees of freedom

    According to `eletype` assigns number of degrees of freedom, number of
    nodes and minimum required number of integration points.

    Parameters
    ----------
    eletype : int
      Type of element. These are:
        1. 4 node bilinear quadrilateral.
        2. 6 node quadratic triangle.
        3. 3 node linear triangle.
        5. 2 node spring.
        6. 2 node truss element.
        7. 2 node beam (3 DOF per node).
        8. 2 node beam with axial force (3 DOF per node).

    Returns
    -------
    ndof : int
      Number of degrees of freedom for the selected element.
    nnodes : int
      Number of nodes for the selected element.
    ngpts : int
      Number of Gauss points for the selected element.

    Raises
    ------
    ValueError
        If an invalid element type is provided.

    """
    elem_id = {
        1: (8, 4, 4),
        2: (12, 6, 7),
        3: (6, 3, 3),
        4: (18, 9, 9),
        5: (4, 2, 3),
        6: (4, 2, 3),
        7: (6, 2, 3),
        8: (6, 2, 3)}
    try:
        return elem_id[eletype]
    except KeyError:
        raise ValueError("You entered an invalid type of element.")


#%% Shape functions and derivatives

# Triangles
def shape_tri3(r: float, s: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Shape functions and derivatives for a linear element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s).
    dNdr : NDArray[np.float64]
        Array with the derivative of the shape functions evaluated at
        the point (r, s).

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (0, 0.5). Thus

    >>> N, _ = shape_tri3(0, 0)
    >>> N_ex = np.array([
    ...    [1, 0, 0],
    ...    [0, 1, 0],
    ...    [0, 0, 1]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N, _ = shape_tri3(0.5, 0.5)
    >>> N_ex = np.array([
    ...    [0, 0.5, 0.5],
    ...    [0, 1, 0],
    ...    [0, 0, 1]])
    >>> np.allclose(N, N_ex)
    True

    """
    N = np.array([1 - r - s, r, s])
    dNdr = np.array([
      [-1, 1, 0],
      [-1, 0, 1]])
    return N, dNdr


def shape_tri6(r: float, s: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Shape functions and derivatives for a quadratic triangular element.

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : NDArray[np.float64]
        Shape functions evaluated at the point (r, s).
    dNdrds : NDArray[np.float64]
        Derivatives of the shape functions with respect to r and s.

    """
    # Convert inputs to NumPy arrays
    r = np.asarray(r)
    s = np.asarray(s)
    
    # Check if inputs are scalar
    scalar_input = False
    if r.ndim == 0 and s.ndim == 0:
        r = r.reshape(1)
        s = s.reshape(1)
        scalar_input = True

    # Compute shape functions N1 to N6
    N = np.vstack((
        (1 - r - s) - 2 * r * (1 - r - s) - 2 * s * (1 - r - s),
        r - 2 * r * (1 - r - s) - 2 * r * s,
        s - 2 * r * s - 2 * s * (1 - r - s),
        4 * r * (1 - r - s),
        4 * r * s,
        4 * s * (1 - r - s)
    )).T  # Shape: (n_points, 6)

    # Compute derivatives with respect to r
    dNdr = np.vstack((
        4 * r + 4 * s - 3,
        4 * r - 1,
        np.zeros_like(r),
        -8 * r - 4 * s + 4,
        4 * s,
        -4 * s
    )).T  # Shape: (n_points, 6)

    # Compute derivatives with respect to s
    dNds = np.vstack((
        4 * r + 4 * s - 3,
        np.zeros_like(r),
        4 * s - 1,
        -4 * r,
        4 * r,
        -4 * r - 8 * s + 4
    )).T  # Shape: (n_points, 6)

    # Stack derivatives into a single array
    dNdrds = np.stack((dNdr, dNds), axis=1)  # Shape: (n_points, 2, 6)

    if scalar_input:
        return N[0], dNdrds[0]
    else:
        return N, dNdrds

# Quadrilaterals
def shape_quad4(r: float, s: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Shape functions and derivatives for a bilinear element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s).
    dNdr : NDArray[np.float64]
        Array with the derivative of the shape functions evaluated at
        the point (r, s).

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (1, 1). Thus

    >>> N, _ = shape_quad4(0, 0)
    >>> N_ex = np.array([
    ...    [0.25, 0.0, 0.25, 0.0, 0.25, 0.0, 0.25, 0.0],
    ...    [0.0, 0.25, 0.0, 0.25, 0.0, 0.25, 0.0, 0.25]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N, _ = shape_quad4(1, 1)
    >>> N_ex = np.array([
    ...    [0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0],
    ...    [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0]])
    >>> np.allclose(N, N_ex)
    True

    """
    N = 0.25*np.array(
        [(1 - r)*(1 - s),
         (1 + r)*(1 - s),
         (1 + r)*(1 + s),
         (1 - r)*(1 + s)])
    dNdr = 0.25*np.array([
        [s - 1, -s + 1, s + 1, -s - 1],
        [r - 1, -r - 1, r + 1, -r + 1]])
    return N, dNdr


def shape_quad9(r: float, s: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Shape functions and derivatives for a biquadratic element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s).
    dNdr : NDArray[np.float64]
        Array with the derivative of the shape functions evaluated at
        the point (r, s).
    """
    N = np.array(
        [0.25*r*s*(r - 1.0)*(s - 1.0),
         0.25*r*s*(r + 1.0)*(s - 1.0),
         0.25*r*s*(r + 1.0)*(s + 1.0),
         0.25*r*s*(r - 1.0)*(s + 1.0),
         0.5*s*(-r**2 + 1.0)*(s - 1.0),
         0.5*r*(r + 1.0)*(-s**2 + 1.0),
         0.5*s*(-r**2 + 1.0)*(s + 1.0),
         0.5*r*(r - 1.0)*(-s**2 + 1.0),
         (-r**2 + 1.0)*(-s**2 + 1.0)])
    dNdr = np.array([
        [0.25*s*(2.0*r - 1.0)*(s - 1.0),
         0.25*s*(2.0*r + 1.0)*(s - 1.0),
         0.25*s*(2.0*r + 1.0)*(s + 1.0),
         0.25*s*(2.0*r - 1.0)*(s + 1.0),
         r*s*(-s + 1.0),
         -0.5*(2.0*r + 1.0)*(s**2 - 1.0),
         -r*s*(s + 1.0),
         0.5*(-2.0*r + 1.0)*(s**2 - 1.0),
         2.0*r*(s**2 - 1.0)],
        [0.25*r*(r - 1.0)*(2.0*s - 1.0),
         0.25*r*(r + 1.0)*(2.0*s - 1.0),
         0.25*r*(r + 1.0)*(2.0*s + 1.0),
         0.25*r*(r - 1.0)*(2.0*s + 1.0),
         0.5*(r**2 - 1.0)*(-2.0*s + 1.0),
         -r*s*(r + 1.0),
         -0.5*(r**2 - 1.0)*(2.0*s + 1.0),
         r*s*(-r + 1.0),
         2.0*s*(r**2 - 1.0)]])
    return N, dNdr


def shape_quad8(r: float, s: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Shape functions and derivatives for an 8-noded serendipity element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s).
    dNdr : NDArray[np.float64]
        Array with the derivative of the shape functions evaluated at
        the point (r, s).
    """

    N = 0.25*np.array(
        [(1.0 - r)*(1.0 - s) - (1.0 - r)*(1.0 - s**2) - (1.0 - s)*(1.0 - r**2),
         (1.0 + r)*(1.0 - s) - (1.0 + r)*(1.0 - s**2) - (1.0 - s)*(1.0 - r**2),
         (1.0 + r)*(1.0 + s) - (1.0 + r)*(1.0 - s**2) - (1.0 + s)*(1.0 - r**2),
         (1.0 - r)*(1.0 + s) - (1.0 - r)*(1.0 - s**2) - (1.0 + s)*(1.0 - r**2),
         2.0*(1.0 - s)*(1.0 - r**2),
         2.0*(1.0 + r)*(1.0 - s**2),
         2.0*(1.0 + s)*(1.0 - r**2),
         2.0*(1.0 - r)*(1.0 - s**2)])
    dNdr = 0.25*np.array([
        [-2.0*r*(s - 1.0) - s**2 + s,
         -2.0*r*(s - 1.0) + s**2 - s,
         -2.0*r*(-s - 1.0) + s**2 + s,
         -2.0*r*(-s - 1.0) - s**2 - s,
         -2.0*r*(2.0 - 2.0*s),
         2.0 - 2.0*s**2,
         -2.0*r*(2.0*s + 2.0),
         2.0*s**2 - 2.0],
        [-r**2 + r - 2.0*s*(r - 1.0),
         -r**2 - r - 2.0*s*(-r - 1.0),
         r**2 + r - 2.0*s*(-r - 1.0),
         r**2 - r - 2.0*s*(r - 1.0),
         2.0*r**2 - 2.0,
         -2.0*s*(2.0*r + 2.0),
         2.0 - 2.0*r**2,
         -2.0*s*(2.0 - 2.0*r)]])
    return N, dNdr



# 3D elements
def shape_tet4(r: float, s: float, t: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Shape functions and derivatives for a linear tetrahedron

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Horizontal coordinate of the evaluation point.
    t : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s, t).
    dNdr : NDArray[np.float64]
        Array with the derivative of the shape functions evaluated at
        the point (r, s, t).
    """
    N = np.array([1 - r - s - t, r, s, t])
    dNdr = np.array([
        [-1, 1, 0, 0],
        [-1, 0, 1, 0],
        [-1, 0, 0, 1]])
    return N, dNdr


def shape_hex8(r: float, s: float, t: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Shape functions and derivatives for a trilinear element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Horizontal coordinate of the evaluation point.
    t : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s, t).
    dNdr : NDArray[np.float64]
        Array with the derivative of the shape functions evaluated at
        the point (r, s, t).
    """
    N = np.array([
        (1 - r)*(1 - s)*(1 - t), (1 - s)*(1 - t)*(r + 1),
        (1 - t)*(r + 1)*(s + 1), (1 - r)*(1 - t)*(s + 1),
        (1 - r)*(1 - s)*(t + 1), (1 - s)*(r + 1)*(t + 1),
        (r + 1)*(s + 1)*(t + 1), (1 - r)*(s + 1)*(t + 1)])
    dNdr = np.array([
        [(1 - t)*(s - 1), (1 - s)*(1 - t),
         (1 - t)*(s + 1), (1 - t)*(-s - 1),
         (1 - s)*(-t - 1), (1 - s)*(t + 1),
         (s + 1)*(t + 1), -(s + 1)*(t + 1)],
        [(1 - t)*(r - 1), (1 - t)*(-r - 1),
         (1 - t)*(r + 1), (1 - r)*(1 - t),
         -(1 - r)*(t + 1), -(r + 1)*(t + 1),
         (r + 1)*(t + 1), (1 - r)*(t + 1)],
        [-(1 - r)*(1 - s), -(1 - s)*(r + 1),
         -(r + 1)*(s + 1), -(1 - r)*(s + 1),
         (1 - r)*(1 - s), (1 - s)*(r + 1),
         (r + 1)*(s + 1), (1 - r)*(s + 1)]])
    return 0.125*N, 0.125*dNdr


#%% Derivative matrices
def elast_diff_2d(
    r: float,
    s: float,
    coord: NDArray[np.float64],
    element: Callable[[float, float], Tuple[NDArray[np.float64], NDArray[np.float64]]]
) -> Tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """
    Interpolation matrices for elements for plane elasticity

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.
    coord : NDArray[np.float64]
        Coordinates of the element.
    element : Callable
        Shape function generator for the element.

    Returns
    -------
    H : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s)
        for each degree of freedom.
    B : NDArray[np.float64]
        Array with the displacement to strain matrix evaluated
        at the point (r, s).
    det : float
        Determinant of the Jacobian.

    """
    N, dNdr = element(r, s)
    det, jaco_inv = jacoper(dNdr, coord)
    dNdx = jaco_inv @ dNdr
    H = np.zeros((2, 2*N.shape[0]))
    B = np.zeros((3, 2*N.shape[0]))
    H[0, 0::2] = N
    H[1, 1::2] = N
    B[0, 0::2] = dNdx[0, :]
    B[1, 1::2] = dNdx[1, :]
    B[2, 0::2] = dNdx[1, :]
    B[2, 1::2] = dNdx[0, :]
    return H, B, det


def elast_diff_axi(
    r: float,
    s: float,
    coord: NDArray[np.float64],
    element: Callable[[float, float], Tuple[NDArray[np.float64], NDArray[np.float64]]]
) -> Tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """
    Interpolation matrices for elements for axisymmetric elasticity

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.
    coord : NDArray[np.float64]
        Coordinates of the element.
    element : Callable
        Shape function generator for the element.

    Returns
    -------
    H : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s)
        for each degree of freedom.
    B : NDArray[np.float64]
        Array with the displacement to strain matrix evaluated
        at the point (r, s).
    det : float
        Determinant of the Jacobian.

    Raises
    ------
    ValueError
        If the radial coordinate `x` is negative.

    """
    N, dNdr = element(r, s)
    x = N.dot(coord[:, 0])
    if x < 0:
        raise ValueError("Horizontal coordinates should be non-negative.")
    det, jaco_inv = jacoper(dNdr, coord)
    dNdx = jaco_inv @ dNdr
    H = np.zeros((2, 2*N.shape[0]))
    B = np.zeros((4, 2*N.shape[0]))
    H[0, 0::2] = N
    H[1, 1::2] = N
    B[0, 0::2] = dNdx[0, :]
    B[1, 1::2] = dNdx[1, :]
    B[2, 0::2] = dNdx[1, :]
    B[2, 1::2] = dNdx[0, :]
    B[2, 1::2] = dNdx[0, :]
    B[3, 0::2] = N/x
    return H, B, det


def elast_diff_3d(
    r: float,
    s: float,
    t: float,
    coord: NDArray[np.float64],
    interp: Callable[[float, float, float], Tuple[NDArray[np.float64], NDArray[np.float64]]] = shape_hex8
) -> Tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """
    Interpolation matrices for a trilinear element for elasticity

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Horizontal coordinate of the evaluation point.
    t : float
        Vertical coordinate of the evaluation point.
    coord : NDArray[np.float64]
        Coordinates of the element.
    interp : Callable, optional
        Shape function generator for the element. Defaults to `shape_hex8`.

    Returns
    -------
    H : NDArray[np.float64]
        Array with the shape functions evaluated at the point (r, s, t)
        for each degree of freedom.
    B : NDArray[np.float64]
        Array with the displacement to strain matrix evaluated
        at the point (r, s, t).
    det : float
        Determinant of the Jacobian.

    """
    N, dNdr = interp(r, s, t)
    det, jaco_inv = jacoper(dNdr, coord)
    dNdx = jaco_inv @ dNdr
    H = np.zeros((3, 3*N.shape[0]))
    B = np.zeros((6, 3*N.shape[0]))
    H[0, 0::3] = N
    H[1, 1::3] = N
    H[2, 2::3] = N
    B[0, 0::3] = dNdx[0, :]
    B[1, 1::3] = dNdx[1, :]
    B[2, 2::3] = dNdx[2, :]
    B[3, 1::3] = dNdx[2, :]
    B[3, 2::3] = dNdx[1, :]
    B[4, 0::3] = dNdx[2, :]
    B[4, 2::3] = dNdx[0, :]
    B[5, 0::3] = dNdx[1, :]
    B[5, 1::3] = dNdx[0, :]
    return H, B, det


#%%
def jacoper(
    dNdr: NDArray[np.float64],
    coord: NDArray[np.float64]
) -> Tuple[float, NDArray[np.float64]]:
    """
    Compute the Jacobian of the transformation evaluated at
    the Gauss point

    Parameters
    ----------
    dNdr : NDArray[np.float64]
      Derivatives of the interpolation function with respect to the
      natural coordinates.
    coord : NDArray[np.float64]
      Coordinates of the nodes of the element (nnodes, ndim).

    Returns
    -------
    det : float
      Determinant of the Jacobian.
    jaco_inv : NDArray[np.float64]
      Inverse of the Jacobian matrix.

    Raises
    ------
    ValueError
        If the Jacobian determinant is close to zero or negative.

    """
    jaco = dNdr @ coord
    det = np.linalg.det(jaco)
    if np.isclose(np.abs(det), 0.0):
        msg = "Jacobian close to zero. Check the shape of your elements!"
        raise ValueError(msg)
    jaco_inv = np.linalg.inv(jaco)
    if det < 0.0:
        msg = "Jacobian is negative. Check your elements orientation!"
        raise ValueError(msg)
    return det, jaco_inv


#%% Material routines
def umat(params: Tuple[float, float]) -> NDArray[np.float64]:
    """2D Elasticity constitutive matrix in plane stress

    For plane strain use effective properties.

    Parameters
    ----------
    params : Tuple[float, float]
        Material parameters in the following order:

            E : float
                Young modulus (>0).
            nu : float
                Poisson coefficient (-1, 0.5).

    Returns
    -------
    C : NDArray[np.float64]
      Constitutive tensor in Voigt notation.

    Examples
    --------

    >>> params = (8/3, 1/3)
    >>> C = umat(params)
    >>> C_ex = np.array([
    ...    [3, 1, 0],
    ...    [1, 3, 0],
    ...    [0, 0, 1]])
    >>> np.allclose(C, C_ex)
    True

    """
    E, nu = params
    C = np.zeros((3, 3))
    enu = E/(1 - nu**2)
    mnu = (1 - nu)/2
    C[0, 0] = enu
    C[0, 1] = nu*enu
    C[1, 0] = C[0, 1]
    C[1, 1] = enu
    C[2, 2] = enu*mnu

    return C


def elast_mat_2d(params: Tuple[float, float]) -> NDArray[np.float64]:
    """2D Elasticity constitutive matrix in plane stress

    For plane strain use effective properties.

    Parameters
    ----------
    params : Tuple[float, float]
        Material parameters in the following order:

            E : float
                Young modulus (>0).
            nu : float
                Poisson coefficient (-1, 0.5).

    Returns
    -------
    C : NDArray[np.float64]
      Constitutive tensor in Voigt notation.

    """
    E, nu = params
    lamda = E*nu/(1 + nu)/(1 - 2*nu)
    mu = 0.5*E/(1 + nu)
    C = np.array([
    [2*mu + lamda, lamda, 0, 0],
    [lamda, 2*mu + lamda, 0, 0],
    [0, 0, mu, 0],
    [0, 0, 0, 2*mu + lamda]])
    return C


def elast_mat_axi(params: Tuple[float, float]) -> NDArray[np.float64]:
    """Elasticity constitutive matrix for axisymmetric problems.

    Parameters
    ----------
    params : Tuple[float, float]
        Material parameters in the following order:

            E : float
                Young modulus (>0).
            nu : float
                Poisson coefficient (-1, 0.5).

    Returns
    -------
    C : NDArray[np.float64]
      Constitutive tensor in Voigt notation.

    """
    E, nu = params
    lamda = E*nu/(1 + nu)/(1 - 2*nu)
    mu = 0.5*E/(1 + nu)
    C = np.array([
        [2*mu + lamda, lamda, 0, 0],
        [lamda, 2*mu + lamda, 0, 0],
        [0, 0, mu, 0],
        [0, 0, 0, 2*mu + lamda]])
    return C


def elast_mat(params: Tuple[float, float]) -> NDArray[np.float64]:
    """3D Elasticity constitutive matrix.

    Parameters
    ----------
    params : Tuple[float, float]
        Material parameters in the following order:

            E : float
                Young modulus (>0).
            nu : float
                Poisson coefficient (-1, 0.5).

    Returns
    -------
    C : NDArray[np.float64]
        Constitutive tensor in Voigt notation.

    """
    E, nu = params
    lamda = E*nu/(1 + nu)/(1 - 2*nu)
    mu = 0.5*E/(1 + nu)
    C = np.array([
    [2*mu + lamda, lamda, lamda, 0, 0, 0],
    [lamda, 2*mu + lamda, lamda, 0, 0, 0],
    [lamda, lamda, 2*mu + lamda, 0, 0, 0],
    [0, 0, 0, mu, 0, 0],
    [0, 0, 0, 0, mu, 0],
    [0, 0, 0, 0, 0, mu]])
    return C


#%% Elemental strains
def str_el3(
    coord: NDArray[np.float64],
    ul: NDArray[np.float64]
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute the strains at each element integration point

    This one is used for 3-noded triangular elements.

    Parameters
    ----------
    coord : NDArray[np.float64]
      Coordinates of the nodes of the element (nn, 2).
    ul : NDArray[np.float64]
      Array with displacements for the element.

    Returns
    -------
    epsGT : NDArray[np.float64]
      Strain components for the Gauss points.
    xl : NDArray[np.float64]
      Configuration of the Gauss points after deformation.

    """
    epsG = np.zeros([3, 3])
    xl = np.zeros([3, 2])
    gpts, _ = gau.gauss_tri(order=1)
    for i in range(gpts.shape[0]):
        ri, si =  gpts[i, :]
        H, B, _ = elast_diff_2d(ri, si, coord, shape_tri3)
        epsG[:, i] = B @ ul
        xl[i, 0] = np.dot(H[0, ::2], coord[:, 0])
        xl[i, 1] = np.dot(H[0, ::2], coord[:, 0])
    return epsG.T, xl


def str_el6(
    coord: NDArray[np.float64],
    ul: NDArray[np.float64]
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute the strains at each element integration point

    This one is used for 6-noded triangular elements.

    Parameters
    ----------
    coord : NDArray[np.float64]
      Coordinates of the nodes of the element (6, 2).
    ul : NDArray[np.float64]
      Array with displacements for the element.

    Returns
    -------
    epsGT : NDArray[np.float64]
      Strain components for the Gauss points.
    xl : NDArray[np.float64]
      Configuration of the Gauss points after deformation.

    """
    epsG = np.zeros([3, 7])
    xl = np.zeros([7, 2])
    gpts, _ = gau.gauss_tri(order=3)
    for i in range(7):
        ri, si = gpts[i, :]
        H, B, _ = elast_diff_2d(ri, si, coord, shape_tri6)
        epsG[:, i] = B @ ul
        xl[i, 0] = np.dot(H[0, ::2], coord[:, 0])
        xl[i, 1] = np.dot(H[0, ::2], coord[:, 0])
    return epsG.T, xl


def str_el4(
    coord: NDArray[np.float64],
    ul: NDArray[np.float64]
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute the strains at each element integration point

    This one is used for 4-noded quadrilateral elements.

    Parameters
    ----------
    coord : NDArray[np.float64]
      Coordinates of the nodes of the element (4, 2).
    ul : NDArray[np.float64]
      Array with displacements for the element.

    Returns
    -------
    epsGT : NDArray[np.float64]
      Strain components for the Gauss points.
    xl : NDArray[np.float64]
      Configuration of the Gauss points after deformation.

    """
    epsG = np.zeros([3, 4])
    xl = np.zeros([4, 2])
    gpts, _ = gau.gauss_nd(2)
    for i in range(gpts.shape[0]):
        ri, si = gpts[i, :]
        H, B, _= elast_diff_2d(ri, si, coord, shape_quad4)
        epsG[:, i] = B @ ul
        xl[i, 0] = np.dot(H[0, ::2], coord[:, 0])
        xl[i, 1] = np.dot(H[0, ::2], coord[:, 0])
    return epsG.T, xl


if __name__ == "__main__":
    import doctest
    doctest.testmod()
