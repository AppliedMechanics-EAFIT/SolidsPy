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


def eletype(eletype):
    """Assigns number to degrees of freedom

    According to iet assigns number of degrees of freedom, number of
    nodes and minimum required number of integration points.

    Parameters
    ----------
    eletype :  int
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
    except:
        raise ValueError("You entered an invalid type of element.")


#%% Shape functions and derivatives

# Triangles
def shape_tri3(r, s):
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
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s).

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (0, 0.5). Thus

    >>> N, _ = shape_tri3(0, 0)
    >>> N_ex = np.array([
    ...    [1, 0, 0, 0, 0, 0],
    ...    [0, 1, 0, 0, 0, 0]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N, _ = shape_tri3(1/2, 1/2)
    >>> N_ex = np.array([
    ...    [0, 0, 1/2, 0, 1/2, 0],
    ...    [0, 0, 0, 1/2, 0, 1/2]])
    >>> np.allclose(N, N_ex)
    True

    """
    N = np.array([1 - r - s, r, s])
    dNdr = np.array([
      [-1, 1, 0],
      [-1, 0, 1]])
    return N, dNdr


def shape_tri6(r, s):
    """
    Shape functions and derivatives for a quadratic element

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.

    Returns
    -------
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s).
    """
    N = np.array(
        [(1 - r - s) - 2*r*(1 - r - s) - 2*s*(1 - r - s),
         r - 2*r*(1 - r - s) - 2*r*s,
         s - 2*r*s - 2*s*(1-r-s),
         4*r*(1 - r - s),
         4*r*s,
         4*s*(1 - r - s)])
    dNdr = np.array([
        [4*r + 4*s - 3, 4*r - 1, 0, -8*r - 4*s + 4, 4*s, -4*s],
        [4*r + 4*s - 3, 0, 4*s - 1, -4*r, 4*r, -4*r - 8*s + 4]])
    return N, dNdr

# Quadrilaterals
def shape_quad4(r, s):
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
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s).
    dNdr : ndarray (float)
        Array with the derivative of the shape functions evaluated at
        the point (r, s).

    Examples
    --------
    We can check evaluating at two different points, namely (0, 0) and
    (1, 1). Thus

    >>> N, _ = shape_quad4(0, 0)
    >>> N_ex = np.array([
    ...    [1/4, 0, 1/4, 0, 1/4, 0, 1/4, 0],
    ...    [0, 1/4, 0, 1/4, 0, 1/4, 0, 1/4]])
    >>> np.allclose(N, N_ex)
    True

    and

    >>> N, _ = shape_quad4(1, 1)
    >>> N_ex = np.array([
    ...    [0, 0, 0, 0, 1, 0, 0, 0],
    ...    [0, 0, 0, 0, 0, 1, 0, 0]])
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


def shape_quad9(r, s):
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
    N : ndarray (float)
        Array with the shape functions evaluated at the point (r, s).
    dNdr : ndarray (float)
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


#%% Derivative matrices
def elast_mat_2d(r, s, coord, element):
    """
    Interpolation matrices for elements for plane elasticity

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.
    coord : ndarray (float)
        Coordinates of the element.

    Returns
    -------
    H : ndarray (float)
        Array with the shape functions evaluated at the point (r, s)
        for each degree of freedom.
    B : ndarray (float)
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


def elast_mat_axi(r, s, coord, element):
    """
    Interpolation matrices for elements for axisymetric elasticity

    Parameters
    ----------
    r : float
        Horizontal coordinate of the evaluation point.
    s : float
        Vertical coordinate of the evaluation point.
    coord : ndarray (float)
        Coordinates of the element.

    Returns
    -------
    H : ndarray (float)
        Array with the shape functions evaluated at the point (r, s)
        for each degree of freedom.
    B : ndarray (float)
        Array with the displacement to strain matrix evaluated
        at the point (r, s).
    det : float
        Determinant of the Jacobian.
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


#%%
def jacoper(dNdr, coord):
    """
    Compute the Jacobian of the transformation evaluated at
    the Gauss point

    Parameters
    ----------
    dNdr : ndarray
      Derivatives of the interpolation function with respect to the
      natural coordinates.
    coord : ndarray
      Coordinates of the nodes of the element (nnodes, ndim).

    Returns
    -------
    jaco_inv : ndarray (ndim, ndim)
      Jacobian of the transformation evaluated at a point.

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


#%% Elemental strains
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
    epsG = np.zeros([3, 3])
    xl = np.zeros([3, 2])
    gpts, _ = gau.gauss_tri(order=1)
    for i in range(gpts.shape[0]):
        ri, si =  gpts[i, :]
        H, B, _ = elast_mat_2d(ri, si, coord, shape_tri3)
        epsG[:, i] = B @ ul
        xl[i, 0] = np.dot(H[0, ::2], coord[:, 0])
        xl[i, 1] = np.dot(H[0, ::2], coord[:, 0])
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
    epsG = np.zeros([3, 7])
    xl = np.zeros([7, 2])
    gpts, _ = gau.gauss_tri(order=3)
    for i in range(7):
        ri, si = gpts[i, :]
        H, B, _ = elast_mat_2d(ri, si, coord, shape_tri6)
        epsG[:, i] = B @ ul
        xl[i, 0] = np.dot(H[0, ::2], coord[:, 0])
        xl[i, 1] = np.dot(H[0, ::2], coord[:, 0])
    return epsG.T, xl


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
    epsG = np.zeros([3, 4])
    xl = np.zeros([4, 2])
    gpts, _ = gau.gauss_nd(2)
    for i in range(gpts.shape[0]):
        ri, si = gpts[i, :]
        H, B, _= elast_mat_2d(ri, si, coord, shape_quad4)
        epsG[:, i] = B @ ul
        xl[i, 0] = np.dot(H[0, ::2], coord[:, 0])
        xl[i, 1] = np.dot(H[0, ::2], coord[:, 0])
    return epsG.T, xl


#%% Material routines
def umat(params):
    """2D Elasticity consitutive matrix in plane stress

    For plane strain use effective properties.

    Parameters
    ----------
    params : tuple
        Material parameters in the following order:

            E : float
                Young modulus (>0).
            nu : float
                Poisson coefficient (-1, 0.5).

    Returns
    -------
    C : ndarray
      Constitutive tensor in Voigt notation.

    Examples
    --------

    >>> params = 8/3, 1/3
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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
