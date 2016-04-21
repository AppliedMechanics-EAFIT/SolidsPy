# -*- coding: utf-8 -*-
"""
Postprocessor subroutines
-------------------------

@author: eafit
"""
from __future__ import division
import numpy as np
import sympy as sym
import femutil as fe
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy.interpolate import griddata, interp2d
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['image.cmap'] = "YlGnBu_r"


def mesh2tri(nodes, elements):
    """Generate a  matplotlib.tri.Triangulation object from the mesh
    
    Parameters
    ----------
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.
    
    Returns
    -------
    tri : Triangulation
        An unstructured triangular grid consisting of npoints points
        and ntri triangles.
    
    """
    x = nodes[:, 1]
    y = nodes[:, 2]
    triangs = []
    for el in elements:
        if el[1]==1:
            triangs.append(el[[3, 4, 5]])
            triangs.append(el[[5, 6, 3]])
        if el[1]==2:
            triangs.append(el[[3, 6, 8]])
            triangs.append(el[[6, 7, 8]])
            triangs.append(el[[6, 4, 7]])
            triangs.append(el[[7, 5, 8]])
        if el[1]==3:
            triangs.append(el[3:])
    
    tri = Triangulation(x, y, np.array(triangs))
    return tri    


def tri_plot(tri, field, title="", figtitle="", levels=12, savefigs=False,
             plt_type="contourf", filename="solution_plot.pdf"):
    
    if plt_type=="pcolor":
        disp_plot = plt.tripcolor
    elif plt_type=="contourf":
        disp_plot = plt.tricontourf

    plt.figure(figtitle)
    disp_plot(tri, field, levels, shading="gouraud")
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.axis("image")
    plt.grid()
    if savefigs:
        plt.savefig(filename)


def plot_disp(UC, nodes, elements, plt_type="contourf", levels=12,
               savefigs=False):
    """Plot the nodal displacement using a triangulation

    Parameters
    ----------
    UC : ndarray (float)
      Array with the displacements.
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    """
    tri = mesh2tri(nodes, elements)
    tri_plot(tri, UC[:, 0], title=r'$u_x$',
             figtitle="Solution: Horizontal displacement",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="ux_sol.pdf")
    tri_plot(tri, UC[:, 1], title=r'$u_y$',
             figtitle="Solution: Vertical displacement",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="uy_sol.pdf")


def plot_strain(E_nodes, nodes, elements, plt_type="contourf", levels=12,
               savefigs=False):
    """Plot the nodal strains using a triangulation
    
    The strains need to be computed at nodes first.

    Parameters
    ----------
    E_nodes : ndarray (float)
      Array with the nodal strains.
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    """
    tri = mesh2tri(nodes, elements)
    tri_plot(tri, E_nodes[:, 0], title=r'$\epsilon_{xx}$',
             figtitle="Solution: epsilon-xx strain",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="epsxx_sol.pdf")
    tri_plot(tri, E_nodes[:, 1], title=r'$\epsilon_{yy}$',
             figtitle="Solution: epsilon-yy strain",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="epsyy_sol.pdf")
    tri_plot(tri, E_nodes[:, 2], title=r'$\gamma_{xy}$',
             figtitle="Solution: gamma-xy strain",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="gammaxy_sol.pdf")


def plot_stress(S_nodes, nodes, elements, plt_type="contourf", levels=12,
               savefigs=False):
    """Plot the nodal stresses using a triangulation
    
    The stresses need to be computed at nodes first.

    Parameters
    ----------
    S_nodes : ndarray (float)
      Array with the nodal stresses.
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    """
    tri = mesh2tri(nodes, elements)
    tri_plot(tri, S_nodes[:, 0], title=r'$\sigma_{xx}$',
             figtitle="Solution: sigma-xx stress",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="sigmaxx_sol.pdf")
    tri_plot(tri, S_nodes[:, 1], title=r'$\sigma_{yy}$',
             figtitle="Solution: sigma-yy stres",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="sigmayy_sol.pdf")
    tri_plot(tri, S_nodes[:, 2], title=r'$\sigma_{xy}$',
             figtitle="Solution: sigma-xy stress",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="sigmaxy_sol.pdf")


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
    """
    Scatter the nodal displacements vector `UG` over the Gauss
    points of each element.

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
      Array with the displacements for each Gauss point in the mesh.

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


def complete_disp(IBC, nodes, UG):
    """
    Fill the displacement vectors with imposed and computed values.
    
    IBC : ndarray (int)
      IBC (Indicator of Boundary Conditions) indicates if the nodes
      has any type of boundary conditions applied to it.
    UG : ndarray (float)
      Array with the computed displacements.
    nodes : ndarray (float)
      Array with number and nodes coordinates:
      
    Returns
    -------
    UC : ndarray (float)
      Array with the displacements.

    """
    nn = nodes.shape[0]
    UC = np.zeros([nn, 2], dtype=np.float)
    for i in range(nn):
        for j in range(2):
            kk = IBC[i, j]
            if kk == -1:
                UC[i, j] = 0.0
            else:
                UC[i, j] = UG[kk]

    return UC


def strain_nodes(IELCON, UU, ne, COORD, elements, mats):
    """Compute averaged strains and stresses at nodes
    
    First, the variable is extrapolated from the Gauss
    point to nodes for each element. Then, these values are averaged
    according to the number of element that share that node. The theory
    for this technique can be found in [1]_.
    
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
    E_nodes : ndarray
        Strains evaluated at the nodes.
      
    
    References
    ----------
    .. [1] O.C. Zienkiewicz and J.Z. Zhu, The Superconvergent patch
        recovery and a posteriori error estimators. Part 1. The
        recovery technique, Int. J. Numer. Methods Eng., 33,
        1331-1364 (1992).
    """
    iet = elements[0, 1]
    ndof, nnodes, ngpts = fe.eletype(iet)

    elcoor = np.zeros([nnodes, 2])
    E_nodes = np.zeros([COORD.shape[0], 3])
    S_nodes = np.zeros([COORD.shape[0], 3])
    el_nodes = np.zeros([COORD.shape[0]], dtype=int)
    ul = np.zeros([ndof])
    for i in range(ne):
        young, poisson = mats[elements[i, 2], :]
        shear = young/(2*(1 + poisson))
        fact1 = young/(1 - poisson**2)
        fact2 = poisson*young/(1 - poisson**2)
        for j in range(nnodes):
            elcoor[j, 0] = COORD[IELCON[i, j], 0]
            elcoor[j, 1] = COORD[IELCON[i, j], 1]
        for j in range(ndof):
            ul[j] = UU[i, j]
        if iet == 1:
            epsG, xl = fe.str_el4(elcoor, ul)
            extrap0 = interp2d(xl[:, 0], xl[:,1], epsG[:, 0])
            extrap1 = interp2d(xl[:, 0], xl[:,1], epsG[:, 1])
            extrap2 = interp2d(xl[:, 0], xl[:,1], epsG[:, 2])
        elif iet == 2:
            epsG, xl = fe.str_el6(elcoor, ul)
            extrap0 = interp2d(xl[:, 0], xl[:,1], epsG[:, 0])
            extrap1 = interp2d(xl[:, 0], xl[:,1], epsG[:, 1])
            extrap2 = interp2d(xl[:, 0], xl[:,1], epsG[:, 2])
        elif iet == 3:
            epsG, xl = fe.str_el3(elcoor, ul)
            extrap0 = lambda x, y: epsG[0, 0]
            extrap1 = lambda x, y: epsG[0, 1]
            extrap2 = lambda x, y: epsG[0, 2]

        for node in IELCON[i, :]:
            x, y = COORD[node, :]
            E_nodes[node, 0] = E_nodes[node, 0] + extrap0(x, y)
            E_nodes[node, 1] = E_nodes[node, 1]  + extrap1(x, y)
            E_nodes[node, 2] = E_nodes[node, 2] + extrap2(x, y)
            S_nodes[node, 0] = S_nodes[node, 0] + fact1*extrap0(x, y) \
                        + fact2*extrap1(x, y)
            S_nodes[node, 1] = S_nodes[node, 1] + fact2*extrap0(x, y) \
                        + fact1*extrap1(x, y)
            S_nodes[node, 2] = S_nodes[node, 2] + shear*extrap2(x, y)
            el_nodes[node] = el_nodes[node] + 1

    E_nodes[:, 0] = E_nodes[:, 0]/el_nodes
    E_nodes[:, 1] = E_nodes[:, 1]/el_nodes
    E_nodes[:, 2] = E_nodes[:, 2]/el_nodes
    S_nodes[:, 0] = S_nodes[:, 0]/el_nodes
    S_nodes[:, 1] = S_nodes[:, 1]/el_nodes
    S_nodes[:, 2] = S_nodes[:, 2]/el_nodes
    return E_nodes, S_nodes


def eigvals(A, tol=1e-6):
    """Eigenvalues and eigenvectors for a 2x2 symmetric matrix/tensor
    
    Parameters
    ----------
    A : ndarray
        Symmetric matrix.
    tol : float (optional)
        Tolerance for considering a matrix diagonal.

    Returns
    -------
    eig1 : float
        First eigenvalue.
    eig2 : float
        Second eigenvalue.
    vec1 : ndarray
        First eigenvector.
    vec2 : ndarray
        Second eigenvector
    
    Examples
    --------
    
    >>> A = np.array([[5, 6],
    ...              [6, 9]])
    >>> eig1, eig2, vec1, vec2 =  eigvals(A)
    >>> np.allclose(eig1, 7 + 2*np.sqrt(10))
    True
    >>> np.allclose(eig2, 7 - 2*np.sqrt(10))
    True
    >>> np.allclose(vec1, np.array([-0.584710284663765, -0.8112421851755609]))
    True
    >>> np.allclose(vec2, np.array([-0.8112421851755609,0.584710284663765]))
    True
    
    """
    if np.abs(A).max() < tol:
        eig1 = 0.0
        eig2 = 0.0
        vec1 = np.array([np.NaN, np.NaN])
        vec2 = np.array([np.NaN, np.NaN])
    elif abs(A[0, 1])/np.abs(A).max() < tol:
        eig1 = A[0, 0]
        eig2 = A[1, 1]
        vec1 = np.array([1, 0])
        vec2 = np.array([0, 1])
    else:
        tr = A[0, 0] + A[1, 1]
        det = A[0, 0]*A[1, 1] - A[0, 1]**2
        eig1 = 0.5*(tr - np.sqrt(tr**2 - 4*det))
        eig2 = 0.5*(tr + np.sqrt(tr**2 - 4*det))
        vec1 = np.array([A[0, 0] - eig2, A[0, 1]])
        vec1 = vec1/np.sqrt(vec1[0]**2 + vec1[1]**2)
        vec2 = np.array([-vec1[1], vec1[0]])
    if abs(eig2) > abs(eig1):
        eig2, eig1 = eig1, eig2
        vec2, vec1 = vec1, vec2

    return eig1, eig2, vec1, vec2


def principal_dirs(field):
    """Compute the principal directions of a tensor field

    Parameters
    ----------
    field : ndarray
        Tensor field. The tensor is written as "vector" using
        Voigt notation.

    Returns
    -------
    eigs1 : ndarray
        Array with the first eigenvalues.
    eigs2 : ndarray
        Array with the second eigenvalues.
    vecs1 : ndarray
        Array with the first eigenvectors.
    vecs2 : ndarray
        Array with the Second eigenvector.

    """
    num = field.shape[0]
    eigs1 = np.empty((num))
    eigs2 = np.empty((num))
    vecs1 = np.empty((num, 2))
    vecs2 = np.empty((num, 2))
    A = np.zeros((2, 2))
    for cont, tensor in enumerate(field):
        A[0, 0] = tensor[0]
        A[1, 1] = tensor[1]
        A[0, 1] = tensor[2]
        eig1, eig2, vec1, vec2 = eigvals(A, tol=1e-6)
        eigs1[cont] = eig1
        eigs2[cont] = eig2
        vecs1[cont, :] = vec1
        vecs2[cont, :] = vec2

    return eigs1, eigs2, vecs1, vecs2
        

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


def gmeshpost(IBC, nn, UG, folder=""):
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
    nomfile1 = folder + 'out.txt'
    np.savetxt(nomfile1, UR, fmt='%.18e', delimiter=' ')


def energy(UG, KG):
    r"""
    Computes the potential energy for the current sln.

    Parameters
    ----------
    UG : ndarray (float)
      Array with the computed displacements.
    KG : ndarray (float)
      Global stiffness matrix.

    Returns
    -------
    EFE : scalar (float)
      Total energy in the system. :math:`-\frac{1}{2} U^T K U`

    """
    f = np.dot(UG.T, KG)
    EFE = -0.5*np.dot(f, UG)

    return EFE


#%%
if __name__ == "__main__":
    import doctest
    doctest.testmod()