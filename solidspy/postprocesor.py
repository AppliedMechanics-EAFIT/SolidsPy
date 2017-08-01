# -*- coding: utf-8 -*-
"""
Postprocessor subroutines

"""
from __future__ import division
import numpy as np
import femutil as fe
import uelutil as uel
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['image.cmap'] = "YlGnBu_r"
rcParams['axes.axisbelow'] = True
rcParams['mathtext.fontset'] = "cm"


#%% Plotting routines
def fields_plot(elements, nodes, UC, E_nodes=None, S_nodes=None):
    """Plot contours for displacements, strains and stresses

    Parameters
    ----------
    nodes : ndarray (float)
        Array with number and nodes coordinates:
         `number coordX coordY BCX BCY`
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each element.
    UC : ndarray (float)
        Array with the displacements.
    E_nodes : ndarray (float)
        Array with strain field in the nodes.
    S_nodes : ndarray (float)
        Array with stress field in the nodes.

    """
    # Check for structural elements in the mesh
    struct_pos = 5 in elements[:, 1] or \
             6 in elements[:, 1] or \
             7 in elements[:, 1]
    if struct_pos:
        # Still not implemented visualization for structural elements
        print(UC)
    else:
        plot_node_field(UC, nodes, elements, title=[r"$u_x$", r"$u_y$"],
                        figtitle=["Horizontal displacement",
                                  "Vertical displacement"])
        if E_nodes is not None:
            plot_node_field(E_nodes, nodes, elements,
                    title=[r"$\epsilon_{xx}$",
                           r"$\epsilon_{yy}$",
                           r"$\gamma_{xy}$",],
                    figtitle=["Strain epsilon-xx",
                              "Strain epsilon-yy",
                              "Strain gamma-xy"])
        if S_nodes is not None:
            plot_node_field(S_nodes, nodes, elements,
                    title=[r"$\sigma_{xx}$",
                           r"$\sigma_{yy}$",
                           r"$\tau_{xy}$",],
                    figtitle=["Stress sigma-xx",
                              "Stress sigma-yy",
                              "Stress tau-xy"])


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
    coord_x = nodes[:, 1]
    coord_y = nodes[:, 2]
    triangs = []
    for el in elements:
        if el[1] == 1:
            triangs.append(el[[3, 4, 5]])
            triangs.append(el[[5, 6, 3]])
        if el[1] == 2:
            triangs.append(el[[3, 6, 8]])
            triangs.append(el[[6, 7, 8]])
            triangs.append(el[[6, 4, 7]])
            triangs.append(el[[7, 5, 8]])
        if el[1] == 3:
            triangs.append(el[3:])

    tri = Triangulation(coord_x, coord_y, np.array(triangs))
    return tri


def tri_plot(tri, field, title="", figtitle="", levels=12, savefigs=False,
             plt_type="contourf", filename="solution_plot.pdf"):
    """Plot contours over triangulation

    Parameters
    ----------
    tri : ndarray (float)
        Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    field : ndarray (float)
        Array with data to be plotted for each node.
    title : string (optional)
        Title of the plot.
    figtitle : string (optional)
        Title for the Figure.
    levels : int (optional)
        Number of levels to be used in ``contourf``.
    savefigs : bool (optional)
        Allow to save the figure.
    plt_type : string (optional)
        Plot the field as one of the options: ``pcolor`` or
        ``contourf``
    filename : string (optional)
        Filename to save the figures.
    """
    if plt_type == "pcolor":
        disp_plot = plt.tripcolor
    elif plt_type == "contourf":
        disp_plot = plt.tricontourf

    plt.figure(figtitle)
    disp_plot(tri, field, levels, shading="gouraud")
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.axis("image")
    plt.grid()
    if savefigs:
        plt.savefig(filename)

def plot_node_field(field, nodes, elements, plt_type="contourf", levels=12,
              savefigs=False, title=None, figtitle=None, filename=None) :
    """Plot the nodal displacement using a triangulation

    Parameters
    ----------
    field : ndarray (float)
          Array with the field to be plotted. The number of columns
          determine the number of plots.
    nodes : ndarray (float)
        Array with number and nodes coordinates
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each  element.
    plt_type : string (optional)
        Plot the field as one of the options: ``pcolor`` or
        ``contourf``.
    levels : int (optional)
        Number of levels to be used in ``contourf``.
    savefigs : bool (optional)
        Allow to save the figure.
    title : Tuple of strings (optional)
        Titles of the plots. If not provided the plots will not have
        a title.
    figtitle : Tuple of strings (optional)
        Titles of the plotting windows. If not provided the
        windows will not have a title.
    filename : Tuple of strings (optional)
        Filenames to save the figures. Only used when `savefigs=True`.
        If not provided the name of the figures would be "outputk.pdf",
        where `k` is the number of the column.
    """
    tri = mesh2tri(nodes, elements)
    _, nfields = field.shape
    if title is None:
        title = ["" for cont in range(nfields)]
    if figtitle is None:
        figtitle = ["" for cont in range(nfields)]
    if filename is None:
        filename = ["output{}.pdf".format(cont) for cont in range(nfields)]
    for cont in range(nfields):
        tri_plot(tri, field[:, cont], title=title[cont],
             figtitle=figtitle[cont], levels=levels,
             plt_type=plt_type, savefigs=savefigs, filename=filename[cont])


#%% Auxiliar variables computation
def complete_disp(IBC, nodes, UG):
    """
    Fill the displacement vectors with imposed and computed values.

    IBC : ndarray (int)
        IBC (Indicator of Boundary Conditions) indicates if the
        nodes has any type of boundary conditions applied to it.
    UG : ndarray (float)
        Array with the computed displacements.
    nodes : ndarray (float)
        Array with number and nodes coordinates

    Returns
    -------
    UC : (nnodes, 2) ndarray (float)
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


def strain_nodes(nodes, elements, mats, UC):
    """Compute averaged strains and stresses at nodes

    First, the variable is extrapolated from the Gauss
    point to nodes for each element. Then, these values are averaged
    according to the number of element that share that node. The theory
    for this technique can be found in [1]_.

    Parameters
    ----------
    nodes : ndarray (float).
        Array with nodes coordinates.
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each element.
    mats : ndarray (float)
        Array with material profiles.
    UC : ndarray (float)
        Array with the displacements. This one contains both, the
        computed and imposed values.

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
    ne = elements.shape[0]
    nn = nodes.shape[0]
    iet = elements[0, 1]
    ndof, nnodes, _ = fe.eletype(iet)

    elcoor = np.zeros([nnodes, 2])
    E_nodes = np.zeros([nn, 3])
    S_nodes = np.zeros([nn, 3])
    el_nodes = np.zeros([nn], dtype=int)
    ul = np.zeros([ndof])
    IELCON = elements[:, 3:]

    for i in range(ne):
        young, poisson = mats[np.int(elements[i, 2]), :]
        shear = young/(2*(1 + poisson))
        fact1 = young/(1 - poisson**2)
        fact2 = poisson*young/(1 - poisson**2)
        for j in range(nnodes):
            elcoor[j, 0] = nodes[IELCON[i, j], 1]
            elcoor[j, 1] = nodes[IELCON[i, j], 2]
            ul[2*j] = UC[IELCON[i, j], 0]
            ul[2*j + 1] = UC[IELCON[i, j], 1]
        if iet == 1:
            epsG, _ = fe.str_el4(elcoor, ul)
        elif iet == 2:
            epsG, _ = fe.str_el6(elcoor, ul)
        elif iet == 3:
            epsG, _ = fe.str_el3(elcoor, ul)

        for cont, node in enumerate(IELCON[i, :]):
            E_nodes[node, 0] = E_nodes[node, 0] + epsG[cont, 0]
            E_nodes[node, 1] = E_nodes[node, 1] + epsG[cont, 1]
            E_nodes[node, 2] = E_nodes[node, 2] + epsG[cont, 2]
            S_nodes[node, 0] = S_nodes[node, 0] + fact1*epsG[cont, 0] \
                        + fact2*epsG[cont, 1]
            S_nodes[node, 1] = S_nodes[node, 1] + fact2*epsG[cont, 0] \
                        + fact1*epsG[cont, 1]
            S_nodes[node, 2] = S_nodes[node, 2] + shear*epsG[cont, 2]
            el_nodes[node] = el_nodes[node] + 1

    E_nodes[:, 0] = E_nodes[:, 0]/el_nodes
    E_nodes[:, 1] = E_nodes[:, 1]/el_nodes
    E_nodes[:, 2] = E_nodes[:, 2]/el_nodes
    S_nodes[:, 0] = S_nodes[:, 0]/el_nodes
    S_nodes[:, 1] = S_nodes[:, 1]/el_nodes
    S_nodes[:, 2] = S_nodes[:, 2]/el_nodes
    return E_nodes, S_nodes


def stress_truss(nodes, elements, mats, disp):
    r"""Compute axial stresses in truss members

    Parameters
    ----------
    nodes : ndarray (float).
        Array with nodes coordinates.
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each element.
    mats : ndarray (float)
        Array with material profiles.
    disp : ndarray (float)
        Array with the displacements. This one contains both, the
        computed and imposed values.

    Returns
    -------
    stress : ndarray
        Stresses for each member on the truss

    Examples
    --------

    The following examples are taken from [1]_. In all the examples
    :math:`A=1`, :math:`E=1`.

    >>> import assemutil as ass
    >>> import solidspy.solutil as sol

    >>> def fem_sol(nodes, elements, mats, loads):
    ...     DME , IBC , neq = ass.DME(nodes, elements)
    ...     KG = ass.assembler(elements, mats, nodes, neq, DME)
    ...     RHSG = ass.loadasem(loads, IBC, neq)
    ...     UG = sol.static_sol(KG, RHSG)
    ...     UC = complete_disp(IBC, nodes, UG)
    ...     return UC

    **Exercise 3.3-18**

    The axial stresses in this example are

    .. math::
        [\sigma] = \left[\frac{1}{2},\frac{\sqrt{3}}{4},\frac{1}{4}\right]


    >>> nodes = np.array([
    ... [0, 0.0,  0.0, 0, -1],
    ... [1, -1.0,  0.0, -1, -1],
    ... [2,  -np.cos(np.pi/6),  -np.sin(np.pi/6),  -1,  -1],
    ... [3,  -np.cos(np.pi/3),  -np.sin(np.pi/3),  -1,  -1]])
    >>> mats = np.array([[1.0, 1.0]])
    >>> elements = np.array([
    ... [0, 6, 0, 1, 0],
    ... [1, 6, 0, 2, 0],
    ... [2, 6, 0, 3, 0]])
    >>> loads = np.array([[0, 1.0, 0]])
    >>> disp = fem_sol(nodes, elements, mats, loads)
    >>> stress = stress_truss(nodes, elements, mats, disp)
    >>> stress_exact = np.array([1/2, np.sqrt(3)/4, 1/4])
    >>> np.allclose(stress_exact, stress)
    True

    **Exercise 3.3-19**

    The axial stresses in this example are

    .. math::

        [\sigma] = \left[\frac{1}{\sqrt{2}+2},
        \frac{\sqrt{2}}{\sqrt{2}+1},
        \frac{1}{\sqrt{2}+2}\right]

    >>> nodes = np.array([
    ...    [0, 0.0,  0.0, 0, 0],
    ...    [1, -1.0,  -1.0, -1, -1],
    ...    [2,  0.0,  -1.0,  -1,  -1],
    ...    [3,  1.0, -1.0,  -1,  -1]])
    >>> mats = np.array([[1.0, 1.0]])
    >>> elements = np.array([
    ...    [0, 6, 0, 1, 0],
    ...    [1, 6, 0, 2, 0],
    ...    [2, 6, 0, 3, 0]])
    >>> loads = np.array([[0, 0, 1]])
    >>> disp = fem_sol(nodes, elements, mats, loads)
    >>> stress = stress_truss(nodes, elements, mats, disp)
    >>> stress_exact = np.array([
    ...     1/(np.sqrt(2) + 2),
    ...     np.sqrt(2)/(np.sqrt(2) + 1),
    ...     1/(np.sqrt(2) + 2)])
    >>> np.allclose(stress_exact, stress)
    True

    **Exercise 3.3-22**

    The axial stresses in this example are

    .. math::

        [\sigma] =\left[\frac{1}{3 \sqrt{2}},\frac{5}{12},
            \frac{1}{2^{\frac{3}{2}}},
            \frac{1}{12},
            -\frac{1}{3 \sqrt{2}}\right]

    >>> cathetus = np.cos(np.pi/4)
    >>> nodes = np.array([
    ...    [0, 0.0,  0.0, 0, 0],
    ...    [1, -1.0,  0.0, -1, -1],
    ...    [2,  -cathetus,  cathetus,  -1,  -1],
    ...    [3,  0.0, 1.0,  -1,  -1],
    ...    [4,  cathetus, cathetus,  -1,  -1],
    ...    [5,  1.0, 0.0,  -1,  -1]])
    >>> mats = np.array([[1.0, 1.0]])
    >>> elements = np.array([
    ...    [0, 6, 0, 1, 0],
    ...    [1, 6, 0, 2, 0],
    ...    [2, 6, 0, 3, 0],
    ...    [3, 6, 0, 4, 0],
    ...    [4, 6, 0, 5, 0]])
    >>> loads = np.array([[0, cathetus, -cathetus]])
    >>> disp = fem_sol(nodes, elements, mats, loads)
    >>> stress = stress_truss(nodes, elements, mats, disp)
    >>> stress_exact = np.array([
    ...     1/(3*np.sqrt(2)),
    ...     5/12,
    ...     1/2**(3/2),
    ...     1/12,
    ...     -1/(3*np.sqrt(2))])
    >>> np.allclose(stress_exact, stress)
    True


    References
    ----------
    .. [1] William Weaver and James Gere. Matrix Analysis
        of Framed Structures. Third Edition, Van Nostrand
        Reinhold, New York (1990).

    """
    neles = elements.shape[0]
    stress = np.zeros((neles))
    for cont, elem in enumerate(elements):
        ini = elements[cont, 3]
        end = elements[cont, 4]
        coords = nodes[[ini, end], 1:3]
        tan_vec = coords[1, :] - coords[0, :]
        length = np.linalg.norm(tan_vec)
        mat_id = elements[cont, 2]
        local_stiff = uel.ueltruss2D(coords, * mats[mat_id, :])
        local_disp = np.hstack((disp[ini, :], disp[end, :]))
        local_forces = local_stiff.dot(local_disp)
        stress[cont] = local_forces[2:].dot(tan_vec)/(length*mats[mat_id, 1])
    return stress


def eigvals(mat, tol=1e-6):
    """Eigenvalues and eigenvectors for a 2x2 symmetric matrix/tensor

    Parameters
    ----------
    mat : ndarray
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

    >>> mat = np.array([[5, 6],
    ...              [6, 9]])
    >>> eig1, eig2, vec1, vec2 =  eigvals(mat)
    >>> np.allclose(eig1, 7 + 2*np.sqrt(10))
    True
    >>> np.allclose(eig2, 7 - 2*np.sqrt(10))
    True
    >>> np.allclose(vec1, np.array([-0.584710284663765, -0.8112421851755609]))
    True
    >>> np.allclose(vec2, np.array([-0.8112421851755609,0.584710284663765]))
    True

    """
    if np.abs(mat).max() < tol:
        eig1 = 0.0
        eig2 = 0.0
        vec1 = np.array([np.NaN, np.NaN])
        vec2 = np.array([np.NaN, np.NaN])
    elif abs(mat[0, 1])/np.abs(mat).max() < tol:
        eig1 = mat[0, 0]
        eig2 = mat[1, 1]
        vec1 = np.array([1, 0])
        vec2 = np.array([0, 1])
    else:
        trace = mat[0, 0] + mat[1, 1]
        det = mat[0, 0]*mat[1, 1] - mat[0, 1]**2
        eig1 = 0.5*(trace - np.sqrt(trace**2 - 4*det))
        eig2 = 0.5*(trace + np.sqrt(trace**2 - 4*det))
        vec1 = np.array([mat[0, 0] - eig2, mat[0, 1]])
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
    mat = np.zeros((2, 2))
    for cont, tensor in enumerate(field):
        mat[0, 0] = tensor[0]
        mat[1, 1] = tensor[1]
        mat[0, 1] = tensor[2]
        eig1, eig2, vec1, vec2 = eigvals(mat, tol=1e-6)
        eigs1[cont] = eig1
        eigs2[cont] = eig2
        vecs1[cont, :] = vec1
        vecs2[cont, :] = vec2

    return eigs1, eigs2, vecs1, vecs2


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
    filename = folder + 'out.txt'
    np.savetxt(filename, UR, fmt='%.18e', delimiter=' ')


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
    f = KG.dot(UG)
    EFE = -0.5*f.dot(UG)

    return EFE


#%% Doc-testing
if __name__ == "__main__":
    import doctest
    doctest.testmod()
