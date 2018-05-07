# -*- coding: utf-8 -*-
"""
Postprocessor subroutines
-------------------------

This module contains functions to postprocess results.

"""
from __future__ import absolute_import, division, print_function
import numpy as np
import solidspy.femutil as fe
import solidspy.uelutil as uel
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# Set plotting defaults
gray = '#757575'
plt.rcParams['image.cmap'] = "YlGnBu_r"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["text.color"] = gray
plt.rcParams["font.size"] = 12
plt.rcParams["xtick.color"] = gray
plt.rcParams["ytick.color"] = gray
plt.rcParams["axes.labelcolor"] = gray
plt.rcParams["axes.edgecolor"] = gray
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False


#%% Plotting routines
def fields_plot(elements, nodes, disp, E_nodes=None, S_nodes=None):
    """Plot contours for displacements, strains and stresses

    Parameters
    ----------
    nodes : ndarray (float)
        Array with number and nodes coordinates:
         `number coordX coordY BCX BCY`
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each element.
    disp : ndarray (float)
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
        print(disp)
    else:
        plot_node_field(disp, nodes, elements, title=[r"$u_x$", r"$u_y$"],
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




def tri_plot(tri, field, title="", levels=12, savefigs=False,
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
    disp_plot(tri, field, levels, shading="gouraud")
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.axis("image")
    if savefigs:
        plt.savefig(filename)


def plot_node_field(field, nodes, elements, plt_type="contourf", levels=12,
                    savefigs=False, title=None, figtitle=None,
                    filename=None):
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
    if len(field.shape) == 1:
        nfields = 1
    else:
        _, nfields = field.shape
    if title is None:
        title = ["" for cont in range(nfields)]
    if figtitle is None:
        figs = plt.get_fignums()
        nfigs = len(figs)
        figtitle = [cont + 1 for cont in range(nfigs, nfigs + nfields)]
    if filename is None:
        filename = ["output{}.pdf".format(cont) for cont in range(nfields)]
    for cont in range(nfields):
        if nfields == 1:
            current_field = field
        else:
            current_field = field[:, cont]
        plt.figure(figtitle[cont])
        tri_plot(tri, current_field, title=title[cont], levels=levels,
                 plt_type=plt_type, savefigs=savefigs,
                 filename=filename[cont])
        if savefigs:
            plt.savefig(filename[cont])


def plot_truss(nodes, elements, mats, stresses=None, max_val=4,
               min_val=0.5, savefigs=False, title=None, figtitle=None,
               filename=None):
    """Plot a truss and encodes the stresses in a colormap

    Parameters
    ----------
    UC : (nnodes, 2) ndarray (float)
      Array with the displacements.
    nodes : ndarray (float)
        Array with number and nodes coordinates
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each  element.
    mats : ndarray (float)
        Array with material profiles.
    loads : ndarray (float)
        Array with loads.
    tol : float (optional)
        Minimum difference between cross-section areas of the members
        to be considered different.
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
    min_area = mats[:, 1].min()
    max_area = mats[:, 1].max()
    areas = mats[:, 1].copy()
    if stresses is None:
        scaled_stress = np.ones_like(elements[:, 0])
    else:
        max_stress = max(-stresses.min(), stresses.max())
        scaled_stress = 0.5*(stresses + max_stress)/max_stress
    if max_area - min_area > 1e-6:
        widths = (max_val - min_val)*(areas - min_area)/(max_area - min_area)\
            + min_val
    else:
        widths = 3*np.ones_like(areas)
    plt.figure(figtitle)
    for elem in elements:
        ini, end = elem[3:]
        color = plt.cm.seismic(scaled_stress[elem[0]])
        plt.plot([nodes[ini, 1], nodes[end, 1]],
                 [nodes[ini, 2], nodes[end, 2]],
                 color=color, lw=widths[elem[2]])

    if title is None:
        title = ''
    if figtitle is None:
        figtitle = ""
    if filename is None:
        filename = "output.pdf"
    plt.title(title)
    plt.axis("image")
    if savefigs:
        plt.savefig(filename)


#%% Auxiliar functions for plotting
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
    for elem in elements:
        if elem[1] == 1:
            triangs.append(elem[[3, 4, 5]])
            triangs.append(elem[[5, 6, 3]])
        if elem[1] == 2:
            triangs.append(elem[[3, 6, 8]])
            triangs.append(elem[[6, 7, 8]])
            triangs.append(elem[[6, 4, 7]])
            triangs.append(elem[[7, 5, 8]])
        if elem[1] == 3:
            triangs.append(elem[3:])

    tri = Triangulation(coord_x, coord_y, np.array(triangs))
    return tri


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
    nnodes = nodes.shape[0]
    UC = np.zeros([nnodes, 2], dtype=np.float)
    for row in range(nnodes):
        for col in range(2):
            cons = IBC[row, col]
            if cons == -1:
                UC[row, col] = 0.0
            else:
                UC[row, col] = UG[cons]

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
    nelems = elements.shape[0]
    nnodes = nodes.shape[0]
    iet = elements[0, 1]
    ndof, nnodes_elem, _ = fe.eletype(iet)

    elcoor = np.zeros([nnodes_elem, 2])
    E_nodes = np.zeros([nnodes, 3])
    S_nodes = np.zeros([nnodes, 3])
    el_nodes = np.zeros([nnodes], dtype=int)
    ul = np.zeros([ndof])
    IELCON = elements[:, 3:]

    for i in range(nelems):
        young, poisson = mats[np.int(elements[i, 2]), :]
        shear = young/(2*(1 + poisson))
        fact1 = young/(1 - poisson**2)
        fact2 = poisson*young/(1 - poisson**2)
        for j in range(nnodes_elem):
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
    for cont in range(neles):
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


def energy(disp, stiff):
    r"""
    Computes the potential energy for the current solution.

    Parameters
    ----------
    disp : ndarray (float)
        Array with the computed displacements.
    stiff : ndarray (float)
        Global stiffness matrix.

    Returns
    -------
    el_energy : scalar (float)
        Total energy in the system. :math:`-\frac{1}{2} U^T K U`

    """
    force = stiff.dot(disp)
    el_energy = -0.5*force.dot(disp)

    return el_energy


#%% Doc-testing
if __name__ == "__main__":
    import doctest
    doctest.testmod()
