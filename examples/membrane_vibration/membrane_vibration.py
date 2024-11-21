# coding: utf-8
"""
# Creation of new elements
SolidsPy supports user elements, that is, you can pass a function that returns
the mass and stiffness matrices for your element to the assembly operator.

This example shows how to create an element for membrane deflection.
## Equation to solve
We want to find the vibration modes of vibration modes of a thin membrane. This
is described by the Helmholtz equation

$$\nabla^2 u + k^2 u = 0\quad \forall \mathbf{x} \in \Omega\, .$$

The weak form for this problem is

$$
\int_\Omega \nabla u \cdot \nabla u \mathrm{d}\Omega
= k^2 \int_\Omega u^2 \mathrm{d}\Omega\, .
$$

And the discretized system is given by

$$[K]\{U\} = k^2[M]\{U\}\, .$$
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
import meshio
import solidspy.assemutil as ass
import solidspy.femutil as fem
import solidspy.gaussutil as gau
import solidspy.postprocesor as pos


def acoust_diff(r, s, coord, element):
    """
    Interpolation matrices for elements for acoustics

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
        Array with the gradient matrix evaluated
        at the point (r, s).
    det : float
        Determinant of the Jacobian.
    """
    N, dNdr = element(r, s)
    N.shape = 1, N.shape[0]
    det, jaco_inv = fem.jacoper(dNdr, coord)
    dNdx = jaco_inv @ dNdr
    return N, dNdx, det


def acoust_tri6(coord, params):
    """
    Triangular element with 6 nodes for acoustics under
    axisymmetric conditions.

    Parameters
    ----------
    coord : coord
        Coordinates of the element.
    params : list
        List with material parameters in the following order:
        [Speed].

    Returns
    -------
    stiff_mat : ndarray (float)
        Local stifness matrix.
    mass_mat : ndarray (float)
        Local mass matrix.
    """
    
    speed = params
    stiff_mat = np.zeros((6, 6))
    mass_mat = np.zeros((6, 6))
    gpts, gwts = gau.gauss_tri(order=3)
    for cont in range(gpts.shape[0]):
        r = gpts[cont, 0]
        s = gpts[cont, 1]
        H, B, det = acoust_diff(r, s, coord, fem.shape_tri6)
        factor = det * gwts[cont]
        stiff_mat += 0.5 * speed**2 * factor * (B.T @ B)
        mass_mat += 0.5 * factor * (H.T @ H)
    return stiff_mat, mass_mat


if __name__ == "__main__":
    mesh = meshio.read("square.msh")
    points = mesh.points
    cells = mesh.cells
    tri6 = cells["triangle6"]
    line3 = cells["line3"]
    npts = points.shape[0]
    nels = tri6.shape[0]

    nodes = np.zeros((npts, 3))
    nodes[:, 1:] = points[:, 0:2] 


    # Constraints
    line_nodes = list(set(line3.flatten()))
    cons = np.zeros((npts, 1), dtype=int)
    cons[line_nodes, :] = -1


    # Elements
    elements = np.zeros((nels, 9), dtype=int)
    elements[:, 1] = 2
    elements[:, 3:] = tri6

    # Materials
    mats = np.array([[1.0]])


    # Assembly
    assem_op, bc_array, neq = ass.DME(cons, elements,
                                    ndof_node=1, ndof_el_max=6)
    stiff_mat, mass_mat = ass.assembler(elements, mats, nodes, neq,
                                        assem_op, uel=acoust_tri6)


    # Solution
    eigvals, eigvecs = eigsh(stiff_mat, M=mass_mat, k=10, which="LM",
                            sigma=1e-6)

    eigvals_exact = np.array(sorted([m**2 + n**2
                                    for m in range(1, 10)
                                    for n in range(1, 10)])[:10])


    plt.figure()
    plt.plot(eigvals_exact, "ko")
    plt.plot(eigvals, "r.")
    plt.xlabel("Eigenvalue number")
    plt.ylabel("Eigenvale")

    # Visualization

    sol = pos.complete_disp(bc_array, nodes, eigvecs[:, 0], ndof_node=1)
    pos.plot_node_field(sol[:, 0], nodes, elements)


    # Exporting to VTK
    for cont in range(10):
        aux = pos.complete_disp(bc_array, nodes, eigvecs[:, cont],
                        ndof_node=1)
        mesh.point_data["mode_%d" % cont] = aux

    meshio.write("membrane.vtk", mesh)

    plt.show()