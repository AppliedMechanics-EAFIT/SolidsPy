# -*- coding: utf-8 -*-
"""
Assembly routines
-----------------

Functions to assemble the system of equations for a finite element
analysis.

"""
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import coo_matrix, csr_matrix, spmatrix
from typing import Callable, Optional, Tuple
import src.uelutil as ue
import src.femutil as fem


def eqcounter(
    cons: NDArray[np.int_], ndof_node: int = 2
) -> Tuple[int, NDArray[np.int_]]:
    """Count active equations

    Creates boundary conditions array bc_array

    Parameters
    ----------
    cons : ndarray.
      Array with constraints for each node.

    ndof_node : int, optional
      Number of degrees of freedom per node. Default is 2.

    Returns
    -------
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    bc_array : ndarray (int)
      Array that maps the nodes with number of equations.

    """
    nnodes = cons.shape[0]
    bc_array = cons.copy().astype(int)
    neq = 0
    for i in range(nnodes):
        for j in range(ndof_node):
            if bc_array[i, j] == 0:
                bc_array[i, j] = neq
                neq += 1

    return neq, bc_array


def DME(
    cons: NDArray[np.int_],
    elements: NDArray[np.int_],
    ndof_node: int = 2,
    ndof_el_max: int = 18,
    ndof_el: Optional[Callable[[int], int]] = None,
) -> Tuple[NDArray[np.int_], NDArray[np.int_], int]:
    """Create assembly array operator

    Count active equations, create boundary conditions array ``bc_array``
    and the assembly operator ``assem_op``.

    Parameters
    ----------
    cons : ndarray.
      Array with constraints for each degree of freedom in each node.
    elements : ndarray
      Array with the number for the nodes in each element.
    ndof_node : int, optional
      Number of degrees of freedom per node. By default it is 2.
    ndof_el_max : int, optional
      Number of maximum degrees of freedom per element. By default it is
      18.
    ndof_el : callable, optional
      Function that returns the number of degrees of freedom for elements. It
      is needed for user elements.

    Returns
    -------
    assem_op : ndarray (int)
      Assembly operator.
    bc_array : ndarray (int)
      Boundary conditions array.
    neq : int
      Number of active equations in the system.

    """
    nels = elements.shape[0]
    assem_op = np.zeros([nels, ndof_el_max], dtype=int)
    neq, bc_array = eqcounter(cons, ndof_node=ndof_node)
    for ele in range(nels):
        iet = elements[ele, 1]
        if ndof_el is None:
            ndof, _, _ = fem.eletype(iet)
        else:
            ndof = ndof_el(iet)
        assem_op[ele, :ndof] = bc_array[elements[ele, 3:]].flatten()
    return assem_op, bc_array, neq


def ele_fun(eletype: int) -> Callable[[NDArray[np.float64], NDArray[np.float64]], Tuple[NDArray[np.float64], NDArray[np.float64]]]:
    """Return the function for the element type given

    Parameters
    ----------
    eletype : int
      Type of element.

    Returns
    -------
    ele_fun : callable
      Function for the element type given.
    """
    elem_id = {
        1: ue.elast_quad4,
        2: ue.elast_tri6,
        3: ue.elast_tri3,
        4: ue.elast_quad9,
        5: ue.spring,
        6: ue.truss2D,
        7: ue.beam2DU,
        8: ue.beam2D,
    }
    try:
        return elem_id[eletype]
    except KeyError:
        raise ValueError("You entered an invalid type of element.")


def retriever(
    elements: NDArray[np.int_],
    mats: NDArray[np.float64],
    nodes: NDArray[np.float64],
    ele: int,
    uel: Optional[
        Callable[[NDArray[np.float64], NDArray[np.float64]], Tuple[NDArray[np.float64], NDArray[np.float64]]]
    ] = None,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Computes the elemental stiffness matrix of element ``ele``

    Parameters
    ----------
    elements : ndarray
      Array with the number for the nodes in each element.
    mats : ndarray.
      Array with the material profiles.
    nodes : ndarray.
      Array with the nodal numbers and coordinates.
    ele : int.
      Identifier of the element to be assembled.
    uel : callable, optional
      User element function that returns stiffness and mass matrices.

    Returns
    -------
    kloc : ndarray (float)
      Array with the local stiffness matrix.
    mloc : ndarray (float)
      Array with the local mass matrix.
    """
    elem_type = elements[ele, 1]
    params = mats[elements[ele, 2], :]
    elcoor = nodes[elements[ele, 3:], 1:]
    if uel is None:
        uel = ele_fun(elem_type)
    kloc, mloc = uel(elcoor, params)
    return kloc, mloc


def assembler(
    elements: NDArray[np.int_],
    mats: NDArray[np.float64],
    nodes: NDArray[np.float64],
    neq: int,
    assem_op: NDArray[np.int_],
    sparse: bool = True,
    uel: Optional[
        Callable[[NDArray[np.float64], NDArray[np.float64]], Tuple[NDArray[np.float64], NDArray[np.float64]]]
    ] = None,
) -> Tuple[spmatrix, spmatrix]:
    """Assembles the global stiffness matrix

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    assem_op : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    sparse : boolean (optional)
      Boolean variable to pick sparse assembler. It is True
      by default.
    uel : callable function (optional)
      Python function that returns the local stiffness matrix.

    Returns
    -------
    kglob : sparse matrix (float)
      Array with the global stiffness matrix. It is sparse in CSR format.
    mglob : sparse matrix (float)
      Array with the global mass matrix. It is sparse in CSR format.

    """
    if sparse:
        kglob, mglob = sparse_assem(elements, mats, nodes, neq, assem_op,
                                    uel=uel)
    else:
        kglob, mglob = dense_assem(elements, mats, nodes, neq, assem_op,
                                   uel=uel)

    return kglob, mglob


def dense_assem(
    elements: NDArray[np.int_],
    mats: NDArray[np.float64],
    nodes: NDArray[np.float64],
    neq: int,
    assem_op: NDArray[np.int_],
    uel: Optional[
        Callable[[NDArray[np.float64], NDArray[np.float64]], Tuple[NDArray[np.float64], NDArray[np.float64]]]
    ] = None,
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Assembles the global stiffness matrix
    using a dense storing scheme

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats : ndarray (float)
      Array with the material profiles.
    nodes : ndarray (float)
      Array with the nodal numbers and coordinates.
    assem_op : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    uel : callable, optional
      Python function that returns the local stiffness matrix.

    Returns
    -------
    kglob : ndarray (float)
      Array with the global stiffness matrix in a dense numpy
      array.
    mglob : ndarray (float)
      Array with the global mass matrix in a dense numpy
      array.

    """
    kglob = np.zeros((neq, neq))
    mglob = np.zeros((neq, neq))
    nels = elements.shape[0]
    for ele in range(nels):
        kloc, mloc = retriever(elements, mats, nodes, ele, uel=uel)
        ndof = kloc.shape[0]
        dme = assem_op[ele, :ndof]
        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        kglob[glob_row, glob_col] += kloc[row, col]
                        mglob[glob_row, glob_col] += mloc[row, col]

    return kglob, mglob


def sparse_assem(
    elements: NDArray[np.int_],
    mats: NDArray[np.float64],
    nodes: NDArray[np.float64],
    neq: int,
    assem_op: NDArray[np.int_],
    uel: Optional[
        Callable[[NDArray[np.float64], NDArray[np.float64]], Tuple[NDArray[np.float64], NDArray[np.float64]]]
    ] = None,
) -> Tuple[csr_matrix, csr_matrix]:
    """
    Assembles the global stiffness matrix
    using a sparse storing scheme

    The scheme used to assemble is COOrdinate list (COO), and
    it is converted to Compressed Sparse Row (CSR) afterward
    for the solution phase [1]_.

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    assem_op : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    uel : callable, optional
      Python function that returns the local stiffness matrix.

    Returns
    -------
    kglob : sparse matrix (float)
      Array with the global stiffness matrix in a sparse
      Compressed Sparse Row (CSR) format.
    mglob : sparse matrix (float)
      Array with the global mass matrix in a sparse
      Compressed Sparse Row (CSR) format.

    References
    ----------
    .. [1] Sparse matrix. (2017, March 8). In Wikipedia,
        The Free Encyclopedia.
        https://en.wikipedia.org/wiki/Sparse_matrix

    """
    rows = []
    cols = []
    stiff_vals = []
    mass_vals = []
    nels = elements.shape[0]
    for ele in range(nels):
        kloc, mloc = retriever(elements, mats, nodes, ele, uel=uel)
        ndof = kloc.shape[0]
        dme = assem_op[ele, :ndof]
        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        rows.append(glob_row)
                        cols.append(glob_col)
                        stiff_vals.append(kloc[row, col])
                        mass_vals.append(mloc[row, col])

    stiff = coo_matrix((stiff_vals, (rows, cols)), shape=(neq, neq)).tocsr()
    mass = coo_matrix((mass_vals, (rows, cols)), shape=(neq, neq)).tocsr()
    return stiff, mass


def loadasem(
    loads: NDArray[np.float64],
    bc_array: NDArray[np.int_],
    neq: int,
    ndof_node: int = 2,
) -> NDArray[np.float64]:
    """Assembles the global Right Hand Side Vector

    Parameters
    ----------
    loads : ndarray
      Array with the loads imposed in the system.
    bc_array : ndarray (int)
      Array that maps the nodes with number of equations.
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    ndof_node : int, optional
      Number of degrees of freedom per node. By default it is 2.

    Returns
    -------
    rhs_vec : ndarray
      Array with the right hand side vector.

    """
    nloads = loads.shape[0]
    rhs_vec = np.zeros([neq])
    for cont in range(nloads):
        node = int(loads[cont, 0])
        for dof in range(ndof_node):
            dof_id = bc_array[node, dof]
            if dof_id != -1:
                rhs_vec[dof_id] = loads[cont, dof + 1]
    return rhs_vec


if __name__ == "__main__":
    import doctest
    doctest.testmod()
