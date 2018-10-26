# -*- coding: utf-8 -*-
"""
Assembly routines
-----------------

Functions to assemble the system of equations for a finite element
analysis.

"""
from __future__ import absolute_import, division, print_function
import numpy as np
from scipy.sparse import coo_matrix
import solidspy.uelutil as ue
import solidspy.femutil as fem


def eqcounter(cons, ndof_node=2):
    """Count active equations

    Creates boundary conditions array bc_array

    Parameters
    ----------
    cons : ndarray.
      Array with constraints for each node.

    Returns
    -------
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    bc_array : ndarray (int)
      Array that maps the nodes with number of equations.

    """
    nnodes = cons.shape[0]
    bc_array = cons.copy()
    neq = 0
    for i in range(nnodes):
        for j in range(ndof_node):
            if bc_array[i, j] == 0:
                bc_array[i, j] = neq
                neq += 1

    return neq, bc_array


def DME(cons, elements, ndof_node=2):
    """Create assembly array operator

    Count active equations, create boundary conditions array ``bc_array``
    and the assembly operator DME.

    Parameters
    ----------
    cons : ndarray.
      Array with constraints for each degree of freedom in each node.
    elements : ndarray
      Array with the number for the nodes in each element.

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
    assem_op = np.zeros([nels, 18], dtype=np.integer)
    neq, bc_array = eqcounter(cons, ndof_node=2)
    for ele in range(nels):
        iet = elements[ele, 1]
        ndof, _, _ = fem.eletype(iet)
        assem_op[ele, :ndof] = bc_array[elements[ele, 3:]].flatten()
    return assem_op, bc_array, neq


def ele_fun(eletype):
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
        1: ue.uel4nquad,
        2: ue.uel6ntrian,
        3: ue.uel3ntrian,
        5: ue.uelspring,
        6: ue.ueltruss2D,
        7: ue.uelbeam2DU}
    try:
        return elem_id[eletype]
    except:
        raise ValueError("You entered an invalid type of element.")



def retriever(elements, mats, nodes, ele, uel=None):
    """Computes the elemental stiffness matrix of element i

    Parameters
    ----------
    elements : ndarray
      Array with the number for the nodes in each element.
    mats : ndarray.
      Array with the material profiles.
    nodes : ndarray.
      Array with the nodal numbers and coordinates.
    i : int.
      Identifier of the element to be assembled.

    Returns
    -------
    kloc : ndarray (float)
      Array with the local stiffness matrix.
    mloc : ndarray (float)
      Array with the local mass matrix.
    """
    elem_type = elements[ele, 1]
    params = mats[elements[ele, 2], :]
    elcoor = nodes[elements[ele, 3:], 1:3]
    if uel is None:
        uel = ele_fun(elem_type)
    kloc, mloc = uel(elcoor, params)
    return kloc, mloc


def assembler(elements, mats, nodes, neq, DME, sparse=True, uel=None):
    """Assembles the global stiffness matrix

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    DME  : ndarray (int)
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
    kglob : ndarray (float)
      Array with the global stiffness matrix. It might be
      dense or sparse, depending on the value of _sparse_

    """
    if sparse:
        kglob = sparse_assem(elements, mats, nodes, neq, DME, uel=uel)
    else:
        kglob = dense_assem(elements, mats, nodes, neq, DME, uel=uel)

    return kglob


def dense_assem(elements, mats, nodes, neq, DME, uel=None):
    """
    Assembles the global stiffness matrix _KG_
    using a dense storing scheme

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    DME  : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    uel : callable function (optional)
      Python function that returns the local stiffness matrix.

    Returns
    -------
    kglob : ndarray (float)
      Array with the global stiffness matrix in a dense numpy
      array.

    """
    kglob = np.zeros((neq, neq))
    nels = elements.shape[0]
    for ele in range(nels):
        kloc, _ = retriever(elements, mats, nodes, ele, uel=uel)
        ndof = kloc.shape[0]
        dme = DME[ele, :ndof]
        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        kglob[glob_row, glob_col] += kloc[row, col]

    return kglob


def sparse_assem(elements, mats, nodes, neq, DME, uel=None):
    """
    Assembles the global stiffness matrix _KG_
    using a sparse storing scheme

    The scheme used to assemble is COOrdinate list (COO), and
    it converted to Compressed Sparse Row (CSR) afterward
    for the solution phase [1]_.

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    DME  : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    uel : callable function (optional)
      Python function that returns the local stiffness matrix.

    Returns
    -------
    KG : ndarray (float)
      Array with the global stiffness matrix in a sparse
      Compressed Sparse Row (CSR) format.

    References
    ----------
    .. [1] Sparse matrix. (2017, March 8). In Wikipedia,
        The Free Encyclopedia.
        https://en.wikipedia.org/wiki/Sparse_matrix

    """
    rows = []
    cols = []
    vals = []
    nels = elements.shape[0]
    for ele in range(nels):
        kloc, _ = retriever(elements, mats, nodes, ele, uel=uel)
        ndof = kloc.shape[0]
        dme = DME[ele, :ndof]

        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        rows.append(glob_row)
                        cols.append(glob_col)
                        vals.append(kloc[row, col])

    return coo_matrix((vals, (rows, cols)), shape=(neq, neq)).tocsr()


def loadasem(loads, IBC, neq):
    """Assembles the global Right Hand Side Vector RHSG

    Parameters
    ----------
    loads : ndarray
      Array with the loads imposed in the system.
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.

    Returns
    -------
    RHSG : ndarray
      Array with the right hand side vector.

    """
    nloads = loads.shape[0]
    RHSG = np.zeros([neq])
    for i in range(nloads):
        il = int(loads[i, 0])
        ilx = IBC[il, 0]
        ily = IBC[il, 1]
        if ilx != -1:
            RHSG[ilx] = loads[i, 1]
        if ily != -1:
            RHSG[ily] = loads[i, 2]

    return RHSG


if __name__ == "__main__":
    import doctest
    doctest.testmod()
