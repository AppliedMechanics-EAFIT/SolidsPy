# -*- coding: utf-8 -*-
"""
assemutil.py
------------

Functions to assemble the system of equations for the Finite Element
Analysis.

"""
from __future__ import division, print_function
import numpy as np
import uelutil as ue
import femutil as fem


def eqcounter(nn, nodes):
    """Counts active equations and creates BCs array IBC

    Parameters
    ----------
    nn : int
      Number of nodes.
    nodes : ndarray
      Array with nodes coordinates and boundary conditions.

    Returns
    -------
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.

    """
    IBC = np.zeros([nn, 2], dtype=np.integer)
    for i in range(nn):
        for k in range(2):
            IBC[i , k] = int(nodes[i , k+3])
    neq = 0
    for i in range(nn):
        for j in range(2):
            if IBC[i, j] == 0:
                IBC[i, j] = neq
                neq = neq + 1

    return neq, IBC


def DME(nn , ne , nodes , elements):
    """Counts active equations, creates BCs array IBC[]
    and the assembly operator DME[]

    Parameters
    ----------
    nn : int
      Number of nodes.
    ne : int
      Number of elements.
    nodes    : ndarray.
      Array with the nodal numbers and coordinates.
    elements : ndarray
      Array with the number for the nodes in each element.

    Returns
    -------
    DME : ndarray (int)
      Assembly operator.
    IBC : ndarray (int)
      Boundary conditions array.
    neq : int
      Number of activ equations in the system.

    """
    IELCON = np.zeros([ne, 9], dtype=np.integer)
    DME = np.zeros([ne, 18], dtype=np.integer)
#
    neq, IBC = eqcounter(nn, nodes)
#
    for i in range(ne):
        iet = elements[i, 1]
        ndof, nnodes, ngpts = fem.eletype(iet)
        for j in range(nnodes):
            IELCON[i, j] = elements[i, j+3]
            kk = IELCON[i, j]
            for l in range(2):
                DME[i, 2*j+l] = IBC[kk, l]

    return DME , IBC , neq


def retriever(elements , mats , nodes , i):
    """Computes the elemental stiffness matrix of element i

    Parameters
    ----------
    elements : ndarray
      Array with the number for the nodes in each element.
    mats    : ndarray.
      Array with the material profiles.
    nodes    : ndarray.
      Array with the nodal numbers and coordinates.
    i    : int.
      Identifier of the element to be assembled.

    Returns
    -------
    kloc : ndarray (float)
      Array with the local stiffness matrix.
    ndof : int.
      Number of degrees of fredom of the current element.
    """
    IELCON = np.zeros([9], dtype=np.integer)
    iet = elements[i , 1]
    ndof, nnodes, ngpts = fem.eletype(iet)
    elcoor = np.zeros([nnodes, 2])
    im = np.int(elements[i, 2])
    par0 = mats[im, 0]
    par1 = mats[im, 1]
    for j in range(nnodes):
        IELCON[j] = elements[i, j+3]
        elcoor[j, 0] = nodes[IELCON[j], 1]
        elcoor[j, 1] = nodes[IELCON[j], 2]
    if iet == 1:
        kloc = ue.uel4nquad(elcoor, par1 , par0)
    elif iet == 2:
        kloc = ue.uel6ntrian(elcoor, par1 , par0)
    elif iet == 3:
        kloc = ue.uel3ntrian(elcoor, par1 , par0)
    elif iet == 5:
        kloc = ue.uelspring(elcoor, par1 , par0)
    elif iet == 6:
        kloc =ue.ueltruss2D(elcoor, par1 , par0)
    elif iet == 7:
        kloc =ue.uelbeam2DU(elcoor, par1 , par0)
    
    return kloc , ndof , iet


def assembler(KG , neq , kloc , ndof , DME , iet , i):
    """Assembles the global stiffness matrix KG[]

    Parameters
    ----------
    KG : ndarray (float)
      Array with the current values of the stiffness matrix.
    neq : int
      Total number of equations in the system.
    kloc : ndarray (float)
      Array with elemental stiffness matrix to be assembled.
      with imposed displacements.
    ndof : int
      Number of degrees of freedom of the element to be assembled.
    DME  : ndarray (int)
      Assembly operator.
    i    : int.
      Identifier of the element to be assembled.

    Returns
    -------
    KGLOB : ndarray (float)
      Array with the global stiffness matrix.

    """
    dme    = np.zeros([ndof], dtype=np.integer)
    if iet == 6:
        dme[0] = DME[i, 0]
        dme[1] = DME[i, 1]
        dme[2] = DME[i, 3]
        dme[3] = DME[i, 4]
    else:
        for ii in range(ndof):
            dme[ii] = DME[i, ii]
#    
    for ii in range(ndof):
        kk = dme[ii]
        if kk != -1:
            for jj in range(ndof):
                ll = dme[jj]
                if ll != -1:
                    KG[kk, ll] = KG[kk, ll] + kloc[ii, jj]
    
    return KG


def coo_assem(rows, cols, vals , neq , kloc , ndof , DME , iet , i):
    """Assembles the global stiffness matrix KG[]

    Parameters
    ----------
    KG : ndarray (float)
      Array with the current values of the stiffness matrix.
    neq : int
      Total number of equations in the system.
    kloc : ndarray (float)
      Array with elemental stiffness matrix to be assembled.
      with imposed displacements.
    ndof : int
      Number of degrees of freedom of the element to be assembled.
    DME  : ndarray (int)
      Assembly operator.
    i    : int.
      Identifier of the element to be assembled.

    Returns
    -------
    KGLOB : ndarray (float)
      Array with the global stiffness matrix.

    """
    dme = np.zeros([ndof], dtype=np.integer)
    if iet == 6:
        dme[0] = DME[i, 0]
        dme[1] = DME[i, 1]
        dme[2] = DME[i, 3]
        dme[3] = DME[i, 4]
    else:
        dme = DME[i, :ndof]

    for row in range(ndof):
        glob_row = dme[row]
        if glob_row != -1:
            for col in range(ndof):
                glob_col = dme[col]
                if glob_col != -1:
                    rows.append(glob_row)
                    cols.append(glob_col)
                    vals.append(kloc[row, col])

    return None


def loadasem(loads, IBC, neq, nl):
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
    nl : int
      Number of loads.

    Returns
    -------
    RHSG : ndarray
      Array with the right hand side vector.

    """
    RHSG = np.zeros([neq])
    for i in range(nl):
        il = int(loads[i, 0])
        ilx = IBC[il, 0]
        ily = IBC[il, 1]
        if ilx != -1:
            RHSG[ilx] = loads[i, 1]
        if ily != -1:
            RHSG[ily] = loads[i, 2]

    return RHSG
