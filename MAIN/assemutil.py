# -*- coding: utf-8 -*-
"""
assemutil.py
------------

Functions to assemble the system of equations for the Finite Element
Analysis.

"""
from __future__ import division
import numpy as np
from sympy import *
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
    neq = 0
    for i in range(nn):
        for j in range(2):
            IBC[i, j] = int(nodes[i, j+3])
            if IBC[i, j] == 0:
                IBC[i, j] = neq
                neq = neq + 1

    return neq, IBC


def DME(IBC, ne, elements):
    """Create the assembly operator

    Create the assembly operator DME and processes the element
    connectivity array IELCON

    Parameters
    ----------
    IBC : ndarray
      Array that maps the nodes with number of equations.
    ne : int
      Number of elements.
    elements : ndarray
      Array with the number for the nodes in each element.

    Returns
    -------
    DME : ndarray (int)
      Assembly operator.
    IELCON : ndarray (int)
      Element connectivity.

    """
    IELCON = np.zeros([ne, 9], dtype=np.integer)
    DME = np.zeros([ne, 18], dtype=np.integer)

    for i in range(ne):
        iet = elements[i, 1]
        ndof, nnodes, ngpts = fem.eletype(iet)
        for j in range(nnodes):
            IELCON[i, j] = elements[i, j+3]
            kk = IELCON[i, j]
            for l in range(2):
                DME[i, 2*j+l] = IBC[kk, l]

    return DME, IELCON


def matassem(IBC, mats, elements, nn, ne, neq, COORD, DME, IELCON):
    """Assembles the global stiffness matrix KG

    Parameters
    ----------
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.
    mats : ndarray
      Material properties.
    elements : ndarray
      Array with the number for the nodes in each element.
    nn : int
      Number of nodes.
    ne : int
      Number of elements.
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    COORD : ndarray
      Coordinates of the nodes.
    DME : ndarray (int)
      Assembly operator.
    IELCON : ndarray (int)
      Element connectivity.

    Returns
    -------
    KG : ndarray
      Global stiffness matrix.

    """
    KG = np.zeros([neq, neq])
    for i in range(ne):
        iet = elements[i, 1]
        ndof, nnodes, ngpts = fem.eletype(iet)
        elcoor = np.zeros([nnodes, 2])
        kloc = np.zeros([ndof, ndof])
        im = elements[i, 2]
        emod = mats[im, 0]
        enu = mats[im, 1]
        dme = np.zeros([ndof], dtype=np.integer)
        for j in range(nnodes):
            elcoor[j, 0] = COORD[IELCON[i, j], 0]
            elcoor[j, 1] = COORD[IELCON[i, j], 1]
        if iet == 1:
            kloc = ue.uel4nquad(elcoor, enu, emod)
        elif iet == 2:
            kloc = ue.uel6ntrian(elcoor, enu, emod)
        elif iet == 3:
            kloc = ue.uel3ntrian(elcoor, enu, emod)
        for ii in range(ndof):
            dme[ii] = DME[i, ii]
        for ii in range(ndof):
            kk = dme[ii]
            if kk != -1:
                for jj in range(ndof):
                    ll = dme[jj]
                    if ll != -1:
                        KG[kk, ll] = KG[kk, ll] + kloc[ii, jj]

    return KG


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
        RHSG[ilx] = loads[i, 1]
        RHSG[ily] = loads[i, 2]

    return RHSG
