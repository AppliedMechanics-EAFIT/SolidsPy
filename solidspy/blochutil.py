#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bloch analysis utilities
------------------------
Routines used in Bloch Analisis.

"""
import numpy as np
from scipy.sparse import coo_matrix


#%% Auxiliar
def reassign_dof_bloch(ndofs, dofs_ima, dofs_ref):
    """Find dof number after applying boundary conditions

    Parameters
    ----------
    ndofs : int
        Number of degrees of freedom.
    dofs_ima : list, int
        List with image degrees of freedom.
    dofs_ref : list, int
        List with reference degrees of freedom.


    Examples
    --------
    
    Let us consider an example for  a square made with spring and masses.
    In that case we have the following dofs mappings.
    
    >>> dofs_ref = [0, 1, 0, 1, 0, 1]
    >>> dofs_ima = [2, 3, 4, 5, 6, 7]
    
    And the new numbering would be
    
    >>> new_num = reassign_dof_bloch(8, dofs_ima, dofs_ref)
    >>> new_num
    [0, 1, 0, 1, 0, 1, 0, 1]

    Returns
    -------
    new_num : list, int
        List with new numbering for the dof after applying Bloch analysis.

    """
    new_num = []
    cont = 0
    for cont_dof in range(ndofs):
        if cont_dof in dofs_ima:
            cont += 1
            ref_dof = dofs_ref[dofs_ima.index(cont_dof)]
            if ref_dof < cont_dof:
                new_num.append(ref_dof)
            else:
                new_num.append(ref_dof - cont)
        else:
            new_num.append(cont_dof - cont)
    return new_num


#%% Dense matrices


#%% Sparse matrices
def bloch_transform_mat(wavevector, coords, ndofs, nodes_ref, nodes_ima,
                        dofs_ref, dofs_ima, nodes_ref_dof, nodes_ima_dof,
                        new_num):
    """Form transformation matrices for Bloch conditions

    Parameters
    ----------
    wavevector: array like
        Wavevector, it should have as many components as coordinates.
    coords : array like
        Coordinates array.
    ndofs : int
        Number of degrees of freedom.
    nodes_ref : list, int
        List with reference nodes.
    nodes_ima : list, int
        List with image nodes.
    dofs_ref : list, int
        List with reference degrees of freedom.        
    dofs_ima : list, int
        List with image degrees of freedom.
    nodes_ref_dof : dictionary
        Dictionary that maps nodes to dof for reference nodes.
    nodes_ima_dof : dictionary
        Dictionary that maps nodes to dof for image nodes.


    Examples
    --------
    
    Let us consider an example for  a square made with spring and masses.
    In that case we have the following dofs mappings.
   
    >>> wavevector = np.array([0, 0])
    >>> coords = np.array([
    ...              [0, 0],
    ...              [1, 0],
    ...              [0, 1],
    ...              [1, 1]])
    >>> nodes_ref = [0, 0, 0]
    >>> nodes_ima = [1, 2, 3]
    >>> dofs_ref = [0, 1, 0, 1, 0, 1]
    >>> dofs_ima = [2, 3, 4, 5, 6, 7]
    >>> nodes_ref_dof = {0: [0, 1]}
    >>> nodes_ima_dof = {1: [2, 3], 2: [4, 5], 3: [6, 7]}    
    >>> new_num = [0, 1, 0, 1, 0, 1, 0, 1]
    >>> T = bloch_transform_mat(wavevector, coords, 8, nodes_ref,
    ...                         nodes_ima, dofs_ref, dofs_ima,
    ...                         nodes_ref_dof, nodes_ima_dof, new_num)
    >>> T_exact = np.array([
    ... [1., 0.],
    ... [0., 1.],
    ... [1., 0.],
    ... [0., 1.],
    ... [1., 0.],
    ... [0., 1.],
    ... [1., 0.],
    ... [0., 1.]])
    >>> np.allclose(T.toarray(), T_exact)
    True

    Now, let us consider an example with a mesh with 2x2 elements with
    one dof per node.
        
    >>> wavevector = 0.5*np.pi*np.array([1, 1])
    >>> coords = np.array([
    ...              [0, 0],
    ...              [1, 0],
    ...              [2, 0],
    ...              [0, 1],
    ...              [1, 1],
    ...              [2, 1],
    ...              [0, 2],
    ...              [1, 2],
    ...              [2, 2]])
    >>> nodes_ref = [0, 0, 0, 1, 3]
    >>> nodes_ima = [2, 6, 8, 7, 5]
    >>> dofs_ref = [0, 0, 0, 1, 3]
    >>> dofs_ima = [2, 6, 8, 7, 5]
    >>> nodes_ref_dof = {0: [0], 1: [1], 3: [3]}
    >>> nodes_ima_dof = {2: [2], 5: [5], 6: [6], 7: [7], 8: [8]}    
    >>> new_num = [0, 1, 0, 2, 3, 2, 0, 1, 0]
    >>> T = bloch_transform_mat(wavevector, coords, 9, nodes_ref,
    ...                         nodes_ima, dofs_ref, dofs_ima,
    ...                         nodes_ref_dof, nodes_ima_dof, new_num)
    >>> T_exact = np.array([
    ... [ 1.,  0.,  0.,  0.],
    ... [ 0.,  1.,  0.,  0.],
    ... [-1.,  0.,  0.,  0.],
    ... [ 0.,  0.,  1.,  0.],
    ... [ 0.,  0.,  0.,  1.],
    ... [ 0.,  0., -1.,  0.],
    ... [-1.,  0.,  0.,  0.],
    ... [ 0., -1.,  0.,  0.],
    ... [ 1.,  0.,  0.,  0.]])
    >>> np.allclose(T.toarray(), T_exact)
    True


    Returns
    -------
    Tmat : sparse matrix
        Transformation matrix to apply Bloch boundary conditions.

    """
    rows = []
    cols = []
    vals = []
        
    # Inner nodes
    for cont_dof in range(ndofs):
        if cont_dof not in dofs_ima:
            rows.append(cont_dof)
            cols.append(new_num[cont_dof])
            vals.append(1)
    
    # Bloch nodes
    for node_ima, node_ref in zip(nodes_ima, nodes_ref):
        pos_diff = coords[node_ima] - coords[node_ref]
        phase_shift = np.exp(1.0j*np.dot(wavevector, pos_diff))
        dof_ref = nodes_ref_dof[node_ref]
        dof_ima = nodes_ima_dof[node_ima]
        for cont_ref, cont_ima in zip(dof_ref, dof_ima):
            rows.append(cont_ima)
            cols.append(new_num[cont_ref])
            vals.append(phase_shift)
    
    nconds = len(dofs_ref)
    return coo_matrix((vals, (rows, cols)),
                      shape=(ndofs, ndofs - nconds)).tocsr()


#%% Doctests
if __name__ == "__main__":    
    import doctest
    doctest.testmod()