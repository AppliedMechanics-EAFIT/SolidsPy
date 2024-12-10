# -*- coding: utf-8 -*-
"""
Test for functions in ``blochutil`` module
"""
import numpy as np
from numpy import array, exp, pi
from solidspy.blochutil import bloch_transform_mat

def test_bloch_transform_mat():

    # Square made with springs and mass in the center
    coords = np.array([
            [0, 0],
            [2, 0],
            [2, 2],
            [0, 2],
            [1, 1]])
    ndofs = 10
    nodes_ref = [0, 0, 0]
    nodes_ima = [1, 2, 3]
    dofs_ima = [2, 3, 4, 5, 6, 7]
    dofs_ref = [0, 1, 0, 1, 0, 1]
    nodes_ref_dof = {0: [0, 1]}
    nodes_ima_dof = {1: [2, 3], 2: [4, 5], 3: [6, 7]}
    new_num = [0, 1, 0, 1, 0, 1, 0, 1, 2, 3]
    nk = 51
    wavenumbers = np.random.uniform(-pi/2, pi/2, (nk, 2))
    for wavenumber in wavenumbers:
        kx, ky = wavenumber
        T = bloch_transform_mat((kx, ky), coords, ndofs, nodes_ref, nodes_ima,
                        dofs_ref, dofs_ima, nodes_ref_dof, nodes_ima_dof,
                        new_num)
        TH = array([
                [1, 0, exp(-2j*kx), 0, exp(-2j*(kx + ky)), 0, exp(-2j*ky), 0, 0, 0],
                [0, 1, 0, exp(-2j*kx), 0, exp(-2j*(kx + ky)), 0, exp(-2j*ky), 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])
        assert np.allclose(T.toarray(), TH.T.conj())
