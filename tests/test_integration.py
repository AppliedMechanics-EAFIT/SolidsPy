# -*- coding: utf-8 -*-
"""
Integration tests for solidspy

"""
import numpy as np
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol


def test_4_elements():
    """2×2 mesh with uniaxial load"""
    nodes = np.array([
            [0, 0, 0, 0, -1],
            [1, 2, 0, 0, -1],
            [2, 2, 2, 0, 0],
            [3, 0, 2, 0, 0],
            [4, 1, 0, -1, -1],
            [5, 2, 1, 0, 0],
            [6, 1, 2, 0, 0],
            [7, 0, 1, 0, 0],
            [8, 1, 1, 0, 0]])
    eles = np.array([
            [0, 1, 0, 0, 4, 8, 7],
            [1, 1, 0, 4, 1, 5, 8],
            [2, 1, 0, 7, 8, 6, 3],
            [3, 1, 0, 8, 5, 2, 6]])
    loads = np.array([
            [3, 0, 1],
            [6, 0, 2],
            [2, 0, 1]])
    mater = np.array([[1.0, 0.3]])
    DME, IBC, neq = ass.DME(nodes, eles)
    KG = ass.assembler(eles, mater, nodes, neq, DME)
    RHSG = ass.loadasem(loads, IBC, neq)
    disp = sol.static_sol(KG, RHSG)
    disp_complete = pos.complete_disp(IBC, nodes, disp)
    disp_anal = np.array([
            [ 0.6, 0.0],
            [-0.6, 0.0],
            [-0.6, 4.0],
            [0.6, 4.0],
            [0.0, 0.0],
            [-0.6, 2.0],
            [0.0, 4.0],
            [0.6, 2.0],
            [0.0, 2.0]])
    assert np.allclose(disp_complete, disp_anal)

    
    
