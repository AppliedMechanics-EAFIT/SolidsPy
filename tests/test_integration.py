# -*- coding: utf-8 -*-
"""
Integration tests for solidspy

"""
import numpy as np
from scipy.sparse.linalg import eigsh
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
    assem_op, bc_array, neq = ass.DME(nodes[:, -2:], eles)
    stiff, _ = ass.assembler(eles, mater, nodes, neq, assem_op)
    load_vec = ass.loadasem(loads, bc_array, neq)
    disp = sol.static_sol(stiff, load_vec)
    disp_complete = pos.complete_disp(bc_array, nodes, disp)
    disp_analytic = np.array([
            [ 0.6, 0.0],
            [-0.6, 0.0],
            [-0.6, 4.0],
            [0.6, 4.0],
            [0.0, 0.0],
            [-0.6, 2.0],
            [0.0, 4.0],
            [0.6, 2.0],
            [0.0, 2.0]])
    assert np.allclose(disp_complete, disp_analytic)


def test_2_elements():
    """2x1 mesh cantilever beam"""
    nodes = np.array([
            [0, 0, 0, -1, -1],
            [1, 1, 0, 0, 0],
            [2, 2, 0, 0, 0],
            [3, 0, 1, -1, -1],
            [4, 1, 1, 0, 0],
            [5, 2, 1, 0, 0]])
    eles = np.array([
            [0, 1, 0, 0, 1, 4, 3],
            [1, 1, 0, 1, 2, 5, 4]])
    loads = np.array([
            [2, 0, -0.5],
            [5, 0, -0.5]])
    mater = np.array([[1.0, 0.3]])
    assem_op, bc_array, neq = ass.DME(nodes[:, -2:], eles)
    stiff, _ = ass.assembler(eles, mater, nodes, neq, assem_op)
    load_vec = ass.loadasem(loads, bc_array, neq)
    disp = sol.static_sol(stiff, load_vec)
    disp_complete = pos.complete_disp(bc_array, nodes, disp)
    disp_analytic = 1/45 * np.array([
            [0, 0],
            [-273, -390],
            [-364, -1144],
            [0, 0],
            [273, -390],
            [364, -1144]])
    
    assert np.allclose(disp_complete, disp_analytic)
    

def test_eigs_truss():
    """Eigenvalues of a bar"""
    nnodes = 513
    
    x = np.linspace(0, np.pi, nnodes)
    nodes = np.zeros((nnodes, 5))
    nodes[:, 0] = range(nnodes)
    nodes[:, 1] = x
    nodes[:, -1] = -1
    nodes[0, -2] = -1 
    nodes[-2, -2] = -1
    mats = np.array([[1.0, 1.0, 1.0]])
    elements = np.zeros((nnodes - 1, 5 ), dtype=int)
    elements[:, 0] = range(nnodes - 1)
    elements[:, 1] = 6
    elements[:, 3] = range(nnodes - 1)
    elements[:, 4] = range(1, nnodes)
    
    assem_op, bc_array, neq = ass.DME(nodes[:, 3:], elements)
    stiff, mass = ass.assembler(elements, mats, nodes, neq, assem_op)
    
    vals, _ = eigsh(stiff, M=mass, which="SM")
    assert np.allclose(vals, np.linspace(1, 6, 6)**2, rtol=1e-2)


def test_eigs_beam():
    """Eigenvalues of a cantilever beam"""
    
    nnodes = 10

    x = np.linspace(0, np.pi, nnodes)
    nodes = np.zeros((nnodes, 6))
    nodes[:, 0] = range(nnodes)
    nodes[:, 1] = x
    nodes[0, 3:] = -1
    nodes[:, -3] = -1
    mats = np.array([[1.0, 1.0, 1.0, 1.0]])
    elements = np.zeros((nnodes - 1, 5 ), dtype=int)
    elements[:, 0] = range(nnodes - 1)
    elements[:, 1] = 7
    elements[:, 3] = range(nnodes - 1)
    elements[:, 4] = range(1, nnodes)
    
    assem_op, bc_array, neq = ass.DME(nodes[:, 3:], elements, ndof_node=3)
    stiff, mass = ass.assembler(elements, mats, nodes, neq, assem_op, sparse=False)
    
    vals, _ = eigsh(stiff, M=mass, which="SM")
    vals_analytic = np.array([0.596864162694467, 1.49417561427335,
                              2.50024694616670,  3.49998931984744,
                              4.50000046151508, 5.49999998005609])

    assert np.allclose(vals**0.25, vals_analytic, rtol=1e-2)
    
