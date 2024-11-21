# -*- coding: utf-8 -*-
"""
Integration tests for solidspy

"""
import numpy as np
from scipy.sparse.linalg import eigsh
import src.postprocesor as pos
import src.assemutil as ass
import src.solutil as sol


def test_4_elements():
    """2×2 mesh with uniaxial load"""
    nodes = np.array([
            [0, 0, 0],
            [1, 2, 0],
            [2, 2, 2],
            [3, 0, 2],
            [4, 1, 0],
            [5, 2, 1],
            [6, 1, 2],
            [7, 0, 1],
            [8, 1, 1]])
    cons = np.array([
            [0, -1],
            [0, -1],
            [0, 0],
            [0, 0],
            [-1, -1],
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0]])
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
    assem_op, bc_array, neq = ass.DME(cons, eles)
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
            [0, 0, 0],
            [1, 1, 0],
            [2, 2, 0],
            [3, 0, 1],
            [4, 1, 1],
            [5, 2, 1]])
    cons = np.array([
            [-1, -1],
            [0, 0],
            [0, 0],
            [-1, -1],
            [0, 0],
            [0, 0]])
    eles = np.array([
            [0, 1, 0, 0, 1, 4, 3],
            [1, 1, 0, 1, 2, 5, 4]])
    loads = np.array([
            [2, 0, -0.5],
            [5, 0, -0.5]])
    mater = np.array([[1.0, 0.3]])
    assem_op, bc_array, neq = ass.DME(cons, eles)
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


def test_beams():
    """Beams with axial force"""
    
    # Analytic problem
    nodes = np.array([
        [0, 0.0, 0.0],
        [1, 0.0, 6.0],
        [2, 4.0, 6.0]])
    cons = np.array([
            [-1, -1, -1],
            [0, 0, 0],
            [-1, -1, -1]])
    mats = np.array([[200e9, 1.33e-4, 0.04]])
    elements = np.array([
        [0, 8, 0, 0, 1],
        [1, 8, 0, 1, 2]])
    loads = np.array([
        [1, -12000, -24000, -6000]])
    assem_op, bc_array, neq = ass.DME(cons, elements, ndof_node=3)
    stiff, _ = ass.assembler(elements, mats, nodes, neq, assem_op,
                             sparse=False)
    load_vec = ass.loadasem(loads, bc_array, neq, ndof_node=3)
    solution = sol.static_sol(stiff, load_vec)
    solution_analytic = np.array([-6.29e-6, -1.695e-5, -0.13e-3])
    assert np.allclose(solution, solution_analytic, rtol=1e-1)


def test_eigs_truss():
    """Eigenvalues of a bar"""
    nnodes = 513
    
    x = np.linspace(0, np.pi, nnodes)
    nodes = np.zeros((nnodes, 3))
    nodes[:, 0] = range(nnodes)
    nodes[:, 1] = x
    cons = np.zeros((nnodes, 2))
    cons[:, 1] = -1
    cons[0, 0] = -1 
    cons[-1, 0] = -1
    mats = np.array([[1.0, 1.0, 1.0]])
    elements = np.zeros((nnodes - 1, 5 ), dtype=int)
    elements[:, 0] = range(nnodes - 1)
    elements[:, 1] = 6
    elements[:, 3] = range(nnodes - 1)
    elements[:, 4] = range(1, nnodes)
    
    assem_op, bc_array, neq = ass.DME(cons, elements)
    stiff, mass = ass.assembler(elements, mats, nodes, neq, assem_op)
    
    vals, _ = eigsh(stiff, M=mass, which="SM")
    assert np.allclose(vals, np.linspace(1, 6, 6)**2, rtol=1e-2)


def test_eigs_beam():
    """Eigenvalues of a cantilever beam"""
    
    nnodes = 10

    x = np.linspace(0, np.pi, nnodes)
    nodes = np.zeros((nnodes, 3))
    nodes[:, 0] = range(nnodes)
    nodes[:, 1] = x
    cons = np.zeros((nnodes, 3))
    cons[0, :] = -1
    cons[:, 0] = -1
    mats = np.array([[1.0, 1.0, 1.0, 1.0]])
    elements = np.zeros((nnodes - 1, 5 ), dtype=int)
    elements[:, 0] = range(nnodes - 1)
    elements[:, 1] = 7
    elements[:, 3] = range(nnodes - 1)
    elements[:, 4] = range(1, nnodes)
    
    assem_op, bc_array, neq = ass.DME(cons, elements, ndof_node=3)
    stiff, mass = ass.assembler(elements, mats, nodes, neq, assem_op)
    
    vals, _ = eigsh(stiff, M=mass, which="SM")
    vals_analytic = np.array([0.596864162694467, 1.49417561427335,
                              2.50024694616670,  3.49998931984744,
                              4.50000046151508, 5.49999998005609])
    assert np.allclose(vals**0.25, vals_analytic, rtol=1e-2)
    
