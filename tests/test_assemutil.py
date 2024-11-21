# -*- coding: utf-8 -*-
"""
Test cases for functions on ``assemutil`` module

"""
import numpy as np
import src.assemutil as ass


def test_sparse_assem():
    """Tests for sparse assembler"""
    # 2 x 2 mesh
    mats = np.array([[16, 1/3]])
    nodes = np.array([
        [0, -1, -1],
        [1, 0, -1],
        [2, 1, -1],
        [3, -1, 0],
        [4, 0, 0],
        [5, 1, 0],
        [6, -1, 1],
        [7, 0, 1],
        [8, 1, 1]])
    elements = np.array([
        [0, 1, 0, 0, 1, 4, 3],
        [1, 1, 0, 1, 2, 5, 4],
        [2, 1, 0, 3, 4, 7, 6],
        [3, 1, 0, 4, 5, 8, 7]])
    assem_op = np.array([
        [0, 1, 2, 3, 8, 9, 6, 7],
        [2, 3, 4, 5, 10, 11, 8, 9],
        [6, 7, 8, 9, 14, 15, 12, 13],
        [8, 9, 10, 11, 16, 17, 14, 15]], dtype=int)

    neq = 18
    K_ass, _ = ass.sparse_assem(elements, mats, nodes, neq, assem_op)
    K_exact = np.array([
        [8, 3, -5, 0, 0, 0, 1, 0, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0],
        [3, 8, 0, 1, 0, 0, 0, -5, -3, -4, 0, 0, 0, 0, 0, 0, 0, 0],
        [-5, 0, 16, 0, -5, 0, -4, 3, 2, 0, -4, -3, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 16, 0, 1, 3, -4, 0, -10, -3, -4, 0, 0, 0, 0, 0, 0],
        [0, 0, -5, 0, 8, -3, 0, 0, -4, 3, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, -3, 8, 0, 0, 3, -4, 0, -5, 0, 0, 0, 0, 0, 0],
        [1, 0, -4, 3, 0, 0, 16, 0, -10, 0, 0, 0, 1, 0, -4, -3, 0, 0],
        [0, -5, 3, -4, 0, 0, 0, 16, 0, 2, 0, 0, 0, -5, -3, -4, 0, 0],
        [-4, -3, 2, 0, -4, 3, -10, 0, 32, 0, -10, 0, -4, 3, 2, 0, -4, -3],
        [-3, -4, 0, -10, 3, -4, 0, 2, 0, 32, 0, 2, 3, -4, 0, -10, -3, -4],
        [0, 0, -4, -3, 1, 0, 0, 0, -10, 0, 16, 0, 0, 0, -4, 3, 1, 0],
        [0, 0, -3, -4, 0, -5, 0, 0, 0, 2, 0, 16, 0, 0, 3, -4, 0, -5],
        [0, 0, 0, 0, 0, 0, 1, 0, -4, 3, 0, 0, 8, -3, -5, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, -5, 3, -4, 0, 0, -3, 8, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, -4, -3, 2, 0, -4, 3, -5, 0, 16, 0, -5, 0],
        [0, 0, 0, 0, 0, 0, -3, -4, 0, -10, 3, -4, 0, 1, 0, 16, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, -4, -3, 1, 0, 0, 0, -5, 0, 8, 3],
        [0, 0, 0, 0, 0, 0, 0, 0, -3, -4, 0, -5, 0, 0, 0, 1, 3, 8]])
    assert np.allclose(K_ass.toarray(), K_exact)


def test_dense_assem():
    """Tests for dense assembler"""

    # 2 x 2 mesh
    mats = np.array([[16, 1/3]])
    nodes = np.array([
        [0, -1, -1],
        [1, 0, -1],
        [2, 1, -1],
        [3, -1, 0],
        [4, 0, 0],
        [5, 1, 0],
        [6, -1, 1],
        [7, 0, 1],
        [8, 1, 1]])
    elements = np.array([
        [0, 1, 0, 0, 1, 4, 3],
        [1, 1, 0, 1, 2, 5, 4],
        [2, 1, 0, 3, 4, 7, 6],
        [3, 1, 0, 4, 5, 8, 7]])
    assem_op = np.array([
        [0, 1, 2, 3, 8, 9, 6, 7],
        [2, 3, 4, 5, 10, 11, 8, 9],
        [6, 7, 8, 9, 14, 15, 12, 13],
        [8, 9, 10, 11, 16, 17, 14, 15]], dtype=int)

    neq = 18
    K_ass, _ = ass.dense_assem(elements, mats, nodes, neq, assem_op)
    K_exact = np.array([
        [8, 3, -5, 0, 0, 0, 1, 0, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0],
        [3, 8, 0, 1, 0, 0, 0, -5, -3, -4, 0, 0, 0, 0, 0, 0, 0, 0],
        [-5, 0, 16, 0, -5, 0, -4, 3, 2, 0, -4, -3, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 16, 0, 1, 3, -4, 0, -10, -3, -4, 0, 0, 0, 0, 0, 0],
        [0, 0, -5, 0, 8, -3, 0, 0, -4, 3, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, -3, 8, 0, 0, 3, -4, 0, -5, 0, 0, 0, 0, 0, 0],
        [1, 0, -4, 3, 0, 0, 16, 0, -10, 0, 0, 0, 1, 0, -4, -3, 0, 0],
        [0, -5, 3, -4, 0, 0, 0, 16, 0, 2, 0, 0, 0, -5, -3, -4, 0, 0],
        [-4, -3, 2, 0, -4, 3, -10, 0, 32, 0, -10, 0, -4, 3, 2, 0, -4, -3],
        [-3, -4, 0, -10, 3, -4, 0, 2, 0, 32, 0, 2, 3, -4, 0, -10, -3, -4],
        [0, 0, -4, -3, 1, 0, 0, 0, -10, 0, 16, 0, 0, 0, -4, 3, 1, 0],
        [0, 0, -3, -4, 0, -5, 0, 0, 0, 2, 0, 16, 0, 0, 3, -4, 0, -5],
        [0, 0, 0, 0, 0, 0, 1, 0, -4, 3, 0, 0, 8, -3, -5, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, -5, 3, -4, 0, 0, -3, 8, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, -4, -3, 2, 0, -4, 3, -5, 0, 16, 0, -5, 0],
        [0, 0, 0, 0, 0, 0, -3, -4, 0, -10, 3, -4, 0, 1, 0, 16, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, -4, -3, 1, 0, 0, 0, -5, 0, 8, 3],
        [0, 0, 0, 0, 0, 0, 0, 0, -3, -4, 0, -5, 0, 0, 0, 1, 3, 8]])
    assert np.allclose(K_ass, K_exact)


    # Test for uel with all ones
    def uel_ones(elcoord, params):
        """Dummy UEL with all ones"""
        return np.ones((8, 8)), np.ones((8, 8))

    nodes = np.zeros((9, 3))
    nodes[:, 0] = range(0, 9)
    nodes[:, 1:3] = np.array([
        [0, 0],
        [1, 0],
        [2, 0],
        [0, 1],
        [1, 1],
        [2, 1],
        [0, 2],
        [1, 2],
        [2, 2]])
    cons = np.zeros((9, 2))
    elements = np.ones((4, 7), dtype=int)
    elements[:, 0] = range(0, 4)
    elements[:, 2] = 0
    elements[:, 3:] = np.array([
        [0, 1, 4, 3],
        [1, 2, 5, 4],
        [3, 4, 7, 6],
        [4, 5, 8, 7]])
    mats = np.array([[1, 0.3]])
    assemp_op, bc_array, neq = ass.DME(cons, elements)
    stiff, _ = ass.assembler(elements, mats, nodes, neq, assem_op,
                             sparse=False, uel=uel_ones)
    stiff_exact = np.array([
        [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0],
        [1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 0, 0, 1, 1, 1, 1, 0, 0],
        [1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 0, 0, 1, 1, 1, 1, 0, 0],
        [1, 1, 2, 2, 1, 1, 2, 2, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1],
        [1, 1, 2, 2, 1, 1, 2, 2, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1],
        [0, 0, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 0, 0, 1, 1, 1, 1],
        [0, 0, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1],
        [0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1]])
    assert np.allclose(stiff, stiff_exact)


    # Test assembly of a truss
    length = 10*np.cos(np.pi/6)
    nodes = np.array([
        [0.0, length, 0.0],
        [1.0, 0.0, 5.0],
        [2.0, 0.0, 0.0]])
    mats = np.array([[1e6, 0.01]])
    elements = np.array([
        [0, 6, 0, 2, 0],
        [1, 6, 0, 1, 0],
        [2, 6, 0, 1, 2]])
    neq = 3
    assem_op = np.array([
            [-1, -1,  0,  1],
            [ 0,  1, -1,  2],
            [-1,  2, -1, -1]])
    stiff, _ = ass.assembler(elements, mats, nodes, neq, assem_op,
                             sparse=False)
    stiff_exact = 250*np.array([
        [8/np.sqrt(3) + 3, -np.sqrt(3), np.sqrt(3)],
        [-np.sqrt(3), 1, -1],
        [np.sqrt(3), -1, 9]])
    assert np.allclose(stiff, stiff_exact)


