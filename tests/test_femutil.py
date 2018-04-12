# -*- coding: utf-8 -*-
"""
Test cases for functions on ``femutil`` module

"""
from __future__ import division, print_function
import numpy as np
import solidspy.femutil as fem


#%% Tests for Shape functions and derivatives
def test_sha4():
    """Tests for 4-nodes quad shape functions"""

    # For point (0, 0)
    N = fem.sha4(0, 0)
    N_ex = 0.25 * np.array([
        [1, 0, 1, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 1, 0, 1]])
    assert np.allclose(N, N_ex)


def test_sha6():
    """Tests for 6-nodes tri shape functions"""

    # For point (0, 0)
    N = fem.sha6(0, 0)
    N_ex = np.array([
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    assert np.allclose(N, N_ex)


def test_stdm4NQ():
    """Tests for 4-nodes quad shape functions derivatives"""

    # For point (1, 1)
    coord = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    det, B = fem.stdm4NQ(1, 1, coord)
    B_ex = 0.5 * np.array([
        [0, 0, 0, 0, 1, 0, -1, 0],
        [0, 0, 0, -1, 0, 1, 0, 0],
        [0, 0, -1, 0, 1, 1, 0, -1]])
    assert np.isclose(det, 1)
    assert np.allclose(B, B_ex)


def test_stdm6NT():
    """Tests for 6-nodes tri shape functions derivatives"""

    # For point (1, 0)
    coord = np.array([
          [0, 0],
          [1, 0],
          [0, 1],
          [0.5, 0],
          [0.5, 0.5],
          [0, 0.5]])
    det, B = fem.stdm6NT(1, 0, coord)
    B_ex = np.array([
        [1, 0, 3, 0, 0, 0, -4, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, -1, 0, -4, 0, 4, 0, 0],
        [1, 1, 0, 3, -1, 0, -4, -4, 4, 0, 0, 0]])
    assert np.isclose(det, 1)
    assert np.allclose(B, B_ex)


def test_stdm3NT():
    """Tests for 3-nodes tri shape functions derivatives"""

    # For point (0, 1)
    coord = np.array([
          [0, 0],
          [1, 0],
          [0, 1]])
    det, B = fem.stdm3NT(0, 1, coord)
    B_ex = np.array([
        [-1, 0, 1, 0, 0, 0],
        [0, -1, 0, 0, 0, 1],
        [-1, -1, 0, 1, 1, 0]])
    assert np.isclose(det, 1)
    assert np.allclose(B, B_ex)


def test_jacoper():
    """Tests for jacobian of the elemental transformation"""

    # Perfect element at (0, 0)
    dhdx = 0.25*np.array([
            [-1, 1, 1, -1],
            [-1, -1, 1, 1]])
    coord = np.array([
        [-1, -1],
        [1, -1],
        [1, 1],
        [-1, 1]])
    det, jaco_inv = fem.jacoper(dhdx, coord)
    jaco_inv_ex = np.eye(2)
    assert np.isclose(det, 1)
    assert np.allclose(jaco_inv, jaco_inv_ex)

    # Shear element at (0, 0)
    dhdx = 0.25*np.array([
            [-1, 1, 1, -1],
            [-1, -1, 1, 1]])
    coord = np.array([
        [-1.5, -1],
        [0.5, -1],
        [1.5, 1],
        [-0.5, 1]])
    det, jaco_inv = fem.jacoper(dhdx, coord)
    jaco_inv_ex = np.eye(2)
    jaco_inv_ex[1, 0] = -0.5
    assert np.isclose(det, 1)
    assert np.allclose(jaco_inv, jaco_inv_ex)
