# -*- coding: utf-8 -*-
"""
Test cases for functions on ``femutil`` module

"""
import numpy as np
import pytest
import solidspy.femutil as fem


#%% Tests for Shape functions and derivatives
def test_shape_tri3():
    # Interpolation condition check
    coords = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0]])
    N, _ = fem.shape_tri3(coords[:, 0], coords[:, 1])
    assert np.allclose(N, np.eye(3))


def test_shape_tri6():
    # Interpolation condition check
    coords = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.0],
        [0.5, 0.5],
        [0.0, 0.5]])
    N, _ = fem.shape_tri6(coords[:, 0], coords[:, 1])
    assert np.allclose(N, np.eye(6))

    # Evaluation at (1/3, 1/3)
    N, dNdr = fem.shape_tri6(1/3, 1/3)
    N_exp = np.array([-1., -1., -1., 4., 4., 4.])/9
    dNdr_exp = np.array([
                [-1.,  1.,  0.,  0.,  4., -4.],
                [-1.,  0.,  1., -4.,  4.,  0.]])/3
    assert np.allclose(N, N_exp)
    assert np.allclose(dNdr, dNdr_exp)


def test_shape_quad4():
    # Interpolation condition check
    coords = np.array([
        [-1.0, -1.0],
        [1.0, -1.0],
        [1.0, 1.0],
        [-1.0, 1.0]])
    N, _ = fem.shape_quad4(coords[:, 0], coords[:, 1])
    assert np.allclose(N, np.eye(4))

    # For point (0, 0)
    N, _ = fem.shape_quad4(0, 0)
    N_ex = 0.25 * np.array([[1, 1, 1, 1]])
    assert np.allclose(N, N_ex)


def test_shape_quad9():
    # Interpolation condition check
    coords = np.array([
        [-1.0, -1.0],
        [ 1.0, -1.0],
        [ 1.0,  1.0],
        [-1.0,  1.0],
        [ 0.0, -1.0],
        [ 1.0,  0.0],
        [ 0.0,  1.0],
        [-1.0,  0.0],
        [ 0.0,  0.0]])
    N, _ = fem.shape_quad9(coords[:, 0], coords[:, 1])
    assert np.allclose(N, np.eye(9))

    # Evaluation at (1/4, 1/4)
    N, dNdr = fem.shape_quad9(0.25, 0.25)
    N_exp = np.array(
        [0.00878906, -0.01464844, 0.02441406, -0.01464844,
        -0.08789062, 0.14648438, 0.14648438, -0.08789062,
        0.87890625])

    dNdr_exp = np.array([
        [0.0234375, -0.0703125, 0.1171875, -0.0390625, 0.046875,
            0.703125, -0.078125, -0.234375, -0.46875],
        [0.0234375, -0.0390625, 0.1171875, -0.0703125, -0.234375,
            -0.078125, 0.703125, 0.046875, -0.46875]])
    assert np.allclose(N, N_exp)
    assert np.allclose(dNdr, dNdr_exp)


#%% Jacobian
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

    # Wrong triangles

    # Repeated node
    dhdx = np.array([
            [-1, 1, 0],
            [-1, 0, 1]])    
    coord = np.array([
            [0, 0],
            [0, 0],
            [0, 1]])

    with pytest.raises(ValueError):
        det, jaco_inv = fem.jacoper(dhdx, coord)

    # Opposite orientation
    coord = np.array([
            [0, 0],
            [0, 1],
            [1, 0]])

    with pytest.raises(ValueError):
        det, jaco_inv = fem.jacoper(dhdx, coord)


    # Wrong quads

    # Opposite orientation
    dhdx = 0.25*np.array([
            [-1, 1, 1, -1],
            [-1, -1, 1, 1]])
    coord = np.array([
            [-1, 1],
            [1, 1],
            [1, -1],
            [-1, -1]])
    with pytest.raises(ValueError):
        det, jaco_inv = fem.jacoper(dhdx, coord)

    # Repeated nodes
    dhdx = 0.25*np.array([
            [-1, 1, 1, -1],
            [-1, -1, 1, 1]])
    coord = np.array([
            [1, -1],
            [1, -1],
            [1, -1],
            [-1, 1]])
    with pytest.raises(ValueError):
        det, jaco_inv = fem.jacoper(dhdx, coord)
