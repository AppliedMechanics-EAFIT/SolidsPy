# -*- coding: utf-8 -*-
"""
Test cases for functions on ``uelutil`` module

"""
from __future__ import division, print_function
import numpy as np
import solidspy.uelutil as uel


def test_uel4nquad():
    """Tests for 4-nodes quad uel"""
    coord = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    stiff = uel.uel4nquad(coord, 1/3, 8/3)
    stiff_ex = 1/6 * np.array([
         [ 8,  3, -5,  0, -4, -3,  1,  0],
         [ 3,  8,  0,  1, -3, -4,  0, -5],
         [-5,  0,  8, -3,  1,  0, -4,  3],
         [ 0,  1, -3,  8,  0, -5,  3, -4],
         [-4, -3,  1,  0,  8,  3, -5,  0],
         [-3, -4,  0, -5,  3,  8,  0,  1],
         [ 1,  0, -4,  3, -5,  0,  8, -3],
         [ 0, -5,  3, -4,  0,  1, -3,  8]])
    assert np.allclose(stiff, stiff_ex)

    
def test_uel6ntrian():
    """Tests for 6-nodes tri uel"""
    coord = np.array([
        [0, 0],
        [1, 0],
        [0, 1],
        [0.5, 0],
        [0.5, 0.5],
        [0, 0.5]])
    stiff = uel.uel6ntrian(coord,1/3, 8/3)
    stiff_ex = 1/6 * np.array([
        [12, 6, 3, 1, 1, 1, -12, -4, 0, 0, -4, -4],
        [6, 12, 1, 1, 1, 3, -4, -4, 0, 0, -4, -12],
        [3, 1, 9, 0, 0, -1, -12, -4, 0, 4, 0, 0],
        [1, 1, 0, 3, -1, 0, -4, -4, 4, 0, 0, 0],
        [1, 1, 0, -1, 3, 0, 0, 0, 0, 4, -4, -4],
        [1, 3, -1, 0, 0, 9, 0, 0, 4, 0, -4, -12],
        [-12, -4, -12, -4, 0, 0, 32, 8, -8, -8, 0, 8],
        [-4, -4, -4, -4, 0, 0, 8, 32, -8, -24, 8, 0],
        [0, 0, 0, 4, 0, 4, -8, -8, 32, 8, -24, -8],
        [0, 0, 4, 0, 4, 0, -8, -24, 8, 32, -8, -8],
        [-4, -4, 0, 0, -4, -4, 0, 8, -24, -8, 32, 8],
        [-4, -12, 0, 0, -4, -12, 8, 0, -8, -8, 8, 32]])
    assert np.allclose(stiff, stiff_ex)


def test_uel3ntrian():
    """Tests for 3-nodes tri uel"""
    coord = np.array([
        [0, 0],
        [1, 0],
        [0, 1]])
    stiff = uel.uel3ntrian(coord, 1/3, 8/3)
    stiff_ex = 1/2 * np.array([
       [4, 2, -3, -1, -1, -1],
       [2, 4, -1, -1, -1, -3],
       [-3, -1, 3, 0, 0, 1],
       [-1, -1, 0, 1, 1, 0],
       [-1, -1, 0, 1, 1, 0],
       [-1, -3, 1, 0, 0, 3]])
    assert np.allclose(stiff, stiff_ex)


def test_uelspring():
    """Tests for 2-nodes springs uel"""
    coord = np.array([
        [0, 0],
        [1, 0]])
    stiff = uel.uelspring(coord, 1/3, 1)
    stiff_ex = np.array([
        [1, 0, -1, 0],
        [0, 0, 0, 0],
        [-1, 0, 1, 0],
        [0, 0, 0, 0]])
    assert np.allclose(stiff, stiff_ex)


def test_ueltruss2D():
    
    coord = np.array([
        [0, 0],
        [1, 0]])
    stiff = uel.ueltruss2D(coord, 1.0 , 1.0)
    stiff_ex =  np.array([
    [1, 0, -1, 0],
    [0, 0, 0, 0],
    [-1, 0, 1, 0],
    [0, 0, 0, 0]])
    np.allclose(stiff, stiff_ex)