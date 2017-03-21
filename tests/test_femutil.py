# -*- coding: utf-8 -*-
"""
Test cases for functions on ``femutil`` module

"""
from __future__ import division, print_function
from os import sys
sys.path.append("../MAIN/")
import numpy as np
import femutil as fem


def test_sha4():

    # For point (0, 0)
    N = fem.sha4(0, 0)
    N_ex = 0.25 * np.array([
        [1, 0, 1, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 1, 0, 1]])
    assert np.allclose(N, N_ex)


def test_sha6():

    # For point (0, 0)
    N = fem.sha6(0, 0)
    N_ex = np.array([
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    assert np.allclose(N, N_ex)


def test_stdm4NQ():

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
