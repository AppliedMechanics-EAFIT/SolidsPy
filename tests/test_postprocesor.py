# -*- coding: utf-8 -*-
"""
Test cases for functions on ``postprocesor`` module

"""
from __future__ import division, print_function
import numpy as np
import solidspy.postprocesor as pos


def test_strain_nodes():
    """Tests for strain/stress calculation at nodes"""

    # 2 x 2 mesh with axial load and rollers on the sides
    mats = np.array([[8/3, 1/3]])
    nodes = np.array([
            [0, -1, -1],
            [1, 0, -1],
            [2, 1, -1],
            [3, -1,  0],
            [4, 0,  0],
            [5, 1,  0],
            [6,-1,  1],
            [7, 0,  1],
            [8, 1,  1]])
    elements =np.array([
            [0, 1, 0, 0, 1, 4, 3],
            [1, 1, 0, 1, 2, 5, 4],
            [2, 1, 0, 3, 4, 7, 6],
            [3, 1, 0, 4, 5, 8, 7]])
    UC = np.array([
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 2],
            [0, 2],
            [0, 2]])
    E_nodes, S_nodes = pos.strain_nodes(nodes , elements, mats, UC)
    E_exact = np.zeros((9, 3))
    E_exact[:, 1] = 1
    S_exact = np.zeros((9, 3))
    S_exact[:, 0] = 1
    S_exact[:, 1] = 3
    assert np.allclose(E_exact, E_nodes)
    assert np.allclose(S_exact, S_nodes)
