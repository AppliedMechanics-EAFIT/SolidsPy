# -*- coding: utf-8 -*-
"""
Test cases for functions on ``preprocesor`` module

"""
import numpy as np
import solidspy.preprocesor as pre


def test_rect_grid():
    """Tests for structured meshes generation"""

    # 2 x 2 mesh
    x, y, els = pre.rect_grid(2, 2, 2, 2)
    assert np.allclose(x, np.array([-1., 0., 1.,
                                    -1., 0., 1.,
                                    -1., 0., 1.]))
    assert np.allclose(y, np.array([-1., -1., -1.,
                                     0., 0., 0.,
                                     1., 1., 1.]))
    assert np.allclose(els, np.array([
            [0, 1, 0, 0, 1, 4, 3],
            [1, 1, 0, 1, 2, 5, 4],
            [2, 1, 0, 3, 4, 7, 6],
            [3, 1, 0, 4, 5, 8, 7]]))
