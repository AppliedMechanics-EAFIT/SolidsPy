# -*- coding: utf-8 -*-
"""
Test cases for functions on ``uelutil`` module

"""
from __future__ import division, print_function
import numpy as np
import solidspy.gaussutil as gauss


def test_gpoints2x2():
    """Test for 2x2 Gauss integration"""
    wts, pts = gauss.gpoints2x2()
    x, y = pts.T
    fun = lambda x, y: x**2 + y**2
    inte = np.sum(fun(x, y) * wts)
    assert np.isclose(inte, 8/3)

    wts, pts = gauss.gpoints2x2()
    x, y = pts.T
    fun = lambda x, y: 3*x**3 - 4*y**2
    inte = np.sum(fun(x, y) * wts)
    assert np.isclose(inte, -16/3)


def test_gauss_nd():
    """Test for ND Gauss integration"""

    npts = 4
    ndim = 2
    pts, wts = gauss.gauss_nd(npts, ndim=ndim)
    fun = lambda x, y: 3*x*y**2 - x**3
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y = pt
        inte += wts[cont] * fun(x, y)

    assert np.isclose(inte, 0)

    npts = 10
    ndim = 2
    pts, wts = gauss.gauss_nd(npts, ndim=ndim)
    fun = lambda x, y: np.exp(-x**2 - y**2)
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y = pt
        inte += wts[cont] * fun(x, y)

    assert np.isclose(inte, 2.23098514140413)
