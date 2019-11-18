# -*- coding: utf-8 -*-
"""
Test cases for functions on ``gaussutil`` module

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


def test_gauss_tri():
    # Zero order
    pts, wts = gauss.gauss_tri(order=1)
    fun = lambda x, y: 1
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y = pt
        inte += 0.5 * wts[cont] * fun(x, y)
    assert np.isclose(inte, 0.5)

    # First order
    pts, wts = gauss.gauss_tri(order=1)
    a, b = np.random.uniform(0, 1, 2)
    fun = lambda x, y: a*x + b*y
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y = pt
        inte += 0.5 * wts[cont] * fun(x, y)
    assert np.isclose(inte, (a + b)/6)

    # Second order
    pts, wts = gauss.gauss_tri(order=2)
    a, b, c, d, e, f = np.random.uniform(0, 1, 6)
    fun = lambda x, y: a*x**2 + b*y**2 + c*x*y + d*x + e*y + f
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y = pt
        inte += 0.5 * wts[cont] * fun(x, y)
    assert np.isclose(inte, (2*a + 2*b + c + 4*d + 4*e + 12*f)/24)

    # Third order
    pts, wts = gauss.gauss_tri(order=3)
    a, b, c, d, e, f = np.random.uniform(0, 1, 6)
    fun = lambda x, y: a*x**2 + b*y**2 + c*x*y + d*x + e*y + f
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y = pt
        inte += 0.5 * wts[cont] * fun(x, y)
    assert np.isclose(inte, (2*a + 2*b + c + 4*d + 4*e + 12*f)/24)

    # Seventh order
    pts, wts = gauss.gauss_tri(order=7)
    a, b, c, d, e, f, g, h = np.random.uniform(0, 1, 8)
    fun = lambda x, y: a*x**7 + b*x**6*y + c*x**5*y**2 + d*x**4*y**3\
                       + e*x**3*y**4 + f*x**2*y**5 + g*x*y**6 + h*y**7   
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y = pt
        inte += 0.5 * wts[cont] * fun(x, y)
    assert np.isclose(inte, (105*a + 15*b + 5*c + 3*d + 3*e + 5*f + 15*g
                             + 105*h)/7560)