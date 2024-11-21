# -*- coding: utf-8 -*-
"""
Test cases for functions on ``gaussutil`` module

"""
import numpy as np
import solidspy.gaussutil as gauss


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


def test_gauss_tet():
    # Zero order
    pts, wts = gauss.gauss_tet(order=1)
    fun = lambda x, y, z: 1
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y, z = pt
        inte += wts[cont] * fun(x, y, z)/6
    assert np.isclose(inte, 1/6)

    # First order
    pts, wts = gauss.gauss_tet(order=1)
    a, b, c = np.random.uniform(0, 1, 3)
    fun = lambda x, y, z : a*x + b*y + c*z
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y, z = pt
        inte += wts[cont] * fun(x, y, z)/6
    assert np.isclose(inte, (a + b + c)/24)

    # Second order
    pts, wts = gauss.gauss_tet(order=2)
    a, b, c, d, e, f = np.random.uniform(0, 1, 6)
    fun = lambda x, y, z: a*x**2 + b*y**2 + c*z**2 + d*x*y + e*y*z + f*z*x
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y, z = pt
        inte += wts[cont] * fun(x, y, z) / 6
    assert np.isclose(inte, (2*a + 2*b + 2*c + d + e + f)/120)

    # Third order
    pts, wts = gauss.gauss_tet(order=3)
    a, b, c, d = np.random.uniform(0, 1, 4)
    fun = lambda x, y, z: a*x**3 + b*y**3 + c*z**3 + d*x*y*z
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y, z = pt
        inte += wts[cont] * fun(x, y, z) / 6
    assert np.isclose(inte, (6*a + 6*b + 6*c + d)/720)

    # Seventh order
    pts, wts = gauss.gauss_tet(order=7)
    a, b, c, d, e, f, g, h = np.random.uniform(0, 1, 8)
    fun = lambda x, y, z: a*x**7 + b*y**7 + c*z**7 + d*x**6 + e*y**6\
          + f*z**6 + g*x*y*z + h  
    inte = 0.0
    for cont, pt in enumerate(pts):
        x, y, z= pt
        inte += wts[cont] * fun(x, y, z) / 6
    assert np.isclose(inte, (7*a + 7*b + 7*c + 10*d + 10*e + 10*f + 7*g\
                            + 840*h)/5040)