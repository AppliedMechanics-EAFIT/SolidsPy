# -*- coding: utf-8 -*-
"""
gaussutil.py
------------
Weights and coordinates for Gauss-Legendre quadrature [1]_. The
values for triangles is presented in section 5.5 of Bathe book [2]_.

References
----------
.. [1] Wikipedia contributors. "Gaussian quadrature." Wikipedia,
  The Free Encyclopedia, 2 Nov.  2015. Web. 25 Dec. 2015.
  url: https://en.wikipedia.org/wiki/Gaussian_quadrature
.. [2] Bathe, Klaus-Jürgen. Finite element procedures. Prentice Hall,
   Pearson Education, 2006.
"""
from __future__ import division, print_function
import numpy as np


def gpoints2x2():
    """Gauss points for a 2 by 2 grid

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([4])
    xp = np.zeros([4, 2])
    xw[:] = 1.0
    xp[0, 0] = -0.577350269189626
    xp[1, 0] = 0.577350269189626
    xp[2, 0] = -0.577350269189626
    xp[3, 0] = 0.577350269189626

    xp[0, 1] = 0.577350269189626
    xp[1, 1] = 0.577350269189626
    xp[2, 1] = -0.577350269189626
    xp[3, 1] = -0.577350269189626

    return xw, xp


def gpoints7():
    """Gauss points for a triangle (7 points)

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([7])
    xp = np.zeros([7, 2])
    xw[0] = 0.1259391805448
    xw[1] = 0.1259391805448
    xw[2] = 0.1259391805448
    xw[3] = 0.1323941527885
    xw[4] = 0.1323941527885
    xw[5] = 0.1323941527885
    xw[6] = 0.225

    xp[0, 0] = 0.1012865073235
    xp[1, 0] = 0.7974269853531
    xp[2, 0] = 0.1012865073235
    xp[3, 0] = 0.4701420641051
    xp[4, 0] = 0.4701420641051
    xp[5, 0] = 0.0597158717898
    xp[6, 0] = 0.3333333333333

    xp[0, 1] = 0.1012865073235
    xp[1, 1] = 0.1012865073235
    xp[2, 1] = 0.7974269853531
    xp[3, 1] = 0.0597158717898
    xp[4, 1] = 0.4701420641051
    xp[5, 1] = 0.4701420641051
    xp[6, 1] = 0.3333333333333

    return xw, xp


def gpoints3():
    """Gauss points for a triangle element (3 points)

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([3])
    xp = np.zeros([3, 2])
    xw[0] = 0.3333333333333
    xw[1] = 0.3333333333333
    xw[2] = 0.3333333333333

    xp[0, 0] = 0.1666666666667
    xp[1, 0] = 0.6666666666667
    xp[2, 0] = 0.1666666666667

    xp[0, 1] = 0.1666666666667
    xp[1, 1] = 0.1666666666667
    xp[2, 1] = 0.6666666666667

    return xw, xp
