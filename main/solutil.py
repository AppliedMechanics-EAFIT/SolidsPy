"""
solutil.py
----------

Utilities for solution of FEM systems

"""
from __future__ import division, print_function
from numpy import ndarray
from numpy.linalg import solve
from scipy.sparse.csr import csr_matrix
from scipy.sparse.linalg import spsolve


def static_sol(mat, rhs):
    """Solve a static problem [mat]{u_sol} = {rhs}
    """
    if type(mat) is csr_matrix:
        u_sol = spsolve(mat, rhs)
    elif type(mat) is ndarray:
        u_sol = solve(mat, rhs)
    else:
        raise Exception("Not supported matrix storage scheme!")

    return u_sol