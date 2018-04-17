"""
Solver routines
---------------

Utilities for solution of FEM systems

"""
from __future__ import absolute_import, division, print_function
from numpy import ndarray
from numpy.linalg import solve
from scipy.sparse.csr import csr_matrix
from scipy.sparse.linalg import spsolve


def static_sol(mat, rhs):
    """Solve a static problem [mat]{u_sol} = {rhs}
    
    Parameters
    ----------
    mat : array
        Array with the system of equations. It can be stored in
        dense or sparse scheme.
    rhs : array
        Array with right-hand-side of the system of equations.
    
    Returns
    -------
    u_sol : array
        Solution of the system of equations.

    Raises
    ------
    
        

    """
    if type(mat) is csr_matrix:
        u_sol = spsolve(mat, rhs)
    elif type(mat) is ndarray:
        u_sol = solve(mat, rhs)
    else:
        raise TypeError("Matrix should be numpy array or csr_matrix.")

    return u_sol
