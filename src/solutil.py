"""
Solver routines
---------------

Utilities for solution of FEM systems

"""
from numpy import ndarray
from numpy.linalg import solve
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from typing import Union


def static_sol(mat: Union[csr_matrix, ndarray], rhs: ndarray) -> ndarray:
    """Solve a static problem [mat]{u_sol} = {rhs}
    
    Parameters
    ----------
    mat : Union[csr_matrix, ndarray]
        Array with the system of equations. It can be stored in
        dense or sparse scheme.
    rhs : ndarray
        Array with right-hand-side of the system of equations.
    
    Returns
    -------
    u_sol : ndarray
        Solution of the system of equations.

    Raises
    ------
    TypeError
        If `mat` is neither a `numpy.ndarray` nor a `scipy.sparse.csr_matrix`.
    """
    if type(mat) is csr_matrix:
        u_sol = spsolve(mat, rhs)
    elif type(mat) is ndarray:
        u_sol = solve(mat, rhs)
    else:
        raise TypeError("Matrix should be numpy array or csr_matrix.")

    return u_sol
