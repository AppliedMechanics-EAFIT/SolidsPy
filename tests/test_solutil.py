# -*- coding: utf-8 -*-
"""
Test cases for functions on ``solutil`` module

"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import rand
import solidspy.solutil as sol


def test_static_solve():
    """Tests for static solver"""
    
    # Identity matrix and rhs = {1, 2, 3}
    mat = np.eye(3)
    rhs = np.array([1, 2, 3])
    u_sol = sol.static_sol(mat, rhs)
    
    row = np.array([0, 1, 2])
    col = np.array([0, 1, 2])
    data = np.array([1, 1, 1])
    mat_sparse = csr_matrix((data, (row, col)), shape=(3, 3))
    u_sol2 = sol.static_sol(mat_sparse, rhs)
    assert np.allclose(u_sol, u_sol2)
    assert np.allclose(u_sol, [1, 2, 3])
    
    # Random matrices and right hand side
    np.random.seed(1)
    ntest = 10
    for cont in range(ntest):
        rhs = np.random.rand(100)
        mat = rand(100, 100, density=0.3)
        mat = 0.5*(mat + mat.transpose())
        u_sol = sol.static_sol(mat.toarray(), rhs)
        u_sol2 = sol.static_sol(mat.tocsr(), rhs)
        print(mat.toarray())
        assert np.allclose(u_sol, u_sol2)
