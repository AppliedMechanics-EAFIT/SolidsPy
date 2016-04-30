# -*- coding: utf-8 -*-
"""

"""
from __future__ import division
from os import sys
sys.path.append("../../MAIN/")
import numpy as np
import postprocesor as pos
import assemutil as ass
import matplotlib.pyplot as plt
from beam import beam_sln


def rect_grid(L, h, nx, ny):
    """
    """
    y, x = np.mgrid[-h/2:h/2:ny*1j, 0:L:nx*1j]
    els = np.zeros(((nx - 1)*(ny - 1), 7), dtype=int)
    els[:, 1] = 1
    for row in range(ny - 1):
        for col in range(nx - 1):
            cont = row*(nx - 1) + col
            els[cont, 0] = cont
            els[cont, 3:7] = [cont + row, cont + row + 1,
                              cont + row + nx + 1, cont + row + nx]
    return x.flatten(), y.flatten(), els


P = -50
E = 1000
nu = 0.3
L = 24
h = 8
I = 42.67
niter = 3
err = np.zeros((niter))
mats = np.array([[E, nu], [E, nu]])
for cont in range(1, niter + 1):
    nx = 3*2**cont + 1
    ny = 2**cont + 1
    x, y, els = rect_grid(L, h, nx, ny)
    nodes = np.zeros((nx*ny, 5))
    nodes[:, 0] = range(nx*ny)
    nodes[:, 1] = x
    nodes[:, 2] = y
    nodes[x==L, 3] = -1
    nodes[nx*(ny//2 + 1) - 1, 4] = -1    
    loads = np.zeros((ny, 3))
    loads[:, 0] = nodes[x==0, 0]
    loads[:, 2] = P/ny

    # Assembly
    nn = nx*ny
    ne = (nx - 1)*(ny - 1)
    neq, IBC = ass.eqcounter(nx*ny, nodes)
    DME, IELCON = ass.DME(IBC, ne, els)
    KG = ass.matassem(IBC, mats, els, nn, ne, neq, nodes[:, 1:3],
                      DME, IELCON)
    RHSG = ass.loadasem(loads, IBC, neq, ny)

    # Solution
    UG = np.linalg.solve(KG, RHSG)
    UC = pos.complete_disp(IBC, nodes, UG)
    u, v, _, _, _  = beam_sln(x, y, nu, P, E, I, L, h)
    U_anal = np.column_stack([u, v])
    aux = np.linalg.norm(U_anal - UC)
    err[cont - 1] = aux/np.linalg.norm(U_anal)

#%% Analysis of error
pos.plot_disp(UC, nodes, els, title="FEM:")
pos.plot_disp(U_anal, nodes, els, title="Exact:")
x = np.linspace(1, niter + 1, niter)
plt.figure()
plt.loglog(1/2**x, err)
plt.xlabel(r"$h$")
plt.ylabel(r"$\frac{\Vert u - u_h \Vert}{\Vert u \Vert}$")
plt.show()