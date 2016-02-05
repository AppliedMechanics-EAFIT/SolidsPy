"""


Created by Juan Gomez and Nicolas Guarin-Zapata as part of the course:

IC0283 COMPUTATIONAL MODELLING
Universidad EAFIT
Departamento de Ingenieria Civil

Last updated January 2016
"""
from os import sys
sys.path.append("../MAIN/")
import numpy as np
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
import matplotlib.pyplot as plt
from scipy.interpolate import SmoothBivariateSpline


folder = "../MESHES/RING/"


#%%  MODEL ASSEMBLY

nodes, mats, elements, loads = pre.readin(folder=folder)
ne, nn, nm, nl, COORD = pre.proini(nodes, mats, elements, loads)
neq, IBC = ass.eqcounter(nn, nodes)
DME, IELCON = ass.DME(IBC, ne, elements)
KG = ass.matassem(IBC, mats, elements, nn, ne, neq, COORD, DME, IELCON)

RHSG = ass.loadasem(loads, IBC, neq, nl)

#%% SYSTEM SOLUTION
UG = np.linalg.solve(KG, RHSG)
if not(np.allclose(np.dot(KG, UG), RHSG)):
    print("The system is not in equilibrium!")



#%% POST-PROCCESSING
UC = pos.complete_disp(IBC, nodes, UG)
pos.plot_disp(UC, nodes, elements)

UU = pos.scatter(DME, UG, ne, neq, elements)
x = nodes[:, 1]
y = nodes[:, 2]
E_gauss, pts_gauss = pos.strainGLO(IELCON, UU, ne, COORD, elements)
E_int1 = SmoothBivariateSpline(pts_gauss[:, 0], pts_gauss[:, 1],
                               E_gauss[:, 0])
E_int2 = SmoothBivariateSpline(pts_gauss[:, 0], pts_gauss[:, 1],
                               E_gauss[:, 1])
E_int3 = SmoothBivariateSpline(pts_gauss[:, 0], pts_gauss[:, 1],
                               E_gauss[:, 2])
E1 = E_int1.ev(x, y)
E2 = E_int2.ev(x, y)
E3 = E_int3.ev(x, y)
E_nodes = np.column_stack([E1, E2, E3])
pos.plot_strain(E_nodes, nodes, elements, plt_type="pcolor")
tri = pos.mesh2tri(nodes, elements)
plt.triplot(tri, color='k', alpha=0.3)

plt.show()

