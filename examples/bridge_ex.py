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



folder = "../MESHES/BRIDGE_FULL/"
#folder = "../MESHES/BRIDGE_HALF/"


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
E_nodes, S_nodes = pos.strain_nodes(IELCON, UU, ne, COORD, elements, mats)
pos.plot_strain(E_nodes, nodes, elements)
tri = pos.mesh2tri(nodes, elements)
eigs1, eigs2, vecs1, vecs2 = pos.principal_dirs(S_nodes)
pos.tri_plot(tri, eigs1)
plt.quiver(x, y, vecs1[:,0], vecs1[:,1], pivot="middle",
           headwidth=1.5, headlength=2.5)

plt.show()

