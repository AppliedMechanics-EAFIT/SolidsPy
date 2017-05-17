"""
PROGRAM SOLIDS
--------------
Computes the displacement solution for a finite element assembly
of 2D solids under point loads using as input easy-to-create
text files containing element, nodal, materials and loads data.
The input files are created out of a Gmsh (.msh) generated file
using python module meshio.py.

Created by Juan Gomez and Nicolas Guarin-Zapata as part of the courses:

IC0283 COMPUTATIONAL MODELLING
IC0602 INTRODUCTION TO THE FINITE ELEMENT METHOD
Universidad EAFIT
Departamento de Ingenieria Civil

"""
from __future__ import division, print_function
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
import solutil as sol

folder, name, echo = pre.initial_params()
start_time = datetime.now()


#%% PRE-PROCESSING
nodes, mats, elements, loads = pre.readin(folder=folder)
if echo:
    pre.echomod(nodes, mats, elements, loads, folder=folder)
DME , IBC , neq = ass.DME(nodes, elements)
print("Number of nodes: {}".format(nodes.shape[0]))
print("Number of elements: {}".format(elements.shape[0]))
print("Number of equations: {}".format(neq))

#%% SYSTEM ASSEMBLY
KG = ass.assembler(elements, mats, nodes, neq, DME)
RHSG = ass.loadasem(loads, IBC, neq)

#%% SYSTEM SOLUTION
UG = sol.static_sol(KG, RHSG)
if not(np.allclose(KG.dot(UG), RHSG)):
    print("The system is not in equilibrium!")
end_time = datetime.now()
del KG
print('Duration for system solution: {}'.format(end_time - start_time))

#%% POST-PROCESSING
start_time = datetime.now()
UC = pos.complete_disp(IBC, nodes, UG)
E_nodes, S_nodes = pos.strain_nodes(nodes , elements, mats, UC, DME)
pos.fields_plot(elements, nodes, UC, E_nodes=E_nodes, S_nodes=S_nodes)
print('Duration for post processing: {}'.format(end_time - start_time))
print('Analysis terminated successfully!')
plt.show()
