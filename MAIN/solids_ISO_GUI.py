"""
PROGRAM SOLIDS
--------------

Computes the displacement solution for a finite element assembly
of 2D solids under point loads using as input easy-to-create
text files containing element, nodal, materials and loads data.
Fortran subroutines mesher.for and contour.for are also available to
write the required input files out of a Gmsh (.msh) generated file
and to convert the results file into Gmsh post-processor files.

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
from scipy.sparse.linalg import spsolve
import preprocesor as pre
import postprocesor as pos
import assemutil as ass

folder, name, echo = pre.initial_params()
start_time = datetime.now()


#%% PRE-PROCESSING
nodes, mats, elements, loads = pre.readin(folder=folder)
if echo:
    pre.echomod(nodes, mats, elements, loads, folder=folder)
ne, nn, nm, nl = pre.proini(nodes, mats, elements, loads)
DME , IBC , neq = ass.DME(nn , ne , nodes , elements)
print("Number of nodes: {}".format(nn))
print("Number of elements: {}".format(ne))
print("Number of equations: {}".format(neq))

#%% SYSTEM ASSEMBLY
KG = ass.assembler(elements, mats, nodes, neq, DME)
RHSG = ass.loadasem(loads, IBC, neq, nl)

#%% SYSTEM SOLUTION
UG = spsolve(KG, RHSG)
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
