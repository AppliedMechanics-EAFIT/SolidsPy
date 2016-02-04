"""
PROGRAM SOLIDS
--------------

Computes the displacement solution for a finite element assembly
of 2D solids under point loads using as input easy-to-create
text files containing element, nodal, materials and loads data.
Fortran subroutines mesher.for and contour.for are also available to
write the required input files out of a Gmesh (.msh) generated file
and to convert the results file into Gmesh post-processor files.

Created by Juan Gomez and Nicolas Guarin-Zapata as part of the course:

IC0283 COMPUTATIONAL MODELLING
Universidad EAFIT
Departamento de Ingenieria Civil

Last updated January 2016
"""
import numpy as np
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
from datetime import datetime
import matplotlib.pyplot as plt
import easygui


folder = easygui.diropenbox(title="Folder for the job") + "/"
name = easygui.enterbox("Enter the job name")
echo = easygui.buttonbox("Do you want to echo files?",
                         choices=["Yes", "No"])                        
start_time = datetime.now()

#%%
#   MODEL ASSEMBLY
#
# Reads the model
nodes, mats, elements, loads = pre.readin(folder=folder)
if echo == "Yes":
    pre.echomod(nodes, mats, elements, loads, folder=folder)
# Retrieves problem parameters
ne, nn, nm, nl, COORD = pre.proini(nodes, mats, elements, loads)
# Counts equations and creates BCs array IBC
neq, IBC = ass.eqcounter(nn, nodes)
# Computes assembly operatorring
DME, IELCON = ass.DME(IBC, ne, elements)
# Assembles Global Stiffness Matrix KG
KG = ass.matassem(IBC, mats, elements, nn, ne, neq, COORD, DME, IELCON)
# Assembles Global Rigth Hand Side Vector RHSG
RHSG = ass.loadasem(loads, IBC, neq, nl)

#%%
#   SYSTEM SOLUTION
#
# Solves the system
UG = np.linalg.solve(KG, RHSG)
if not(np.allclose(np.dot(KG, UG), RHSG)):
    print("The system is not in equilibrium!")
end_time = datetime.now()
print('Duration for system solution: {}'.format(end_time - start_time))

#%%
#   POST-PROCCESSING
#
start_time = datetime.now()
UC = pos.complete_disp(IBC, nodes, UG)
pos.plot_disp(UC, nodes, elements)

# Scatter displacements over the elements
UU = pos.scatter(DME, UG, ne, neq, elements)
# Generates points inside the elements and computes strain solution
E_nodes, S_nodes = pos.strain_nodes(IELCON, UU, ne, COORD, elements, mats)
pos.plot_strain(E_nodes, nodes, elements)
end_time = datetime.now()
print('Duration for post processing: {}'.format(end_time - start_time))
print('Program terminated succesfully!')
plt.show()
