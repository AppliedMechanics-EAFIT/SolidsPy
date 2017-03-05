"""
PROGRAM SOLIDS
--------------

Computes the displacement solution for a finite element assembly
of 2D solids under point loads using as input easy-to-create
text files containing element, nodal, materials and loads data.
Fortran subroutines mesher.for and contour.for are also available to
write the required input files out of a Gmesh (.msh) generated file
and to convert the results file into Gmesh post-processor files.

Created by Juan Gomez and Nicolas Guarin-Zapata as part of the courses:

IC0283 COMPUTATIONAL MODELLING
IC0602 INTRODUCTION TO THE FINITE ELEMENT METHOD
Universidad EAFIT
Departamento de Ingenieria Civil

Last updated February 2017
"""
import sys
import numpy as np
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
from datetime import datetime
import matplotlib.pyplot as plt

version = sys.version_info.major

if version == 3:
    raw_input = input
elif version == 2:
    pass
else:
    raise ValueError("You should use Python 2.x at least!")

try:
    import easygui
    folder = easygui.diropenbox(title="Folder for the job") + "/"
    name = easygui.enterbox("Enter the job name")
    echo = easygui.buttonbox("Do you want to echo files?",
                             choices=["Yes", "No"])  
except:
    folder = raw_input('Enter folder (empty for the current one): ')
    name   = raw_input('Enter the job name: ')
    echo   = raw_input('Do you want to echo files? (y/N):')

start_time = datetime.now()
"""
   PRE-PROCESSING
"""
nodes, mats, elements, loads = pre.readin(folder=folder)
if echo.capitalize() in ["YES", "Y"]:
    pre.echomod(nodes, mats, elements, loads, folder=folder)
ne, nn, nm, nl = pre.proini(nodes, mats, elements, loads)
DME , IBC , neq = ass.DME(nn , ne , nodes , elements)
KG = np.zeros([neq, neq])
"""
   SYSTEM ASSEMBLY
"""
for i in range(ne):
    kloc , ndof  = ass.retriever(elements , mats  , nodes , i)
    KG = ass.assembler(KG , neq , kloc , ndof , DME , i)
RHSG = ass.loadasem(loads, IBC, neq, nl)
"""
   SYSTEM SOLUTION
"""
UG = np.linalg.solve(KG, RHSG)
if not(np.allclose(np.dot(KG, UG), RHSG)):
    print("The system is not in equilibrium!")
end_time = datetime.now()
print('Duration for system solution: {}'.format(end_time - start_time))

"""
  POST-PROCCESSING
"""
start_time = datetime.now()
UC = pos.complete_disp(IBC, nodes, UG)
pos.plot_disp(UC, nodes, elements)
# Scatter displacements over the elements
UU = pos.scatter(DME, UG, ne, neq, elements)
pos.gmeshpost(IBC, nn, UG, folder=folder)
# Generates points inside the elements and computes strain solution
E_nodes, S_nodes = pos.strain_nodes(nodes , UU , ne , nn , elements , mats)
pos.plot_strain(E_nodes, nodes, elements)
pos.plot_stress(S_nodes, nodes, elements, plt_type="pcolor")
eigs1, eigs2, vecs1, vecs2 = pos.principal_dirs(S_nodes)
tri = pos.mesh2tri(nodes, elements)
pos.tri_plot(tri, eigs1, title=r"Maximum stress: $\sigma_1$")
plt.quiver(nodes[:, 1], nodes[:, 2], vecs1[:,0], vecs1[:,1], pivot="middle",
           headwidth=1.5, headlength=2.5)
end_time = datetime.now()
print('Duration for post processing: {}'.format(end_time - start_time))
print('Program terminated succesfully!')
plt.show()
