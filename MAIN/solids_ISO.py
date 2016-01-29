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
import shutil as shu
from datetime import datetime
import matplotlib.pyplot as plt


start_time = datetime.now()
name = raw_input('Enter the job name: ')
folder = raw_input('Enter folder (empty for the current one): ')
echo = raw_input('Do you want to echo files? (y/N): ')


## Some Examples

## Ring with radial load
#name = "ring"
#folder = "../MESHES/RING/"
#echo = ""

## Arch bridge
#name = "Limyra bridge arch"
#folder = "../MESHES/Limyra bridge arch/"
#echo = ""


#%%
#   MODEL ASSEMBLY
#
# Reads the model
nodes, mats, elements, loads = pre.readin(folder=folder)
if echo.upper() in ["Y", "YES"]:
    pre.echomod(nodes, mats, elements, loads, folder=folder)
# Retrieves problem parameters
ne, nn, nm, nl, COORD = pre.proini(nodes, mats, elements, loads)
# Counts equations and creates BCs array IBC
neq, IBC = ass.eqcounter(nn, nodes)
# Computes assembly operator
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
print(np.allclose(np.dot(KG, UG), RHSG))
end_time = datetime.now()
print('Duration for system solution: {}'.format(end_time - start_time))

#%%
#   POST-PROCCESSING
#
start_time = datetime.now()
# Sets axis for visualization window
xmin, xmax, ymin, ymax = pos.axisscale(COORD, nn)
# Plot displacement solution
pos.plot_disp2(IBC, UG, nodes, elements)
# Generates displacement solution to be post-processed via Gmesh
pos.gmeshpost(IBC, nn, UG, folder=folder)
nomfile1 = folder + name + '.msh'
nomfileH = folder + name + 'H.msh'
nomfileV = folder + name + 'V.msh'
nomfileF = folder + name + 'F.msh'
shu.copy(nomfile1, nomfileH)
shu.copy(nomfile1, nomfileV)
shu.copy(nomfile1, nomfileF)

# Scatters nodal displacements over the elements
# and plots strain solution. (Activate if required.)
#
# Scatter displacements over the elements
UU = pos.scatter(DME, UG, ne, neq, elements)
# Generates points inside the elements and computes strain solution
E_nodes = pos.strain_nodes(IELCON, UU, ne, COORD, elements)
pos.plot_strain2(E_nodes, nodes, elements)
end_time = datetime.now()
print('Duration for post processing: {}'.format(end_time - start_time))
print('Program terminated succesfully!')
plt.show()
