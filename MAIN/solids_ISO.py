############################################################################
#-------------------------PROGRAM SOLIDS---------------------------------- #
# Computes the displacement solution for a finite element assembly         #
# of finite elements under point loads using as input easy-to-create       #
# text files containing element, nodal, materials and loads data.          #
# Fortran subroutines mesher.for and contour.for are also availbale to     #
# write the required input files out of a Gmesh (.msh) generated file and  #
# and to convert the results file into Gmesh post-processor files.         #
#                                                                          #
# Created by Juan Gomez as part of the course:                             #
# IC0283 COMPUTATIONAL MODELLING                                           #
# Universidad EAFIT                                                        #
# Departamento de Ingenieria Civil                                         #   
#                                                                          #
#  (Last updated December 2015)                                            #
############################################################################
"""
"""
import numpy as np
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
import shutil as shu
import os
from datetime import datetime
start_time = datetime.now()
print 'Enter the job name?'
nombre = raw_input ()
#
#     MODEL ASSEMBLY
#
nodes, mats, elements, loads =pre.readin()                                # Reads the model
#pre.echomod(nodes,mats,elements,loads)                                   # Activate to generate echo files
ne,nn,nm,nl,COORD =pre.proini(nodes,mats,elements,loads)                  # Retrieves problem parameters
neq, IBC = ass.eqcounter(nn,nodes)                                        # Counts equations and creates BCs array IBC
DME, IELCON =ass.DME(IBC,ne,elements)                                     # Computes assembly operator     
KG = ass.matassem(IBC, mats,elements,nn,ne,neq,COORD, DME, IELCON)        # Assembles Global Stiffness Matrix KG
RHSG=ass.loadasem(loads,IBC,neq,nl)                                       # Assembles Global Rigth Hand Side Vector RHSG
#
#    SYSTEM SOLUTION
#
UG = np.linalg.solve(KG, RHSG)                                            # Solves the system
l=np.allclose(np.dot(KG, UG), RHSG)
print l
end_time = datetime.now()
print('Duration for system solution: {}'.format(end_time - start_time))
#
#    POST-PROCCESSING
#
start_time = datetime.now()
xmin,xmax,ymin,ymax=pos.axisscale(COORD,nn)                               # Sets axis for visualization window.
pos.plotdis(IBC,UG,nodes,nn,xmin,xmax,ymin,ymax)                          # Plot displacement solution
#
pos.gmeshpost(IBC,nn,UG)                                                  # Generates displacement solution to be post-processed via Gmesh.
nomfile1='../MESHUTIL/'+nombre+'.msh'
nomfileH='../MESHUTIL/'+nombre+'H.msh'
nomfileV='../MESHUTIL/'+nombre+'V.msh'
nomfileF='../MESHUTIL/'+nombre+'F.msh'
shu.copy(nomfile1,nomfileH)
shu.copy(nomfile1,nomfileV)
shu.copy(nomfile1,nomfileF)
#
# Scatters nodal displacements over the elements
# and plots strain solution. (Activate if required.)
#
UU=pos.scatter(DME,UG,ne,neq,elements)                                   # Scatter displacements over the elements
EG,XS=pos.strainGLO(IELCON,UU,ne,COORD,elements)                         # Generates points inside the elements and computes strain solution.                                 
pos.plotstrain(EG,XS,ne,xmin,xmax,ymin,ymax)                             # Plot strain solution
#
end_time = datetime.now()
print('Duration for post processing: {}'.format(end_time - start_time))                      
print 
print('Program terminated succesfuly')
