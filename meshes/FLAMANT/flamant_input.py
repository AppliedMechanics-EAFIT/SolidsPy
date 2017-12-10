# -*- coding: utf-8 -*-
"""
Solucion de Flamant
@author: Juan Gomez
Ver Timoshenko and Young: Theory of elasticity
"""
from __future__ import division, print_function
import numpy as np
import meshio
import solidspy.preprocesor as msh
#
points, cells, point_data, cell_data, field_data = \
    meshio.read("flamant.msh")
    
    
nodes_array = msh.node_writer(points, point_data)
nf , els_array = msh.ele_writer(cells, cell_data, "triangle", 10000 , 3 , 0 , 0)
nodes_array = msh.boundary_conditions(cells, cell_data,  100 , nodes_array, -1, -1)
#
# Datos de los elementos
#
#eles = cells["triangle"]
#els_array = np.zeros([eles.shape[0], 6], dtype=int)
#els_array[:, 0] = range(eles.shape[0])              # Asigna numero de elemento
#els_array[:, 1] = 3                                 # Asigna tipo de elemento
#els_array[:, 3::] = eles                            # Asigna nudos de los elementos
##
## Datos de los nudos
##
#nodes_array = np.zeros([points.shape[0], 5])
#nodes_array[:, 0] = range(points.shape[0])          # Asigna numeros de nudos
#nodes_array[:, 1:3] = points[:, :2]                 # Asigna coordenadas
##
#lines = cells["line"]
#bounds = cell_data["line"]["physical"]              # Bounds contiene la infromacion de la linea fisica
#id_frontera = [cont for cont in range(len(bounds)) if bounds[cont] == 7]
#nodes_frontera = lines[id_frontera]
#nodes_frontera.shape = nodes_frontera.shape[0] * 2
#nodes_frontera = list(set(nodes_frontera))
#nodes_array[nodes_frontera, 3:5] = -1
#
# Datos de los materiales
#
#
# Escribe los archivos 
#
np.savetxt("eles.txt" , els_array   , fmt="%d")
np.savetxt("nodes.txt", nodes_array , fmt=("%d", "%.4f", "%.4f", "%d", "%d"))