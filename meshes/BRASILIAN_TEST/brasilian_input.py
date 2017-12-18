# -*- coding: utf-8 -*-
"""
Genera archivos de entrada para el programa de elementos finitos
FEM_iso para la prueba brasilera, usando 2 simetrías en el modelo.

@author: Nicolas Guarin-Zapata
@date: Mayo 18, 2017
"""
from __future__ import division, print_function
import numpy as np
import meshio
import solidspy.preprocesor as msh 


points, cells, point_data, cell_data, field_data = \
    meshio.read("Prueba_brasilera.msh")

nodes_array = msh.node_writer(points, point_data)
nf , els_array = msh.ele_writer(cells, cell_data, "triangle", 1000 , 3 , 0 , 0)
nodes_array = msh.boundary_conditions(cells, cell_data, 100, nodes_array, -1, 0)
nodes_array = msh.boundary_conditions(cells, cell_data, 200, nodes_array, 0, -1)
np.savetxt("eles.txt" , els_array   , fmt="%d")
np.savetxt("nodes.txt", nodes_array , fmt=("%d", "%.4f", "%.4f", "%d", "%d"))