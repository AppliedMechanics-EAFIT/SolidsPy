# -*- coding: utf-8 -*-
"""
Template to generate the input files for the FEM code solids_ISO.
The script uses module meshio.py to read a GMSH mesh and produce
text files nodes.txt, eles.txt , mater.txt and loads.txt

@authors: Juan Gomez
          Nicolas Guarin-Zapata
"""
from __future__ import division, print_function
import numpy as np
import meshio
import solidspy.preprocesor as msh


# Read the mesh file
points, cells, point_data, cell_data, field_data = \
    meshio.read("bridge.msh")

# Create elements
nodes_array = msh.node_writer(points, point_data)
nf , els1_array = msh.ele_writer(cells, cell_data, "triangle", 10000 , 3 , 0 , 0)
nf , els2_array = msh.ele_writer(cells, cell_data, "triangle", 20000 , 3 , 1, nf)
els_array =np.append(els1_array, els2_array, axis=0)

# Create nodes
nodes_array = msh.boundary_conditions(cells, cell_data, 100, nodes_array, 0 ,-1)
nodes_array = msh.boundary_conditions(cells, cell_data, 200, nodes_array, -1, 0)
cargas = msh.loading(cells, cell_data, 300, 0.0, -1.0e6)


# Save files
np.savetxt("eles.txt" , els_array   , fmt="%d")
np.savetxt("loads.txt", cargas, fmt=("%d", "%.6f", "%.6f"))
np.savetxt("nodes.txt", nodes_array , fmt=("%d", "%.4f", "%.4f", "%d", "%d"))
