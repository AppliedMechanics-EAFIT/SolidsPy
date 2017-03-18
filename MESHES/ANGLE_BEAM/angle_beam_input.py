"""
Generate input files from a Gmsh mesh file for an angle beam

It imposes encastre constraints for the left side (x==0)
and shear loading for the right side (x==10).
"""
from __future__ import division
import meshio
import numpy as np 


points, cells, point_data, cell_data, field_data = \
    meshio.read("angle_beam.msh")

# Elements data
elements = cells["quad"]
els_array = np.ones([elements.shape[0], 7], dtype=int)
els_array[:, 0] = range(elements.shape[0])
els_array[:, 3::] = elements
# Nodes data
nodes_array = np.zeros([points.shape[0], 5])
nodes_array[:, 0] = range(points.shape[0])
nodes_array[:, 1:3] = points[:, :2]
nodes_array[nodes_array[:, 1]==0, 3:] = -1
# Loads data
nloads = points[points[:, 0] == 10, 0].shape[0]
loads_array = np.zeros((nloads, 3))
loads_array[:, 0] = nodes_array[points[:, 0]==10, 0]
loads_array[:, -1] = -1/nloads
# Material data
mater_array = np.array([[1e3, 1/3], [1e3, 1/3]])

np.savetxt("eles.txt", els_array, fmt="%d")
np.savetxt("nodes.txt", nodes_array, fmt=("%d", "%.4f", "%.4f", "%d", "%d"))
np.savetxt("loads.txt", loads_array, fmt=("%d", "%.6f", "%.6f"))
np.savetxt("mater.txt", mater_array, fmt="%.6f")