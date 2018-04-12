"""
Generate input files from a Gmsh mesh file for a ring with inner pressure
 with symmetry with respect to the horizontal and vertical axws.

It imposes roller constraints for the left side (x==0), and bottom (y==0).
The loading is 1 in the radial direction.
"""
from __future__ import division
import meshio
import numpy as np 


points, cells, point_data, cell_data, field_data = \
    meshio.read("ring.msh")


# Elements data
elements = cells["triangle"]
els_array = np.zeros([elements.shape[0], 6], dtype=int)
els_array[:, 0] = range(elements.shape[0])
els_array[:, 1] = 3
els_array[:, 3::] = elements
# Nodes data
nodes_array = np.zeros([points.shape[0], 5])
nodes_array[:, 0] = range(points.shape[0])
nodes_array[:, 1:3] = points[:, :2]
nodes_array[nodes_array[:, 1]==0, 3] = -1
nodes_array[nodes_array[:, 2]==0, 4] = -1
# Loads data
radius = np.sqrt(points[:, 0]**2 + points[:, 1]**2)
nloads = points[np.abs(radius - 1.5) <= 1e-6, 0].shape[0]
loads_array = np.zeros((nloads, 3))
loads_array[:, 0] = nodes_array[np.abs(radius - 1.5) <= 1e-6, 0]
loads_array[:, 1] = 2.0*(points[np.abs(radius - 1.5) <= 1e-6, 0]/radius[np.abs(radius - 1.5) <= 1e-6])
loads_array[:, 2] = 2.0*(points[np.abs(radius - 1.5) <= 1e-6, 1]/radius[np.abs(radius - 1.5) <= 1e-6])
## Material data
mater_array = np.array([[1e3, 1/3]])

np.savetxt("eles.txt", els_array, fmt="%d")
np.savetxt("nodes.txt", nodes_array, fmt=("%d", "%.4f", "%.4f", "%d", "%d"))
np.savetxt("loads.txt", loads_array, fmt=("%d", "%.6f", "%.6f"))
np.savetxt("mater.txt", mater_array, fmt="%.6f")