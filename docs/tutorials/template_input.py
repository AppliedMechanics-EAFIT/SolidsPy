# -*- coding: utf-8 -*-
"""
Template to generate the input files for the FEM code SolidsPy.
The script uses meshio to read a GMSH mesh and produce
text files nodes.txt, eles.txt , mater.txt and loads.txt

@authors: Juan Gomez
         Nicolas Guarin-Zapata
"""
import meshio
import numpy as np
#
# Read the GMSH file using meshio
#
mesh = meshio.read("template.msh")
points = mesh.points
cells = mesh.cells
point_data = mesh.point_data
cell_data = mesh.cell_data


# Process element data. In this case we have used 3-noded triangles.
# (In SolidsPy 1 stands for quad; 2 for triangle6; 3 for triangle)
eles = cells["triangle"]
els_array = np.zeros([eles.shape[0], 6], dtype=int)
# Asigns element id
els_array[:, 0] = range(eles.shape[0])
# Assigns element type according to SolidsPy
els_array[:, 1] = 3
# Assign element connectivities
els_array[:, 3::] = eles
# Creates the materials array and assigns a profile (either 0 or 1) to
# each element according to the physical surface
mater_array = np.array([[1.0, 0.30],
                        [5.0, 0.30]])
maters = cell_data["triangle"]["gmsh:physical"]
els_array[:, 2]  = [0 if mater == 100 else 1 for mater in maters]

# Process nodal data
nodes_array = np.zeros([points.shape[0], 5])
nodes_array[:, 0] = range(points.shape[0])  # Assigns nodal id
nodes_array[:, 1:3] = points[:, :2]  # Assigns space coordinates

# Process physical lines to assign displacement boundary conditions
# and nodal point loads
#
# Displacement boundary conditions
lines = cells["line"]
# Bounds contains data corresponding to the physical line.
bounds = cell_data["line"]["gmsh:physical"]
# (-1 means a restraint degree of freedom)
id_frontera_lat = [cont for cont in range(len(bounds)) if bounds[cont] == 300]
nodes_frontera_lat = lines[id_frontera_lat]
nodes_frontera_lat = nodes_frontera_lat.flatten()
nodes_frontera_lat = list(set(nodes_frontera_lat))
nodes_array[nodes_frontera_lat, 3] = -1
id_frontera_abajo  = [cont for cont in range(len(bounds)) if bounds[cont] == 400]
nodes_frontera_abajo = lines[id_frontera_abajo]
nodes_frontera_abajo = nodes_frontera_abajo.flatten()
nodes_frontera_abajo = list(set(nodes_frontera_abajo))
nodes_array[nodes_frontera_abajo, 4] = -1
# Nodal point loads
id_carga = [cont for cont in range(len(bounds)) if bounds[cont] == 500]
nodes_carga = lines[id_carga]
nodes_carga = nodes_carga.flatten()
nodes_carga = list(set(nodes_carga))
ncargas = len(nodes_carga)
carga_total = -2.0
cargas = np.zeros((ncargas, 3))
cargas[:, 0] = nodes_carga
cargas[:, 2] = carga_total/ncargas

# Write the model .txt files
np.savetxt("eles.txt" , els_array   , fmt="%d")
np.savetxt("loads.txt", cargas, fmt=("%d", "%.6f", "%.6f"))
np.savetxt("nodes.txt", nodes_array , fmt=("%d", "%.4f", "%.4f", "%d", "%d"))
np.savetxt("mater.txt", mater_array , fmt="%.6f")
