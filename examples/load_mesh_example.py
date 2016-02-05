# -*- coding: utf-8 -*-
"""
Load a GMSH mesh and plot contours for its height.
"""
import meshio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['image.cmap'] = "YlGnBu_r"

points, cells, point_data, cell_data, field_data = \
    meshio.read("../MESHES/DAM/dam.msh")

x = points[:, 0]
y = points[:, 1]
tri = Triangulation(x, y, cells['triangle'])
plt.tricontourf(tri, y, 12, shading="gourad")
plt.axis("image")
plt.show()