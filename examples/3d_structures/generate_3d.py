import numpy as np
from typing import Tuple
from numpy.typing import NDArray

def structure_3d(
    L: float = 10.0,
    H: float = 10.0,
    W: float = 10.0,
    E: float = 206.8e9,
    v: float = 0.28,
    nx: int = 10,
    ny: int = 10,
    nz: int = 10,
    dirs: NDArray[np.float64] = np.array([]),
    positions: NDArray[np.float64] = np.array([])
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.int_], NDArray[np.float64], NDArray[np.int_]]:
    """
    Make the mesh for a cubic model with cantilever constraints in 3D.

    Parameters
    ----------
    L : float, optional
        Length of the structure (x-direction), by default 10.0
    H : float, optional
        Height of the structure (y-direction), by default 10.0
    W : float, optional
        Width of the structure (z-direction), by default 10.0
    E : float, optional
        Young's modulus, by default 206.8e9
    v : float, optional
        Poisson's ratio, by default 0.28
    nx : int, optional
        Number of elements in the x direction, by default 10
    ny : int, optional
        Number of elements in the y direction, by default 10
    nz : int, optional
        Number of elements in the z direction, by default 10
    dirs : NDArray[np.float64], optional
        Array of load directions, shape (n_loads, 3).
        Example: [[0,1,0],[1,0,0],[0,0,-1]]
        By default an empty array.
    positions : NDArray[np.float64], optional
        Array of load positions (in grid indexing), shape (n_loads, 3).
        Example: [[5,5,9], [1,1,9], [8,8,9]]
        By default an empty array.

    Returns
    -------
    nodes : NDArray[np.float64]
        Node array of shape ((nx+1)*(ny+1)*(nz+1), 7).
        Columns: [node_id, x, y, z, BC_x, BC_y, BC_z].
    mats : NDArray[np.float64]
        Material properties array of shape (number_of_elements, 3).
        Each row: [E, v, placeholder].
    els : NDArray[np.int_]
        Element connectivity array of shape (nx*ny*nz, 11).
        For example: [element_id, material_id, 0, n0, n2, n3, n1, n4, n6, n7, n5].
    loads : NDArray[np.float64]
        Loads array of shape (dirs.shape[0], 4).
        Columns: [node_id, Fx, Fy, Fz].
    idx_BC : NDArray[np.int_]
        Indices of nodes on the constrained face (z = 0).
    """
    # Generate 3D grid of nodes
    x = np.linspace(-L / 2, L / 2, nx + 1)
    y = np.linspace(-H / 2, H / 2, ny + 1)
    z = np.linspace(0, W, nz + 1)
    xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')

    nodes = np.zeros(((nx + 1) * (ny + 1) * (nz + 1), 7))
    nodes[:, 0] = np.arange(nodes.shape[0])   # Node IDs
    nodes[:, 1] = xv.ravel()
    nodes[:, 2] = yv.ravel()
    nodes[:, 3] = zv.ravel()

    # Build hexahedral elements
    els_list = []
    count = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                n0 = k + (nz + 1) * (j + (ny + 1) * i)
                n1 = n0 + 1
                n2 = n0 + (nz + 1)
                n3 = n2 + 1
                n4 = n0 + (nz + 1) * (ny + 1)
                n5 = n4 + 1
                n6 = n4 + (nz + 1)
                n7 = n6 + 1
                # [element_id, material_id, 0, node0, node2, node3, node1, node4, node6, node7, node5]
                els_list.append([count, 9, 0, n0, n2, n3, n1, n4, n6, n7, n5])
                count += 1
    els = np.array(els_list, dtype=int)

    # Assign material properties
    mats = np.zeros((els.shape[0], 3))
    mats[:, 0] = E
    mats[:, 1] = v
    mats[:, 2] = 1

    # Apply constraints on face z = 0
    mask_z0 = (nodes[:, 3] == 0.0)
    nodes[mask_z0, 4] = -1  # BC_x
    nodes[mask_z0, 5] = -1  # BC_y
    nodes[mask_z0, 6] = -1  # BC_z
    idx_BC = nodes[mask_z0, 0].astype(int)

    # Define loads on face z = W
    loads = np.zeros((dirs.shape[0], 4))
    mask_load = (nodes[:, 3] == W)
    node_indices = np.where(mask_load)[0]

    # The user supplies positions in a "grid" sense: [x_index, y_index, z_index].
    # We slice out the correct nodes from node_indices. Implementation is up to the user.
    # Example: node_indices_selected = ...
    # For demonstration, the code uses the provided positions as if
    # positions[:,0] is x_index, positions[:,1] is y_index,
    # and we assume z_index = positions[:,2] (already at z=W).
    # This line is a placeholder; adjust indexing logic as needed.
    node_indices_selected = node_indices[
        positions[:, 0] * (ny + 1) + positions[:, 1]
    ]

    loads[:, 0] = node_indices_selected
    loads[:, 1] = dirs[:, 0]  # Fx
    loads[:, 2] = dirs[:, 1]  # Fy
    loads[:, 3] = dirs[:, 2]  # Fz

    return nodes, mats, els, loads, idx_BC