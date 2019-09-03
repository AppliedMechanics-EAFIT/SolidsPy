# -*- coding: utf-8 -*-
"""
solids_GUI: simple interface
----------------------------

Computes the displacement solution for a finite element assembly
of 2D solids under point loads using as input easy-to-create
text files containing element, nodal, materials and loads data.
The input files are created out of a Gmsh (.msh) generated file
using the Python module ``meshio``.

Created by Juan Gomez and Nicolas Guarin-Zapata as part of the courses:

- IC0283 Computational Modeling
- IC0602 Introduction to the Finite Element Method

Which are part of the Civil Engineering Department at Universidad
EAFIT.

"""
from __future__ import absolute_import, division, print_function
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import solidspy.preprocesor as pre
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol


def solids_GUI(plot_contours=True, compute_strains=False, folder=None):
    """
    Run a complete workflow for a Finite Element Analysis

    Parameters
    ----------
    plot_contours : Bool (optional)
        Boolean variable to plot contours of the computed variables.
        By default it is True.
    compute_strains : Bool (optional)
        Boolean variable to compute Strains and Stresses at nodes.
        By default it is False.
    folder : string (optional)
        String with the path to the input files. If not provided
        it would ask for it in a pop-up window.

    Returns
    -------
    UC : ndarray (nnodes, 2)
        Displacements at nodes.
    E_nodes : ndarray (nnodes, 3), optional
        Strains at nodes. It is returned when `compute_strains` is True.
    S_nodes : ndarray (nnodes, 3), optional
        Stresses at nodes. It is returned when `compute_strains` is True.

    """
    if folder is None:
        folder = pre.initial_params()
    start_time = datetime.now()
    echo = False

    # Pre-processing
    nodes, mats, elements, loads = pre.readin(folder=folder)
    if echo:
        pre.echomod(nodes, mats, elements, loads, folder=folder)
    assem_op, bc_array, neq = ass.DME(nodes[:, -2:], elements)
    print("Number of nodes: {}".format(nodes.shape[0]))
    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))

    # System assembly
    stiff_mat, _ = ass.assembler(elements, mats, nodes[:, :3], neq, assem_op)
    rhs_vec = ass.loadasem(loads, bc_array, neq)

    # System solution
    disp = sol.static_sol(stiff_mat, rhs_vec)
    if not np.allclose(stiff_mat.dot(disp)/stiff_mat.max(),
                       rhs_vec/stiff_mat.max()):
        print("The system is not in equilibrium!")
    end_time = datetime.now()
    print('Duration for system solution: {}'.format(end_time - start_time))

    # Post-processing
    start_time = datetime.now()
    disp_complete = pos.complete_disp(bc_array, nodes, disp)
    strain_nodes, stress_nodes = None, None
    if compute_strains:
        strain_nodes, stress_nodes = pos.strain_nodes(nodes, elements, mats,
                                            disp_complete)
    if plot_contours:
        pos.fields_plot(elements, nodes, disp_complete, E_nodes=strain_nodes,
                        S_nodes=stress_nodes)
    end_time = datetime.now()
    print('Duration for post processing: {}'.format(end_time - start_time))
    print('Analysis terminated successfully!')
    if compute_strains:
        return (disp_complete, strain_nodes, stress_nodes)
    else:
        return disp_complete



def solids_auto(data, plot_contours=True, compute_strains=False):
    """
    Run a complete workflow for a Finite Element Analysis

    Parameters
    ----------
    data : dict
        Simulation data composed of nodes, constrains, elements,
        materials and loads.
    plot_contours : Bool (optional)
        Boolean variable to plot contours of the computed variables.
        By default it is True.
    compute_strains : Bool (optional)
        Boolean variable to compute Strains and Stresses at nodes.
        By default it is False.

    Returns
    -------
    UC : ndarray (nnodes, 2)
        Displacements at nodes.
    E_nodes : ndarray (nnodes, 3), optional
        Strains at nodes. It is returned when `compute_strains` is True.
    S_nodes : ndarray (nnodes, 3), optional
        Stresses at nodes. It is returned when `compute_strains` is True.

    """
    # Retrieving data
    nodes = data["nodes"]
    cons = data["cons"]
    elements = data["elements"]
    mats = data["mats"]
    loads = data["loads"]


    # Pre-processing
    assem_op, bc_array, neq = ass.DME(cons, elements)
    print("Number of nodes: {}".format(nodes.shape[0]))
    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))

    # System assembly
    stiff_mat, _ = ass.assembler(elements, mats, nodes, neq, assem_op)
    rhs_vec = ass.loadasem(loads, bc_array, neq)

    # System solution
    start_time = datetime.now()
    disp = sol.static_sol(stiff_mat, rhs_vec)
    if not np.allclose(stiff_mat.dot(disp)/stiff_mat.max(),
                       rhs_vec/stiff_mat.max()):
        print("The system is not in equilibrium!")
    end_time = datetime.now()
    print('Duration for system solution: {}'.format(end_time - start_time))

    # Post-processing
    start_time = datetime.now()
    disp_complete = pos.complete_disp(bc_array, nodes, disp)
    strain_nodes, stress_nodes = None, None
    if compute_strains:
        strain_nodes, stress_nodes = pos.strain_nodes(nodes, elements, mats,
                                                      disp_complete)
    if plot_contours:
        pos.fields_plot(elements, nodes, disp_complete, E_nodes=strain_nodes,
                        S_nodes=stress_nodes)
    end_time = datetime.now()
    print('Duration for post processing: {}'.format(end_time - start_time))
    print('Analysis terminated successfully!')
    if compute_strains:
        return (disp_complete, strain_nodes, stress_nodes)
    else:
        return disp_complete


if __name__ == '__main__':
    displacement = solids_GUI()
    plt.show()
