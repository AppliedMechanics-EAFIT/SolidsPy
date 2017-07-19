# -*- coding: utf-8 -*-
"""
Computes the displacement solution for a finite element assembly
of 2D solids under point loads using as input easy-to-create
text files containing element, nodal, materials and loads data.
The input files are created out of a Gmsh (.msh) generated file
using python module meshio.py.

Created by Juan Gomez and Nicolas Guarin-Zapata as part of the courses:

IC0283 COMPUTATIONAL MODELLING
IC0602 INTRODUCTION TO THE FINITE ELEMENT METHOD
Universidad EAFIT
Departamento de Ingenieria Civil

"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import solidspy.preprocesor as pre
import solidspy.postprocesor as pos
import solidspy.assemutil as ass
import solidspy.solutil as sol

def solids_GUI(compute_strains=False, plot_contours=True):
    """
    Run a complete workflow for a Finite Element Analysis
    
    Parameters
    ----------
    compute_strains : Bool (optional)
      Boolean variable to compute Strains and Stresses at nodes.
      By default it is False.
    plot_contours : Bool (optional)
      Boolean variable to plot contours of the computed variables.
      By default it is True.

    Returns
    -------
    UC : ndarray (nnodes, 2)
      Displacements at nodes.
    E_nodes : ndarray (nnodes, 3), optional
      Strains at nodes. It is returned when `compute_strains` is True.
    S_nodes : ndarray (nnodes, 3), optional
      Stresses at nodes. It is returned when `compute_strains` is True.

    """    
    folder = pre.initial_params()
    start_time = datetime.now()
    echo = False
    
    #%% PRE-PROCESSING
    nodes, mats, elements, loads = pre.readin(folder=folder)
    if echo:
        pre.echomod(nodes, mats, elements, loads, folder=folder)
    DME , IBC , neq = ass.DME(nodes, elements)
    print("Number of nodes: {}".format(nodes.shape[0]))
    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))
    
    #%% SYSTEM ASSEMBLY
    KG = ass.assembler(elements, mats, nodes, neq, DME)
    RHSG = ass.loadasem(loads, IBC, neq)
    
    #%% SYSTEM SOLUTION
    UG = sol.static_sol(KG, RHSG)
    if not(np.allclose(KG.dot(UG)/KG.max(), RHSG/KG.max())):
        print("The system is not in equilibrium!")
    end_time = datetime.now()
    print('Duration for system solution: {}'.format(end_time - start_time))
    
    #%% POST-PROCESSING
    start_time = datetime.now()
    UC = pos.complete_disp(IBC, nodes, UG)
    E_nodes, S_nodes = None, None
    if compute_strains:
        E_nodes, S_nodes = pos.strain_nodes(nodes , elements, mats, UC)
    if plot_contours:
        pos.fields_plot(elements, nodes, UC, E_nodes=E_nodes, S_nodes=S_nodes)
    end_time = datetime.now()
    print('Duration for post processing: {}'.format(end_time - start_time))
    print('Analysis terminated successfully!')
    return UC, E_nodes, S_nodes if compute_strains else UC

if __name__ == '__main__':
    solids_GUI()
    plt.show()