# -*- coding: utf-8 -*-
"""
This module contains functions to preprocess the input files to compute
a Finite Element Analysis.

"""
from __future__ import division, print_function
import sys
import numpy as np


def readin(folder=""):
    """Read the input files"""
    nodes = np.loadtxt(folder + 'nodes.txt', ndmin=2)
    mats = np.loadtxt(folder + 'mater.txt', ndmin=2)
    elements = np.loadtxt(folder + 'eles.txt', ndmin=2, dtype=np.int)
    loads = np.loadtxt(folder + 'loads.txt', ndmin=2)

    return nodes, mats, elements, loads


def echomod(nodes, mats, elements, loads, folder=""):
    """Create echoes of the model input files"""
    np.savetxt(folder + "KNODES.txt", nodes, fmt='%5.2f', delimiter=' ')
    np.savetxt(folder + "KMATES.txt", mats, fmt='%5.2f', delimiter=' ')
    np.savetxt(folder + "KELEMS.txt", elements, fmt='%d', delimiter=' ')
    np.savetxt(folder + "KLOADS.txt", loads, fmt='%5.2f', delimiter=' ')


def initial_params():
    """Read initial parameters for the simulation
    
    The parameters to be read are:

    - folder: location of the input files.
    - name: name for the output files (if echo is True).
    - echo: echo output files.
    """
    # Check Python version
    version = sys.version_info.major
    if version == 3:
        global raw_input
        raw_input = input
    elif version == 2:
        pass
    else:
        raise ValueError("You should use Python 2.x at least!")

    # Try to run with easygui
    try:
        import easygui
        folder = easygui.diropenbox(title="Folder for the job") + "/"
    except:
        folder = raw_input('Enter folder (empty for the current one): ')

    return folder


def ele_writer(cells , cell_data , ele_tag , phy_sur , ele_type , mat_tag , nini):
    """
    Extracts a subset of elements from a complete mesh according to the physical surface
    phy_sur and writes down the proper fields into an elements array.
    INPUT PARAMTERS:
    ---------------
        cell and cell_data: Are the dictionaries creatd by meshio.
        ele_tag : String defining the element type according to meshio (e.g., quad9 , line3, etc).
        phy_sur : Integer defining the physical surface for the subset.
        ele_type: Integer defining the element type according to WAVES.
        mat_tag : Integer defining the material profile for the subset.
        ndof    : Integer defining the number of degrees of freedom for the elements.
        nnode   : Integer defining the number of nodes for the element.
        nini   : Integer defining the element id for the first element in the set.
    OUTPUT PARAMTERS:
    ----------------
        nf        : Integer defining the element id for the last element in the set
        els_array : Integer array with the elemental data according to WAVES.
    """
    eles = cells[ele_tag]           # Element connectivities (adds 1 to fortranize)
    dict_nnode = {'triangle': 3 , 'triangle6':6 }
    nnode = dict_nnode[ele_tag]
    phy_surface = cell_data[ele_tag]['physical']
    ele_id = [cont for cont, _ in enumerate(phy_surface[:]) if phy_surface[cont] == phy_sur]    
    els_array = np.zeros([len(ele_id) , 3 + nnode], dtype=int)
    els_array[: , 0] = range(nini , len(ele_id) + nini )
    els_array[: , 1] = ele_type
    els_array[: , 2] = mat_tag
    els_array[: , 3::] = eles[ele_id, :]
    nf = nini + len(ele_id)
    
    return nf , els_array

def node_writer(points , point_data):
    """
    Writes down the nodal data as required by SOLIDSpy.
    
    INPUT PARAMETERS
    ----------------    
    points : Dictionary
        Python dictionary storing the nodal points 
    point_data : Dictionary.
        Python dictionary with physical data associatted to the nodes.
        
    OUTPUT PARAMTERS
    ----------------
    nodes_array : ndarray (int)
        Integer array with the nodal data according to SOLIDSpy.
        
    """
#
    nodes_array = np.zeros([points.shape[0], 5])
    nodes_array[:, 0] = range(points.shape[0])
    nodes_array[:, 1:3] = points[:, :2]    
    
    return nodes_array
#
def boundary_conditions(cells , cell_data , phy_lin , nodes_array , bc_x , bc_y):
    """
   Imposes the nodal point boundary conditions as required by SOLIDSpy.
        INPUT PARAMTERS:
    ---------------
        cell and cell_data: Are the dictionaries creatd by meshio.
        phy_lin     : Integer defining the physical line where BCs are to be imposed.
        nodes_array : Integer array with the nodal data and to be modified by BCs.
        bc_x, bc_y  : Boundary condition flag along the x and y direction (-1: restraind; 0:free)
    OUTPUT PARAMTERS: 
    ----------------
        nodes_array : Integer array with the nodal data after imposing BCs according to SOLIDSpy.
    """
#
    lines = cells["line"]
    phy_line = cell_data["line"]["physical"]              # Bounds contains data corresponding to the physical line.
    id_frontera = [cont for cont in range(len(phy_line)) if phy_line[cont] == phy_lin]
    nodes_frontera = lines[id_frontera]
    nodes_frontera = nodes_frontera.flatten()
    nodes_frontera = list(set(nodes_frontera))
    nodes_array[nodes_frontera, 3] = bc_x
    nodes_array[nodes_frontera, 4] = bc_y 
    
    return nodes_array
#
def loading(cells , cell_data , phy_lin , P_x , P_y):
    """
   Imposes the nodal point boundary conditions as required by SOLIDSpy.
        INPUT PARAMTERS:
    ---------------
        cell and cell_data: Are the dictionaries creatd by meshio.
        phy_lin     : Integer defining the physical line where BCs are to be imposed.
        nodes_array : Integer array with the nodal data and to be modified by BCs.
        bc_x, bc_y  : Boundary condition flag along the x and y direction (-1: restraind; 0:free)
    OUTPUT PARAMTERS: 
    ----------------
        nodes_array : Integer array with the nodal data after imposing BCs according to SOLIDSpy.
    """
#
    lines = cells["line"]
    phy_line = cell_data["line"]["physical"]              # Bounds contains data corresponding to the physical line.
    id_carga = [cont for cont in range(len(phy_line)) if phy_line[cont] == 500]
    nodes_carga = lines[id_carga]
    nodes_carga = nodes_carga.flatten()
    nodes_carga = list(set(nodes_carga))
    ncargas = len(nodes_carga)
    cargas = np.zeros((ncargas, 3))
    cargas[:, 0] = nodes_carga
    cargas[:, 1] = P_x/ncargas
    cargas[:, 2] = P_y/ncargas   
#
    
    return cargas