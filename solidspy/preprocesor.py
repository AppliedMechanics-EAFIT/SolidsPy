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