# -*- coding: utf-8 -*-
"""
This module contains functions to preprocess the input files to compute
a Finite Element Analysis.

"""
from __future__ import division, print_function
import numpy as np


def readin(folder=""):
    """Read the input file"""
    nodes = np.loadtxt(folder + 'nodes.txt')
    mats = np.loadtxt(folder + 'mater.txt')
    elements = np.loadtxt(folder + 'eles.txt')
    loads = np.loadtxt(folder + 'loads.txt')

    return nodes, mats, elements, loads

def proini(nodes, mats, elements, loads):
    """Extract problem parameters and nodal coordinates"""
    ne = len(elements[:, 0])
    nn = len(nodes[:, 0])
    nm = len(mats)
    nl = len(loads[:, 0])

    return ne, nn, nm, nl


def echomod(nodes, mats, elements, loads, folder=""):
    """Create echoes of the model input files"""
    np.savetxt(folder + "KNODES.txt", nodes, fmt='%5.2f', delimiter=' ')
    np.savetxt(folder + "KMATES.txt", mats, fmt='%5.2f', delimiter=' ')
    np.savetxt(folder + "KELEMS.txt", elements, fmt='%5.2f', delimiter=' ')
    np.savetxt(folder + "KLOADS.txt", loads, fmt='%5.2f', delimiter=' ')
