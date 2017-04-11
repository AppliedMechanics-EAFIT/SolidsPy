"""
Example that assemble the stiffness matrix using a user element
function

"""
from __future__ import division, print_function
from os import sys
sys.path.append("../MAIN/")
import numpy as np
import matplotlib.pyplot as plt
import assemutil as ass


def uel_ones(elcoord, par1, par0):
    return np.ones((8, 8)), 8, 1 


nodes = np.zeros((9, 5))
nodes[:, 0] = range(0, 9)
nodes[:, 1:3] = np.array([
        [0, 0],
        [1, 0],
        [2, 0],
        [0, 1],
        [1, 1],
        [2, 1],
        [0, 2],
        [1, 2],
        [2, 2]])
elements = np.ones((4, 7), dtype=np.int)
elements[:, 0] = range(0, 4)
elements[:, 2] = 0
elements[:, 3:] = np.array([
        [0, 1, 4, 3],
        [1, 2, 5, 4],
        [3, 4, 7, 6],
        [4, 5, 8, 7]])
mats = np.array([[1, 0.3]])
DME , IBC , neq = ass.DME(nodes, elements)
KG = ass.assembler(elements, mats, nodes, neq, DME, sparse=False,
                   uel=uel_ones)
plt.spy(KG)
plt.show()