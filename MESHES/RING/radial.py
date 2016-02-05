"""
Computes radial loads for the ring problem 

Last updated December 2015

"""
import numpy as np

nodes = np.loadtxt('nodloads.txt')
nn = len(nodes[:,0])
F = 1.0
F_comp = np.zeros([nn, 2])
x = nodes[:, 1]
y = nodes[:, 2]
theta = np.atan2(y, x)
F_comp[:, 0] = F*np.cos(theta)
F_comp[:, 1] = F*np.sin(theta)
np.savetxt("KNODES.txt", F_comp, fmt='%5.2f', delimiter=' ')