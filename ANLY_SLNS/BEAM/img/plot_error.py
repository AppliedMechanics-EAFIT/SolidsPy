from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['image.cmap'] = "YlGnBu_r"

data = np.loadtxt("error_vs_h.txt", skiprows=1, delimiter=",")
els = data[:, 1]
h = data[:, 2]
err = data[:, 3]
fig = plt.figure(figsize=(8, 5))
ax1 = fig.add_subplot(111)
ax1.loglog(h, err, '-bo')
plt.xlabel(r"Element size: $h$", fontsize=14)
plt.ylabel(r"Relative error: $\frac{\Vert u - u_h \Vert}{\Vert u \Vert}$",
           fontsize=14)
xticks, xlabels = plt.xticks()
ax2 = ax1.twiny()
ax2.set_xscale("log")
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(xticks[1:-1])
labels = 3/xticks[1:-1]**2
ax2.set_xticklabels(["%g"%label for label in labels])
plt.xlabel("Number of elements")
plt.savefig("Beam_convergence.pdf", bbox_inches="tight")