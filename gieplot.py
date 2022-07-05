from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

etap = 0.5
lamb = 1

def n2(etap, lamb, alpha, eps):
    func = -(etap/(2*alpha*lamb))*(1-np.sqrt(1-4*alpha*(1-alpha)))*eps
    return func

nn2 = n2(etap, lamb, 0.5, np.linspace(10e-4, 10e8, 10000))

rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(num=None, figsize=(5,10), dpi=100, facecolor='w', edgecolor='k')
# ax1 = fig.add_subplot(1,1,1)

plt.plot(nn2)
plt.xlim([10e-4,10e8])
plt.ylim([-0.5,0])