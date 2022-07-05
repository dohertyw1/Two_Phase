from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

mesh = Mesh('wedge20.xml')
fig = plt.figure(num=None, figsize=(12,16), dpi=100, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(2,1,1)
circle3 = plt.Circle((0, 0), 1, color='r', fill=False, lw=5)
ax1.add_patch(circle3)
plot(mesh, color='blue', linewidth=0.25)
ax1.set_xticks([-20,0,20])
ax1.set_xlim([-20,20])
ax1.set_ylim([0,2])
ax1.set_title('Single-Phase Viscoelastic Flow Schematic')



ax1 = fig.add_subplot(4,1,2)
circle3 = plt.Circle((0, 0), 0.98, color='r', fill=False, lw=3)
ax1.add_patch(circle3)
plot(mesh, color='blue', linewidth=0.25)
ax1.text(-3.9,0.9, r'$\Gamma_{L}$', fontsize = 26,color='red')
ax1.text(3.6,0.9, r'$\Gamma_{R}$', fontsize = 26,color='red')
ax1.text(-0.1,1.7, r'$\Gamma_{T}$', fontsize = 26,color='red')
ax1.text(-0.1,1.05, r'$\Gamma_{C}$', fontsize = 26,color='red')
ax1.text(-2.1,0.15, r'$\Gamma_{B}$', fontsize = 26,color='red')
ax1.text(1.9,0.15, r'$\Gamma_{B}$', fontsize = 26,color='red')

ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1.tick_params(axis='both', which='major', labelsize=13)
ax1.set_xlim([-4,4])
ax1.set_ylim([0,2])
ax1.set_xticks([-4,-2,0,2,4])
ax1.set_yticks([0,1,2])
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
fig.tight_layout()
