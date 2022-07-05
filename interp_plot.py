from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def point(point1,point2):
    x = [point1[0], point2[0]]
    y = [point1[1], point2[1]]
    return x, y

rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 22
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(num=None, figsize=(12, 6), dpi=100,
                 facecolor='w', edgecolor='k') 
ax1 = fig.add_subplot(1,2,1)
x=np.linspace(-1,1,3)
y = np.zeros(3)
xb=np.linspace(-1,1,101)
yb = np.zeros(101)

x1, y1 = point([-1,0], [0, 1])
x2, y2 = point([0,1], [1, 0])
x3, y3 = point([-1,1], [0, 0])
x4, y4 = point([0,0], [1, 1])

ax1.plot([10],[10], color='blue')
ax1.plot([10],[10], color='red',ls='-.')

ax1.plot(x1, y1,color='blue')
ax1.plot(x2, y2,color='blue')
ax1.plot(x3, y3,color='blue')
ax1.plot(x4, y4,color='blue')

ax1.plot(xb,(1+xb)*(1-xb), color='red',ls='-.')
ax1.plot(xb,(-xb/2)*(1-xb), color='red',ls='-.')
ax1.plot(xb,(xb/2)*(1+xb), color='red',ls='-.')

ax1.plot(point([-1,0], [-1, 1])[0],point([-1,0], [-1, 1])[1], ls=':', color='black')
ax1.plot(point([0,0], [0, 1])[0],point([-1,0], [-1, 1])[1], ls=':', color='black')
ax1.plot(point([1,0], [1, 1])[0],point([-1,0], [-1, 1])[1], ls=':', color='black')



ax1.plot(x,y, marker = 'o', ms = 7, mec = 'red', mfc = 'yellow', color='black')

ax1.set_xlabel(r'$x$')
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1.tick_params(axis='both', which='major', labelsize=13)
ax1.set_xticks([-1, 0, 1])
ax1.set_yticks([0,1])
ax1.set_ylabel(r'$y$')
ax1.set_ylim([-0.2, 1.1])
ax1.legend([r'$\textbf{CG1}$',r'$\textbf{CG2}$'], loc='upper right',prop={'size': 16})
ax1.set_xlim([-1.1,1.1])
# ax1.grid()
ax1.set_title(r'a)')
ax1 = fig.add_subplot(1,2,2)

ax1.plot(10,10, marker = 's', ms = 12,  mfc = 'red', mec='None', color='white')
ax1.plot(10,10, marker = 'o', ms = 10,  mfc = 'blue', mec='None', color='white')
ax1.plot(10,10, marker = 'D', ms = 10,  mfc = 'green', mec='None', color='white')


ax1.plot(point([0.25,0.25], [0.75, 0.25])[0],point([0.25,0.25], [0.75, 0.25])[1],color='black',linewidth=3)
ax1.plot(point([0.25,0.25], [0.25, 0.75])[0],point([0.25,0.25], [0.25, 0.75])[1],color='black',linewidth=3)
ax1.plot(point([0.25,0.75], [0.75, 0.25])[0],point([0.25,0.75], [0.75, 0.25])[1],color='black',linewidth=3)

ax1.plot(0.25,0.25, marker = 's', ms = 20,  mfc = 'red', mec='None')
ax1.plot(0.25,0.75, marker = 's', ms = 20,  mfc = 'red', mec='None')
ax1.plot(0.75,0.25, marker = 's', ms = 20,  mfc = 'red', mec='None')

ax1.plot(0.5,0.5, marker = 'o', ms = 15,  mfc = 'blue', mec='None')
ax1.plot(0.5,0.25, marker = 'o', ms = 15,  mfc = 'blue', mec='None')
ax1.plot(0.25,0.5, marker = 'o', ms = 15,  mfc = 'blue', mec='None')

ax1.plot(0.25,0.25, marker = 'o', ms = 15,  mfc = 'blue', mec='None')
ax1.plot(0.25,0.75, marker = 'o', ms = 15,  mfc = 'blue', mec='None')
ax1.plot(0.75,0.25, marker = 'o', ms = 15,  mfc = 'blue', mec='None')

ax1.plot(5/12,5/12, marker = 'D', ms = 15,  mfc = 'green', mec='None')
ax1.legend([r'$\textbf{CG1}$',r'$\textbf{CG2}$',r'$\textbf{DG0}$'], loc='upper right',prop={'size': 16})
ax1.set_xlim([0.23,0.77])
ax1.set_ylim([0.23,0.77])
ax1.set_title(r'b)')
ax1.set_yticks([])
ax1.set_xticks([])
plt.axis('off')
