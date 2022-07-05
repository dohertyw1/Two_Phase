from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def txx(alpha, eps):
    func = (-(1-2*eps)+np.sqrt( ( (1-2*eps)**2) + 8*alpha*eps ))/(4*alpha)
    return func

def tyy(alpha, eps):
    func = (-(1+eps)+np.sqrt(((1+eps)**2)-4*alpha*eps))/(4*alpha)
    return func

def ext(alpha, eps):
    return 3/2+((txx(alpha,eps)-tyy(alpha,eps))/eps)

def obext(eps):
    func = 3/2+(3/2)/((1-2*eps)*(1+eps))
    return func

a=10e-3
b=10e2

rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 22
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(num=None, figsize=(10,10), dpi=100, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(1,1,1)
ax1.plot(np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5]))),obext(np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5])))),lw=4,linestyle='dashed')
ax1.plot(np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5]))),ext(0.1,np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5])))),lw=4)
ax1.plot(np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5]))),ext(0.25,np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5])))),lw=4)
ax1.plot(np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5]))),ext(0.5,np.concatenate((np.linspace(a,b,1000000),np.array([10e3,10e4,10e5])))),lw=4)
ax1.grid()
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1.tick_params(axis='both', which='major', labelsize=13)
ax1.set_ylim([-2,12])
ax1.set_xscale('log')
ax1.set_xlim([10e-2,10e5])
ax1.set_yticks([0,2,4,6,8,10,12])
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

ax1.set_ylabel(r'$\eta_E$')
ax1.set_xlabel(r'$\dot{\epsilon}$')
ax1.legend([r'$\alpha = 0\:\:(\textrm{OB})$',r'$\alpha = 0.1$',r'$\alpha = 0.25$',r'$\alpha = 0.5$'],loc='upper right',bbox_to_anchor=(1,0.9),frameon=True,prop={'size': 16})
ax1.set_title('Extensional Viscosity')