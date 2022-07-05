from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

with open('Cons_Viscoelastic_CG_NonDim_g0.1/Cons_Viscoelastic_CG_NonDim_g0.1_data.csv', newline='') as csvfile: #70
    datanr1601 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('Cons_Viscoelastic_CG_NonDim_g0.5/Cons_Viscoelastic_CG_NonDim_g0.5_data.csv', newline='') as csvfile: #70
    datanr1602 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('Cons_Viscoelastic_CG_NonDim_g1/Cons_Viscoelastic_CG_NonDim_g1_data.csv', newline='') as csvfile: #70
    datanr1603 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('Cons_Viscoelastic_CG_NonDim_g1/Cons_Viscoelastic_CG_NonDim_g1_data.csv', newline='') as csvfile: #70
    datanr1604 = np.array(list(csv.reader(csvfile, delimiter='\t')))


timescale1601 = []
area1601 = []
xcom1601 = []
ycom1601 = []
circ1601=[]
urise1601=[]
vrise1601=[]

for i in range(np.shape(datanr1601)[0]):
    timescale1601.append(float(datanr1601[i][0]))

for i in range(np.shape(datanr1601)[0]):
    area1601.append(float(datanr1601[i][1]))

for i in range(np.shape(datanr1601)[0]):
    xcom1601.append(float(datanr1601[i][2]))

for i in range(np.shape(datanr1601)[0]):
    ycom1601.append(float(datanr1601[i][3]))

for i in range(np.shape(datanr1601)[0]):
    circ1601.append(float(datanr1601[i][4]))

for i in range(np.shape(datanr1601)[0]):
    urise1601.append(float(datanr1601[i][5]))

for i in range(np.shape(datanr1601)[0]):
    vrise1601.append(float(datanr1601[i][6]))

timescale1602 = []
area1602 = []
xcom1602 = []
ycom1602 = []
circ1602=[]
urise1602=[]
vrise1602=[]

for i in range(np.shape(datanr1602)[0]):
    timescale1602.append(float(datanr1602[i][0]))

for i in range(np.shape(datanr1602)[0]):
    area1602.append(float(datanr1602[i][1]))

for i in range(np.shape(datanr1602)[0]):
    xcom1602.append(float(datanr1602[i][2]))

for i in range(np.shape(datanr1602)[0]):
    ycom1602.append(float(datanr1602[i][3]))

for i in range(np.shape(datanr1602)[0]):
    circ1602.append(float(datanr1602[i][4]))

for i in range(np.shape(datanr1602)[0]):
    urise1602.append(float(datanr1602[i][5]))

for i in range(np.shape(datanr1602)[0]):
    vrise1602.append(float(datanr1602[i][6]))

timescale1603 = []
area1603 = []
xcom1603 = []
ycom1603 = []
circ1603=[]
urise1603=[]
vrise1603=[]

timescale1604 = []
area1604 = []
xcom1604 = []
ycom1604 = []
circ1604=[]
urise1604=[]
vrise1604=[]

for i in range(np.shape(datanr1604)[0]):
    timescale1604.append(float(datanr1604[i][0]))

for i in range(np.shape(datanr1604)[0]):
    area1604.append(float(datanr1604[i][1]))

for i in range(np.shape(datanr1604)[0]):
    xcom1604.append(float(datanr1604[i][2]))

for i in range(np.shape(datanr1604)[0]):
    ycom1604.append(float(datanr1604[i][3]))

for i in range(np.shape(datanr1604)[0]):
    circ1604.append(float(datanr1604[i][4]))

for i in range(np.shape(datanr1604)[0]):
    urise1604.append(float(datanr1604[i][5]))

for i in range(np.shape(datanr1604)[0]):
    vrise1604.append(float(datanr1604[i][6]))


for i in range(np.shape(datanr1603)[0]):
    timescale1603.append(float(datanr1603[i][0]))

for i in range(np.shape(datanr1603)[0]):
    area1603.append(float(datanr1603[i][1]))

for i in range(np.shape(datanr1603)[0]):
    xcom1603.append(float(datanr1603[i][2]))

for i in range(np.shape(datanr1603)[0]):
    ycom1603.append(float(datanr1603[i][3]))

for i in range(np.shape(datanr1603)[0]):
    circ1603.append(float(datanr1603[i][4]))

for i in range(np.shape(datanr1603)[0]):
    urise1603.append(float(datanr1603[i][5]))

for i in range(np.shape(datanr1603)[0]):
    vrise1603.append(float(datanr1603[i][6]))

mesh = RectangleMesh(Point(0, 0), Point(2, 4), 80, 160, 'crossed')
V = FunctionSpace(mesh, 'CG', 2)
phi_651 = Function(V)
phi_652 = Function(V)
phi_653 = Function(V)

with XDMFFile('Cons_Viscoelastic_CG_NonDim_g0.1/frames/frame_10.1036/phi_read_10.1036.xdmf') as infile:
    infile.read_checkpoint(phi_651, "phi")
with XDMFFile('Cons_Viscoelastic_CG_NonDim_g0.5/frames/frame_10.1036/phi_read_10.1036.xdmf') as infile:
    infile.read_checkpoint(phi_652, "phi")
with XDMFFile('Cons_Viscoelastic_CG_NonDim_g1/frames/frame_10.1036/phi_read_10.1036.xdmf') as infile:
    infile.read_checkpoint(phi_653, "phi")

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]


fig = plt.figure(num=None, figsize=(12,10), dpi=100, facecolor='w', edgecolor='k')

ax1= fig.add_subplot(3,2,1)

ax1.plot(timescale1601, ycom1601)
ax1.plot(timescale1602, ycom1602)
ax1.plot(timescale1603, ycom1603)
# ax1.plot(timescale1604, ycom1604)

# ax1.set_xlim([0,3])
# ax1.set_ylim([0.4,1.2])
ax1.set_xlabel('t',**{'fontname':'Times', 'size':'18'})
ax1.set_ylabel('Centre of Mass',**{'fontname':'Times', 'size':'18'})
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1.tick_params(axis='both', which='major', labelsize=13)

ax1.grid()
# ax1.set_title('Centre of Mass')
ax1.legend([r'$\alpha = 0.1$',r'$\alpha = 0.5$',r'$\alpha = 1$'],loc='lower right',frameon=True,prop={'size': 12})

ax5= fig.add_subplot(3,2,3)

ax5.plot(timescale1601, vrise1601)
ax5.plot(timescale1602, vrise1602)
ax5.plot(timescale1603, vrise1603)
# ax5.plot(timescale1604, vrise1604)

# ax5.set_xlim([0,3])
# ax5.set_ylim([0,0.3])
ax5.set_xlabel('t',**{'fontname':'Times', 'size':'18'})
ax5.set_ylabel('Rise Velocity',**{'fontname':'Times', 'size':'18'})
ax5.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax5.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax5.tick_params(axis='both', which='major', labelsize=13)

ax5.grid()
# ax5.set_title('Rise Velocity')
ax5.legend([r'$\alpha = 0.1$',r'$\alpha = 0.5$',r'$\alpha = 1$'],loc='upper right',frameon=True,prop={'size': 12})

axc= fig.add_subplot(3,2,5)

axc.plot(timescale1601, circ1601)
axc.plot(timescale1602, circ1602)
axc.plot(timescale1603, circ1603)
# axc.plot(timescale1604, circ1604)

# axc.set_xlim([0,3])
axc.set_xlabel('t',**{'fontname':'Times', 'size':'18'})
axc.set_ylabel('Circularity',**{'fontname':'Times', 'size':'18'})
axc.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
axc.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
axc.tick_params(axis='both', which='major', labelsize=13)

axc.grid()
# axc.set_title('Circularity')
axc.legend([r'$\alpha = 0.1$',r'$\alpha = 0.5$',r'$\alpha = 1$'],loc='lower left',frameon=True,prop={'size': 12})

ax1= fig.add_subplot(1,2,2)
plot(phi_651, levels=[0.5], linewidths=3,mode='contour', colors = '#1f77b4')
plot(phi_652, levels=[0.5], linewidths=3,mode='contour', colors = '#ff7f0e')
plot(phi_653, levels=[0.5], linewidths=3,mode='contour', colors = '#2ca02c')
plt.plot([1],[4],color='#1f77b4')
plt.plot([1],[5],color='#ff7f0e')
plt.plot([1],[6],color='#2ca02c')
ax1.set_xlim([0.5,1.5])
ax1.set_ylim([1,2.75])
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Bubble Shape')
ax1.grid()
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1.tick_params(axis='both', which='major', labelsize=13)
ax1.legend([r'$\alpha = 0.1$',r'$\alpha = 0.5$',r'$\alpha = 1$'],loc='lower left',frameon=True,prop={'size': 15})
