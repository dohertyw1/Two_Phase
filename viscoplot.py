from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# with open('/media/wdoherty1/Samsung_T5/mass_frac/Cons_Viscoelastic_CG_Dim_c=6.5/Cons_Viscoelastic_CG_Dim_c=6.5_data.csv', newline='') as csvfile: #70
#     datanr1601 = np.array(list(csv.reader(csvfile, delimiter='\t')))

# with open('/media/wdoherty1/Samsung_T5/mass_frac/Cons_Viscoelastic_CG_Dim_c=13/Cons_Viscoelastic_CG_Dim_c=6.5_data.csv', newline='') as csvfile: #70
#     datanr1602 = np.array(list(csv.reader(csvfile, delimiter='\t')))

# with open('/media/wdoherty1/Samsung_T5/mass_frac/Cons_Viscoelastic_CG_Dim_c=26/Cons_Viscoelastic_CG_Dim_c=6.5_data.csv', newline='') as csvfile: #70
#     datanr1603 = np.array(list(csv.reader(csvfile, delimiter='\t')))


# timescale1601 = []
# area1601 = []
# xcom1601 = []
# ycom1601 = []
# circ1601=[]
# urise1601=[]
# vrise1601=[]

# for i in range(np.shape(datanr1601)[0]):
#     timescale1601.append(float(datanr1601[i][0]))

# for i in range(np.shape(datanr1601)[0]):
#     area1601.append(float(datanr1601[i][1]))

# for i in range(np.shape(datanr1601)[0]):
#     xcom1601.append(float(datanr1601[i][2]))

# for i in range(np.shape(datanr1601)[0]):
#     ycom1601.append(float(datanr1601[i][3]))

# for i in range(np.shape(datanr1601)[0]):
#     circ1601.append(float(datanr1601[i][4]))

# for i in range(np.shape(datanr1601)[0]):
#     urise1601.append(float(datanr1601[i][5]))

# for i in range(np.shape(datanr1601)[0]):
#     vrise1601.append(float(datanr1601[i][6]))

# timescale1602 = []
# area1602 = []
# xcom1602 = []
# ycom1602 = []
# circ1602=[]
# urise1602=[]
# vrise1602=[]

# for i in range(np.shape(datanr1602)[0]):
#     timescale1602.append(float(datanr1602[i][0]))

# for i in range(np.shape(datanr1602)[0]):
#     area1602.append(float(datanr1602[i][1]))

# for i in range(np.shape(datanr1602)[0]):
#     xcom1602.append(float(datanr1602[i][2]))

# for i in range(np.shape(datanr1602)[0]):
#     ycom1602.append(float(datanr1602[i][3]))

# for i in range(np.shape(datanr1602)[0]):
#     circ1602.append(float(datanr1602[i][4]))

# for i in range(np.shape(datanr1602)[0]):
#     urise1602.append(float(datanr1602[i][5]))

# for i in range(np.shape(datanr1602)[0]):
#     vrise1602.append(float(datanr1602[i][6]))

# timescale1603 = []
# area1603 = []
# xcom1603 = []
# ycom1603 = []
# circ1603=[]
# urise1603=[]
# vrise1603=[]

# for i in range(np.shape(datanr1603)[0]):
#     timescale1603.append(float(datanr1603[i][0]))

# for i in range(np.shape(datanr1603)[0]):
#     area1603.append(float(datanr1603[i][1]))

# for i in range(np.shape(datanr1603)[0]):
#     xcom1603.append(float(datanr1603[i][2]))

# for i in range(np.shape(datanr1603)[0]):
#     ycom1603.append(float(datanr1603[i][3]))

# for i in range(np.shape(datanr1603)[0]):
#     circ1603.append(float(datanr1603[i][4]))

# for i in range(np.shape(datanr1603)[0]):
#     urise1603.append(float(datanr1603[i][5]))

# for i in range(np.shape(datanr1603)[0]):
#     vrise1603.append(float(datanr1603[i][6]))

mesh = RectangleMesh(Point(0, 0), Point(2, 4), 200, 400, 'crossed')
V = FunctionSpace(mesh, 'CG', 2)
Vvec = VectorFunctionSpace(mesh, 'CG', 2)
Vten = TensorFunctionSpace(mesh, 'DG', 0)
phi_65 = Function(V)
u_65=Function(Vvec)
tau_65 = Function(Vten)
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=6.5_T0.13/frames/frame_0.13/phi_read_0.13.xdmf') as infile:
    infile.read_checkpoint(phi_65, "phi")
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=6.5_T0.13/frames/frame_0.13/u_read_0.13.xdmf') as infile:
    infile.read_checkpoint(u_65, "u")
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=6.5_T0.13/frames/frame_0.13/tau_read_0.13.xdmf') as infile:
    infile.read_checkpoint(tau_65, "tau")

phi_13 = Function(V)
u_13=Function(Vvec)
tau_13 = Function(Vten)
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=13_T0.13/frames/frame_0.13/phi_read_0.13.xdmf') as infile:
    infile.read_checkpoint(phi_13, "phi")
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=13_T0.13/frames/frame_0.13/u_read_0.13.xdmf') as infile:
    infile.read_checkpoint(u_13, "u")
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=13_T0.13/frames/frame_0.13/tau_read_0.13.xdmf') as infile:
    infile.read_checkpoint(tau_13, "tau")

phi_26 = Function(V)
u_26=Function(Vvec)
tau_26 = Function(Vten)
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=26_T0.13/frames/frame_0.13/phi_read_0.13.xdmf') as infile:
    infile.read_checkpoint(phi_26, "phi")
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=26_T0.13/frames/frame_0.13/u_read_0.13.xdmf') as infile:
    infile.read_checkpoint(u_26, "u")
with XDMFFile('Cons_Viscoelastic_CG_Dim_c=26_T0.13/frames/frame_0.13/tau_read_0.13.xdmf') as infile:
    infile.read_checkpoint(tau_26, "tau")

# u_65_wake = Function(Vvec)
# u_26_wake = Function(Vvec)
# with XDMFFile('Cons_Viscoelastic_CG_Dim_c=6.5/frames/frame_0.201/u_read_0.201.xdmf') as infile:
#     infile.read_checkpoint(u_65_wake, "u")
# with XDMFFile('Cons_Viscoelastic_CG_Dim_c=26/frames/frame_0.201/u_read_0.201.xdmf') as infile:
#     infile.read_checkpoint(u_26_wake, "u")


# mesh = RectangleMesh(Point(0, 0), Point(2, 4), 50, 100, 'crossed')
# V = FunctionSpace(mesh, 'CG', 2)
# phi200 = Function(V)
# with XDMFFile('phi_read.xdmf') as infile:
#     infile.read_checkpoint(phi200, "phi")
# mesh = RectangleMesh(Point(0, 0), Point(2, 4), 25, 50, 'crossed')
# V = FunctionSpace(mesh, 'CG', 2)
# phi300 = Function(V)
# with XDMFFile('Cons_Viscoelastic_CG_Dim_c=26/phi_read.xdmf') as infile:
#     infile.read_checkpoint(phi300, "phi")
# mesh = RectangleMesh(Point(0.5,0),Point(1,2),160,640,'crossed')
# V = FunctionSpace(mesh, 'CG', 2)
# phi400 = Function(V)
# with XDMFFile('Cons_Newtonian_CG_Dim_test2400/phi_read.xdmf') as infile:
#     infile.read_checkpoint(phi400, "phi")
# mesh = RectangleMesh(Point(0,0),Point(2,4),50,100,'crossed')
# Vu = VectorFunctionSpace(mesh, 'CG', 2)
# u200 = Function(Vu)
# with XDMFFile('Cons_Viscoelastic_CG_Dim_c=6.5/u_read.xdmf') as infile:
#     infile.read_checkpoint(u200, "u")

# mesh = RectangleMesh(Point(0,0),Point(2,4),200,400,'crossed')
# V = FunctionSpace(mesh, 'CG', 2)
# phi300 = Function(V)
# with XDMFFile('/media/wdoherty1/Samsung_T5/mass_frac/Cons_Viscoelastic_CG_Dim_c=13/phi_read.xdmf') as infile:
#     infile.read_checkpoint(phi300, "phi")

# mesh = RectangleMesh(Point(0,0),Point(2,4),200,400,'crossed')
# V = FunctionSpace(mesh, 'CG', 2)
# phi400 = Function(V)
# with XDMFFile('/media/wdoherty1/Samsung_T5/mass_frac/Cons_Viscoelastic_CG_Dim_c=26/phi_read.xdmf') as infile:
#     infile.read_checkpoint(phi400, "phi")


# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(num=None, figsize=(20, 20), dpi=100,
                 facecolor='w', edgecolor='k') 

ax1 = fig.add_subplot(3, 4, 1)
ax1.set_title(r'$1\textrm{A}:\:\textbf{u}\:(6.5\:wt.\%)$')
plot(phi_65, levels=[0.5], linewidths=3,mode='contour', colors = 'black')
plot(sqrt(u_65**2), cmap='jet')
plt.colorbar(plot(sqrt(u_65**2), cmap='jet', vmin = 0, vmax = 3.6),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax2 = fig.add_subplot(3, 4, 2)
ax2.set_title(r'$1\textrm{B}:\:\boldsymbol{\tau}_{xx}\:(6.5\:wt.\%)$')
plot(phi_65, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
plot(tau_65[0,0], cmap='jet')
plt.colorbar(plot(tau_65[0,0], cmap='jet', vmin = 0, vmax = 480),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax3 = fig.add_subplot(3, 4, 3)
ax3.set_title(r'$1\textrm{C}:\:\boldsymbol{\tau}_{xy}\:(6.5\:wt.\%)$')
plot(phi_65, levels=[0.5], linewidths=3,mode='contour', colors = 'black')
plot(tau_65[1,0], cmap='jet')
plt.colorbar(plot(tau_65[1,0], cmap='jet', vmin = -200, vmax = 200),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax4 = fig.add_subplot(3, 4, 4)
ax4.set_title(r'$1\textrm{D}:\:\boldsymbol{\tau}_{yy}\:(6.5\:wt.\%)$')
plot(phi_65, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
plot(tau_65[1,1], cmap='jet')
cbar = plt.colorbar(plot(tau_65[1,1], cmap='jet', vmin = 0, vmax = 640),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax5 = fig.add_subplot(3, 4, 5)
ax5.set_title(r'$2\textrm{A}:\:\textbf{u}\:(13\:wt.\%)$')
plot(phi_13, levels=[0.5], linewidths=3,mode='contour', colors = 'black')
plot(sqrt(u_13**2), cmap='jet')
plt.colorbar(plot(sqrt(u_13**2), cmap='jet', vmin = 0, vmax = 3.6),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax6 = fig.add_subplot(3, 4, 6)
ax6.set_title(r'$2\textrm{B}:\:\boldsymbol{\tau}_{xx}\:(13\:wt.\%)$')
plot(phi_13, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
plot(tau_13[0,0], cmap='jet')
plt.colorbar(plot(tau_13[0,0], cmap='jet', vmin = 0, vmax = 480),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax7 = fig.add_subplot(3, 4, 7)
ax7.set_title(r'$2\textrm{C}:\:\boldsymbol{\tau}_{xy}\:(13\:wt.\%)$')
plot(phi_13, levels=[0.5], linewidths=3,mode='contour', colors = 'black')
plt.colorbar(plot(tau_13[1,0], cmap='jet',vmin = -200, vmax = 200),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax8 = fig.add_subplot(3, 4, 8)
ax8.set_title(r'$2\textrm{D}:\:\boldsymbol{\tau}_{yy}\:(13\:wt.\%)$')
plot(phi_13, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
cbar = plt.colorbar(plot(tau_13[1,1], cmap='jet', vmin = 0, vmax = 640),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax9 = fig.add_subplot(3, 4, 9)
ax9.set_title(r'$3\textrm{A}:\:\textbf{u}\:(26\:wt.\%)$')
plot(phi_26, levels=[0.5], linewidths=3,mode='contour', colors = 'black')
plot(sqrt(u_26**2), cmap='jet')
plt.colorbar(plot(sqrt(u_26**2), cmap='jet', vmin = 0, vmax = 3.6),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax10 = fig.add_subplot(3, 4, 10)
ax10.set_title(r'$3\textrm{B}:\:\boldsymbol{\tau}_{xx}\:(26\:wt.\%)$')
plot(phi_26, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
plot(tau_26[0,0], cmap='jet')
plt.colorbar(plot(tau_26[0,0], cmap='jet', vmin = 0, vmax = 480),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax11 = fig.add_subplot(3, 4, 11)
ax11.set_title(r'$3\textrm{C}:\:\boldsymbol{\tau}_{xy}\:(26\:wt.\%)$')
plot(phi_26, levels=[0.5], linewidths=3,mode='contour', colors = 'black')
plot(tau_26[1,0], cmap='jet')
plt.colorbar(plot(tau_26[1,0], cmap='jet',vmin = -200, vmax = 200),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

ax12 = fig.add_subplot(3, 4, 12)
ax12.set_title(r'$3\textrm{D}:\:\boldsymbol{\tau}_{yy}\:(26\:wt.\%)$')
plot(phi_26, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
plot(tau_26[1,1], cmap='jet')
cbar = plt.colorbar(plot(tau_26[1,1], cmap='jet', vmin = 0, vmax = 640),fraction=0.09, pad=0.04)
plt.xlim([0.5,1.5])
plt.ylim([0.5,2.5])

plt.subplots_adjust(wspace=0)
plt.subplots_adjust(hspace=0)
fig.tight_layout()



# ax1 = fig.add_subplot(1,5,5)
# plot(phi100, mode='contour', levels=[0.5],colors='#1f77b4')
# plot(phi200, mode='contour', levels=[0.5],colors='#ff7f0e')
# plot(phi300, mode='contour', levels=[0.5],colors='#2ca02c')
# plot(phi400, mode='contour', levels=[0.5],colors='#d62728')
# ax1.legend(['dt=0.01','dt=0.005','dt=0.0025','dt=0.00124'],loc='upper right')
# ax1.set_xlim([0.5,1])
# ax1.set_ylim([0.5, 1.5])

# ax7.grid()

# ax7.set_xlim([0.5,1.5])
# ax7.set_ylim([1.25, 3.25])
# ax7.set_xlabel('x',**{'fontname':'Times', 'size':'18'})
# ax7.set_ylabel('y',**{'fontname':'Times', 'size':'18'})
# ax7.tick_params(axis='both', which='major', labelsize=13)

# ax7.set_title('Bubble Shape',**{'fontname':'Times', 'size':'18'})


# ax_ins = inset_axes(ax7, width="100%", height="100%", loc='upper left',
#                    bbox_to_anchor=(0.44,0.68,.45,.45), bbox_transform=ax7.transAxes)

# # plot(phi100, mode='contour',levels=[0.5],colors='#1f77b4')
# plot(phi200, mode='contour', levels=[0.5],colors='#1f77b4')
# plot(phi300, mode='contour', levels=[0.5],colors='#ff7f0e')
# plot(phi400, mode='contour', levels=[0.5],colors='#2ca02c')
# ax_ins.scatter(fMooNMD_shapex,fMooNMD_shapey,s=3, color='red')
# ax_ins.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey,s=3, color='purple')

# ax_ins.grid()
# ax_ins.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
# ax_ins.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')

# ax_ins.set_xlim([0.55,0.6])
# ax_ins.set_ylim([0.99, 1.025])
# ax_ins.set_xticks([0.55,0.575,0.6])
# ax_ins.tick_params(axis='both', which='major', labelsize=14)

# ax_ins2 = inset_axes(ax7, width="100%", height="100%", loc='upper left',
#                    bbox_to_anchor=(0.09,0.09,.35,.35), bbox_transform=ax7.transAxes)

# # plot(phi100, mode='contour',levels=[0.5],colors='#1f77b4')
# plot(phi200, mode='contour', levels=[0.5],colors='#1f77b4')
# plot(phi300, mode='contour', levels=[0.5],colors='#ff7f0e')
# plot(phi400, mode='contour', levels=[0.5],colors='#2ca02c')
# ax_ins2.scatter(fMooNMD_shapex,fMooNMD_shapey,s=3, color='red')
# ax_ins2.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey,s=3, color='purple')
# ax_ins2.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
# ax_ins2.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
# ax_ins2.grid()

# ax_ins2.set_xlim([0.76,0.85])
# ax_ins2.set_ylim([0.57, 0.77])
# plt.subplots_adjust(wspace=0)

# ax_ins2.set_xticks([0.55,0.6])
