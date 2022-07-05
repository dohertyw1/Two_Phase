from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

mesh = RectangleMesh(Point(0, 0), Point(2, 4), 200, 400, 'crossed')
V = FunctionSpace(mesh, 'CG', 2)
Vvec = VectorFunctionSpace(mesh, 'CG', 2)
Vten = TensorFunctionSpace(mesh, 'DG', 0)
phi_651 = Function(V)
u_651=Function(Vvec)
tau_651 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.0505/phi_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(phi_651, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.0505/u_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(u_651, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.0505/tau_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(tau_651, "tau")

phi_131 = Function(V)
u_131=Function(Vvec)
tau_131 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.0505/phi_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(phi_131, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.0505/u_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(u_131, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.0505/tau_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(tau_131, "tau")

phi_261 = Function(V)
u_261=Function(Vvec)
tau_261 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.0505/phi_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(phi_261, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.0505/u_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(u_261, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.0505/tau_read_0.0505.xdmf') as infile:
    infile.read_checkpoint(tau_261, "tau")

mesh = RectangleMesh(Point(0, 0), Point(2, 4), 200, 400, 'crossed')
V = FunctionSpace(mesh, 'CG', 2)
Vvec = VectorFunctionSpace(mesh, 'CG', 2)
Vten = TensorFunctionSpace(mesh, 'DG', 0)
phi_65 = Function(V)
u_65=Function(Vvec)
tau_65 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.1005/phi_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(phi_65, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.1005/u_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(u_65, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.1005/tau_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(tau_65, "tau")

phi_13 = Function(V)
u_13=Function(Vvec)
tau_13 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.1005/phi_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(phi_13, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.1005/u_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(u_13, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.1005/tau_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(tau_13, "tau")

phi_26 = Function(V)
u_26=Function(Vvec)
tau_26 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.1005/phi_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(phi_26, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.1005/u_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(u_26, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.1005/tau_read_0.1005.xdmf') as infile:
    infile.read_checkpoint(tau_26, "tau")

mesh = RectangleMesh(Point(0, 0), Point(2, 4), 200, 400, 'crossed')
V = FunctionSpace(mesh, 'CG', 2)
Vvec = VectorFunctionSpace(mesh, 'CG', 2)
Vten = TensorFunctionSpace(mesh, 'DG', 0)
phi_652 = Function(V)
u_652=Function(Vvec)
tau_652 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.2005/phi_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(phi_652, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.2005/u_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(u_652, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=6.5_OB_NS/frames/frame_0.2005/tau_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(tau_652, "tau")

phi_132 = Function(V)
u_132=Function(Vvec)
tau_132 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.2005/phi_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(phi_132, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.2005/u_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(u_132, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=13_OB_NS/frames/frame_0.2005/tau_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(tau_132, "tau")

phi_262 = Function(V)
u_262=Function(Vvec)
tau_262 = Function(Vten)
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.2005/phi_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(phi_262, "phi")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.2005/u_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(u_262, "u")
with XDMFFile('/media/wdoherty1/Samsung_T5/data_for_paper/Cons_Viscoelastic_CG_Dim_c=26_OB_NS/frames/frame_0.2005/tau_read_0.2005.xdmf') as infile:
    infile.read_checkpoint(tau_262, "tau")



# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(num=None, figsize=(15, 15), dpi=100,
                 facecolor='w', edgecolor='k') 

ax1 = fig.add_subplot(1, 3, 1)
plot(phi_651, levels=[0.5], linewidths=3,mode='contour', colors = '#1f77b4')
plot(phi_131, levels=[0.5], linewidths=3,mode='contour', colors = '#ff7f0e')
plot(phi_261, levels=[0.5], linewidths=3,mode='contour', colors = '#2ca02c')
plt.scatter([1],[4],color='#1f77b4')
plt.scatter([1],[5],color='#ff7f0e')
plt.scatter([1],[6],color='#2ca02c')
ax1.set_xlim([0.5,1.5])
ax1.set_ylim([0.5,2.5])
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title(r'$T=0.05$')
ax1.legend([r'$c = 6.5\:wt.\%$',r'$c = 13\:wt.\%$',r'$c = 26\:wt.\%$'],loc='lower left',frameon=True)

ax1.grid()
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1 = fig.add_subplot(1, 3, 2)
plot(phi_65, levels=[0.5], linewidths=3,mode='contour', colors = '#1f77b4')
plot(phi_13, levels=[0.5], linewidths=3,mode='contour', colors = '#ff7f0e')
plot(phi_26, levels=[0.5], linewidths=3,mode='contour', colors = '#2ca02c')
plt.scatter([1],[4],color='#1f77b4')
plt.scatter([1],[5],color='#ff7f0e')
plt.scatter([1],[6],color='#2ca02c')
ax1.set_xlim([0.5,1.5])
ax1.set_ylim([0.5,2.5])
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title(r'$T=0.1$')
ax1.legend([r'$c = 6.5\:wt.\%$',r'$c = 13\:wt.\%$',r'$c = 26\:wt.\%$'],loc='lower left',frameon=True)
ax1.grid()
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')

ax1 = fig.add_subplot(1, 3, 3)
plot(phi_652, levels=[0.5], linewidths=3,mode='contour', colors = '#1f77b4')
plot(phi_132, levels=[0.5], linewidths=3,mode='contour', colors = '#ff7f0e')
plot(phi_262, levels=[0.5], linewidths=3,mode='contour', colors = '#2ca02c')
plt.scatter([1],[4],color='#1f77b4')
plt.scatter([1],[5],color='#ff7f0e')
plt.scatter([1],[6],color='#2ca02c')
ax1.set_xlim([0.5,1.5])
ax1.set_ylim([0.5,2.5])
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title(r'$T=0.2$')
ax1.legend([r'$c = 6.5\:wt.\%$',r'$c = 13\:wt.\%$',r'$c = 26\:wt.\%$'],loc='lower left',frameon=True)
ax1.grid()
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
plt.subplots_adjust(wspace=0)


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

fig.tight_layout()

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
