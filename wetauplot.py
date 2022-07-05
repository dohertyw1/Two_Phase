import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

mesh = RectangleMesh(Point(0, 0), Point(2, 4), 160, 320, 'crossed')
V = FunctionSpace(mesh, 'CG', 2)
Vten = TensorFunctionSpace(mesh, 'CG', 2)
phi_4 = Function(V)
tau_4 = Function(Vten)
phi_8 = Function(V)
tau_8 = Function(Vten)
phi_16 = Function(V)
tau_16 = Function(Vten)
phi_32 = Function(V)
tau_32 = Function(Vten)

# with XDMFFile('Cons_Viscoelastic_CG_NonDim_w4/frames/frame_12.1244/phi_read_12.1244.xdmf') as infile:
#     infile.read_checkpoint(phi_4, "phi")
# with XDMFFile('Cons_Viscoelastic_CG_NonDim_w8/frames/frame_12.1244/phi_read_12.1244.xdmf') as infile:
#     infile.read_checkpoint(phi_8, "phi")
# with XDMFFile('Cons_Viscoelastic_CG_NonDim_w16/frames/frame_12.1244/phi_read_12.1244.xdmf') as infile:
#     infile.read_checkpoint(phi_16, "phi")
# with XDMFFile('Cons_Viscoelastic_CG_NonDim_w32/frames/frame_12.1244/phi_read_12.1244.xdmf') as infile:
#     infile.read_checkpoint(phi_32, "phi")

with XDMFFile('Cons_Viscoelastic_CG_NonDim_w4_supg/frames/frame_12.1244/tau_read_12.1244.xdmf') as infile:
    infile.read_checkpoint(tau_4, "tau")
# with XDMFFile('Cons_Viscoelastic_CG_NonDim_w8/frames/frame_12.1244/tau_read_12.1244.xdmf') as infile:
#     infile.read_checkpoint(tau_8, "tau")
# with XDMFFile('Cons_Viscoelastic_CG_NonDim_w16/frames/frame_12.1244/tau_read_12.1244.xdmf') as infile:
#     infile.read_checkpoint(tau_16, "tau")
# with XDMFFile('Cons_Viscoelastic_CG_NonDim_w32/frames/frame_12.1244/tau_read_12.1244.xdmf') as infile:
#     infile.read_checkpoint(tau_32, "tau")

fig = plt.figure(num=None, figsize=(24, 12), dpi=100,
                 facecolor='w', edgecolor='k')

ax1 = fig.add_subplot(1, 1, 1)
plot(phi_4, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
plot(tau_4[1,1], cmap='jet')
# plt.colorbar(plot(tau_4[1,1], cmap='jet'),fraction=0.09, pad=0.04)
# ax1.set_xlim([0.5,1.5])
# ax1.set_ylim([0.5,3.5])
# ax1 = fig.add_subplot(1, 4, 2)
# plot(phi_8, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
# plot(tau_8[1,1], mode='contour', cmap='jet')
# # plt.colorbar(plot(tau_8[1,1], cmap='jet'),fraction=0.09, pad=0.04)
# ax1.set_xlim([0.5,1.5])
# ax1.set_ylim([0.5,3.5])
# ax1 = fig.add_subplot(1, 4, 3)
# plot(phi_16, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
# plot(tau_16[1,1], mode='contour', cmap='jet')
# # plt.colorbar(plot(tau_16[1,1], cmap='jet'),fraction=0.09, pad=0.04)
# ax1.set_xlim([0.5,1.5])
# ax1.set_ylim([0.5,3.5])
# ax1 = fig.add_subplot(1, 4, 4)
# plot(phi_32, levels=[0.5], linewidths=3,mode='contour', colors = 'white')
# plot(tau_32[1,1], mode='contour', cmap='jet')
# # plt.colorbar(plot(tau_32[1,1], cmap='jet'),fraction=0.09, pad=0.04)

# ax1.set_xlim([0.5,1.5])
# ax1.set_ylim([0.5,3.5])