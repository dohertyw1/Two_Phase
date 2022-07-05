import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

with open('bubble_benchmarks/FreeLIFE_case2_shape.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_shapex = list(zip(*reader))[1]

with open('bubble_benchmarks/FreeLIFE_case2_shape.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_shapey = list(zip(*reader))[2]

fFreeLIFE_shapex=[]
for i in range(np.size(FreeLIFE_shapex)):
    fFreeLIFE_shapex.append(float(FreeLIFE_shapex[i]))

fFreeLIFE_shapey=[]
for i in range(np.size(FreeLIFE_shapey)):
    fFreeLIFE_shapey.append(float(FreeLIFE_shapey[i]))

with open('bubble_benchmarks/MooNMD_case2_shape.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_shapex = list(zip(*reader))[1]

with open('bubble_benchmarks/MooNMD_case2_shape.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_shapey = list(zip(*reader))[2]

fMooNMD_shapex=[]
for i in range(np.size(MooNMD_shapex)):
    fMooNMD_shapex.append(float(MooNMD_shapex[i]))

fMooNMD_shapey=[]
for i in range(np.size(MooNMD_shapey)):
    fMooNMD_shapey.append(float(MooNMD_shapey[i]))


mesh = RectangleMesh(Point(0.5,0),Point(1,2),100,400,'crossed')
V = FunctionSpace(mesh, 'CG', 2) 
phi100 = Function(V)
with XDMFFile('Cons_Newtonian_CG_Dim_rein1/frames/frame_3.0/phi_read_3.0.xdmf') as infile:
    infile.read_checkpoint(phi100, "phi")

phi200 = Function(V)
with XDMFFile('Cons_Newtonian_CG_Dim_rein3/frames/frame_3.0/phi_read_3.0.xdmf') as infile:
    infile.read_checkpoint(phi200, "phi")

phi300 = Function(V)
with XDMFFile('Cons_Newtonian_CG_Dim_rein9/frames/frame_3.0/phi_read_3.0.xdmf') as infile:
    infile.read_checkpoint(phi300, "phi")

phi400 = Function(V)
with XDMFFile('Cons_Newtonian_CG_Dim_rein27/frames/frame_3.0/phi_read_3.0.xdmf') as infile:
    infile.read_checkpoint(phi400, "phi")

import matplotlib as mpl

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2

fig = plt.figure(num=None, figsize=(13,11), dpi=100, facecolor='w', edgecolor='k')

plt.subplot(151)
plot(phi100, cmap='jet')
plt.xlim([0.775,0.85])
plt.ylim([0.575,1])
plt.title('Rein 1')

plt.subplots_adjust(wspace=0)

plt.subplot(152)
plot(phi200, cmap='jet')
plt.xlim([0.775,0.85])
plt.ylim([0.575,1])
plt.title('Rein 3')

plt.subplots_adjust(wspace=0)

plt.subplot(153)
plot(phi300, cmap='jet')
plt.xlim([0.775,0.85])
plt.ylim([0.575,1])
plt.title('Rein 9')

plt.subplots_adjust(wspace=0)

plt.subplot(154)
plot(phi400, cmap='jet')
plt.xlim([0.775,0.85])
plt.ylim([0.575,1])
plt.title('Rein 27')

plt.subplot(155)
plot(phi100, mode='contour', levels=[0.5],colors='#1f77b4')
plot(phi200, mode='contour', levels=[0.5],colors='#ff7f0e')
plot(phi300, mode='contour', levels=[0.5],colors='#2ca02c')
plot(phi400, mode='contour', levels=[0.5],colors='#d62728')
plt.scatter([1],[3],color='#1f77b4')
plt.scatter([1],[4],color='#ff7f0e')
plt.scatter([1],[5],color='#2ca02c')
plt.scatter([1],[6],color='#d62728')
plt.scatter(fMooNMD_shapex,fMooNMD_shapey,s=3, color='red')
plt.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey,s=3, color='purple')

plt.xlim([0.775,0.85])
plt.ylim([0.575,1])
plt.title('Comparison')
plt.legend([1,3,9,27,'MooNMD','FreeLIFE'],loc='upper right')



plt.subplots_adjust(wspace=0)














# mesh = RectangleMesh(Point(0, 0), Point(2, 4), 40, 80, 'crossed')
# dist = Expression('sqrt( (pow((x[0]-A),2)) + (pow((x[1]-B),2)) )-r',
#                 degree=2, A=1,B=1,r=0.5)
# eps = Constant(1.5*mesh.hmin())
# dtau = Constant(0.2*mesh.hmin())
# dist2 = Expression('(1/(1+exp((dist/eps))))',degree=2, eps=eps, dist=dist)
# Q = FiniteElement('CG', mesh.ufl_cell(), 2)
# Qs = FunctionSpace(mesh, Q)
# phi0 = interpolate(dist2, Qs)
# from matplotlib import rc
# #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)
# plt.rcParams['font.size'] = 22
# plt.rcParams['axes.linewidth'] = 2
# fig = plt.figure(num=None, figsize=(10,10), dpi=100,facecolor='w', edgecolor='k')
# plot(phi0, mode='contour',levels=[0.5],colors='red',linewidth=10)
# plot(mesh,color='blue',linewidth=0.25)
# plt.xlim([0.35,1.65])
# plt.ylim([0.35,1.65])
# plt.xlabel('x')
# plt.ylabel('y')
# plt.tick_params(axis = 'both', which='major', size=5, width=2, direction='in')
# plt.title('M1 mesh with overset level-set function')