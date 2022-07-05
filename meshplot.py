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

mesh00 = RectangleMesh(Point(0,0),Point(0.5,2),20,80,'crossed')

mesh0 = RectangleMesh(Point(0,0),Point(0.5,2),160,640,'crossed')
mesh1 = RectangleMesh(Point(0.5,0),Point(1,2),160,640,'crossed')
mesh2 = RectangleMesh(Point(0,0),Point(1,1),320,320,'crossed')

dist = Expression('x[1]-0.5-0.1*sin(8*x[0])',degree=2)

# dist2 = Expression('dist > eps ? p1 : dist < -eps ? p2 : (dist >= -eps && dist <= eps ? p3 : 0)',
#                 p1=Expression('0', degree=0, domain=mesh),
#                 p2=Expression('1', degree=0, domain=mesh),
#                 p3=Expression('(dist+eps)/(2*eps)+(sin(pi*dist/eps))/(2*pi)',
#                 degree = 2, domain=mesh, dist = dist, pi=np.pi, eps=eps),
#                 degree = 2, domain=mesh, dist = dist, pi=np.pi, eps=eps)


eps = 0.015 #Constant(1.5*mesh1.hmin())
dist2 = Expression('(1/(1+exp((dist/eps))))',degree=2, eps=eps, dist=dist)

Q = FiniteElement('CG', mesh0.ufl_cell(), 2)
Qs = FunctionSpace(mesh0, Q)
phi0 = interpolate(dist2, Qs)

Q = FiniteElement('CG', mesh1.ufl_cell(), 2)
Qs = FunctionSpace(mesh2, Q)
phi1 = interpolate(dist2, Qs)

# mesh = RectangleMesh(Point(0.5,0),Point(1,2),160,640,'crossed')
# V = VectorFunctionSpace(mesh, 'CG', 2) 
# u400 = Function(V)
# with XDMFFile('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_160x640test_case2/u_read.xdmf') as infile:
#     infile.read_checkpoint(u400, "u")



rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 2
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(num=None, figsize=(6,6), dpi=250, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(1,1,1)
# plot(mesh00, color='blue', linewidth=0.25)

plot(phi1, cmap='jet')
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1.tick_params(axis='both', which='major', labelsize=13)

# plot(phi400, mode='contour',levels=[0.5], colors='white',linewidths=0.5)
# plot(phi400, cmap='jet')
plt.colorbar(plot(phi1, cmap='jet', vmin = 0, vmax = 1),fraction=0.088, pad=0.04)
# a = plot(phi400, mode='contour',levels=[0.5], colors='red',linewidths=0.5).allsegs[0][0]
# q=0
# for i in a[:,0]:
#     a[q,0]=-i
#     q+=1
# b=a[:,0]
# c = a[:,1]
# ax1.plot(b+1,c, color='red')

ax1.text(0.6, 0.5, r'$\textbf{n}_{\Gamma}$', fontsize = 22,color='white')
ax1.annotate("", xy=(0.57, 0.58), xytext=(0.57, 0.41),
            arrowprops=dict(width=0.75,color='white'))
ax1.text(0.075, 0.175, r'$\textbf{g}$', fontsize = 22,color='white')
ax1.annotate("", xy=(0.05, 0.05), xytext=(0.05, 0.25),
            arrowprops=dict(width=0.75,color='white'))
# ax1.annotate("", xy=(0.5, 0.5), xytext=(0.75, 0.5),
#             arrowprops=dict(arrowstyle='<->',color='white'))
# ax1.text(0.61, 0.52, r'$r$', fontsize = 22,color='white')
# ax1.text(0.05, 0.7, r'$\Gamma_N$', fontsize = 22,color='white')
# ax1.text(0.88, 0.7, r'$\Gamma_N$', fontsize = 22,color='white')
# ax1.text(0.5, 0.05, r'$\Gamma_D$', fontsize = 22,color='white')
# ax1.text(0.5, 0.88, r'$\Gamma_D$', fontsize = 22,color='white')
plot(phi1, mode='contour',levels=[0.5], colors='black', linewidths=2)

# ax1.text(0.52, 0.78, r'$\Gamma_I$', fontsize = 22,color='white')
# ax1.text(0.58, 0.76, r'$\phi_{0}$', fontsize = 22,color='white')
ax1.text(0.65, 0.8, r'$\Omega_{1}$', fontsize = 22,color='white')
ax1.text(0.35, 0.25, r'$\Omega_{2}$', fontsize = 22,color='white')

ax1.set_title('Two-Phase Flow Schematic')
plt.xticks([])
plt.yticks([])

# ax1.set_xlabel('x')
# ax1.set_ylabel('y')
ax1.set_xlim([0,1])
ax1.set_ylim([0,1])