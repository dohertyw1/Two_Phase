import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


with open('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_20x80test_case2/Cons_Newtonian_CG_Dim_20x80test_case2_data.csv', newline='') as csvfile: #70
    datanr1601 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_40x160test_case2/Cons_Newtonian_CG_Dim_40x160test_case2_data.csv', newline='') as csvfile: #70
    datanr1602 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_80x320test_case2/Cons_Newtonian_CG_Dim_80x320test_case2_data.csv', newline='') as csvfile: #70
    datanr1603 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_160x640test_case2/Cons_Newtonian_CG_Dim_160x640_case2_data.csv', newline='') as csvfile: #70  
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


with open('bubble_benchmarks/FreeLIFE_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_timescale = list(zip(*reader))[1]

fFreeLIFE_timescale=[]
for i in range(np.size(FreeLIFE_timescale)):
    fFreeLIFE_timescale.append(float(FreeLIFE_timescale[i]))

with open('bubble_benchmarks/FreeLIFE_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_area = list(zip(*reader))[2]

fFreeLIFE_area=[]
for i in range(np.size(FreeLIFE_area)):
    fFreeLIFE_area.append(float(FreeLIFE_area[i]))

with open('bubble_benchmarks/FreeLIFE_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_circ = list(zip(*reader))[3]

fFreeLIFE_circ=[]
for i in range(np.size(FreeLIFE_circ)):
    fFreeLIFE_circ.append(float(FreeLIFE_circ[i]))

with open('bubble_benchmarks/FreeLIFE_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_ycom = list(zip(*reader))[4]

fFreeLIFE_ycom=[]
for i in range(np.size(FreeLIFE_ycom)):
    fFreeLIFE_ycom.append(float(FreeLIFE_ycom[i]))

with open('bubble_benchmarks/FreeLIFE_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_vrise = list(zip(*reader))[5]

fFreeLIFE_vrise=[]
for i in range(np.size(FreeLIFE_vrise)):
    fFreeLIFE_vrise.append(float(FreeLIFE_vrise[i]))


with open('bubble_benchmarks/MooNMD_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_timescale = list(zip(*reader))[1]

fMooNMD_timescale=[]
for i in range(np.size(MooNMD_timescale)):
    fMooNMD_timescale.append(float(MooNMD_timescale[i]))


with open('bubble_benchmarks/MooNMD_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_area = list(zip(*reader))[2]

fMooNMD_area=[]
for i in range(np.size(MooNMD_area)):
    fMooNMD_area.append(float(MooNMD_area[i]))


with open('bubble_benchmarks/MooNMD_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_circ = list(zip(*reader))[3]

fMooNMD_circ=[]
for i in range(np.size(MooNMD_circ)):
    fMooNMD_circ.append(float(MooNMD_circ[i]))

with open('bubble_benchmarks/MooNMD_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_ycom = list(zip(*reader))[4]

fMooNMD_ycom=[]
for i in range(np.size(MooNMD_ycom)):
    fMooNMD_ycom.append(float(MooNMD_ycom[i]))

with open('bubble_benchmarks/MooNMD_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_vrise = list(zip(*reader))[5]

fMooNMD_vrise=[]
for i in range(np.size(MooNMD_vrise)):
    fMooNMD_vrise.append(float(MooNMD_vrise[i]))


with open('bubble_benchmarks/TP2D_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    TP2D_timescale = list(zip(*reader))[1]

fTP2D_timescale=[]
for i in range(np.size(TP2D_timescale)):
    fTP2D_timescale.append(float(TP2D_timescale[i]))


with open('bubble_benchmarks/TP2D_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    TP2D_area = list(zip(*reader))[2]

fTP2D_area=[]
for i in range(np.size(TP2D_area)):
    fTP2D_area.append(float(TP2D_area[i]))


with open('bubble_benchmarks/TP2D_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    TP2D_circ = list(zip(*reader))[3]

fTP2D_circ=[]
for i in range(np.size(TP2D_circ)):
    fTP2D_circ.append(float(TP2D_circ[i]))

with open('bubble_benchmarks/TP2D_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    TP2D_ycom = list(zip(*reader))[4]

fTP2D_ycom=[]
for i in range(np.size(TP2D_ycom)):
    fTP2D_ycom.append(float(TP2D_ycom[i]))

with open('bubble_benchmarks/TP2D_case2.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    TP2D_vrise = list(zip(*reader))[5]

fTP2D_vrise=[]
for i in range(np.size(TP2D_vrise)):
    fTP2D_vrise.append(float(TP2D_vrise[i]))

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

mesh = RectangleMesh(Point(0.5,0),Point(1,2),20,80,'crossed')
V = FunctionSpace(mesh, 'CG', 2) 
phi100 = Function(V)
with XDMFFile('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_20x80test_case2/phi_read.xdmf') as infile:
    infile.read_checkpoint(phi100, "phi")

mesh = RectangleMesh(Point(0.5,0),Point(1,2),40,160,'crossed')
V = FunctionSpace(mesh, 'CG', 2) 
phi200 = Function(V)
with XDMFFile('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_40x160test_case2/phi_read.xdmf') as infile:
    infile.read_checkpoint(phi200, "phi")

mesh = RectangleMesh(Point(0.5,0),Point(1,2),80,320,'crossed')
V = FunctionSpace(mesh, 'CG', 2) 
phi300 = Function(V)
with XDMFFile('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_80x320test_case2/phi_read.xdmf') as infile:
    infile.read_checkpoint(phi300, "phi")

mesh = RectangleMesh(Point(0.5,0),Point(1,2),160,640,'crossed')
V = FunctionSpace(mesh, 'CG', 2) 
phi400 = Function(V)
with XDMFFile('/media/wdoherty1/Samsung_T5/case1case2/Cons_Newtonian_CG_Dim_160x640test_case2/phi_read.xdmf') as infile:
    infile.read_checkpoint(phi400, "phi")

import matplotlib as mpl

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2

fig = plt.figure(num=None, figsize=(13,9), dpi=100, facecolor='w', edgecolor='k')

ax1= fig.add_subplot(3,2,1)

# ax1.plot(timescale1601, ycom1601)
ax1.plot(timescale1602, ycom1602)
ax1.plot(timescale1603, ycom1603)
ax1.plot(timescale1604, ycom1604)
ax1.plot(fMooNMD_timescale,fMooNMD_ycom,color='red',linestyle=':')
ax1.plot(fFreeLIFE_timescale,fFreeLIFE_ycom,color='blue',linestyle='-.')

ax1.set_xlim([0,3])
ax1.set_ylim([0.4,1.2])
ax1.set_xlabel('t',**{'fontname':'Times', 'size':'18'})
ax1.set_ylabel('Centre of Mass',**{'fontname':'Times', 'size':'18'})
ax1.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax1.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax1.tick_params(axis='both', which='major', labelsize=13)

ax1.grid()
# ax1.set_title('Centre of Mass')
ax1.legend(['M1','M2','M3','MooNMD','FreeLIFE'],loc='upper left',frameon=True)

# ax3= fig.add_subplot(2,3,2)

# area1 = [abs(x-area1601[0]) for x in area1601]
# area2 = [abs(x-area1602[0]) for x in area1602]
# area3 = [abs(x-area1603[0]) for x in area1603]
# area4 = [abs(x-area1604[0]) for x in area1604]

# ax3.plot(timescale1601, area1)
# ax3.plot(timescale1602, area2)
# ax3.plot(timescale1603, area3)
# ax3.plot(timescale1604, area4)

# moo_area = [x-fMooNMD_area[0] for x in fMooNMD_area]
# free_area = [x-fFreeLIFE_area[0] for x in fFreeLIFE_area]

# ax3.plot(fMooNMD_timescale,moo_area,color='red',linestyle=':')
# ax3.plot(fFreeLIFE_timescale,free_area,color='blue',linestyle='-.')
# ax3.set_xlim([0,3])
# ax3.set_xlabel('t')
# ax3.set_ylabel('A')
# ax3.grid()
# ax3.set_title('Area')
# ax3.legend(['M1','M2','M3','M4','MooNMD','FreeLIFE'],loc='center left')

ax5= fig.add_subplot(3,2,3)

# ax5.plot(timescale1601, vrise1601)
ax5.plot(timescale1602, vrise1602)
ax5.plot(timescale1603, vrise1603)
ax5.plot(timescale1604, vrise1604)

ax5.plot(fMooNMD_timescale,fMooNMD_vrise,color='red',linestyle=':')
ax5.plot(fFreeLIFE_timescale,fFreeLIFE_vrise,color='blue',linestyle='-.')
ax5.set_xlim([0,3])
ax5.set_ylim([0,0.3])
ax5.set_xlabel('t',**{'fontname':'Times', 'size':'18'})
ax5.set_ylabel('Rise Velocity',**{'fontname':'Times', 'size':'18'})
ax5.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax5.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax5.tick_params(axis='both', which='major', labelsize=13)

ax5.grid()
# ax5.set_title('Rise Velocity')
ax5.legend(['M1','M2','M3','MooNMD','FreeLIFE'],loc='lower center',frameon=True)

axc= fig.add_subplot(3,2,5)

# axc.plot(timescale1601, circ1601)
axc.plot(timescale1602, circ1602)
axc.plot(timescale1603, circ1603)
axc.plot(timescale1604, circ1604)

axc.plot(fMooNMD_timescale,fMooNMD_circ,color='red',linestyle=':')
axc.plot(fFreeLIFE_timescale,fFreeLIFE_circ,color='blue',linestyle='-.')
axc.set_xlim([0,3])
axc.set_xlabel('t',**{'fontname':'Times', 'size':'18'})
axc.set_ylabel('Circularity',**{'fontname':'Times', 'size':'18'})
axc.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
axc.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
axc.tick_params(axis='both', which='major', labelsize=13)

axc.grid()
# axc.set_title('Circularity')
axc.legend(['M1','M2','M3','MooNMD','FreeLIFE'],loc='upper right',frameon=True)

ax7= fig.add_subplot(1,2,2)

# plot(phi100, mode='contour',levels=[0.5],colors='#1f77b4')
plot(phi200, mode='contour', levels=[0.5],colors='#1f77b4')
plot(phi300, mode='contour', levels=[0.5],colors='#ff7f0e')
plot(phi400, mode='contour', levels=[0.5],colors='#2ca02c')
# plot(phi400, cmap='jet')

ax7.plot([1.5,2.2],[1.5,2.4],color='#1f77b4')
ax7.plot([1.5,2.2],[1.5,2.4],color='#ff7f0e')
ax7.plot([1.5,2.2],[1.5,2.4],color='#2ca02c')
# ax7.plot([1.5,2.2],[1.5,2.4],color='#d62728')

ax7.scatter(fMooNMD_shapex,fMooNMD_shapey,s=1.5, color='red')
ax7.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey,s=1.5, color='purple')

ax7.legend(['M1','M2','M3','MooNMD','FreeLIFE'],loc='lower right')

ax7.grid()

ax7.set_xlim([0.5,1.025])
ax7.set_ylim([0.5, 1.5])
ax7.set_xlabel('x',**{'fontname':'Times', 'size':'18'})
ax7.set_ylabel('y',**{'fontname':'Times', 'size':'18'})
ax7.tick_params(axis='both', which='major', labelsize=13)


ax7.set_title('Bubble Shape',**{'fontname':'Times', 'size':'18'})

fig.tight_layout()

ax_ins = inset_axes(ax7, width="100%", height="100%", loc='upper left',
                   bbox_to_anchor=(0.44,0.68,.45,.45), bbox_transform=ax7.transAxes)

# plot(phi100, mode='contour',levels=[0.5],colors='#1f77b4')
plot(phi200, mode='contour', levels=[0.5],colors='#1f77b4')
plot(phi300, mode='contour', levels=[0.5],colors='#ff7f0e')
plot(phi400, mode='contour', levels=[0.5],colors='#2ca02c')
ax_ins.scatter(fMooNMD_shapex,fMooNMD_shapey,s=3, color='red')
ax_ins.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey,s=3, color='purple')

ax_ins.grid()
ax_ins.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax_ins.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')

ax_ins.set_xlim([0.55,0.6])
ax_ins.set_ylim([0.99, 1.025])
ax_ins.set_xticks([0.55,0.575,0.6])
ax_ins.tick_params(axis='both', which='major', labelsize=14)

ax_ins2 = inset_axes(ax7, width="100%", height="100%", loc='upper left',
                   bbox_to_anchor=(0.09,0.09,.35,.35), bbox_transform=ax7.transAxes)

# plot(phi100, mode='contour',levels=[0.5],colors='#1f77b4')
plot(phi200, mode='contour', levels=[0.5],colors='#1f77b4')
plot(phi300, mode='contour', levels=[0.5],colors='#ff7f0e')
plot(phi400, mode='contour', levels=[0.5],colors='#2ca02c')
ax_ins2.scatter(fMooNMD_shapex,fMooNMD_shapey,s=3, color='red')
ax_ins2.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey,s=3, color='purple')
ax_ins2.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
ax_ins2.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
ax_ins2.grid()

ax_ins2.set_xlim([0.76,0.85])
ax_ins2.set_ylim([0.57, 0.77])
plt.subplots_adjust(wspace=0)

# ax_ins2.set_xticks([0.55,0.6])

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