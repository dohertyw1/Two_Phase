import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors

with open('Cons_Newtonian/Cons_Newtonian_data.csv', newline='') as csvfile:
    datanr1601 = np.array(list(csv.reader(csvfile, delimiter='\t')))

# with open('eggs_oot_ncls.csv', newline='') as csvfile:
#     datanr1602 = np.array(list(csv.reader(csvfile, delimiter='\t')))

# datanr1601 = datanr1601[2825:]
datanr1602=datanr1601

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

with open('bubble_benchmarks/FreeLIFE_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_timescale = list(zip(*reader))[1]

fFreeLIFE_timescale=[]
for i in range(np.size(FreeLIFE_timescale)):
    fFreeLIFE_timescale.append(float(FreeLIFE_timescale[i]))

with open('bubble_benchmarks/FreeLIFE_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_area = list(zip(*reader))[2]

fFreeLIFE_area=[]
for i in range(np.size(FreeLIFE_area)):
    fFreeLIFE_area.append(float(FreeLIFE_area[i]))

with open('bubble_benchmarks/FreeLIFE_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_circ = list(zip(*reader))[3]

fFreeLIFE_circ=[]
for i in range(np.size(FreeLIFE_circ)):
    fFreeLIFE_circ.append(float(FreeLIFE_circ[i]))

with open('bubble_benchmarks/FreeLIFE_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_ycom = list(zip(*reader))[4]

fFreeLIFE_ycom=[]
for i in range(np.size(FreeLIFE_ycom)):
    fFreeLIFE_ycom.append(float(FreeLIFE_ycom[i]))

with open('bubble_benchmarks/FreeLIFE_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_vrise = list(zip(*reader))[5]

fFreeLIFE_vrise=[]
for i in range(np.size(FreeLIFE_vrise)):
    fFreeLIFE_vrise.append(float(FreeLIFE_vrise[i]))


with open('bubble_benchmarks/MooNMD_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_timescale = list(zip(*reader))[1]

fMooNMD_timescale=[]
for i in range(np.size(MooNMD_timescale)):
    fMooNMD_timescale.append(float(MooNMD_timescale[i]))


with open('bubble_benchmarks/MooNMD_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_area = list(zip(*reader))[2]

fMooNMD_area=[]
for i in range(np.size(MooNMD_area)):
    fMooNMD_area.append(float(MooNMD_area[i]))


with open('bubble_benchmarks/MooNMD_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_circ = list(zip(*reader))[3]

fMooNMD_circ=[]
for i in range(np.size(MooNMD_circ)):
    fMooNMD_circ.append(float(MooNMD_circ[i]))

with open('bubble_benchmarks/MooNMD_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_ycom = list(zip(*reader))[4]

fMooNMD_ycom=[]
for i in range(np.size(MooNMD_ycom)):
    fMooNMD_ycom.append(float(MooNMD_ycom[i]))

with open('bubble_benchmarks/MooNMD_case1.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    MooNMD_vrise = list(zip(*reader))[5]

fMooNMD_vrise=[]
for i in range(np.size(MooNMD_vrise)):
    fMooNMD_vrise.append(float(MooNMD_vrise[i]))


with open('bubble_benchmarks/FreeLIFE_case1_shape.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_shapex = list(zip(*reader))[1]

with open('bubble_benchmarks/FreeLIFE_case1_shape.txt') as inf:
    reader = csv.reader(inf, delimiter=" ")
    FreeLIFE_shapey = list(zip(*reader))[2]

fFreeLIFE_shapex=[]
for i in range(np.size(FreeLIFE_shapex)):
    fFreeLIFE_shapex.append(float(FreeLIFE_shapex[i]))

fFreeLIFE_shapey=[]
for i in range(np.size(FreeLIFE_shapey)):
    fFreeLIFE_shapey.append(float(FreeLIFE_shapey[i]))

# with open('bubble_benchmarks/MooNMD_case1_shape.txt') as inf:
#     reader = csv.reader(inf, delimiter=" ")
#     MooNMD_shapex = list(zip(*reader))[1]

# with open('bubble_benchmarks/MooNMD_case1_shape.txt') as inf:
#     reader = csv.reader(inf, delimiter=" ")
#     MooNMD_shapey = list(zip(*reader))[2]

# fMooNMD_shapex=[]
# for i in range(np.size(MooNMD_shapex)):
#     fMooNMD_shapex.append(float(MooNMD_shapex[i]))

# fMooNMD_shapey=[]
# for i in range(np.size(MooNMD_shapey)):
#     fMooNMD_shapey.append(float(MooNMD_shapey[i]))

nx=50
mesh =  RectangleMesh(Point(0,0),Point(1,2),nx,2*nx) 
V = FunctionSpace(mesh, 'CG', 2)
phi2 = Function(V)
with XDMFFile("Cons_Newtonian/phi_read.xdmf") as infile:
    infile.read_checkpoint(phi2, "phi")

# nxi=100
# meshi =  RectangleMesh(Point(0,0),Point(1,2),nxi,2*nxi) 
# Vi = FunctionSpace(meshi, 'CG', 1)
# phi2i = Function(Vi)
# with XDMFFile("oot_ncls/phi_read.xdmf") as infile:
#     infile.read_checkpoint(phi2i, "phi")

# mesh = Mesh()
# with XDMFFile("dune_bubble/phi_read.xdmf") as infile:
#     infile.read(mesh)
# V = FunctionSpace(mesh, 'CG', 2) 
# phi0 = Function(V)
# with XDMFFile("dune_bubble/phi_read.xdmf") as infile:
#     infile.read_checkpoint(phi0, "phi")

plt.figure(num=None, figsize=(10, 10), dpi=100, facecolor='w', edgecolor='k')
plt.subplot(221)
# plt.plot(timescale1601, ycom1601,color='black')
# # plt.plot(timescale1602, ycom1602,color='black')
# plt.plot(fMooNMD_timescale,fMooNMD_ycom,color='red',linestyle=':')
# plt.plot(fFreeLIFE_timescale,fFreeLIFE_ycom,color='blue',linestyle='-.')

plt.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey, s=1, color='blue')
# plt.scatter(fMooNMD_shapex,fMooNMD_shapey, s=1, color='red')
plot(phi2, mode='contour', levels=[0.5],color='blue')
# # plot(phi2i, mode='contour', levels=[0],color='green'    )
plt.legend(['MooNMD Benchmark','FreeLIFE Benchmark'])

# plt.xlim([0.1,0.9])
# plt.ylim([0.9, 1.3])
plt.grid()
plt.title('Bubble shape')
plt.legend(['MooNMD Benchmark','FreeLIFE Benchmark'])

plt.subplot(222)
plt.plot(timescale1602, area1602,color='green')
plt.plot(timescale1601, area1601,color='black')
plt.plot(fMooNMD_timescale,fMooNMD_area,color='red',linestyle=':')
plt.plot(fFreeLIFE_timescale,fFreeLIFE_area,color='blue',linestyle='-.')
plt.xlim([0,3])
plt.grid()
plt.title('Area')
plt.legend(['100 CELLS O1','500 CELLS O2','MooNMD Benchmark','FreeLIFE Benchmark'])

plt.subplot(223) 
plt.plot(timescale1602, circ1602,color='green')
plt.plot(timescale1601, circ1601,color='black')
plt.plot(fMooNMD_timescale,fMooNMD_circ,color='red',linestyle=':')
plt.plot(fFreeLIFE_timescale,fFreeLIFE_circ,color='blue',linestyle='-.')
plt.xlim([0,3])
# plt.ylim([0.8,1.1])
plt.grid()
plt.title('Degree of circularity')
plt.legend(['100 CELLS O1','500 CELLS O2','MooNMD Benchmark','FreeLIFE Benchmark'])

plt.subplot(224)
plt.plot(timescale1602, vrise1602,color='green')
plt.plot(timescale1601, vrise1601,color='black')
plt.plot(fMooNMD_timescale,fMooNMD_vrise,color='red',linestyle=':')
plt.plot(fFreeLIFE_timescale,fFreeLIFE_vrise,color='blue',linestyle='-.')
plt.xlim([0,3])
plt.grid()
plt.title('Y-rise velocity')
plt.legend(['100 CELLS O1','500 CELLS O2','MooNMD Benchmark','FreeLIFE Benchmark'])
# plt.suptitle('TEST CASE 1 CONSERVATIVE LEVEL SET')
plt.tight_layout()
# plt.savefig('CONcase1.png')
plt.show()
