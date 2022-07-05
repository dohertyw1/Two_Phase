import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors

with open('50x100/data.csv', newline='') as csvfile:  # 70
    datanr1601 = np.array(list(csv.reader(csvfile, delimiter='\t')))
with open('100x200/data.csv', newline='') as csvfile:  # 70
    datanr1602 = np.array(list(csv.reader(csvfile, delimiter='\t')))


with open('200x400/data.csv', newline='') as csvfile:  # 70
    datanr1603 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('400x800eff/data.csv', newline='') as csvfile:  # 70
    datanr1604 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('data.csv', newline='') as csvfile:  # 70
    datanr1605 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('DATAR/data.csv', newline='') as csvfile:  # 70
    datanr1606 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('100x200/data.csv', newline='') as csvfile:  # 70
    datanr1607 = np.array(list(csv.reader(csvfile, delimiter='\t')))


timescale1601 = []
area1601 = []
xcom1601 = []
ycom1601 = []
circ1601 = []
urise1601 = []
vrise1601 = []

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
circ1602 = []
urise1602 = []
vrise1602 = []

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
circ1603 = []
urise1603 = []
vrise1603 = []

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
circ1604 = []
urise1604 = []
vrise1604 = []

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

timescale1605 = []
area1605 = []
xcom1605 = []
ycom1605 = []
circ1605 = []
urise1605 = []
vrise1605 = []

for i in range(np.shape(datanr1605)[0]):
    timescale1605.append(float(datanr1605[i][0]))

for i in range(np.shape(datanr1605)[0]):
    area1605.append(float(datanr1605[i][1]))

for i in range(np.shape(datanr1605)[0]):
    xcom1605.append(float(datanr1605[i][2]))

for i in range(np.shape(datanr1605)[0]):
    ycom1605.append(float(datanr1605[i][3]))

for i in range(np.shape(datanr1605)[0]):
    circ1605.append(float(datanr1605[i][4]))

for i in range(np.shape(datanr1605)[0]):
    urise1605.append(float(datanr1605[i][5]))

for i in range(np.shape(datanr1605)[0]):
    vrise1605.append(float(datanr1605[i][6]))

timescale1606 = []
area1606 = []
xcom1606 = []
ycom1606 = []
circ1606 = []
urise1606 = []
vrise1606 = []

for i in range(np.shape(datanr1606)[0]):
    timescale1606.append(float(datanr1606[i][0]))

for i in range(np.shape(datanr1606)[0]):
    area1606.append(float(datanr1606[i][1]))

for i in range(np.shape(datanr1606)[0]):
    xcom1606.append(float(datanr1606[i][2]))

for i in range(np.shape(datanr1606)[0]):
    ycom1606.append(float(datanr1606[i][3]))

for i in range(np.shape(datanr1606)[0]):
    circ1606.append(float(datanr1606[i][4]))

for i in range(np.shape(datanr1606)[0]):
    urise1606.append(float(datanr1606[i][5]))

for i in range(np.shape(datanr1606)[0]):
    vrise1606.append(float(datanr1606[i][6]))

timescale1607 = []
area1607 = []
xcom1607 = []
ycom1607 = []
circ1607 = []
urise1607 = []
vrise1607 = []

for i in range(np.shape(datanr1607)[0]):
    timescale1607.append(float(datanr1607[i][0]))

for i in range(np.shape(datanr1607)[0]):
    area1607.append(float(datanr1607[i][1]))

for i in range(np.shape(datanr1607)[0]):
    xcom1607.append(float(datanr1607[i][2]))

for i in range(np.shape(datanr1607)[0]):
    ycom1607.append(float(datanr1607[i][3]))

for i in range(np.shape(datanr1607)[0]):
    circ1607.append(float(datanr1607[i][4]))

for i in range(np.shape(datanr1607)[0]):
    urise1607.append(float(datanr1607[i][5]))

for i in range(np.shape(datanr1607)[0]):
    vrise1607.append(float(datanr1607[i][6]))

plt.figure(num=None, figsize=(5, 5), dpi=250, facecolor='w', edgecolor='k')
# # plt.subplot(224)
plt.plot(timescale1601, vrise1601)
# plt.plot(timescale1602, vrise1602)
# plt.plot(timescale1603, vrise1603)
# plt.plot(timescale1604, vrise1604)
plt.plot(timescale1605, vrise1605)
plt.plot(timescale1606, vrise1606)
# plt.plot(timescale1607, vrise1607)
# plt.plot(timescale1602[0:840], vrise1602[0:840])

plt.grid()
plt.title('Y-rise velocity')
# plt.legend(['v=0.15', 'v=0.3', 'v=0.6', 'r=0.2',
#            'r=0.25', 'r=0.3', 'r=0.35', 'r=0.4'])
# plt.suptitle('TEST CASE 1 CONSERVATIVE LEVEL SET')
plt.tight_layout()
plt.savefig('incp.png')
# plt.show()
