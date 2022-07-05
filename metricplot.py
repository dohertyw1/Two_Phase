from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

with open('Cons_Viscoelastic_CG_NonDim_dt4/Cons_Viscoelastic_CG_NonDim_dt4_data.csv', newline='') as csvfile: #70
    datanr1601 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('Cons_Viscoelastic_CG_NonDim_dt8/Cons_Viscoelastic_CG_NonDim_dt8_data.csv', newline='') as csvfile: #70
    datanr1602 = np.array(list(csv.reader(csvfile, delimiter='\t')))


with open('Cons_Viscoelastic_CG_NonDim_d16/Cons_Viscoelastic_CG_NonDim_d16_data.csv', newline='') as csvfile: #70
    datanr1603 = np.array(list(csv.reader(csvfile, delimiter='\t')))

with open('Cons_Viscoelastic_CG_NonDim_d32/Cons_Viscoelastic_CG_NonDim_d32_data.csv', newline='') as csvfile: #70
    datanr1604 = np.array(list(csv.reader(csvfile, delimiter='\t')))
drag1 = []

for i in range(np.shape(datanr1601)[0]):
    drag1.append(float(datanr1601[i][7]))
drag2 = []

for i in range(np.shape(datanr1602)[0]):
    drag2.append(float(datanr1602[i][7]))

drag3 = []

for i in range(np.shape(datanr1603)[0]):
    drag3.append(float(datanr1603[i][7]))

drag4 = []

for i in range(np.shape(datanr1604)[0]):
    drag4.append(float(datanr1604[i][7]))

plt.plot(drag1)
plt.plot(drag2)
plt.plot(drag3)
plt.plot(drag4)