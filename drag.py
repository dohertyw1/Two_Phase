from fenics import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import meshio
import csv
from mpi4py import MPI
from matplotlib import rc
import matplotlib as mpl
import csv

import matplotlib.pyplot as plt
import numpy as np
from fenics import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
with open('Cons_Viscoelastic_CG_NonDim_w4/Cons_Viscoelastic_CG_NonDim_w4_data.csv', newline='') as csvfile:
    d1 = np.array(list(csv.reader(csvfile, delimiter='\t')))

drag01 = []
for i in range(np.shape(d1)[0]):
    drag01.append(float(d1[i][7]))

with open('Cons_Viscoelastic_CG_NonDim_w8/Cons_Viscoelastic_CG_NonDim_w8_data.csv', newline='') as csvfile:
    d2 = np.array(list(csv.reader(csvfile, delimiter='\t')))

drag02 = []
for i in range(np.shape(d2)[0]):
    drag02.append(float(d2[i][7]))

with open('Cons_Viscoelastic_CG_NonDim_w16/Cons_Viscoelastic_CG_NonDim_w16_data.csv', newline='') as csvfile:
    d3 = np.array(list(csv.reader(csvfile, delimiter='\t')))

drag03 = []
for i in range(np.shape(d3)[0]):
    drag03.append(float(d3[i][7]))

with open('Cons_Viscoelastic_CG_NonDim_w32/Cons_Viscoelastic_CG_NonDim_w32_data.csv', newline='') as csvfile:
    d4 = np.array(list(csv.reader(csvfile, delimiter='\t')))

drag04 = []
for i in range(np.shape(d4)[0]):
    drag04.append(float(d4[i][7]))
