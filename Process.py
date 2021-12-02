from fenics import *
import csv
import matplotlib.pyplot as plt
import numpy as np

# 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

class Process():

    def __init__(self):

        pass
        # self.method = 'Cons'
        # self.fluid = 'Newtonian'

    """Load in data from the Multiphase.py solver."""
    def load_data(self, suffix):

        with open(f'{self.method}_{self.fluid}_{self.element}/{self.method}_{self.fluid}_{self.element}_{suffix}.csv', newline='') as csvfile:

            self.raw_data = np.array(list(csv.reader(csvfile, delimiter='\t')))

    """Sort raw data into respective lists of data."""
    def categorise_data(self):

        self.timescale = []
        self.area = []
        self.xcom = []
        self.ycom = []
        self.circ = []
        self.urise = []
        self.vrise = []

        variables = [self.timescale, self.area, self.xcom, self.ycom, self.circ, self.urise, self.vrise]

        for variable in variables:

            for i in range(np.shape(self.raw_data)[0]):

                variable.append(float(self.raw_data[i][variables.index(variable)]))

    """Load in the Newtonian benchmark data from bubble_benchmarks file."""
    def load_benchmark_data(self, case):

        raw_FreeLIFE_timescale = []
        raw_FreeLIFE_area  = []
        raw_FreeLIFE_circ  = []
        raw_FreeLIFE_ycom  = []
        raw_FreeLIFE_vrise  = []

        raw_FreeLIFE = [raw_FreeLIFE_timescale, raw_FreeLIFE_area, raw_FreeLIFE_circ, raw_FreeLIFE_ycom, raw_FreeLIFE_vrise]

        self.FreeLIFE_timescale = []
        self.FreeLIFE_area  = []
        self.FreeLIFE_circ  = []
        self.FreeLIFE_ycom  = []
        self.FreeLIFE_vrise  = []

        self.FreeLIFE = [self.FreeLIFE_timescale, self.FreeLIFE_area, self.FreeLIFE_circ, self.FreeLIFE_ycom, self.FreeLIFE_vrise]

        for variable in raw_FreeLIFE:

            with open(f'bubble_benchmarks/FreeLIFE_{case}.txt') as inf:

                reader = csv.reader(inf, delimiter=" ")

                variable.append(list(zip(*reader))[raw_FreeLIFE.index(variable)+1])

        for variable in self.FreeLIFE:

            for i in range(np.size(raw_FreeLIFE[self.FreeLIFE.index(variable)])):

                variable.append(float(raw_FreeLIFE[self.FreeLIFE.index(variable)][0][i]))

        raw_MooNMD_timescale = []
        raw_MooNMD_area  = []
        raw_MooNMD_circ  = []
        raw_MooNMD_ycom  = []
        raw_MooNMD_vrise  = []

        raw_MooNMD = [raw_MooNMD_timescale, raw_MooNMD_area, raw_MooNMD_circ, raw_MooNMD_ycom, raw_MooNMD_vrise]

        self.MooNMD_timescale = []
        self.MooNMD_area  = []
        self.MooNMD_circ  = []
        self.MooNMD_ycom  = []
        self.MooNMD_vrise  = []

        self.MooNMD = [self.MooNMD_timescale, self.MooNMD_area, self.MooNMD_circ, self.MooNMD_ycom, self.MooNMD_vrise]

        for variable in raw_MooNMD:

            with open(f'bubble_benchmarks/MooNMD_{case}.txt') as inf:

                reader = csv.reader(inf, delimiter=" ")

                variable.append(list(zip(*reader))[raw_MooNMD.index(variable)+1])

        for variable in self.MooNMD:

            for i in range(np.size(raw_MooNMD[self.MooNMD.index(variable)])):

                variable.append(float(raw_MooNMD[self.MooNMD.index(variable)][0][i]))

    """Load in the saved bubble shape from the last temporal iteration."""
    def load_bubble_shape(self, mesh_saved):

        if (mesh_saved):

            mesh = Mesh()

            with XDMFFile(f'{self.method}_{self.fluid}_{self.element}/phi_read.xdmf') as infile:

                infile.read(mesh)

        else:

            mesh = RectangleMesh(Point(0,0), Point(self.height, self.length), self.nx, self.ny)

        V = FunctionSpace(mesh, 'CG', self.ls_order) 
        self.level_set = Function(V)

        with XDMFFile(f'{self.method}_{self.fluid}_{self.element}/phi_read.xdmf') as infile:

            infile.read_checkpoint(self.level_set, "phi")

    """Plot for various quantities."""
    def individual_plotter(self, variable):

        if (variable == 'Area'):

            current_var = self.area
            MooNMD_var = self.MooNMD_area
            FreeLIFE_var = self.FreeLIFE_area

        elif (variable == 'Circularity'):

            current_var = self.circ
            MooNMD_var = self.MooNMD_circ
            FreeLIFE_var = self.FreeLIFE_circ

        elif (variable == 'Centre of mass'):

            current_var = self.ycom
            MooNMD_var = self.MooNMD_ycom
            FreeLIFE_var = self.FreeLIFE_ycom

        elif (variable == 'Rise Velocity'):

            current_var = self.vrise
            MooNMD_var = self.MooNMD_vrise
            FreeLIFE_var = self.FreeLIFE_vrise

        plt.plot(self.timescale, current_var, color='black')
        # plt.plot(self.MooNMD_timescale, MooNMD_var,color='red',linestyle=':')
        # plt.plot(self.FreeLIFE_timescale,FreeLIFE_var,color='blue',linestyle='-.')
        # plt.xlim([0,3])
        plt.grid()
        plt.title(f'{variable}')
        plt.legend(['Current study', 'MooNMD Benchmark','FreeLIFE Benchmark'])
        plt.tight_layout()

    """Run the post processing."""
    def post_process(self):

        self.load_data('data')

        self.categorise_data()

        self.load_benchmark_data('case1')

        # self.load_bubble_shape('Cons', 'Newtonian')

        if (self.rank == 0):

            plt.figure(num=None, figsize=(10, 10), dpi=100, facecolor='w', edgecolor='k')

            plt.subplot(221)
            self.individual_plotter('Area')
            plt.subplot(222)
            self.individual_plotter('Circularity')
            plt.subplot(223)
            self.individual_plotter('Centre of mass')
            plt.subplot(224)
            self.individual_plotter('Rise Velocity')
            plt.suptitle(f'{self.method}_{self.fluid}_{self.element}')
            plt.savefig(f'test{self.element}')
            plt.show()
