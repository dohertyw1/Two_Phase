from fenics import *
import csv
import matplotlib.pyplot as plt
import numpy as np

# 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

class Process():

    def __init__(self, simulations):

        self.simulations = simulations
        self.timescale = []
        self.area = []
        self.xcom = []
        self.ycom = []
        self.circ = []
        self.urise = []
        self.vrise = []

    def load_data(self, suffix, method, fluid):

        with open(f'{method}_{fluid}/{method}_{fluid}_{suffix}.csv', newline='') as csvfile:

            self.raw_data = np.array(list(csv.reader(csvfile, delimiter='\t')))

    def categorise_data(self):

        variables = [self.timescale, self.area, self.xcom, self.ycom, self.circ, self.urise, self.vrise]

        for variable in variables:

            for i in range(np.shape(self.raw_data[0])[0]):

                variable.append(float(self.raw_data[i][variables.index(variable)]))

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

    def load_bubble_shape(self, method, fluid, mesh_saved):

        if (mesh_saved):

            mesh = Mesh()

            with XDMFFile(f'{method}_{fluid}/phi_read.xdmf') as infile:

                infile.read(mesh)

        else:

            mesh = RectangleMesh(Point(0,0), Point(self.height, self.length), self.nx, self.ny)

        V = FunctionSpace(mesh, 'CG', self.ls_order) 
        self.level_set = Function(V)

        with XDMFFile(f'{method}_{fluid}/phi_read.xdmf') as infile:

            infile.read_checkpoint(self.level_set, "phi")

    def individual_plotter(self, variable):

        plt.figure(num=None, figsize=(10, 10), dpi=100, facecolor='w', edgecolor='k')
        plt.plot(self.timescale, variable, color='black')
        plt.plot(self.MooNMD_timescale, self.MooNMD_circ,color='red',linestyle=':')
        plt.plot(self.FreeLIFE_timescale,self.FreeLIFE_circ,color='blue',linestyle='-.')
        plt.xlim([0,3])
        plt.ylim([0.8,1.1])
        plt.grid()
        plt.title('Degree of circularity')
        plt.legend(['Current study', 'MooNMD Benchmark','FreeLIFE Benchmark'])

    def run(self):

        self.load_data('data', 'Cons', 'Newtonian')

        self.categorise_data()

        self.load_benchmark_data('case1')

        # self.load_bubble_shape('Cons', 'Newtonian')

        self.individual_plotter(self.circ)

case = Process(1)
case.run()



#
# plt.subplot(221)
# # plt.plot(timescale1601, ycom1601,color='black')
# # # plt.plot(timescale1602, ycom1602,color='black')
# # plt.plot(fMooNMD_timescale,fMooNMD_ycom,color='red',linestyle=':')
# # plt.plot(fFreeLIFE_timescale,fFreeLIFE_ycom,color='blue',linestyle='-.')

# plt.scatter(fFreeLIFE_shapex,fFreeLIFE_shapey, s=1, color='blue')
# # plt.scatter(fMooNMD_shapex,fMooNMD_shapey, s=1, color='red')
# plot(phi2, mode='contour', levels=[0.5],color='blue')
# # # plot(phi2i, mode='contour', levels=[0],color='green'    )
# # plt.legend(['MooNMD Benchmark','FreeLIFE Benchmark'])

# # plt.xlim([0.1,0.9])
# # plt.ylim([0.9, 1.3])
# plt.grid()
# plt.title('Bubble shape')
# plt.legend(['MooNMD Benchmark','FreeLIFE Benchmark'])

# plt.subplot(222)
# plt.plot(timescale1602, area1602,color='green')
# plt.plot(timescale1601, area1601,color='black')
# plt.plot(fMooNMD_timescale,fMooNMD_area,color='red',linestyle=':')
# plt.plot(fFreeLIFE_timescale,fFreeLIFE_area,color='blue',linestyle='-.')
# # plt.xlim([0,5.25])
# plt.grid()
# plt.title('Area')
# plt.legend(['100 CELLS O1','500 CELLS O2','MooNMD Benchmark','FreeLIFE Benchmark'])

# plt.subplot(223) 
# plt.plot(timescale1602, circ1602,color='green')
# plt.plot(timescale1601, circ1601,color='black')
# plt.plot(fMooNMD_timescale,fMooNMD_circ,color='red',linestyle=':')
# plt.plot(fFreeLIFE_timescale,fFreeLIFE_circ,color='blue',linestyle='-.')
# # plt.xlim([0,3])
# # plt.ylim([0.8,1.1])
# plt.grid()
# plt.title('Degree of circularity')
# plt.legend(['100 CELLS O1','500 CELLS O2','MooNMD Benchmark','FreeLIFE Benchmark'])

# plt.subplot(224)
# plt.plot(timescale1602, vrise1602,color='green')
# plt.plot(timescale1601, vrise1601,color='black')
# plt.plot(fMooNMD_timescale,fMooNMD_vrise,color='red',linestyle=':')
# plt.plot(fFreeLIFE_timescale,fFreeLIFE_vrise,color='blue',linestyle='-.')
# # plt.xlim([0,3])
# plt.grid()
# plt.title('Y-rise velocity')
# plt.legend(['100 CELLS O1','500 CELLS O2','MooNMD Benchmark','FreeLIFE Benchmark'])
# # plt.suptitle('TEST CASE 1 CONSERVATIVE LEVEL SET')
# plt.tight_layout()
# # plt.savefig('CONcase1.png')
# plt.show()
