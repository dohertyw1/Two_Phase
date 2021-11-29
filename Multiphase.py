from fenics import *
import numpy as np
from tqdm import tqdm
import csv
import os
from mpi4py import MPI
from Flow import *
from Funcs import *
from Process import *

class Multiphase (Flow, Process):

    """"Initialise material parameters, mesh and test case data"""
    def __init__(self, test_case, method, fluid, constitutive_equation, dimensional, stability, ls_scheme, lsr_method,
                 rein_steps, size, location, dimension, ls_order):

        super().__init__()

        self.test_case = test_case
        self.method = method
        self.fluid = fluid
        self.constitutive_equation = constitutive_equation
        self.dimensional = dimensional
        self.stability = stability
        self.ls_scheme = ls_scheme
        self.lsr_method = lsr_method
        self.size = size
        self.location = location
        self.dimension = dimension
        self.ls_order = ls_order

        if (dimension == '2D'):

            self.mesh = RectangleMesh(Point(0,0),Point(1,2),self.size,2*self.size) 
            self.g = Constant((0,0.98))

        elif (dimension == '3D'):

            self.mesh = BoxMesh(Point(0,0,0), Point(1,2,1), self.size, 2*self.size, self.size)
            self.g = Constant((0,0.98,0))
        
        self.amin = RectangleMesh(Point(0,0),Point(1,2),300,600).hmin()
        self.hmin = self.mesh.hmin()

        self.pli_in = None
        self.pli_out = 2

        self.T = 3
        self.num_steps = int(self.T/self.dt)
        
        self.rein_div = 1
        self.rein_steps = rein_steps

        self.comm = MPI.comm_world
        self.rank = MPI.rank(self.comm)
        self.Id = Identity(int(self.dimension[0]))

        self.facet_normal = FacetNormal(self.mesh)

        if (self.method == 'NCons'):

            self.alph = Constant(0.0625*self.hmin)
            self.eps1 = Constant(self.hmin)

        if (self.test_case == 1):

            self.height = 2
            self.length = 1
            self.mu1 = [1, "Inside"]
            self.rho1 = [100, "Inside"]
            self.mu2 = [10, "Outside"]
            self.rho2 = [1000, "Outside"]
            self.sigma = 24.5
            self.d = 0.10
            self.eps = 0.01 #Constant(0.50*self.hmin**(1-self.d)) 
            self.dtau = Constant(0.50*self.hmin**(1+self.d)) 
            self.centre = [0.5, 0.5]
            self.radius = 0.25
            self.dt = (1/self.size)/2

        elif (self.test_case == 2):
            
            self.height = 2
            self.length = 1
            self.mu1 = [0.1, "Inside"]
            self.rho1 = [1, "Inside"]
            self.mu2 = [10, "Outside"]
            self.rho2 = [1000, "Outside"]
            self.sigma = 1.96
            self.d = 0.10
            self.eps = Constant(0.50*self.hmin**(1-self.d)) 
            self.dtau = Constant(0.50*self.hmin**(1+self.d)) 
            self.centre = [0.5, 0.5]
            self.radius = 0.25
            self.dt = (1/self.size)/2

        if (self.fluid == 'Viscoelastic'):

            self.eta_s_in = 1.025
            self.eta_p_in = 0
            self.eta_s_out = 0.7175
            self.eta_p_out = 9.5325
            self.lamb_in = 0
            self.lamb_out = 0.2
            self.rho1 = [0.1, "Inside"]
            self.rho2 = [1, "Outside"]
            self.sigma = 10
            self.alpha = 0.1

            self.beta1 = 1
            self.beta2 = 0.07 
            self.Re2 = 1.419436759
            self.Re_eps = 10
            self.Wi1 = 0
            self.Wi2 = 8.082903769
            self.Eo = 35.28

            self.dtau = Constant(0.5*self.hmin**1.10)

            if (self.dimensional == 'Dim'):
                
                self.height = 4
                self.length = 2
                self.mesh = RectangleMesh(Point(0,0),Point(self.length,self.height),self.size,2*self.size) 
                self.g = Constant((0,980))
                self.centre = [1, 1]
                self.radius = 0.3
                self.eps = 0.02
                self.curvature = Constant(self.sigma)
                
                self.T = 0.13
                self.dt = 0.001
                self.num_steps = int(self.T/self.dt)

            if (self.dimensional == 'NonDim'):

                self.height = 6.666
                self.length = 3.333
                self.mesh = RectangleMesh(Point(0,0),Point(self.length,self.height),self.size,2*self.size) 
                self.g = Constant((0,1))
                self.centre = [1/0.6, 1/0.6]
                self.radius = 0.3/0.6
                self.eps = 0.0333333
                self.curvature = Constant(1/self.Eo)

                self.T = 0.3*(np.sqrt(980*0.6))/0.6 
                self.dt = (0.001*(np.sqrt(980*0.6))/0.6)
                self.num_steps = int(self.T/self.dt)

    """"Define signed distance function (ncls and cls) and heaviside function (cls)"""
    def level_set(self):

        if (self.dimension == '2D'):

            self.sdf = Expression('sqrt((pow((x[0]-A),2))+(pow((x[1]-B),2)))-r', degree=2, A=self.centre[0], B=self.centre[1], r=self.radius)

        elif (self.dimension == '3D'):

            self.sdf = Expression('sqrt((pow((x[0]-A),2))+(pow((x[1]-B),2))+(pow((x[2]-C),2)))-r', degree=2, A=0.5, B=0.5, C=0.5, r=0.25)

        if (self.method == 'Cons'):

            self.hdf = Expression('(1/(1+exp((dist/eps))))',degree=2, eps=self.eps, dist=self.sdf)

    """Initialise FEM spaces and functions"""
    def fem_data(self):
        
        """ Level Set spaces/functions """
        self.Q = FunctionSpace(self.mesh, 'CG', self.ls_order)
        self.Vnorm = VectorFunctionSpace(self.mesh, 'CG', 1)

        self.phi = TrialFunction(self.Q)
        self.psi = TestFunction(self.Q)
        self.phi_rein = Function(self.Q)

        if (self.method == 'NCons'):

            self.phi0 = interpolate(self.sdf,self.Q)
            self.phi00 = interpolate(self.sdf,self.Q)
            self.phiic = interpolate(self.sdf,self.Q)
            self.sign = sgn(self.phi0, self.eps1)
            self.phin = Function(self.Vnorm)
            self.psin = TestFunction(self.Vnorm)
        
        elif (self.method == 'Cons'):

            self.phi0 = interpolate(self.hdf,self.Q)
            self.phi00 = interpolate(self.hdf,self.Q)
            self.phiic = interpolate(self.hdf,self.Q)
            self.phin = Function(self.Vnorm)
            self.psin = TestFunction(self.Vnorm)

        """ Pressure spaces/functions """
        self.P = FunctionSpace(self.mesh, 'CG', 1)

        self.p = TrialFunction(self.P)
        self.q = TestFunction(self.P)
        self.p0 = Function(self.P)
        self.p_  = Function(self.P)

        """ Velocity spaces/functions """
        if (self.stability == None):

            self.V = VectorFunctionSpace(self.mesh, 'CG', 2)

            self.u = TrialFunction(self.V) 
            self.v = TestFunction(self.V)
            self.u0 = Function(self.V) 
            self.u_ = Function(self.V)

        """ Stress spaces/functions """
        if (self.stability == None):

            self.T = TensorFunctionSpace(self.mesh, 'CG', 2)

            self.tau = TrialFunction(self.T)
            self.zeta = TestFunction(self.T)
            self.tau0 = Function(self.T)

        """ DEVSSG-DG spaces/functions"""

        if (self.stability == 'DEVSSG-DG'):

            self.T = TensorFunctionSpace(self.mesh, 'DG', 0)
            self.V_element = VectorElement('CG', self.mesh.ufl_cell(), 2) 
            self.D_element = TensorElement('DG', self.mesh.ufl_cell(), 2)
            self.VD = self.V_element*self.D_element
            self.W = FunctionSpace(self.mesh, self.VD)

            self.u, self.G = TrialFunctions(self.W)
            self.v, self.R = TestFunctions(self.W)

            self.w0 = Function(self.W)
            (self.u0, self.G0) = split(self.w0)

            self.w_ = Function(self.W)
            (self.u_, self.G_) = split(self.w_)

    """Construct boundary conditions"""
    def bcs(self):
        
        self.bcs = []

        walls   = f'near(x[1], 0) || near(x[1], {self.height})'
        fswalls = f'near(x[0], 0) || near(x[0], {self.length})'

        fswalls_x = 'near(x[0], 0) || near(x[0], 1)'
        fswalls_z = 'near(x[2], 0) || near(x[2], 1)'

        if (self.dimension == '2D'):

            bcu_noslip  = DirichletBC(self.V, Constant((0, 0)), walls)
            bcu_fslip  = DirichletBC(self.V.sub(0), Constant(0), fswalls)

            self.bc_ns = [bcu_noslip, bcu_fslip]

        elif (self.dimension == '3D'):

            bcu_noslip  = DirichletBC(self.V, Constant((0, 0, 0)), walls)
            bcu_fslip_x = DirichletBC(self.V.sub(0), Constant(0), fswalls_x)
            bcu_fslip_z = DirichletBC(self.V.sub(2), Constant(0), fswalls_z)

            self.bc_ns = [bcu_noslip, bcu_fslip_x, bcu_fslip_z]

    """Calculate density/viscosity functions depending on level set and type of fluid."""
    def mat_pams(self):

        if (self.method == 'NCons'):

            self.rho = phasify_noncon(self.phi, self.rho1[0], self.rho2[0], self.eps)
            self.rho0 = phasify_noncon(self.phi00, self.rho1[0], self.rho2[0], self.eps)

            if (self.fluid == 'Newtonian'):

                self.mu = phasify_noncon(self.phi, self.mu1[0], self.mu2[0], self.eps)

            elif (self.fluid == 'Viscoelastic'):

                pass

        elif (self.method == 'Cons'):

            self.rho = phasify_con(self.phi, self.rho1[0], self.rho2[0])
            self.rho0 = phasify_con(self.phi00, self.rho1[0], self.rho2[0])

            if (self.fluid == 'Newtonian'):

                self.mu = phasify_con(self.phi, self.mu1[0], self.mu2[0])
            
            elif (self.fluid == 'Viscoelastic'):

                if (self.constitutive_equation == 'OB' or self.constitutive_equation == 'Giesekus') and (self.dimensional == 'Dim'):

                    self.eta_s = phasify_con(self.phi, self.eta_s_in, self.eta_s_out)

                    self.eta_p = phasify_con(self.phi, self.eta_p_in, self.eta_p_out)

                    self.lamb = phasify_con(self.phi, self.lamb_in, self.lamb_out)

                if (self.constitutive_equation == 'OB' or self.constitutive_equation == 'Giesekus') and (self.dimensional == 'NonDim'):

                    self.beta = phasify_con(self.phi, self.beta1, self.beta2)

                    self.Re = phasify_con(self.phi, self.Re2*self.Re_eps, self.Re2)

                    self.Wi = phasify_con(self.phi, 0, self.Wi2)

    """Write the solution to file."""
    def write_to_file(self):

        self.xdmf_file_phi.write(self.phi0, self.t)
        self.xdmf_file_u.write(self.u0, self.t)
        self.xdmf_file_p.write(self.p0, self.t)

        if (self.fluid == 'Viscoelastic'):

            self.xdmf_file_tau.write(self.tau0, self.t)

    """Calculate benchmark data."""
    def process_data(self):

        if (self.method == 'NCons' and self.dimension == '2D'):

            area = assemble(conditional(lt(self.phi0, 0), 1.0, 0.0)*dx)
            x_com = assemble(Expression("x[0]", degree = 1)*(conditional(lt(self.phi0, 0), 1.0, 0.0))*dx)/area
            y_com = assemble(Expression("x[1]", degree = 1)*(conditional(lt(self.phi0, 0), 1.0, 0.0))*dx)/area

            Pa = 2.0*sqrt(np.pi*area)
            Pb = assemble(mgrad(self.phi0)*delta_func(self.phi0/sqrt(dot(grad(self.phi0),grad(self.phi0))), self.eps)*dx)
            circ = Pa/Pb

            u_rise = assemble(self.u0[0]*(conditional(lt(self.phi0, 0), 1.0, 0.0))*dx)/area
            v_rise = assemble(self.u0[1]*(conditional(lt(self.phi0, 0), 1.0, 0.0))*dx)/area

        if (self.method == 'Cons' and self.dimension == '2D'):

            area = assemble(self.phi0*dx)
            x_com = assemble(Expression("x[0]", degree = 1)*self.phi0*dx)/area
            y_com = assemble(Expression("x[1]", degree = 1)*self.phi0*dx)/area

            Pa = 2.0*sqrt(np.pi*area)
            Pb = assemble(mgrad(self.phi0)*dx)
            circ = Pa/Pb

            u_rise = assemble(self.u0[0]*self.phi0*dx)/area
            v_rise = assemble(self.u0[1]*self.phi0*dx)/area

        self.timeseries = [self.t, area, x_com, y_com, circ, u_rise, v_rise]

        if self.rank == 0:
            with open((f"{self.method}_{self.fluid}/{self.method}_{self.fluid}_data.csv"), 'a') as csvfile:
                f = csv.writer(csvfile, delimiter='\t',lineterminator='\n',)
                f.writerow(self.timeseries)

    """Write final level-set function to separate file in order to retrieve the bubble shape"""
    def process_shape(self):

        if (self.method == 'NCons'):

            with XDMFFile(f"{self.method}_{self.fluid}/phi_read.xdmf") as outfile:

                outfile.write_checkpoint(self.phi0, "phi", 0, append=True)

        elif (self.method == 'Cons'):

            with XDMFFile(f"{self.method}_{self.fluid}/phi_read.xdmf") as outfile:

                outfile.write_checkpoint(self.phi0, "phi", 0, append=True)

    """Set up files to save data and delete old data."""
    def set_up_files(self):

        self.xdmf_file_phi = XDMFFile(f'{self.method}_{self.fluid}/phi.xdmf')
        self.xdmf_file_u = XDMFFile(f'{self.method}_{self.fluid}/u.xdmf')
        self.xdmf_file_p = XDMFFile(f'{self.method}_{self.fluid}/p.xdmf')

        if (self.fluid == 'Viscoelastic'):

            self.xdmf_file_tau = XDMFFile(f'{self.method}_{self.fluid}/tau.xdmf')
            self.xdmf_file_tau.parameters['flush_output'] = True 

        if (self.rank == 0 and os.path.isfile(f'{self.method}_{self.fluid}/{self.method}_{self.fluid}_data.csv') == True):

            os.remove(f'{self.method}_{self.fluid}/{self.method}_{self.fluid}_data.csv')
          
        self.xdmf_file_phi.parameters['flush_output'] = True
        self.xdmf_file_u.parameters['flush_output'] = True
        self.xdmf_file_p.parameters['flush_output'] = True  

    """Run the solver"""
    def run(self):

        set_log_active(False)

        self.set_up_files()

        self.level_set()

        self.fem_data()

        self.bcs()

        self.ls_form(self.phi, self.phi0, self.psi, self.u0)

        if (self.method == 'NCons'):

            self.nclsr_form(self.phi_rein, self.phi0, self.psi, self.sign)

        elif (self.method == 'Cons'):

            self.clsr_form(self.phi_rein, self.phi0, self.psi, self.phin)

        self.phi = Function(self.Q)
        
        self.mat_pams()

        self.ns_form(self.u, self.u0, self.u_, self.p, self.p0, self.p_, self.phi, self.v, self.q, self.tau0)

        if (self.fluid == 'Viscoelastic'):

            self.constitutive_equation_form(self.tau, self.tau0, self.zeta, self.u0)
            self.tau = Function(self.T)

        self.t = 0
        self.n = 0
        self.q = 0

        for self.n in tqdm(range(self.num_steps)):

            self.t += self.dt

            self.ls_solve()

            if (self.n % self.rein_div == 0):

                if (self.method == 'NCons'):

                    self.nclsr_solve(self.phi0, self.phi00, self.phi, self.phi_rein)

                elif (self.method == 'Cons'):

                    self.clsr_solve(self.phi0, self.phi00, self.phi, self.phi_rein, self.phin, self.psin)
            else:
                
                self.phi00.assign(self.phi0)

            self.mat_pams()

            if (self.fluid == 'Viscoelastic'):

                self.constitutive_equation_solve()

            self.ns_solve(self.bc_ns, self.u_, self.p_)

            self.u0.assign(self.u_)
            self.p0.assign(self.p_)
            self.phi0.assign(self.phi)

            if (self.fluid == 'Viscoelastic'):

                self.tau0.assign(self.tau)

            if (self.q % 1 == 0):

                self.write_to_file()

            self.process_data()

            self.q += 1
        
        self.process_shape()

        self.post_process()


case1 = Multiphase(test_case = None,
                   method = 'Cons',
                   fluid = 'Viscoelastic',
                   constitutive_equation = 'OB',
                   dimensional = 'NonDim',
                   stability = None,
                   ls_scheme = 'Euler',
                   lsr_method = 'old_method',
                   rein_steps = 1,
                   size = 50,
                   location = 'HAWK',
                   dimension = '2D',
                   ls_order = 2)

case2 = Multiphase(test_case = 1,
                   method = 'Cons',
                   fluid = 'Newtonian',
                   constitutive_equation = None,
                   dimensional = 'Dim',
                   stability = None,
                   ls_scheme = 'Euler',
                   lsr_method = 'old_method',
                   rein_steps = 3,
                   size = 50,
                   location = 'HAWK',
                   dimension = '2D',
                   ls_order = 2)

case1.run()

# follows
# follows
# false
# for b d


"""good results were eps: 0.9 dtau: 1.10 size: 200 euler scheme old method order 2"""

            # elif (self.fluid == 'GNF_PL'):

            #     self.mu = phasify_con(self.phi, self.mu1[0], powerlaw(self.u0, self.pli_out))

            # elif (self.fluid == 'GNF_C'):

            #     self.mu = phasify_con(self.phi, 
            #                     carreau(self.u0, self.pli_in, self.eta0_in, self.etainf_in, self.lamb_in), 
            #                     carreau(self.u0, self.pli_out, self.eta0_out, self.etainf_out, self.lamb_out))

# if u dot n > 0 => u dot n = abs(u dot n)

# <=> (u dot n) / abs(u dot n) = 1 since abs(u dot n) > 0

# <=> (u / abs(u dot n)) dot n = 1

# => u / abs(u dot n) = n

# => u = abs(u dot n) * n 