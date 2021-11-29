import matplotlib.pyplot as plt
from fenics import *
import numpy as np
from tqdm import tqdm
import math
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.rank
set_log_active(False)

def ngamma(phi):
    return grad(phi)/sqrt(dot(grad(phi),grad(phi)))

if comm.size > 1 :
    set_log_active(False)

parameters["std_out_all_processes"] = False;
#parameters['krylov_solver']['nonzero_initial_guess'] = True
# parameters['ghost_mode'] = 'shared_facet' 
parameters["mesh_partitioner"] = "ParMETIS"

T=1
num_steps = 1000
dt=T/num_steps
nx = 100
ny=nx
rein_div = 1
num_steps_rein = 5

original_mesh = UnitSquareMesh(nx,nx) 
mesh = UnitSquareMesh(nx,nx) 
V = FunctionSpace(mesh, 'CG', 2)

center = Point(0.5, 0.75)
radius = 0.15

eps = 0.01 #Constant(0.5*mesh.hmin()**0.9) #### !!!!!!!!!!!!!!!!!!! ####
dtau = Constant(0.5*mesh.hmin()**1.1)

dist = Expression('sqrt((x[0]-A)*(x[0]-A) + (x[1]-B)*(x[1]-B))-r',
            degree=2, A=center[0],    B=center[1],r=radius)

dist2 = Expression('(1/(1+exp((dist/eps))))',
                    degree=2, eps=eps, dist=dist)

phi0=interpolate(dist2, V)

t=0

u0 = Expression(('2*sin(2*p*x[1])*sin(p*x[0])*sin(p*x[0])*cos(p*t)',
                 '-2*sin(2*p*x[0])*sin(p*x[1])*sin(p*x[1])*cos(p*t)'), degree=2,p=np.pi,t=t)

Vel = VectorFunctionSpace(mesh, 'CG', 2)

uic=interpolate(u0,Vel)
u=interpolate(u0,Vel)

file_phi = XDMFFile('vortex_adapt/phi.xdmf')
file_phi.parameters['flush_output'] = True

phi = TrialFunction (V )
w = TestFunction (V)
phi_rein = Function(V)
phi00 = Function(V)

phi_h = Function (V)

Vnormal = VectorFunctionSpace(mesh, 'CG', 1)

phigrad = Function(Vnormal)
vnorm = TestFunction(Vnormal)

F = (phi/dt)*w*dx - (phi0/dt)*w*dx + dot(u, grad(phi0))*w*dx

n = FacetNormal(mesh)
flux = (dot(u,n) + abs(dot(u,n)))/2

def dg_form():

    F = (1/dt)*(phi-phi0)*w*dx \
    - dot ( grad ( w ), u0 * phi ) * dx \
    + dot ( jump ( w ), flux('+') * phi('+') - flux('-') * phi('-') ) * dS \
    + dot ( w, flux * phi ) * ds

    return F

F_grad = inner((phigrad-ngamma(phi0)),vnorm)*dx

F_rein = (phi_rein - phi0)/dtau*w*dx \
        - phi_rein*(1.0 - phi_rein)*inner(grad(w), phigrad)*dx \
        + eps*inner(grad(phi_rein), grad(w))*dx

dt = (0.1*mesh.hmin())/assemble(sqrt(dot(u,u))*dx)

t=0
q=0

while t < T:

    t+=dt

    solve(lhs(F) == rhs(F), phi_h)

    phi0.assign(phi_h)

    solve(F_grad == 0, phigrad)

    for n in range(num_steps_rein):

        solve(F_rein == 0, phi_rein ,solver_parameters={"newton_solver": {"linear_solver": 'gmres', "preconditioner": 'default',\
                                    "maximum_iterations": 20, "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                                    form_compiler_parameters={"optimize": True})

        phi0.assign(phi_rein)

    phi_h.assign(phi_rein)

    file_phi.write(phi_rein, t)
    
    u.t = t
    q+=1
