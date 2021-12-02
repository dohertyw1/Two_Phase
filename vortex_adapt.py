import matplotlib.pyplot as plt
from fenics import *
import numpy as np
from tqdm import tqdm
import math
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.rank
set_log_active(False)
import time
start = time.time()
def ngamma(phi):
    return grad(phi)/sqrt(dot(grad(phi),grad(phi)))

if comm.size > 1 :
    set_log_active(False)

parameters["std_out_all_processes"] = False;
#parameters['krylov_solver']['nonzero_initial_guess'] = True
parameters["mesh_partitioner"] = "ParMETIS"

element = 'DG'

if element == 'DG':
    parameters['ghost_mode'] = 'shared_facet' 

file_phi = XDMFFile('vortex_adapt/test1.xdmf')
file_phi.parameters['flush_output'] = True

T=0.5
nx = 50

rein_div = 1
num_steps_rein = 1

original_mesh = UnitSquareMesh(nx,nx) 
mesh = UnitSquareMesh(nx,nx) 
V = FunctionSpace(mesh, element, 2)

center = Point(0.5, 0.75)
radius = 0.15

eps = Constant(0.01) #Constant(0.5*mesh.hmin()**0.9) #### !!!!!!!!!!!!!!!!!!! ####
dtau = Constant(0.5*mesh.hmin()**1.1)

dist = Expression('sqrt((x[0]-A)*(x[0]-A) + (x[1]-B)*(x[1]-B))-r',
            degree=2, A=center[0],    B=center[1],r=radius)

dist2 = Expression('(1/(1+exp((dist/eps))))',
                    degree=2, eps=eps, dist=dist)

phi0=interpolate(dist2, V)

t=0

u0 = Expression(('2*sin(2*p*x[1])*sin(p*x[0])*sin(p*x[0])*cos(p*t)',
                 '-2*sin(2*p*x[0])*sin(p*x[1])*sin(p*x[1])*cos(p*t)'), degree=2,p=np.pi,t=t)

Vel = VectorFunctionSpace(mesh, 'DG', 2)

uic=interpolate(u0,Vel)
u=interpolate(u0,Vel)


phi = TrialFunction (V )
w = TestFunction (V)
phi_rein = Function(V)
phi00 = Function(V)
phi_h = Function (V)

dt = (0.1*mesh.hmin())/assemble(sqrt(dot(u,u))*dx)

if element == 'CG':

    Vnormal = VectorFunctionSpace(mesh, 'CG' , 1)

    phigrad = Function(Vnormal)
    vnorm = TestFunction(Vnormal)

    F_grad = inner((phigrad-ngamma(phi0)),vnorm)*dx

    F = (phi/dt)*w*dx - (phi0/dt)*w*dx + dot(u, grad(phi0))*w*dx

    F_reincg = (phi_rein - phi0)/dtau*w*dx \
            - phi_rein*(1.0 - phi_rein)*inner(grad(w), phigrad)*dx \
            + eps*inner(grad(phi_rein), grad(w))*dx  
        # - 0.5*dot(grad(w),(phi_rein+phi0)*phigrad)*dx \
        # + dot(phi_rein*phi0*phigrad,grad(w))*dx \
        # + (eps/2)*(dot(grad(phi_rein)+grad(phi0),phigrad)*dot(grad(w),phigrad))*dx

if element == 'DG':

    phigrad = ngamma(phi0)

    n = FacetNormal(mesh)
    un = abs(dot(u('+'), n('+')))
    nn0 = abs(dot(phigrad('+'), n('+')))

    F = (1/dt)*(phi-phi0)*w*dx \
    - dot ( grad ( w ), u0 * phi ) * dx \
    + (dot(u('+'), jump(w, n))*avg(phi) \
    + 0.5*un*dot(jump(phi, n), jump(w, n)))*dS \

    F_reindg = (phi_rein - phi0)/dtau*w*dx \

    alpha = 50    

    F_reindg_compr = - dot ( grad ( w ), phigrad * phi_rein*(1-phi_rein) ) * dx \
                    + (dot(phigrad('+'), jump(w, n))*avg(phi_rein) \
                    + 0.5*nn0*dot(jump(phi_rein, n), jump(w, n)))*dS \
                    - (dot(phigrad('+'), jump(w, n))*avg(phi_rein*phi_rein) \
                    + 0.5*nn0*dot(jump(phi_rein*phi_rein, n), jump(w, n)))*dS \

    F_reindg_diff = eps*(inner(grad(phi_rein), grad(w))*dx \
                + (alpha/mesh.hmin())*dot(jump(w, n), jump(phi_rein, n))*dS \
                - dot(avg(grad(w)), jump(phi_rein, n))*dS \
                - dot(jump(w, n), avg(grad(phi_rein)))*dS)

    F_reindg += F_reindg_compr
    F_reindg += F_reindg_diff

t=0
q=0

while t < T:
# for n in tqdm(range(num_steps)):

    t+=dt

    solve(lhs(F) == rhs(F), phi_h)

    phi0.assign(phi_h)

    if element == 'CG':

        solve(F_grad == 0, phigrad)
    
    if element == 'DG':

        phigrad = ngamma(phi0)

    for n in range(num_steps_rein):

        if element == 'CG':

            solve(F_reincg == 0, phi_rein ,solver_parameters={"newton_solver": {"linear_solver": 'gmres', "preconditioner": 'default',\
                            "maximum_iterations": 20, "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                            form_compiler_parameters={"optimize": True})

        elif element == 'DG':

            solve(F_reindg == 0, phi_rein ,solver_parameters={"newton_solver": {"linear_solver": 'gmres', "preconditioner": 'default',\
                                        "maximum_iterations": 20, "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                                        form_compiler_parameters={"optimize": True})

        phi0.assign(phi_rein)

    phi_h.assign(phi_rein)

    file_phi.write(phi_h, t)
    
    u.t = t
    q+=1


print('It took', time.time()-start, 'seconds.')

