import matplotlib.pyplot as plt
from fenics import *
import numpy as np
from tqdm import tqdm
import math
import ufl as uf
import csv
import os
#os.remove('eggs.csv') 
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

element = 'DG'
parameters["std_out_all_processes"] = False;
#parameters['krylov_solver']['nonzero_initial_guess'] = True
parameters['ghost_mode'] = 'shared_facet' 
parameters["mesh_partitioner"] = "ParMETIS"

def case1():

    mu1 = 1 # INSIDE
    mu2 = 10 # OUTSIDE
    rho1 = 100 # INSIDE
    rho2 = 1000 # OUTSIDE
    sig = 24.5

    return mu1, mu2, rho1, rho2, sig

mu1, mu2, rho1, rho2, sig = case1()

def magnitude(u):
    return sqrt(u**2)

def epsilon(u): 
    return sym(grad(u))

def trans(u):
    return u.T

def mgrad(phi):

    return sqrt(dot(grad(phi),grad(phi)))

def ngamma(phi):
    
    return grad(phi)/sqrt(dot(grad(phi),grad(phi)))

def rho(phi):

    return rho1*phi + rho2*(1-phi)

def mu(phi):
    
    return mu1*phi + mu2*(1-phi)

def sigma(u, p, phi, Re):

    return ((2*mu(phi))/Re)*epsilon(u) - p*Identity(len(u))

T=3

if rank == 0:
    if (os.path.isfile('eggs_test.csv') == True):
        os.remove('eggs_test.csv') 
    if (os.path.isfile('eggs_bubs.csv') == True):
        os.remove('eggs_bubs.csv') 
set_log_active(False)

num_steps_rein = 1
rein_div=1

nx=50
ny=2*nx
beta = 2
d=0.1
radius = 0.25

grav = 0.98

Fr = 1
g = Constant((0, grav))

dt= 0.01
num_steps =int( T / dt)


def curvature(phi, sig, u0, v):

    return (sig)*mgrad(phi)*inner((Identity(len(u0)) - outer(ngamma(phi), ngamma(phi))), epsilon(v))*dx

mesh = RectangleMesh(Point(0,0),Point(1,2),nx,ny) 

V = FunctionSpace(mesh, element, 1)

P = FunctionSpace(mesh, 'CG', 1)
Vel = VectorFunctionSpace(mesh, 'CG', 2)

eps = Constant(0.01)
dtau = Constant(0.5*mesh.hmin()**1.1)

dist = Expression('sqrt( (pow((x[0]-A),2)) + (pow((x[1]-B),2)) )-r',
                degree=2, A=0.5,B=0.5,r=0.25)

dist2 = Expression('(1/(1+exp((dist/eps))))',
                    degree=2, eps=eps, dist=dist)

walls   = 'near(x[1], 0) || near(x[1], 2)'
fswalls = 'near(x[0], 0) || near(x[0], 1)'

bcu_noslip  = DirichletBC(Vel, Constant((0, 0)), walls)
bcu_fslip  = DirichletBC(Vel.sub(0), Constant(0), fswalls)

bcu = [bcu_noslip, bcu_fslip]

file_phi = XDMFFile('cons_test/phi50.xdmf')
file_u = XDMFFile('cons_test/u50.xdmf')
file_p = XDMFFile('cons_test/p50.xdmf')

file_phi.parameters['flush_output'] = True
file_u.parameters['flush_output'] = True
file_p.parameters['flush_output'] = True

phi = TrialFunction (V )
phi0 = Function (V)
phi00 = Function (V)
phins = Function (V)
w = TestFunction (V)

phiic = interpolate(dist2, V)
phi0 = interpolate(dist2, V)

u = TrialFunction(Vel) 
v = TestFunction(Vel)
p = TrialFunction(P)
q = TestFunction(P)

u0 = Function(Vel) 
u_  = Function(Vel) 
p0 = Function(P)
p_  = Function(P)

phi_rein = Function(V)

if element == 'CG':

    Vnormal = VectorFunctionSpace(mesh, 'CG' , 1)

    phigrad = Function(Vnormal)
    vnorm = TestFunction(Vnormal)

    F_grad = inner((phigrad-ngamma(phi0)),vnorm)*dx

    F = (phi/dt)*w*dx - (phi0/dt)*w*dx + dot(u0, grad(phi0))*w*dx

    F_reincg = (phi_rein - phi0)/dtau*w*dx \
            - phi_rein*(1.0 - phi_rein)*inner(grad(w), phigrad)*dx \
            + eps*inner(grad(phi_rein), grad(w))*dx  
        # - 0.5*dot(grad(w),(phi_rein+phi0)*phigrad)*dx \
        # + dot(phi_rein*phi0*phigrad,grad(w))*dx \
        # + (eps/2)*(dot(grad(phi_rein)+grad(phi0),phigrad)*dot(grad(w),phigrad))*dx

if element == 'DG':

    # phigrad = ngamma(phi0)

    Vnormal = VectorFunctionSpace(mesh, 'DG' , 1)

    phigrad = Function(Vnormal)
    vnorm = TestFunction(Vnormal)

    F_grad = inner((phigrad-ngamma(phi0)),vnorm)*dx

    n = FacetNormal(mesh)
    un = abs(dot(u0('+'), n('+')))
    nn0 = abs(dot(phigrad('+'), n('+')))

    F = (1/dt)*(phi-phi0)*w*dx \
    - dot ( grad ( w ), u0 * phi ) * dx \
    + (dot(u0('+'), jump(w, n))*avg(phi) \
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

phi = Function (V)

ns1 = (rho(phi)/dt)*inner(u - u0, v)*dx \
    + rho(phi)*inner(dot(u0, nabla_grad(u)), v)*dx \
    + (mu(phi))*inner(grad(u),grad(v))*dx \
    - p0*div(v)*dx \
    + rho(phi)*inner(g,v)*dx \
    + curvature(phi, sig, u0, v)

ns2 = (1/rho(phi))*(dot(grad(p),grad(q)) - dot(grad(p0),grad(q)))*dx \
    + div(u_)*q*dx \

ns3 = inner(u,v)*dx - inner(u_,v)*dx \
    + (1/rho(phi))*inner(grad(p_-p0),v)*dx

bcs=[]

tol = 1.0e-4

t=0

# while t < T:
for n in tqdm(range(num_steps)):

    t+=dt

    a = assemble(lhs(F))
    L = assemble(rhs(F))

    solve(a, phi.vector(), L, 'gmres', 'default') 

    phi0.assign(phi)

    if element == 'CG':

        solve(F_grad == 0, phigrad)
    
    if element == 'DG':

        # phigrad = ngamma(phi0)
        solve(F_grad == 0, phigrad)


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

    phi.assign(phi_rein)

    file_phi.write(phi, t)

    A1=assemble(lhs(ns1))
    M1=assemble(rhs(ns1))
    [bc.apply(A1) for bc in bcu]
    [bc.apply(M1) for bc in bcu]
    solve(A1, u_.vector(), M1, 'gmres', 'default') 

    M2=assemble(rhs(ns2))
    A2=assemble(lhs(ns2))
    solve(A2, p_.vector(), M2, 'gmres', 'default') 
        
    A3=assemble(lhs(ns3))
    M3=assemble(rhs(ns3))
    [bc.apply(A3) for bc in bcu]
    [bc.apply(M3) for bc in bcu]
    solve(A3, u_.vector(), M3, 'gmres', 'default') 

    file_u.write(u_, t)
    file_p.write(p_, t)

    u0.assign(u_)
    p0.assign(p_)

    phi0.assign(phi)

