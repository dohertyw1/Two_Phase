import matplotlib.pyplot as plt
from fenics import *
import numpy as np
from tqdm import tqdm
import math
import ufl as uf
import csv
import os
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if (rank == 0 and os.path.isfile(f'droplet_ve401/data.csv') == True):

    os.remove('droplet_ve401/data.csv')

set_log_active(False)

parameters["std_out_all_processes"] = False;
parameters['ghost_mode'] = 'shared_facet' 
parameters["mesh_partitioner"] = "ParMETIS"

file_phi = XDMFFile('droplet_ve401/phidim.xdmf')
file_u = XDMFFile('droplet_ve401/udim.xdmf')
file_p = XDMFFile('droplet_ve401/pdim.xdmf')
file_tau = XDMFFile('droplet_ve401/taudim.xdmf')

file_phi.parameters['flush_output'] = True
file_u.parameters['flush_output'] = True
file_p.parameters['flush_output'] = True
file_tau.parameters['flush_output'] = True

rho_in = 1.204
eta_s_in = 0.00001825
eta_p_in = 0 

rho_out = 1000.9
eta_s_out = 0.031
eta_p_out = 1.511

lamb_in = 0
lamb_out = 0.207
sig = 0.076
grav = 9.81

# alpha = 0.1

nx=100
ny=200

T=0.3
num_steps = 1500
qdiv=15
num_steps_rein = 3
rein_div=3

dt= T/num_steps

g = Constant((0, grav))

scaling = 0.1
toler = 1.0e-3
centrex=0.05
centrey=0.01
radius = 0.00229
x0 = 0.05
y0 = 0
x1 = 0.1
y1 = 0.1

"""widen refined region"""
mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),nx,ny,'left') 


# cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
# cell_markers.set_all(False)
# for cell in cells(mesh):
#     for facet in facets(cell): 
#         for vertex in vertices(facet):
#             if (0.05 <= vertex.point().array()[0] <= 0.075) and (0 <= vertex.point().array()[1] <= 0.05):
#                 cell_markers[cell] = True

# mesh = refine(mesh, cell_markers)

# cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
# cell_markers.set_all(False)
# for cell in cells(mesh):
#     for facet in facets(cell): 
#         for vertex in vertices(facet):
#             if (0.05 <= vertex.point().array()[0] <= 0.075) and (0 <= vertex.point().array()[1] <= 0.05):
#                 cell_markers[cell] = True

# mesh = refine(mesh, cell_markers)

# cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
# cell_markers.set_all(False)
# for cell in cells(mesh):
#     for facet in facets(cell): 
#         for vertex in vertices(facet):
#             if (0.05 <= vertex.point().array()[0] <= 0.075) and (0 <= vertex.point().array()[1] <= 0.05):
#                 cell_markers[cell] = True

# mesh = refine(mesh, cell_markers)

eps = Constant(1.5*mesh.hmin())
dtau = Constant(0.2*mesh.hmin())

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

def CDelta(phi, eps):

    return conditional(lt(abs(phi),eps), 1.0/(2.0*eps)*(1.0 + uf.cos(np.pi*phi/eps)), 0.0)

def CHeaviside(phi, eps):

    return conditional(lt(abs(phi),eps), 0.5*(1.0 + phi/eps + 1/np.pi*uf.sin(np.pi*phi/eps)), (uf.sign(phi) + 1)/2.0)

def rho(phi):

    return rho_in*phi + rho_out*(1-phi)

def lamb(phi):

    return lamb_in*phi + lamb_out*(1-phi)

def eta_s(phi):

    return eta_s_in*phi + eta_s_out*(1-phi)

def eta_p(phi):

    return eta_p_in*phi + eta_p_out*(1-phi)

def c(phi):

    c_in = 0
    c_out = 1

    return c_in*phi + c_out*(1-phi)


def epsilon_axi(u,x):

    return sym(as_tensor([[u[0].dx(0), 0, u[0].dx(1)],
                          [0, u[0] / x[0], 0],
                          [u[1].dx(0), 0, u[1].dx(1)]]))

def div_axi(u, x):

    return (1/x[0])*(x[0]*u[0]).dx(0) + u[1].dx(1)

def nabla_grad_axi(u,x):

    return as_tensor([[u[0].dx(0), 0, u[0].dx(1)],
                          [0, u[0] / x[0], 0],
                          [u[1].dx(0), 0, u[1].dx(1)]])



dist = Expression('sqrt((pow((x[0]-A),2))+(pow((x[1]-B),2)))-r', degree=2, A=centrex, B=centrey, r=radius)
dist2 = Expression('(1/(1+exp((dist/eps))))',degree=2, eps=eps, dist=dist)


""" LEVEL SET """
Qs = FunctionSpace(mesh, 'CG', 2)
phi = TrialFunction(Qs)
psi = TestFunction(Qs)
phi00 = Function(Qs)
phi0 = interpolate(dist2, Qs)

phi_rein = Function(Qs)

""" LEVEL SET REIN """
Vnormal = VectorFunctionSpace(mesh, 'CG', 2)
phigrad = Function(Vnormal)
vnorm = TestFunction(Vnormal)

""" PRESSURE """
Ps = FunctionSpace(mesh, 'CG', 1)
p = TrialFunction(Ps)
q = TestFunction(Ps)
p0 = Function(Ps)
p_  = Function(Ps)

""" VELOCITY """
Vs = VectorFunctionSpace(mesh, 'CG', 2)
u = TrialFunction(Vs) 
v = TestFunction(Vs)
u0 = Function(Vs) 
u_  = Function(Vs) 

""" MIXED """
# V = VectorElement('CG', mesh.ufl_cell(), 2)
# D = TensorElement('DG', mesh.ufl_cell(), 2)
# Dtt = FiniteElement('DG', mesh.ufl_cell(), 2)

# mixed = FunctionSpace(mesh, MixedElement([V, D, Dtt]))

# u, G, G_tt = TrialFunctions(mixed)
# v, R, R_tt = TestFunctions(mixed)

# w1 = Function(mixed)
# u1, G1, G1_tt = split(w1)

# w0 = Function(mixed)
# u0, G0, G0_tt = split(w0)

# w_ = Function(mixed)
# u_, G_, G__tt = split(w_)

""" STRESS """
Ts = TensorFunctionSpace(mesh, 'DG', 0)
tau = TrialFunction(Ts)
S = TestFunction(Ts)
tau0 = Function(Ts)

""" STRESS tt"""
Ttt = FunctionSpace(mesh, 'DG', 0)
tautt = TrialFunction(Ttt)
zetatt = TestFunction(Ttt)
tau0tt = Function(Ttt)

""" DEVSSG """
T_devssg = TensorFunctionSpace(mesh, 'DG', 2)
G = TrialFunction(T_devssg)
R = TestFunction(T_devssg)
G1 = Function(T_devssg)

""" DEVSSG """
T_tt_devssg = FunctionSpace(mesh, 'DG', 2)
G_tt = TrialFunction(T_tt_devssg)
R_tt = TestFunction(T_tt_devssg)
G1_tt = Function(T_tt_devssg)

walls   = f'near(x[1], {y0}) || near(x[1], {y1})'
fswalls = f'near(x[0], {x0}) || near(x[0], {x1})'

bcu_noslip  = DirichletBC(Vs, Constant((0, 0)), walls)
bcu_fslip  = DirichletBC(Vs.sub(0), Constant(0), fswalls)

def flux(ui, ni, taui):
    return (dot(ui,ni)*taui + abs(dot(ui,ni))*taui)/2

n=FacetNormal(mesh)
x = SpatialCoordinate(mesh)
r=x[0]

bcu = [bcu_noslip, bcu_fslip]

F_levelset = (phi/dt)*psi*r*dx - (phi0/dt)*psi*r*dx + inner(u0, grad(phi))*psi*r*dx \
           + ((phi - phi0)/dt + inner(u0, grad(phi))) \
           * scaling*mesh.hmin()/uf.Max(2.0*sqrt(inner(u0, u0)),toler/mesh.hmin())*inner(u0, grad(psi))*r*dx

F_rein = (phi_rein - phi0)/dtau*psi*r*dx \
        - 4*phi_rein*(1.0 - phi_rein)*inner(grad(psi), phigrad)*r*dx \
        + eps*inner(grad(phi_rein), grad(psi))*r*dx

F_grad = inner((phigrad-ngamma(phi0)),vnorm)*r*dx

phi = Function(Qs)

h = 0.5*CellDiameter(mesh)
supg = (h/magnitude(u0)+0.000001)

F_devssg = inner(G - grad(u0),R)*r*dx

F_tt_devssg = inner(G_tt - u0[0]/r,R_tt)*r*dx
# G0 = grad(u0)
# G0_tt = u0[0]/r

F_stress = (1/dt)*inner(lamb(phi)*(tau-tau0),S)*r*dx \
         + 0.5*((inner(dot(lamb(phi)*u0,nabla_grad(tau)) \
         - dot(lamb(phi)*tau, trans(G1)) \
         - dot(lamb(phi)*G1, tau),S)*r*dx) \
         + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tau,S)*r*dx) \
         + 0.5*((inner(dot(lamb(phi)*u0,nabla_grad(tau0)) \
         - dot(lamb(phi)*tau0, trans(G1)) \
         - dot(lamb(phi)*G1, tau0),S)*r*dx) \
         + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tau0,S)*r*dx) \
         - eta_p(phi)*inner(G1+trans(G1),S)*r*dx \
         + 0.5*0.1*inner(flux(u0,n,tau0)('+') - flux(u0,n,tau0)('-'),jump(S))*r*dS \
         + 0.5*0.1*inner(flux(u0,n,tau)('+') - flux(u0,n,tau)('-'),jump(S))*r*dS


        #  + 0.5*lamb(phi)*inner(dot(u0, nabla_grad(S)), supg*dot(u0,nabla_grad(tau)))*r*dx \
        #  + 0.5*lamb(phi)*inner(dot(u0, nabla_grad(S)), supg*dot(u0,nabla_grad(tau0)))*r*dx \

        #  + (1/dt)*inner(lamb(phi)*(tau-tau0),supg*dot(u0, nabla_grad(S)))*r*dx \
        #  + 0.5*((inner(dot(lamb(phi)*u0,nabla_grad(tau)) \
        #  - dot(lamb(phi)*tau, trans(G1)) \
        #  - dot(lamb(phi)*G1, tau),supg*dot(u0, nabla_grad(S)))*r*dx) \
        #  + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tau,supg*dot(u0, nabla_grad(S)))*r*dx) \
        #  + 0.5*((inner(dot(lamb(phi)*u0,nabla_grad(tau0)) \
        #  - dot(lamb(phi)*tau0, trans(G1)) \
        #  - dot(lamb(phi)*G1, tau0),supg*dot(u0, nabla_grad(S)))*r*dx) \
        #  + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tau0,supg*dot(u0, nabla_grad(S)))*r*dx) \
        # - eta_p(phi)*inner(G1+trans(G1),supg*dot(u0, nabla_grad(S)))*r*dx \
         



Ftt = (1/dt)*lamb(phi)*(tautt-tau0tt)*zetatt*r*dx \
    + 0.5*((lamb(phi)*u0[1]*tautt.dx(1) + lamb(phi)*u0[0]*tautt.dx(0) \
    - 2*lamb(phi)*tautt*G1_tt)*zetatt*r*dx \
    + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tautt,zetatt)*r*dx) \
    + 0.5*((lamb(phi)*u0[1]*tau0tt.dx(1) + lamb(phi)*u0[0]*tau0tt.dx(0) \
    - 2*lamb(phi)*tau0tt*G1_tt)*zetatt*r*dx \
    + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tau0tt,zetatt)*r*dx ) \
    - 2*eta_p(phi)*G1_tt*zetatt*r*dx \
    + 0.5*0.1*inner(flux(u0,n,tautt)('+') - flux(u0,n,tautt)('-'),jump(zetatt))*r*dS \
    + 0.5*0.1*inner(flux(u0,n,tau0tt)('+') - flux(u0,n,tau0tt)('-'),jump(zetatt))*r*dS


    # + 0.5*lamb(phi)*inner((u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0)), supg*(u0[1]*tau0tt.dx(1) + u0[0]*tau0tt.dx(0)))*r*dx \
    # + 0.5*lamb(phi)*inner((u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0)), supg*(u0[1]*tautt.dx(1) + u0[0]*tautt.dx(0)))*r*dx

    # + ((1/dt)*lamb(phi)*(tautt-tau0tt)*(supg*(u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0)))*r*dx \
    # + 0.5*((lamb(phi)*u0[1]*tautt.dx(1) + lamb(phi)*u0[0]*tautt.dx(0) \
    # - 2*lamb(phi)*tautt*G1_tt)*(supg*(u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0)))*r*dx \
    # + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tautt,(supg*(u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0))))*r*dx) \
    # + 0.5*((lamb(phi)*u0[1]*tau0tt.dx(1) + lamb(phi)*u0[0]*tau0tt.dx(0) \
    # - 2*lamb(phi)*tau0tt*G1_tt)*(supg*(u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0)))*r*dx \
    # + inner((exp((0.05*lamb(phi)/eta_p(phi))*(tau0[0,0]+tau0[1,1]+tau0tt)))*tau0tt,(supg*(u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0))))*r*dx ) \
    # - 2*eta_p(phi)*G1_tt*(supg*(u0[1]*zetatt.dx(1) + u0[0]*zetatt.dx(0)))*r*dx) \




outn = as_vector([ngamma(phi0)[0],0,ngamma(phi0)[1]])

out = outer(outn, outn)
"""""""""""""""" maybe this below? """""""""""""""
G3D = as_tensor([[G1[0,0],0, G1[1,0]],
                    [0,G1_tt,0],
                    [G1[0,1],0,G1[1,1]]]) 

ns1 = (1/dt)*inner((rho(phi)*u - rho(phi)*u0), v)*r*dx \
    + inner(dot(rho(phi)*u0, nabla_grad(u)), v)*r*dx \
    + 2*inner(eta_s(phi)*epsilon_axi(u,x),epsilon_axi(v,x))*r*dx \
    - p0*div_axi(v,x)*r*dx \
    + inner(rho(phi)*g,v)*r*dx \
    + sig*mgrad(phi0)*inner((Identity(3) - out), epsilon_axi(v,x))*r*dx \
    + 2*inner((1-eta_s(phi))*epsilon_axi(u0,x),epsilon_axi(v,x))*r*dx \
    - inner((1-eta_s(phi))*(G3D+trans(G3D)),epsilon_axi(v,x))*r*dx \
    #     + sig*div_axi(ngamma(phi0),x)*mgrad(phi0)*inner(ngamma(phi0), v)*r*dx \
    # + inner((G - grad(u0)),R)*r*dx \
    # + inner((G_tt - u0[0]/r),R_tt)*r*dx \

tau_axi = as_tensor([[tau0[0,0], 0, tau0[1,0]],
                        [0, tau0tt, 0],
                        [tau0[1,0], 0, tau0[1,1]]])

# tau_axi = as_vector([(1/r)*(r*tau0[0,0]).dx(0)+(tau0[1,0]).dx(1)-tau0tt/r,
#                     (1/r)*(r*tau0[1,0]).dx(0)+(tau0[1,1]).dx(1)])

ns1 += c(phi)*inner(tau_axi, nabla_grad_axi(v, x))*r*dx \

ns2 = (1/rho(phi))*(dot(nabla_grad(p),nabla_grad(q)) - dot(nabla_grad(p0),nabla_grad(q)))*r*dx \
    + (1/dt)*div_axi(u_,x)*q*r*dx \

ns3 = inner(u,v)*r*dx - inner(u_,v)*r*dx \
    + (dt/rho(phi))*inner(nabla_grad(p_-p0),v)*r*dx \
    # + inner((G - grad(u_)),R)*r*dx \
    # + inner((G_tt - u_[0]/r),R_tt)*r*dx

bcs=[]

tol = 1.0e-4

tau =Function(Ts)
tautt =Function(Ttt)

t=0
q=0

for k in range(int(8)):

    a = assemble(lhs(F_levelset))
    L = assemble(rhs(F_levelset))

    solve(a, phi.vector(), L, 'gmres', 'default') 

    phi0.assign(phi)

for n in tqdm(range(num_steps)):

    t+=dt
    q+=1

    a = assemble(lhs(F_levelset))
    L = assemble(rhs(F_levelset))

    solve(a, phi.vector(), L, 'gmres', 'default') 

    if n % rein_div == 0:

        phi0.assign(phi)

        solve(F_grad == 0, phigrad)

        for n in range(num_steps_rein):

            solve(F_rein == 0, phi_rein ,solver_parameters={"newton_solver": {"linear_solver": 'gmres', "preconditioner": 'default',\
                                     "maximum_iterations": 20, "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                                        form_compiler_parameters={"optimize": True})

            phi00.assign(phi0)
            phi0.assign(phi_rein)

        phi.assign(phi_rein)

    else:

        phi0.assign(phi)

    outn = as_vector([ngamma(phi0)[0],0,ngamma(phi0)[1]])

    out = outer(outn, outn)

    G3D = as_tensor([[G1[0,0],0, G1[1,0]],
                    [0,G1_tt,0],
                    [G1[0,1],0,G1[1,1]]]) 

    tau_axi = as_tensor([[tau0[0,0], 0, tau0[1,0]],
                        [0, tau0tt, 0],
                        [tau0[1,0], 0, tau0[1,1]]])

    # tau_axi = as_vector([(1/r)*(r*tau0[0,0]).dx(0)+(tau0[1,0]).dx(1)-tau0tt/r,
    #                 (1/r)*(r*tau0[1,0]).dx(0)+(tau0[1,1]).dx(1)])

    A1=assemble(lhs(ns1))
    M1=assemble(rhs(ns1))
    [bc.apply(A1) for bc in bcu]
    [bc.apply(M1) for bc in bcu]
    solve(A1, u_.vector(), M1, 'gmres', 'default') 

    # u_, G_, G__tt = split(w_)

    M2=assemble(rhs(ns2))
    A2=assemble(lhs(ns2))
    solve(A2, p_.vector(), M2, 'gmres', 'default') 
        
    A3=assemble(lhs(ns3))
    M3=assemble(rhs(ns3))
    [bc.apply(A3) for bc in bcu]
    [bc.apply(M3) for bc in bcu]
    solve(A3, u_.vector(), M3, 'gmres', 'default') 

    u0.assign(u_)
    p0.assign(p_)

    A_devssg = assemble(lhs(F_devssg))
    M_devssg = assemble(rhs(F_devssg))

    solve(A_devssg, G1.vector(), M_devssg, 'gmres', 'default')

    A_tt_devssg = assemble(lhs(F_tt_devssg))
    M_tt_devssg = assemble(rhs(F_tt_devssg))

    solve(A_tt_devssg, G1_tt.vector(), M_tt_devssg, 'gmres', 'default')

    # u0, Gb, Gb_tt = split(w0)

    A = assemble(lhs(F_stress))
    M = assemble(rhs(F_stress))

    solve(A, tau.vector(), M, 'gmres', 'default')

    Att = assemble(lhs(Ftt))
    Mtt = assemble(rhs(Ftt))

    solve(Att, tautt.vector(), Mtt, 'gmres', 'default')

    tau0.assign(tau)
    tau0tt.assign(tautt)



    if q % qdiv == 0:

        # ua, Ga, Ga_tt = w0.split(deepcopy=True)

        file_phi.write(phi, t)
        file_u.write(u_, t)
        # file_p.write(p_, t)
        file_tau.write(tau, t)



    phi0.assign(phi)

    area = assemble(phi0*dx)
    area_x2 = 2*area

    x_com = assemble(Expression("x[0]", degree = 1)*phi0*dx)/area
    y_com = assemble(Expression("x[1]", degree = 1)*phi0*dx)/area

    Pa = 2.0*sqrt(np.pi*area_x2)
    Pb = 2*assemble(mgrad(phi0)*dx)
    circ = Pa/Pb

    u_rise = assemble(u0[0]*phi0*dx)/area
    v_rise = assemble(u0[1]*phi0*dx)/area


    timeseries = [t, area_x2, x_com, y_com, circ, u_rise, v_rise]

    if rank == 0:
        with open((f"droplet_ve401/data.csv"), 'a') as csvfile:
            f = csv.writer(csvfile, delimiter='\t',lineterminator='\n',)
            f.writerow(timeseries)

# if rank == 0:

#     plt.figure(num=None, figsize=(8, 8), dpi=100, facecolor='w', edgecolor='k')
#     plot(phi,mode = 'contour',levels=[0],linestyles='dashed',colors='red')
#     plot(phi,cmap='jet')
#     plt.title('T = 0')
#     plt.savefig('/media/lamb(phi)ll/Samsung_T5/bubble_para/levelset_results/phi.png')

#     plt.figure(num=None, figsize=(8, 8), dpi=100, facecolor='w', edgecolor='k')
#     plt.subplot(121)
#     plot(u_)
#     plt.subplot(122)
#     plot(p_)
#     plt.savefig('/media/lamb(phi)ll/Samsung_T5/bubble_para/levelset_results/u_p.png')


# print('Level set ic: ', assemble(sqrt (phiic.dx (0)**2 + phiic.dx (1)**2)*r*dx))
# print('Level set fin: ', assemble(sqrt (phi.dx (0)**2 + phi.dx (1)**2)*r*dx))