from fenics import *
import ufl as uf
from Funcs import *

class Flow():

    """Choose linear solvers for solving equations"""
    def __init__(self):
        
        self.ls_linear_solver = 'gmres'
        self.ls_preconditioner = 'default'
        self.lsr_linear_solver = 'gmres'
        self.lsr_preconditioner = 'default'
        self.ns_linear_solver = 'gmres'
        self.ns_preconditioner = 'default'
        self.oldroyd_b_linear_solver = 'gmres'
        self.oldroyd_b_preconditioner = 'default'
        
    """Construct weak form for level set equation"""
    def ls_form(self, phi, phi0, psi, u0):

        if (self.ls_scheme == 'Euler'):

            F = (phi/self.dt)*psi*dx \
              - (phi0/self.dt)*psi*dx \
              + inner(u0, grad(phi))*psi*dx \

        elif (self.ls_scheme == 'CN'):

            F = (phi/self.dt)*psi*dx \
              - (phi0/self.dt)*psi*dx \
              + 0.5*inner(u0, grad(phi))*psi*dx \
              + 0.5*inner(u0, grad(phi0))*psi*dx \

        elif (self.ls_scheme == 'CN_SUPG'):

            scaling = 0.1
            toler = 1.0e-3

            F = (phi/self.dt)*psi*dx \
              - (phi0/self.dt)*psi*dx \
              + inner(u0, grad(phi))*psi*dx \
              + ((phi - phi0)/self.dt \
              + inner(u0, grad(phi))) \
              * scaling*self.hmin \
              / uf.Max(2.0*sqrt(inner(u0, u0)),toler/self.hmin) \
              * inner(u0, grad(psi))*dx
        
        self.a_ls = lhs(F)
        self.m_ls = rhs(F)

        self.A_ls = PETScMatrix()
        self.M_ls = PETScVector()

    """Construct weak form for non-conservative level set reinitialisation"""
    def nclsr_form(self, phi_rein, phi0, psi, sign):

        F = (phi_rein/self.dtau)*psi*dx - (phi0 /self.dtau)*psi*dx \
            - sign*(1-sqrt(dot(grad(phi0),grad(phi0))))*psi*dx \
            + self.alph*inner(grad(phi0),grad(psi))*dx

        self.F_nclsr = F

    """Construct weak form for conservative level set reinitialisation"""
    def clsr_form(self, phi_rein, phi0, psi, phin):

        F = (phi_rein - phi0)/self.dtau*psi*dx

        if (self.lsr_method == 'old_method'):

            terms = - phi_rein*(1.0 - phi_rein)*inner(grad(psi), phin)*dx \
                    + self.eps*inner(grad(phi_rein), grad(psi))*dx  

        elif (self.lsr_method == 'new_method'):
            
            terms = - 0.5*dot(grad(psi),(phi_rein+phi0)*phin)*dx \
                    + dot(phi_rein*phi0*phin,grad(psi))*dx \
                    + (self.eps/2)*(dot(grad(phi_rein)+grad(phi0),phin)*dot(grad(psi),phin))*dx

        F += terms

        self.F_clsr = F

    """Construct the weak form for the Oldroyd-B viscoelastic constitutive equation."""
    def oldroyd_b_form(self, tau, tau0, zeta, u0):

        F = (1/self.dt)*inner(tau-tau0,zeta)*dx \
          + 0.5*(inner(dot(u0,nabla_grad(tau)) \
          - dot(tau, trans(grad(u0))) \
          - dot(grad(u0), tau),zeta)*dx \
          + (1/self.lamb1)*inner(tau-self.Id,zeta)*dx) \
          + 0.5*(inner(dot(u0,nabla_grad(tau0)) \
          - dot(tau0, trans(grad(u0))) \
          - dot(grad(u0), tau0),zeta)*dx \
          + (1/self.lamb1)*inner(tau0-self.Id,zeta)*dx)

        self.a_ob = lhs(F)
        self.m_ob = rhs(F)

        self.A_ob = PETScMatrix()
        self.M_ob = PETScVector()

    """Construct weak form for navier stokes equations (IPCS scheme)"""
    def ns_form(self, rho, rho0, mu, u, u0, u_, 
                p, p0, p_, phi, v, q, phig, phic,
                tau0, eta_s, eta_p):

        if (self.method == 'NCons'):

            curv_term  = Constant(self.sigma)*mgrad(phi)*inner((self.Id \
                       - outer(ngamma(phi), ngamma(phi))), epsilon(v)) \
                       * delta_func(phi,self.eps)*dx

        elif (self.method == 'Cons'):

            curv_term = Constant(self.sigma)*mgrad(phi)*inner((self.Id \
                       - outer(self.phin, self.phin)), epsilon(v))*dx
            
                        # - Constant(self.sigma)*inner(phig*phic,v)*dx

        if (self.fluid != 'Viscoelastic'):

            ns1 = (1/self.dt)*inner(rho*u - rho0*u0, v)*dx \
                + inner(dot(rho*u0, nabla_grad(u)), v)*dx \
                + Constant(2.0)*inner(mu*epsilon(u), epsilon(v))*dx \
                - p0*div(v)*dx \
                + inner(rho*self.g,v)*dx \

        elif (self.fluid == 'Viscoelastic'):

            ns1 = (1/self.dt)*inner(rho*u - rho0*u0, v)*dx \
                + inner(dot(rho*u0, nabla_grad(u)), v)*dx \
                + inner(eta_s*grad(u), grad(v))*dx \
                + (eta_p/self.lamb1)*inner(tau0-self.Id, grad(v))*dx \
                - p0*div(v)*dx \
                + inner(rho*self.g,v)*dx \

        ns2 = (1/rho)*(dot(grad(p),grad(q)) \
            - dot(grad(p0),grad(q)))*dx \
            + (1/self.dt)*div(u_)*q*dx \

        ns3 = inner(u,v)*dx - inner(u_,v)*dx \
            + (self.dt/rho)*inner(grad(p_-p0),v)*dx

        ns1 += curv_term

        self.a_ns1 = lhs(ns1)
        self.m_ns1 = rhs(ns1)

        self.A_ns1 = PETScMatrix()
        self.M_ns1 = PETScVector()

        self.a_ns2 = lhs(ns2)
        self.m_ns2 = rhs(ns2)

        self.A_ns2 = PETScMatrix()
        self.M_ns2 = PETScVector()

        self.a_ns3 = lhs(ns3)
        self.m_ns3 = rhs(ns3)

        self.A_ns3 = PETScMatrix()
        self.M_ns3 = PETScVector()

    """Solve the level set equation"""
    def ls_solve(self):

        assemble(self.a_ls, tensor = self.A_ls)
        assemble(self.m_ls, tensor = self.M_ls)

        solve(self.A_ls, self.phi.vector(), self.M_ls, self.ls_linear_solver, self.ls_preconditioner)

    """Solve the non-conservative level set reinitialisation equation"""
    def nclsr_solve(self, phi0, phi00, phi, phi_rein):

        phi0.assign(phi)

        self.sign = sgn(phi, self.eps1)

        self.E = 0
        self.tol = 1.0e-4
        
        for self.n in range(self.rein_steps):
            
            solve(self.F_nclsr == 0, phi_rein, solver_parameters={"newton_solver": {"linear_solver": 'gmres', \
                                               "preconditioner": 'default', "maximum_iterations": 20, \
                                               "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                                               form_compiler_parameters={"optimize": True})

            phi00.assign(phi0)

            error = (((phi_rein - phi0)/self.dtau)**2)*dx
            self.E = sqrt(assemble(error))

            if self.E < self.tol:
                break

            phi0.assign(phi_rein)

        phi.assign(phi_rein)

    """Solve the conservative level set reinitialisation equation"""
    def clsr_solve(self, phi0, phi00, phi, phi_rein, phin, psin):
        
        phi0.assign(phi)

        F_norm = inner((phin-new_norm(phi0)),psin)*dx
        
        solve(F_norm == 0, phin, solver_parameters={"newton_solver": {"linear_solver": 'gmres', \
                                               "preconditioner": 'default', "maximum_iterations": 20, \
                                               "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                                               form_compiler_parameters={"optimize": True})

        self.E = 0
        self.tol = 1.0e-4
        
        for self.n in range(self.rein_steps):
            
            solve(self.F_clsr == 0, phi_rein, solver_parameters={"newton_solver": {"linear_solver": 'gmres', \
                                               "preconditioner": 'default', "maximum_iterations": 20, \
                                               "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                                               form_compiler_parameters={"optimize": True})

            phi00.assign(phi0)

            error = (((phi_rein - phi0)/self.dtau)**2)*dx
            self.E = sqrt(assemble(error))

            if self.E < self.tol:
                break

            phi0.assign(phi_rein)

        phi.assign(phi_rein)

    """Solve the Oldroyd-b equation."""
    def oldroydb_solve(self):

        assemble(self.a_ob, tensor = self.A_ob)
        assemble(self.m_ob, tensor = self.M_ob)

        solve(self.A_ob, self.tau.vector(), self.M_ob, self.oldroyd_b_linear_solver, self.oldroyd_b_preconditioner)

    """Solve the navier stokes equations"""
    def ns_solve(self, bc_ns, u_, p_, phig, phic, psig, psic, phi):

        # F_grad = inner((phig-grad(phi)),psig)*dx

        # solve(F_grad == 0, phig)

        # F_curv = inner(phic,psic)*dx - inner(phig/(sqrt(phig[0]**2+phig[1]**2)),grad(psic))*dx + inner(grad(phic),grad(psic))*dx

        # solve(F_curv == 0, phic)

        assemble(self.a_ns1, tensor = self.A_ns1)
        assemble(self.m_ns1, tensor = self.M_ns1)

        for bc in bc_ns:
            bc.apply(self.A_ns1)
            bc.apply(self.M_ns1)
            
        solve(self.A_ns1, u_.vector(), self.M_ns1, self.ns_linear_solver, self.ns_preconditioner)

        assemble(self.a_ns2, tensor = self.A_ns2)
        assemble(self.m_ns2, tensor = self.M_ns2)

        solve(self.A_ns2, p_.vector(), self.M_ns2, self.ns_linear_solver, self.ns_preconditioner)

        assemble(self.a_ns3, tensor = self.A_ns3)
        assemble(self.m_ns3, tensor = self.M_ns3)

        solve(self.A_ns3, u_.vector(), self.M_ns3, self.ns_linear_solver, self.ns_preconditioner)
