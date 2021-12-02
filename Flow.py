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
        self.constitutive_equation_linear_solver = 'gmres'
        self.constitutive_equation_preconditioner = 'default'
        
    """Construct weak form for level set equation"""
    def ls_form(self, phi, phi0, psi, u0, fn):

        if (self.ls_scheme == 'Euler'):

            if (self.element == 'CG'):

                F = (phi/self.dt)*psi*dx \
                - (phi0/self.dt)*psi*dx \
                + inner(u0, grad(phi))*psi*dx \

            elif (self.element == 'DG'):

                un = abs(dot(u0('+'), fn('+')))

                F = (1/self.dt)*(phi-phi0)*psi*dx \
                  - dot ( grad ( psi ), u0 * phi ) * dx \
                  + (dot(u0('+'), jump(psi, fn))*avg(phi) \
                  + 0.5*un*dot(jump(phi, fn), jump(psi, fn)))*dS \

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
    def clsr_form(self, phi_rein, phi0, psi, phin, fn):

        F = (phi_rein - phi0)/self.dtau*psi*dx

        if (self.element == 'CG'):

            if (self.lsr_method == 'old_method'):

                terms = - phi_rein*(1.0 - phi_rein)*inner(grad(psi), phin)*dx \
                        + self.eps*inner(grad(phi_rein), grad(psi))*dx  

            elif (self.lsr_method == 'new_method'):
                
                terms = - 0.5*dot(grad(psi),(phi_rein+phi0)*phin)*dx \
                        + dot(phi_rein*phi0*phin,grad(psi))*dx \
                        + (self.eps/2)*(dot(grad(phi_rein)+grad(phi0),phin)*dot(grad(psi),phin))*dx
        
            F += terms

        elif (self.element == 'DG'):

            nn0 = abs(dot(phin('+'), fn('+')))
            self.penalty = 50

            compressive_terms = - dot ( grad ( psi ), phin * phi_rein*(1-phi_rein) ) * dx \
                              + (dot(phin('+'), jump(psi, fn))*avg(phi_rein) \
                              + 0.5*nn0*dot(jump(phi_rein, fn), jump(psi, fn)))*dS \
                              - (dot(phin('+'), jump(psi, fn))*avg(phi_rein*phi_rein) \
                              + 0.5*nn0*dot(jump(phi_rein*phi_rein, fn), jump(psi, fn)))*dS \

            diffusive_terms = self.eps*(inner(grad(phi_rein), grad(psi))*dx \
                            + (self.penalty/self.hmin)*dot(jump(psi, fn), jump(phi_rein, fn))*dS \
                            - dot(avg(grad(psi)), jump(phi_rein, fn))*dS \
                            - dot(jump(psi, fn), avg(grad(phi_rein)))*dS)

            F += compressive_terms
            F += diffusive_terms

        self.F_clsr = F

    """Construct the weak form for the Oldroyd-B viscoelastic constitutive equation."""
    def constitutive_equation_form(self, tau, tau0, zeta, u0):

        if (self.constitutive_equation == 'OB') or (self.constitutive_equation == 'Giesekus'):

            F = inner(tau, zeta)*dx \

            if (self.dimensional == 'Dim'):

                F += (1/self.dt)*inner(self.lamb*(tau-tau0),zeta)*dx \
                  + inner(dot(self.lamb*u0,nabla_grad(tau)) \
                  - dot(self.lamb*tau, trans(grad(u0))) \
                  - dot(self.lamb*grad(u0), tau),zeta)*dx \
                  - self.eta_p*inner(grad(u0)+grad(u0).T,zeta)*dx
            
                if (self.constitutive_equation == 'Giesekus'):

                    F += (self.gmf*self.lamb/self.eta_p)*inner(tau0*tau0,zeta)*dx \

            elif (self.dimensional == 'NonDim'):

                F += (1/self.dt)*inner(self.Wi*(tau-tau0),zeta)*dx \
                  + inner(dot(self.Wi*u0,nabla_grad(tau)) \
                  - dot(self.Wi*tau, trans(grad(u0))) \
                  - dot(self.Wi*grad(u0), tau),zeta)*dx \
                  - (1-self.beta)*inner(grad(u0)+grad(u0).T,zeta)*dx

                if (self.constitutive_equation == 'Giesekus'):

                    + (self.gmf*self.Wi/(1-self.beta))*inner(tau0*tau0,zeta)*dx \

        self.a_ce = lhs(F)
        self.m_ce = rhs(F)

        self.A_ce = PETScMatrix()
        self.M_ce = PETScVector()

    """Construct weak form for navier stokes equations (IPCS scheme)"""
    def ns_form(self, u, u0, u_, p, p0, p_, phi, v, q, tau0):

        if (self.method == 'NCons'):

            curv_term  = self.curvature*mgrad(phi)*inner((self.Id \
                       - outer(ngamma(phi), ngamma(phi))), epsilon(v)) \
                       * delta_func(phi,self.eps)*dx

        elif (self.method == 'Cons'):

            curv_term = self.curvature*mgrad(phi)*inner((self.Id \
                       - outer(self.phin, self.phin)), epsilon(v))*dx
            
        if (self.fluid == 'Newtonian' or self.fluid == 'GNF'):

            ns1 = (1/self.dt)*inner(self.rho*u - self.rho0*u0, v)*dx \
                + inner(dot(self.rho*u0, nabla_grad(u)), v)*dx \
                + Constant(2.0)*inner(self.mu*epsilon(u), epsilon(v))*dx \
                - p0*div(v)*dx \
                + inner(self.rho*self.g,v)*dx \

        elif (self.fluid == 'Viscoelastic'):

            if (self.dimensional == 'Dim'):

                ns1 = (1/self.dt)*inner(self.rho*u - self.rho0*u0, v)*dx \
                    + inner(dot(self.rho*u0, nabla_grad(u)), v)*dx \
                    + self.eta_s*inner(grad(u), grad(v))*dx \
                    - p0*div(v)*dx \
                    + inner(self.rho*self.g,v)*dx \

                if (self.constitutive_equation == 'OB'):

                    ns1 += inner(tau0, grad(v))*dx \

            elif (self.dimensional == 'NonDim'):

                ns1 = (1/self.dt)*inner(self.rho*u - self.rho0*u0, v)*dx \
                    + inner(dot(self.rho*u0, nabla_grad(u)), v)*dx \
                    + (self.beta/self.Re)*inner(grad(u), grad(v))*dx \
                    - p0*div(v)*dx \
                    + inner(self.rho*self.g,v)*dx \

                if (self.constitutive_equation == 'OB'):

                    ns1 += (1/self.Re2)*inner(tau0, grad(v))*dx \

        ns2 = (1/self.rho)*(dot(grad(p),grad(q)) \
            - dot(grad(p0),grad(q)))*dx \
            + (1/self.dt)*div(u_)*q*dx \

        ns3 = inner(u,v)*dx - inner(u_,v)*dx \
            + (self.dt/self.rho)*inner(grad(p_-p0),v)*dx

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

        if (self.element == 'CG'):

            F_norm = inner((phin-new_norm(phi0)),psin)*dx

            solve(F_norm == 0, phin, solver_parameters={"newton_solver": {"linear_solver": 'gmres', \
                                        "preconditioner": 'default', "maximum_iterations": 20, \
                                        "absolute_tolerance": 1e-8, "relative_tolerance": 1e-6}}, \
                                        form_compiler_parameters={"optimize": True})

        if (self.element == 'DG'):

            phin = ngamma(phi0)
        
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
    def constitutive_equation_solve(self):

        assemble(self.a_ce, tensor = self.A_ce)
        assemble(self.m_ce, tensor = self.M_ce)

        solve(self.A_ce, self.tau.vector(), self.M_ce, self.constitutive_equation_linear_solver, self.constitutive_equation_preconditioner)

    """Solve the navier stokes equations"""
    def ns_solve(self, bc_ns, u_, p_):

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

    # if (self.constitutive_equation == 'OB-Conf'):

    #     F = (1/self.dt)*inner(tau-tau0,zeta)*dx \
    #       + inner(dot(u0,nabla_grad(tau)) \
    #       - dot(tau, trans(grad(u0))) \
    #       - dot(grad(u0), tau),zeta)*dx \

    #     if (self.dimensional == 'Dim'):

    #         F += (1/self.lamb)*inner(tau - self.Id, zeta)*dx

        # elif (self.dimensional == 'NonDim'):

        #     F += (1/self.Wi)*inner(tau - self.Id, zeta)*dx

    # if (self.constitutive_equation == 'Gie'):

    #     if (self.dimensional == 'Dim'):

    #         pass

    #     elif (self.dimensional == 'NonDim'):

    #         pass

    # if (self.constitutive_equation == 'Gie-Conf'):

    #     if (self.dimensional == 'Dim'):

    #         pass

    #     elif (self.dimensional == 'NonDim'):

    #         pass

    # elif (self.constitutive_equation == 'Giesekus'):

    #     if (self.stability == 'DEVSSG-DG'):

    #         pass

    #     elif (self.stability == None):

    #         F = (1/self.dt)*inner(tau-tau0,zeta)*dx \
    #         + (inner(dot(u0,nabla_grad(tau0)) \
    #         - dot(tau0, trans(grad(u0))) \
    #         - dot(grad(u0), tau0),zeta)*dx \
    #         + (1/Wi)*inner((tau0-self.Id) \
    #         + gmf*(tau0-self.Id)*(tau0-self.Id), zeta)*dx) \
    #         + inner(poly_flux(u0,self.facet_normal,tau0)('+') - poly_flux(u0,self.facet_normal,tau0)('-'),jump(zeta))*dS

