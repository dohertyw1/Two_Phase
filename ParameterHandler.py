class bcolors():

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class ParameterHandler():

    def __init__(self):

        pass

    def print_parameters(self):

        print("Inner fluid density is [{}] and outer fluid density is [{}].".format(self.rho1, self.rho2))
        print("Inner fluid viscosity is [{}] and outer fluid viscosity is [{}].".format(self.mu1, self.mu2))
        print("Acceleration due to gravity is: {}.".format(self.grav))
        print("Surface tension coefficient: {}.".format(self.sigma))

    def one_by_one(self):

            try:

                self.rho1, self.rho2 = [float(x) for x in input("Enter the density of the inner fluid then the outer fluid: ").split()]
                self.mu1, self.mu2 = [float(x) for x in input("Enter the viscosity of the inner fluid then the outer fluid: ").split()]

            except ValueError:
                
                print(f"{bcolors.FAIL}Error: Please enter two values with a space in between them.{bcolors.ENDC}")
                quit()

            self.grav = [float(x) for x in input("Enter the aceleration due to gravity: ").split()]
            self.sigma = [float(x) for x in input("Enter the surface tension coefficient: ").split()]

            self.print_parameters()

    def read_in(self, file):

        self.parameters = {}

        try:

            with open(f'{file}.txt') as File:

                for line in File:
                    
                    index = line.find(' = ')
                    self.parameters[line[0:index]] = float(line[index+2:].replace('\n', ''))
        
        except FileNotFoundError:

            print(f"{bcolors.FAIL}Error: Parameter file not found. Please ensure it is located in the correct directory.{bcolors.ENDC}")
            quit()

        self.T = self.parameters['T']
        self.cox1 = self.parameters['cox1']
        self.coy1 = self.parameters['coy1']
        self.cox2 = self.parameters['cox2']
        self.coy2 = self.parameters['coy2']
        self.sizex = int(self.parameters['sizex'])
        self.sizey = int(self.parameters['sizey'])
        self.centrex = self.parameters['centrex']
        self.centrey = self.parameters['centrey']
        self.radius = self.parameters['radius']
        self.rho1 = self.parameters['rho_inner']
        self.rho2 = self.parameters['rho_outer']
        self.grav = self.parameters['gravity']
        self.curvature = self.parameters['surface_tension_coefficient']

        if (self.fluid == 'Newtonian'):

            self.mu1 = self.parameters['mu_inner']
            self.mu2 = self.parameters['mu_outer']

        elif (self.fluid == 'Viscoelastic' and self.dimensional == 'Dim'):
            
            self.etas_in = self.parameters['etas_inner']
            self.etas_out = self.parameters['etas_outer']
            self.etap_out = self.parameters['etas_outer']
            self.lamb_in = self.parameters['relaxation_time_in']
            self.lamb_out = self.parameters['relaxation_time_out']
            self.gmf = self.parameters['giesekus_mobility_factor']

        if (self.rank == 0):

            for pair in self.parameters.items():

                print(pair)

    def begin(self):

        if (self.rank == 0):

            print('Hello. Welcome to my rising bubble solver.\nWould you like to input your parameters one by one [0], or read them in from the Parameters.txt file [1]?')

            choice = int(input())

            if (choice == 0):

                print('Enter one-by-one.')

            elif (choice == 1):

                print('Read from file.')

        else:

            choice = None

        choice = self.comm.bcast(choice, root = 0)

        if (choice == 0):

            self.one_by_one()

        elif (choice == 1):

            self.read_in('Newtonian')