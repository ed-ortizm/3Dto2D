#!/usr/bin/python
import numpy as np
# arguments for the script
#https://www.tutorialspoint.com/python/python_command_line_arguments.htm
import sys
# sys.arg[0] is the name of the script itself (careful)
################################ code starts ###################################
# raising an assetion error when the script is not provided with the right number of n_arguments or the right format
#https://stackoverflow.com/questions/1964126/what-exception-to-raise-if-wrong-number-of-arguments-passed-in-to-args

n_arguments = len(sys.argv)
if n_arguments == 3:
    cube_name = sys.argv[1]
    assert (cube_name[-5:] == ".fits"), "The format of the cube is .fits: cube_name.fits (lower case)"
    filter_name= sys.argv[2]
    assert (filter_name[-4:] == ".dat"), "The format of the filter is .fits: filter_name.dat (lower case)"
#assert (filter_name[-5:] == ".fits"), "Please remember that the format for the filter name is fits: fiter_name.fits, (lower case)"
# or wavelength range
elif n_arguments == 4:
    cube_name = sys.argv[1]
    assert (cube_name[-5:] == ".fits"), "The format of the cube is .fits: cube_name.fits (lower case)"
# ADD CONDITIONALS FOR FORBIDEN VALUES
    lambda1 = float(sys.argv[2])
    lambda2 = float(sys.argv[3])
    assert (0 < lambda1) & (0 < lambda2), "wavelength must be positive!!"
    assert (lambda1 < lambda2) , "First wavelength must be smaller than the second one"
else:
    assert False, "# of args must be 2 (cube and filter name) OR 3 (cube name, lower wavelength and higher wavelength)"
 # Taking a filter file and converting A to nm

class filter_handler():
    def __init__(self,file):
        self.file = np.loadtxt(file)
    def wavelength(self):
        lambdas = self.file[:,0] * 0.1
        return lambdas
    def energy(self):
        # Following the instructions, I will just multiply by the wavelength
        photons = self.file[:,1]
        energies = self.wavelength() *photons
        norm = np.trapz(energies,self.wavelength())
        n_energies = energies/norm
        return n_energies


#working
filter = filter_handler(filter_name)
E = filter.energy()
x = filter.wavelength()
print(E)
print(x)
print(np.trapz(E,x))
