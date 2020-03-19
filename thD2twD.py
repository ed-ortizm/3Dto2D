#!/usr/bin/python
import numpy as np
from scipy import interpolate
from astropy.table import Table
# arguments for the script
#https://www.tutorialspoint.com/python/python_command_line_arguments.htm
import sys
# sys.arg[0] is the name of the script itself (careful)
################################ code starts ###################################
# raising an assetion error when the script is not provided with the right number of n_arguments or the right format
#https://stackoverflow.com/questions/1964126/what-exception-to-raise-if-wrong-number-of-arguments-passed-in-to-args
# I used this as a reference for the keywords of the header:
# https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html
from astropy.io import fits
from astropy import units as u
#from spectral_cube import SpectralCube
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

class Cube_handler():
    def __init__(self,cube):
        self.hdul = fits.open(cube)
        self.unit = u.Unit(self.hdul[1].header['BUNIT'])
        self.hdul.close()
    # convert the units in the fits file to 'J/(nm m2 s)'
    def u_convert(self):
        req_unit = self.unit.to(u.Unit('J/(nm m2 s)'))
        return req_unit
    # return the data in the cube with fluxes converted
    # to 'J/(nm m2 s)'
    def cube(self):
        return self.hdul[1].data * self.u_convert()

class Filter_handler():
    def __init__(self,filter):
        self.filter = np.loadtxt(filter)
    def wavelength(self):
        #converting filter wavelength to nm
        lambdas = self.filter[:,0] * 0.1
        return lambdas
    # phtons to energy and normalizing the integral of the filter.
    def energy(self):
        # Multiply by the wavelength.
        photons = self.filter[:,1]
        energies = self.wavelength() * photons
        norm = np.trapz(energies,self.wavelength())
        n_energies = energies/norm
        return n_energies

    #def wave_array(self):
        #hdul = fits.open("../edgar.fits")

# Function to compute the array of wavelengths (5) for the spectra
def lamb_s(cube):
    hdul = fits.open(cube)
    # this is the number of slices present in the cube
    i_max = hdul[1].data.shape[0]
    #Following instructions
    CRVAL3 = hdul[1].header['CRVAL3']
    CD3_3 = hdul[1].header['CD3_3']
    CRPIX3 = hdul[1].header['CRPIX3']
    lamb = np.array([CRVAL3 + CD3_3*(i-CRPIX3) for i in range(i_max)])
    hdul.close()
    return lamb*0.1

def lamb_inter(arr_1,arr_2):
    stack = np.concatenate((arr_1,arr_2))
    # np.unique eliminates the duplicates and returns the array sorted :)
    return np.unique(stack)
#working
#cube = Cube_handler(cube_name)
#print(cube.unit)
#type(print(cube.u_convert()))
filter = Filter_handler(filter_name)
#E = filter.energy()
x = filter.wavelength()
#print(E)
#print(x)
#print(np.trapz(E,x))
print(lamb_s(cube_name))
print(x)
print(lamb_inter(lamb_s(cube_name),x))
