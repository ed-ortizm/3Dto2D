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
# Function to compute the array of wavelengths (5) for the spectra
    def lamb_s(self):
        # this is the number of slices present in the cube
        i_max = self.hdul[1].data.shape[0]
        #Following instructions
        CRVAL3 = self.hdul[1].header['CRVAL3']
        CD3_3 = self.hdul[1].header['CD3_3']
        CRPIX3 = self.hdul[1].header['CRPIX3']
        lamb = np.array([CRVAL3 + CD3_3*(i-CRPIX3) for i in range(i_max)])
        hdul.close()
        return lamb*0.1

    def cube_interpolate(self,interval):
        # axis = 0 since this is the one containing the slices of the cube
        f = interpolate.interp1d(self.lamb_s(),self.cube(),axis=0,fill_value='extrapolate')
        return f(interval)

class Filter_handler():
    def __init__(self,filter):
        self.filter = np.loadtxt(filter)
    def lamb_f(self):
        #converting filter wavelength to nm
        lambdas = self.filter[:,0] * 0.1
        return lambdas
    # phtons to energy and normalizing the integral of the filter.
    def energy(self):
        # Multiply by the wavelength.
        photons = self.filter[:,1]
        energies = self.lamb_f() * photons
        norm = np.trapz(energies,self.lamb_f())
        n_energies = energies/norm
        return n_energies

    def filter_interpolate(self,interval):
        f = interpolate.interp1d(self.lamb_f(),self.energy(),fill_value='extrapolate')
        return f(interval)



def lamb_inter(arr_1,arr_2):
    stack = np.concatenate((arr_1,arr_2))
    # np.unique eliminates the duplicates and returns the array sorted :)
    return np.unique(stack)

def image(lambdas,filter,cube):
    # Computing the image
    Tf = cube.T*filter_name
    Tf = Tf.T
    flux_filter = np.trapz(Tf,lambdas,axis=0)
    # Creating the fits file for the image
    # https://python4astronomers.github.io/astropy/fits.html
    # Encapsulating the data
    hdu = fits.PrimaryHDU()
    hdu.data = flux_filter
    # header keywords
    hdu.header['OBSERVER'] = 'Edgar Ortiz Test'
    hdu.header['NAXIS']  = 2 #/ number of data axes'
    hdu.header['NAXIS1'] = flux_filter.shape[0]# '/ length of data axis 1'
    hdu.header['NAXIS2'] = flux_filter.shape[1]#'/ length of data axis 2'
    hdu.header['BUNIT']  = 'J/(nm m2 s)'
    # I just took the values from the data cube
    hdu.header['CRPIX1']  = 201.678490290293# / Pixel coordinate of reference point'
    hdu.header['CRPIX2']  = 211.082586048144# / Pixel coordinate of reference point'
    hdu.header['CD1_1']  = -5.55555555555556E-05# / Coordinate transformation matrix element'
    hdu.header['CD1_2']  = 0.# / Coordinate transformation matrix element'
    hdu.header['CUNIT1']  = 'deg     '#           / Units of coordinate increment and value'
    hdu.header['CUNIT2']  = 'deg     '#           / Units of coordinate increment and value'
    hdu.header['CTYPE1']  = 'RA---TAN'#           / Right ascension, gnomonic projection'
    hdu.header['CTYPE2']  = 'DEC--TAN'#           / Declination, gnomonic projection'
    hdu.header['CRVAL1']  = 136.957401
    hdu.header['CRVAL2']  = 1.02961
#    hdu.header['']  = ''
    # Writing the image
    hdu.writeto('../image.fits')
    #
#working
#cube = Cube_handler(cube_name)
#print(cube.unit)
#type(print(cube.u_convert()))
filter = Filter_handler(filter_name)
#E = filter.energy()
x = filter.lamb_f()
#print(E)
#print(x)
#print(np.trapz(E,x))
print(lamb_s(cube_name))
print(x)
print(lamb_inter(lamb_s(cube_name),x))
