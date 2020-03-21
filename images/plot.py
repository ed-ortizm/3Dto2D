#!/usr/bin/python
#https://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as plt
import sys
folder = sys.argv[1]

for i in range(9):
    image_data = fits.getdata('./' + folder + '/image_'+ folder +'_'+str(i+1)+ '.fits', ext=0)
    plt.figure()
    plt.imshow(image_data)
    plt.colorbar()
#plt.show()
    plt.savefig('./' + folder + '/image_'+str(i+1))
