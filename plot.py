#!/usr/bin/python
#https://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as plt
import sys
n= sys.argv[1]

image_data = fits.getdata('../image' + n +'.fits', ext=0)

plt.figure()
plt.imshow(image_data)
plt.colorbar()
#plt.show()
plt.savefig('image' + n)
