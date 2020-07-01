import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt

F475 = pf.open('Horse_475.fits')[1].data
F814 = pf.open('Horse_814.fits')[1].data
F606 = pf.open('Horse_606.fits')[1].data

n1,n2 = np.shape(F475)
print(np.shape(F475), np.shape(F814), np.shape(F606))

cube = np.zeros((3,600,600))
cube[0,:,:] = F814[3000:3600,1800:2400]
cube[1,:,:] = F606[3000:3600,1800:2400]
cube[2,:,:] = F475[3000:3600,1800:2400]

hdus = pf.PrimaryHDU(cube)
lists = pf.HDUList([hdus])
lists.writeto('Horseshoe.fits', clobber=True)