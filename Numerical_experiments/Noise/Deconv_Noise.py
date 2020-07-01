import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import gaussian as gs
import SLIT
from scipy import signal

def convolve(X,Y):
    return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X*Y))))
def Dconvolve(X,Y):
    return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X/ Y))))

minion = pf.open('minion1.fits')[0].data
minion = minion[::-1,:]
minion = minion[:850,:850]
minion = 255-minion
n1,n2 = np.shape(minion)
kernel = gs.gaussian(n1,n2,n1/2.,n2/2.,1,50/(2.*(2*np.log(2))**0.5),50/(2.*(2*np.log(2))**0.5),0)
kernel = kernel / np.sum(kernel)

fft_kernel = np.fft.fftshift(np.fft.fftn(kernel))
fft_minion = np.fft.fftshift(np.fft.fftn(minion))


Image = convolve(fft_minion , fft_kernel)
fft_Image = np.fft.fftshift(np.fft.fftn(Image))


dec_minion = Dconvolve(fft_Image , fft_kernel)

plt.figure(0)
plt.imshow(minion, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(1)
plt.imshow(Image, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(2)
plt.imshow(dec_minion, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(3)
plt.imshow(kernel, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.show()

####################Now with Noise############################
Noisy = pf.open('stuart.fits')[0].data #(From https://vignette.wikia.nocookie.net/despicableme/images/2/2b/Stuart.png/revision/latest/scale-to-width-down/520?cb=20161108162855)
Noisy = Noisy[::-1,:]
Noisy = Noisy[50:400, 225:575]
Noisy = (255-Noisy)
Noisy = Noisy/np.sum(Noisy)
nn1,nn2 = np.shape(Noisy)
kernel = gs.gaussian(nn1,nn2,nn1/2.,nn2/2.,1,20/(2.*(2*np.log(2))**0.5),20/(2.*(2*np.log(2))**0.5),0)
kernel = kernel / np.sum(kernel**2)


fft_kernel = np.fft.fftshift(np.fft.fftn(kernel))
fft_kernelT = np.fft.fftshift(np.fft.fftn(np.rot90(kernel)))
fft_Noisy = np.fft.fftshift(np.fft.fftn(Noisy))

sigma = 0.00001
Image_Noisy = convolve(fft_Noisy , fft_kernel)+np.random.randn(nn1,nn2)*sigma
fft_Image_Noisy = np.fft.fftshift(np.fft.fftn(Image_Noisy))


dec_minion_Noisy = Dconvolve(fft_Image_Noisy , fft_kernel)

plt.figure(4)
plt.imshow(Noisy, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(5)
plt.imshow(Image_Noisy, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(6)
plt.imshow(dec_minion_Noisy, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.show()



def conv(x):
    x =  np.fft.fftshift(np.fft.fftn(x))
    return convolve(x,fft_kernel)#signal.fftconvolve(x, kernel, mode = 'same')#
def convT(x):
    x =  np.fft.fftshift(np.fft.fftn(x))
    return convolve(x, fft_kernelT)#signal.fftconvolve(x,kernel.T, mode = 'same')#

L = SLIT.Solve.spectralNorm(nn1,nn2,10,10e-10,conv, convT)
niter =200
gamma =0.0001
print(L,gamma)
X = np.copy(Noisy)*0
Xreg = np.copy(X)
for i in range(niter):
    print(i, np.sqrt(np.sum((Image_Noisy-conv( X))**2)))
    X = X+gamma*convT(Image_Noisy-conv(X) )
    Xreg = Xreg+gamma*convT(Image_Noisy-conv(Xreg) )
    Xreg[Xreg<0] = 0
    #Xreg = Xreg-0.1*Xreg




plt.figure(7)
plt.imshow(X, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(8)
plt.imshow(Noisy, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(9)
plt.imshow(Image_Noisy, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.figure(10)
plt.imshow(Xreg, cmap = 'hsv', interpolation = None)#; plt.colorbar()
plt.axis('off')
plt.show()








    
